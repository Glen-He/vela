"""实验结合态局部恢复控制的输入、任务和不可变计划。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import (
    BoundStateDefinition,
    LocalRecoveryControl,
    ValidationError,
)
from vela.validation.readiness import assess_validation_readiness
from vela.validation.rosetta import verify_flexpepdock_tool

PLAN_NAME = "qualification_plan.json"


@dataclass(frozen=True, slots=True)
class ControlInput:
    """一个经过重新核对的标准环肽控制输入。"""

    definition: BoundStateDefinition
    complex_path: Path
    disulfide_path: Path
    fixed_histidine_pose_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ControlTask:
    """一个控制对象的独立随机恢复任务。"""

    task_id: str
    control: LocalRecoveryControl
    control_input: ControlInput
    batch_id: str
    seed: int


@dataclass(frozen=True, slots=True)
class QualificationPlan:
    """已经写入磁盘的阶段三资格运行计划。"""

    run_id: str
    run_dir: Path
    tasks: tuple[ControlTask, ...]


def _standard_sequence(chain: gemmi.Chain) -> str:
    sequence: list[str] = []
    for residue in chain:
        info = gemmi.find_tabulated_residue(residue.name)
        if not info.is_amino_acid() or len(info.one_letter_code) != 1:
            raise ValidationError(
                f"control chain contains non-standard residue: {residue.name}"
            )
        sequence.append(info.one_letter_code)
    return "".join(sequence)


def inspect_control_input(
    *, config: AppConfig, definition: BoundStateDefinition
) -> ControlInput:
    """从通用结合态登记核对局部控制文件; 不依赖特定 PDB 或残基数。"""
    state_dir = (
        config.paths.data_dir / "validation" / "bound_states" / definition.state_id
    )
    complex_path = state_dir / "flexpepdock_control.pdb"
    disulfide_path = state_dir / "fix_disulfide.txt"
    if not complex_path.is_file() or not disulfide_path.is_file():
        raise ValidationError(
            f"local control assets are missing for {definition.state_id}"
        )
    try:
        structure = gemmi.read_structure(str(complex_path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid local control PDB: {complex_path}") from exc
    if len(structure) != 1:
        raise ValidationError("local control PDB must contain one model")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError(
            "local control PDB must contain only receptor A and peptide P"
        )
    receptor = chains[RECEPTOR_CHAIN]
    peptide = chains[PEPTIDE_CHAIN]
    if definition.ligand_sequence is None:
        raise ValidationError(
            f"local control sequence is unresolved: {definition.state_id}"
        )
    sequence = _standard_sequence(peptide)
    if sequence != definition.ligand_sequence:
        raise ValidationError(
            f"local control sequence mismatch: expected {definition.ligand_sequence}, got {sequence}"
        )
    receptor_count = len(receptor)
    expected_disulfides = "".join(
        f"{receptor_count + bond.first} {receptor_count + bond.second}\n"
        for bond in definition.disulfide_bonds
    )
    if disulfide_path.read_text(encoding="utf-8") != expected_disulfides:
        raise ValidationError(
            f"local control disulfide indices are stale: {definition.state_id}"
        )
    return ControlInput(
        definition=definition,
        complex_path=complex_path,
        disulfide_path=disulfide_path,
        fixed_histidine_pose_indices=tuple(
            receptor_count + item.position for item in definition.histidines
        ),
    )


def build_control_tasks(config: AppConfig) -> tuple[ControlTask, ...]:
    """由配置中的任意标准环肽控制登记展开独立 seed 任务。"""
    readiness = assess_validation_readiness(config)
    if not readiness.setup_ready:
        raise ValidationError(
            "Stage 3 setup is not ready: "
            + "; ".join(f"{item.code}: {item.message}" for item in readiness.issues)
        )
    bound_state_by_id = {
        definition.state_id: definition for definition in config.validation.bound_states
    }
    inputs: dict[str, ControlInput] = {}
    tasks: list[ControlTask] = []
    for control in config.validation.local_controls:
        definition = bound_state_by_id[control.bound_state_id]
        control_input = inputs.setdefault(
            definition.state_id,
            inspect_control_input(config=config, definition=definition),
        )
        for batch_index, batch in enumerate(control.seed_batches, 1):
            batch_id = f"batch_{batch_index:02d}"
            for seed in batch:
                tasks.append(
                    ControlTask(
                        task_id=f"{control.control_id}__{batch_id}__seed_{seed}",
                        control=control,
                        control_input=control_input,
                        batch_id=batch_id,
                        seed=seed,
                    )
                )
    return tuple(tasks)


def rosetta_parameters(config: AppConfig) -> dict[str, JsonValue]:
    """记录资格运行真正使用的 Rosetta 协议参数。"""
    settings = config.validation.rosetta
    return {
        "parallel_tasks": settings.parallel_tasks,
        "decoys_per_seed": settings.decoys_per_seed,
        "score_function": settings.score_function,
        "lowres_preoptimize": settings.lowres_preoptimize,
        "receptor_chain": RECEPTOR_CHAIN,
        "peptide_chain": PEPTIDE_CHAIN,
        "disulfide_source": "explicit_pose_indices",
        "histidine_tautomer_policy": "configured_and_fixed",
    }


def _control_parameters(control: LocalRecoveryControl) -> dict[str, JsonValue]:
    return {
        "control_id": control.control_id,
        "bound_state_id": control.bound_state_id,
        "prepack_seed": control.prepack_seed,
        "seed_batches": [list(batch) for batch in control.seed_batches],
        "random_translation_A": control.random_translation_A,
        "random_rotation_degrees": control.random_rotation_degrees,
        "ranking_score": control.ranking_score,
        "recovery_rmsd_score": control.recovery_rmsd_score,
        "top_clusters": control.top_clusters,
        "max_recovery_rmsd_A": control.max_recovery_rmsd_A,
        "max_cluster_backbone_rmsd_A": control.max_cluster_backbone_rmsd_A,
        "min_cluster_seed_support": control.min_cluster_seed_support,
        "max_batch_pose_rmsd_A": control.max_batch_pose_rmsd_A,
    }


def write_qualification_plan(*, config: AppConfig, run_id: str) -> QualificationPlan:
    """冻结当前控制对象、输入哈希、扰动和验收条件。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    run_dir = config.paths.outputs_dir / "validation" / "controls" / run_id
    if run_dir.exists():
        raise ValidationError(f"qualification run directory already exists: {run_dir}")
    tasks = build_control_tasks(config)
    tool = verify_flexpepdock_tool(config.validation.rosetta)
    snapshot_path = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot_path, config.source_snapshot_text)
    preparation_manifest = (
        config.paths.data_dir
        / "validation"
        / "bound_states"
        / "preparation_manifest.json"
    )
    controls: list[dict[str, JsonValue]] = []
    inputs_by_control: dict[str, ControlInput] = {}
    for task in tasks:
        inputs_by_control.setdefault(task.control.control_id, task.control_input)
    for control in config.validation.local_controls:
        control_input = inputs_by_control[control.control_id]
        controls.append(
            {
                **_control_parameters(control),
                "input": {
                    "complex_path": control_input.complex_path.relative_to(
                        config.paths.data_dir
                    ).as_posix(),
                    "complex_sha256": sha256_file(control_input.complex_path),
                    "disulfide_path": control_input.disulfide_path.relative_to(
                        config.paths.data_dir
                    ).as_posix(),
                    "disulfide_sha256": sha256_file(control_input.disulfide_path),
                    "fixed_histidine_pose_indices": list(
                        control_input.fixed_histidine_pose_indices
                    ),
                },
            }
        )
    atomic_write_json(
        run_dir / PLAN_NAME,
        {
            "schema": "vela.validation-qualification-plan/1",
            "stage": "validation_qualification",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "method_id": config.validation.method_id,
            "evidence_category": "method_positive_control",
            "ligand_candidate_evidence": False,
            "software": {
                **vela_software_identity(),
                "rosetta_version": tool.version,
                "flexpepdock_sha256": tool.executable_sha256,
            },
            "inputs": {
                "config_snapshot": {
                    "path": snapshot_path.name,
                    "sha256": sha256_file(snapshot_path),
                },
                "bound_state_preparation_manifest": {
                    "path": preparation_manifest.relative_to(
                        config.paths.data_dir
                    ).as_posix(),
                    "sha256": sha256_file(preparation_manifest),
                },
            },
            "rosetta_parameters": rosetta_parameters(config),
            "controls": controls,
            "tasks": [
                {
                    "task_id": task.task_id,
                    "control_id": task.control.control_id,
                    "bound_state_id": task.control.bound_state_id,
                    "batch_id": task.batch_id,
                    "seed": task.seed,
                    "status": "planned",
                }
                for task in tasks
            ],
        },
    )
    return QualificationPlan(run_id, run_dir, tasks)
