"""实验标准环肽骨架到配体局部精修起点的显式交付。"""

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
    is_current_vela_software,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_list, object_mapping
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import GuidedTemplate, ValidationError
from vela.validation.records import (
    file_record,
    read_document,
    resume_completed_result,
    validate_record,
)
from vela.validation.refinement.reconstruction import (
    validate_flexpepdock_input,
    write_chemistry_protocol,
    write_disulfide_indices,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    run_rosetta_command,
    single_rosetta_pdb_output,
    verify_rosetta_scripts_tool,
)
from vela.validation.scores import read_rosetta_scorefile

GUIDED_PLAN_NAME = "guided_plan.json"
GUIDED_EVIDENCE = "guided_site_compatibility_handoff"
BACKBONE_ATOMS = frozenset({"N", "CA", "C", "O"})


@dataclass(frozen=True, slots=True)
class GuidedTask:
    """一个由配置声明的实验骨架线程化任务。"""

    task_id: str
    template: GuidedTemplate
    receptor_id: str
    target: str
    source_path: Path


@dataclass(frozen=True, slots=True)
class GuidedPlan:
    """已冻结的 guided 起点计划。"""

    run_dir: Path
    tasks: tuple[GuidedTask, ...]


@dataclass(frozen=True, slots=True)
class GuidedOutcome:
    """guided 全原子交付清单。"""

    manifest_path: Path
    task_count: int


def build_guided_tasks(config: AppConfig) -> tuple[GuidedTask, ...]:
    """从项目登记生成任务; 不在代码中指定 PDB 或配体。"""
    states = {item.state_id: item for item in config.validation.bound_states}
    receptors = {item.receptor_id: item for item in config.receptors}
    root = config.paths.data_dir / "validation" / "bound_states"
    tasks: list[GuidedTask] = []
    for index, template in enumerate(config.validation.guided_templates, 1):
        state = states[template.bound_state_id]
        receptor = receptors.get(state.receptor_id)
        if receptor is None:
            raise ValidationError(
                f"guided template receptor is not registered: {state.receptor_id}"
            )
        source = root / state.state_id / "flexpepdock_control.pdb"
        if not source.is_file():
            raise ValidationError(f"guided template asset is missing: {source}")
        tasks.append(
            GuidedTask(
                task_id=f"guided_{index:03d}",
                template=template,
                receptor_id=receptor.receptor_id,
                target=receptor.target,
                source_path=source,
            )
        )
    if not tasks:
        raise ValidationError("no guided templates are configured")
    return tuple(tasks)


def _parameters(config: AppConfig) -> dict[str, JsonValue]:
    return {
        "chemistry_seed": config.validation.handoff.chemistry_seed,
        "score_function": config.validation.rosetta.score_function,
        "sequence": config.chemistry.sequence,
        "disulfide_bonds": [
            [bond.first, bond.second] for bond in config.chemistry.disulfide_bonds
        ],
    }


def write_guided_plan(*, config: AppConfig, run_id: str) -> GuidedPlan:
    """冻结实验骨架、逐位映射和目标配体化学身份。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    tasks = build_guided_tasks(config)
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    run_dir = config.paths.outputs_dir / "validation" / "guided" / run_id
    if run_dir.exists():
        raise ValidationError(f"guided run directory already exists: {run_dir}")
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    atomic_write_json(
        run_dir / GUIDED_PLAN_NAME,
        {
            "schema": "vela.validation-guided-plan/1",
            "stage": "validation_guided_handoff",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": GUIDED_EVIDENCE,
            "known_site_information_used": True,
            "software": {
                **vela_software_identity(),
                "rosetta_version": tool.version,
                "rosetta_scripts_sha256": tool.executable_sha256,
            },
            "config_snapshot": file_record(snapshot, root=run_dir),
            "parameters": _parameters(config),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "template_id": task.template.template_id,
                    "bound_state_id": task.template.bound_state_id,
                    "receptor_id": task.receptor_id,
                    "target": task.target,
                    "ligand_positions": list(task.template.ligand_positions),
                    "selection_reason": task.template.selection_reason,
                    "source": {
                        "path": task.source_path.relative_to(
                            config.paths.data_dir
                        ).as_posix(),
                        "sha256": sha256_file(task.source_path),
                    },
                    "status": "planned",
                }
                for task in tasks
            ],
        },
    )
    return GuidedPlan(run_dir, tasks)


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], tuple[GuidedTask, ...]]:
    plan = read_document(run_dir / GUIDED_PLAN_NAME, name="guided plan")
    if (
        plan.get("schema") != "vela.validation-guided-plan/1"
        or plan.get("stage") != "validation_guided_handoff"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
        or plan.get("evidence_category") != GUIDED_EVIDENCE
        or plan.get("known_site_information_used") is not True
        or plan.get("parameters") != _parameters(config)
    ):
        raise ValidationError("guided plan identity or parameters are invalid")
    snapshot, _ = validate_record(
        root=run_dir,
        raw=plan.get("config_snapshot"),
        name="guided config snapshot",
    )
    if sha256_file(snapshot) != config.source_snapshot_sha256:
        raise ValidationError("current project config differs from the guided plan")
    tasks = build_guided_tasks(config)
    try:
        rows = object_list(plan.get("tasks"), name="guided tasks")
    except TypeError as exc:
        raise ValidationError("guided plan tasks are invalid") from exc
    if len(rows) != len(tasks):
        raise ValidationError("guided task count differs from the frozen plan")
    for task, raw in zip(tasks, rows, strict=True):
        try:
            row = object_mapping(raw, name="guided task")
            source = object_mapping(row.get("source"), name="guided task source")
        except TypeError as exc:
            raise ValidationError("guided task record is invalid") from exc
        expected_path = task.source_path.relative_to(config.paths.data_dir).as_posix()
        if (
            row.get("task_id") != task.task_id
            or row.get("template_id") != task.template.template_id
            or row.get("status") != "planned"
            or source.get("path") != expected_path
            or source.get("sha256") != sha256_file(task.source_path)
        ):
            raise ValidationError("guided task differs from the frozen plan")
    verify_rosetta_scripts_tool(config.validation.rosetta)
    return plan, tasks


def _amino_acids(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(
        residue
        for residue in chain
        if gemmi.find_tabulated_residue(residue.name).is_amino_acid()
    )


def _threaded_residue(
    *, source: gemmi.Residue, target_code: str, target_position: int
) -> gemmi.Residue:
    target_name = gemmi.expand_one_letter(target_code, gemmi.ResidueKind.AA)
    if source.name == target_name:
        residue = source.clone()
    else:
        residue = gemmi.Residue()
        residue.name = target_name
        residue.entity_type = gemmi.EntityType.Polymer
        for atom in source:
            if atom.name in BACKBONE_ATOMS:
                residue.add_atom(atom.clone())
    residue.seqid = gemmi.SeqId(target_position, " ")
    return residue


def _write_threaded_template(
    *, config: AppConfig, task: GuidedTask, destination: Path
) -> int:
    try:
        structure = gemmi.read_structure(str(task.source_path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid guided template: {task.source_path}") from exc
    if len(structure) != 1:
        raise ValidationError("guided template must contain one model")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError("guided template must contain only A/P chains")
    receptor = chains[RECEPTOR_CHAIN].clone()
    ligand = _amino_acids(chains[PEPTIDE_CHAIN])
    peptide = gemmi.Chain(PEPTIDE_CHAIN)
    for target_position, (target_code, source_position) in enumerate(
        zip(
            config.chemistry.sequence,
            task.template.ligand_positions,
            strict=True,
        ),
        1,
    ):
        try:
            source = ligand[source_position - 1]
        except IndexError as exc:
            raise ValidationError(
                f"guided template position is outside ligand: {source_position}"
            ) from exc
        peptide.add_residue(
            _threaded_residue(
                source=source,
                target_code=target_code,
                target_position=target_position,
            )
        )
    output = gemmi.Structure()
    output.name = task.template.template_id
    output.add_model(gemmi.Model(1))
    output[0].add_chain(receptor)
    output[0].add_chain(peptide)
    output.setup_entities()
    ssbond = "".join(
        f"SSBOND  {index:2d} CYS {PEPTIDE_CHAIN} {bond.first:5d}   "
        f"CYS {PEPTIDE_CHAIN} {bond.second:5d}\n"
        for index, bond in enumerate(config.chemistry.disulfide_bonds, 1)
    )
    atomic_write_text(destination, ssbond + output.make_pdb_string())
    return len(_amino_acids(receptor))


def _run_task(
    *, config: AppConfig, task: GuidedTask, run_dir: Path, plan_hash: str
) -> tuple[Path, Path]:
    task_dir = run_dir / "tasks" / task.task_id
    result_path = task_dir / "task_result.json"
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="guided task result",
        schema="vela.validation-guided-task-result/1",
        identity={"task_id": task.task_id},
        plan_hash_key="guided_plan_sha256",
        plan_hash=plan_hash,
        records={"flexpepdock_input": "guided FlexPepDock input"},
        stale_label="guided task result",
    )
    if resumed is not None:
        return resumed.path, resumed.files["flexpepdock_input"]
    if task_dir.exists():
        raise ValidationError(f"incomplete guided task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    threaded = task_dir / "threaded_template.pdb"
    receptor_count = _write_threaded_template(
        config=config, task=task, destination=threaded
    )
    disulfide = task_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide,
        receptor_residue_count=receptor_count,
        chemistry=config.chemistry,
    )
    protocol = task_dir / "restore_chemistry.xml"
    write_chemistry_protocol(
        destination=protocol,
        receptor_residue_count=receptor_count,
        chemistry=config.chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    chemistry_dir = task_dir / "chemistry"
    chemistry_dir.mkdir()
    command = build_chemistry_command(
        settings=config.validation.rosetta,
        input_path=threaded,
        protocol_path=protocol,
        disulfide_path=disulfide,
        output_dir=chemistry_dir,
        seed=config.validation.handoff.chemistry_seed,
    )
    log = task_dir / "restore_chemistry.log"
    run_rosetta_command(command=command, log_path=log)
    read_rosetta_scorefile(chemistry_dir / "chemistry.sc")
    source_output = single_rosetta_pdb_output(chemistry_dir)
    output = task_dir / "flexpepdock_input.pdb"
    atomic_write_text(output, source_output.read_text(encoding="utf-8"))
    confirmed_count, histidines = validate_flexpepdock_input(
        path=output,
        chemistry=config.chemistry,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    if confirmed_count != receptor_count:
        raise ValidationError("guided receptor numbering changed during handoff")
    atomic_write_json(
        result_path,
        {
            "schema": "vela.validation-guided-task-result/1",
            "status": "completed",
            "task_id": task.task_id,
            "template_id": task.template.template_id,
            "bound_state_id": task.template.bound_state_id,
            "receptor_id": task.receptor_id,
            "target": task.target,
            "guided_plan_sha256": plan_hash,
            "completed_at": utc_now(),
            "command": list(command),
            "qc": {
                "fixed_histidine_pose_indices": list(histidines),
                "mapped_ligand_positions": list(task.template.ligand_positions),
            },
            "threaded_template": file_record(threaded, root=task_dir),
            "restore_chemistry_protocol": file_record(protocol, root=task_dir),
            "fix_disulfide": file_record(disulfide, root=task_dir),
            "log": file_record(log, root=task_dir),
            "flexpepdock_input": file_record(output, root=task_dir),
        },
    )
    return result_path, output


def run_guided_handoff(*, config: AppConfig, run_dir: Path) -> GuidedOutcome:
    """执行或恢复全部实验骨架线程化和配体化学恢复任务。"""
    manifest_path = run_dir / "guided_manifest.json"
    if manifest_path.exists():
        raise ValidationError(f"guided manifest already exists: {manifest_path}")
    _, tasks = _verify_plan(config=config, run_dir=run_dir)
    plan_path = run_dir / GUIDED_PLAN_NAME
    plan_hash = sha256_file(plan_path)
    results = tuple(
        _run_task(config=config, task=task, run_dir=run_dir, plan_hash=plan_hash)
        for task in tasks
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-guided-manifest/1",
            "stage": "validation_guided_handoff",
            "status": "completed",
            "completed_at": utc_now(),
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": GUIDED_EVIDENCE,
            "known_site_information_used": True,
            "guided_plan": file_record(plan_path, root=run_dir),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "candidate_id": f"guided__{task.template.template_id}",
                    "receptor_site_id": (
                        f"experimental__{task.template.bound_state_id}"
                    ),
                    "pose_id": task.template.template_id,
                    "receptor_id": task.receptor_id,
                    "target": task.target,
                    "source_seed": None,
                    "task_result": file_record(result, root=run_dir),
                    "flexpepdock_input": file_record(output, root=run_dir),
                }
                for task, (result, output) in zip(tasks, results, strict=True)
            ],
        },
    )
    return GuidedOutcome(manifest_path, len(tasks))
