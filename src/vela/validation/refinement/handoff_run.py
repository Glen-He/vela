"""阶段三 candidate pose 全原子交接的执行与恢复。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_list, object_mapping
from vela.validation.models import ValidationError
from vela.validation.records import (
    file_record,
    read_document,
    resume_completed_result,
    validate_record,
)
from vela.validation.refinement.handoff_plan import (
    HANDOFF_PLAN_NAME,
    CandidateHandoffTask,
    build_handoff_tasks,
    handoff_task_records,
)
from vela.validation.refinement.reconstruction import (
    assess_reconstruction,
    build_cg2all_command,
    validate_flexpepdock_input,
    verify_cg2all_tool,
    write_cg2all_input,
    write_chemistry_protocol,
    write_disulfide_indices,
    write_reference_receptor_complex,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    run_rosetta_command,
    single_rosetta_pdb_output,
    verify_rosetta_scripts_tool,
)
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class HandoffTaskResult:
    """一个已经通过重建和化学核对的 FlexPepDock 起点。"""

    task: CandidateHandoffTask
    result_path: Path
    flexpepdock_input_path: Path


@dataclass(frozen=True, slots=True)
class HandoffOutcome:
    """完整 candidate 交接的最终 manifest。"""

    manifest_path: Path
    task_count: int


def _candidate_ids(plan: dict[str, object]) -> tuple[str, ...]:
    try:
        selection = object_mapping(plan.get("selection"), name="handoff selection")
        values = object_list(
            selection.get("requested_candidate_ids"), name="requested candidate IDs"
        )
    except TypeError as exc:
        raise ValidationError("handoff requested candidate IDs are invalid") from exc
    if any(not isinstance(value, str) for value in values):
        raise ValidationError("handoff requested candidate IDs must be strings")
    return tuple(value for value in values if isinstance(value, str))


def _discovery_run_dir(*, config: AppConfig, plan: dict[str, object]) -> Path:
    try:
        inputs = object_mapping(plan.get("inputs"), name="handoff inputs")
        discovery = object_mapping(inputs.get("discovery_run"), name="discovery run")
        manifests = object_list(discovery.get("manifests"), name="discovery manifests")
    except TypeError as exc:
        raise ValidationError("handoff discovery input is invalid") from exc
    relative = discovery.get("path")
    if not isinstance(relative, str):
        raise ValidationError("handoff discovery run path is invalid")
    path = (config.paths.outputs_dir / relative).resolve()
    if not path.is_relative_to((config.paths.outputs_dir / "runs").resolve()):
        raise ValidationError("handoff source is outside discovery runs")
    for raw in manifests:
        try:
            record = object_mapping(raw, name="discovery manifest")
        except TypeError as exc:
            raise ValidationError("handoff discovery manifest is invalid") from exc
        validate_record(root=path, raw=record, name="discovery manifest")
    return path


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], tuple[CandidateHandoffTask, ...]]:
    plan = read_document(run_dir / HANDOFF_PLAN_NAME, name="handoff plan")
    if (
        plan.get("schema") != "vela.validation-handoff-plan/3"
        or plan.get("stage") != "validation_candidate_handoff"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
    ):
        raise ValidationError("handoff plan identity is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="handoff inputs")
        snapshot = object_mapping(inputs.get("config_snapshot"), name="config snapshot")
        recorded_tasks = object_list(plan.get("tasks"), name="handoff tasks")
    except TypeError as exc:
        raise ValidationError("handoff plan structure is invalid") from exc
    snapshot_path, _ = validate_record(
        root=run_dir, raw=snapshot, name="config snapshot"
    )
    if sha256_file(snapshot_path) != config.source_snapshot_sha256:
        raise ValidationError("current project config differs from the handoff plan")
    discovery_run = _discovery_run_dir(config=config, plan=plan)
    tasks = build_handoff_tasks(
        config=config,
        discovery_run_dir=discovery_run,
        candidate_ids=_candidate_ids(plan),
    )
    for raw in recorded_tasks:
        try:
            row = object_mapping(raw, name="handoff task")
        except TypeError as exc:
            raise ValidationError("handoff task is invalid") from exc
        task_id = row.get("task_id")
        if not isinstance(task_id, str) or row.get("status") != "planned":
            raise ValidationError("handoff tasks must be planned and well-formed")
    expected_tasks = handoff_task_records(
        tasks=tasks,
        discovery_run_dir=discovery_run,
        data_dir=config.paths.data_dir,
    )
    if recorded_tasks != expected_tasks:
        raise ValidationError("current handoff tasks differ from the frozen plan")
    verify_cg2all_tool(config.validation.cg2all)
    verify_rosetta_scripts_tool(config.validation.rosetta)
    return plan, tasks


def _resume_task(
    *, task: CandidateHandoffTask, task_dir: Path, plan_sha256: str
) -> HandoffTaskResult | None:
    result = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="handoff task result",
        schema="vela.validation-handoff-task-result/2",
        identity={"task_id": task.task_id},
        plan_hash_key="handoff_plan_sha256",
        plan_hash=plan_sha256,
        records={"flexpepdock_input": "FlexPepDock input"},
        stale_label="handoff task result",
    )
    if result is None:
        return None
    return HandoffTaskResult(task, result.path, result.files["flexpepdock_input"])


def _run_task(
    *,
    config: AppConfig,
    task: CandidateHandoffTask,
    run_dir: Path,
    plan_sha256: str,
) -> HandoffTaskResult:
    task_dir = run_dir / "tasks" / task.task_id
    resumed = _resume_task(task=task, task_dir=task_dir, plan_sha256=plan_sha256)
    if resumed is not None:
        return resumed
    if task_dir.exists():
        raise ValidationError(f"incomplete handoff task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    started_at = utc_now()
    cg_input_path = task_dir / "cg2all_input.pdb"
    cg_input = write_cg2all_input(
        source_path=task.pose.model_path,
        model_index=task.pose.model_index,
        destination=cg_input_path,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
    )
    cg_output_path = task_dir / "cg2all_raw.pdb"
    cg_command = build_cg2all_command(
        settings=config.validation.cg2all,
        input_path=cg_input_path,
        output_path=cg_output_path,
    )
    cg_log_path = task_dir / "cg2all.log"
    run_rosetta_command(
        command=cg_command,
        log_path=cg_log_path,
        thread_count=config.validation.cg2all.processes,
    )
    grafted_path = task_dir / "all_atom_reference_receptor.pdb"
    write_reference_receptor_complex(
        coarse_pose_path=cg_input_path,
        reconstructed_path=cg_output_path,
        reference_receptor_path=task.reference_receptor_path,
        destination=grafted_path,
        chemistry=config.chemistry,
    )
    metrics = assess_reconstruction(
        input_path=cg_input_path,
        output_path=grafted_path,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
    )
    disulfide_path = task_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=cg_input.receptor_residue_count,
        chemistry=config.chemistry,
    )
    protocol_path = task_dir / "restore_chemistry.xml"
    write_chemistry_protocol(
        destination=protocol_path,
        receptor_residue_count=cg_input.receptor_residue_count,
        chemistry=config.chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    chemistry_dir = task_dir / "chemistry"
    chemistry_dir.mkdir()
    chemistry_command = build_chemistry_command(
        settings=config.validation.rosetta,
        input_path=grafted_path,
        protocol_path=protocol_path,
        disulfide_path=disulfide_path,
        output_dir=chemistry_dir,
        seed=config.validation.handoff.chemistry_seed,
    )
    chemistry_log_path = task_dir / "restore_chemistry.log"
    run_rosetta_command(command=chemistry_command, log_path=chemistry_log_path)
    read_rosetta_scorefile(chemistry_dir / "chemistry.sc")
    rosetta_output = single_rosetta_pdb_output(chemistry_dir)
    final_path = task_dir / "flexpepdock_input.pdb"
    atomic_write_text(final_path, rosetta_output.read_text(encoding="utf-8"))
    receptor_count, fixed_histidines = validate_flexpepdock_input(
        path=final_path,
        chemistry=config.chemistry,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    if receptor_count != cg_input.receptor_residue_count or fixed_histidines != (
        cg_input.fixed_histidine_pose_indices
    ):
        raise ValidationError("FlexPepDock input pose numbering changed during handoff")
    result_path = task_dir / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": "vela.validation-handoff-task-result/2",
            "status": "completed",
            "task_id": task.task_id,
            "candidate_id": task.candidate_id,
            "receptor_site_id": task.receptor_site_id,
            "pose_id": task.pose.pose_id,
            "receptor_id": task.pose.receptor_id,
            "target": task.pose.target,
            "source_seed": task.pose.seed,
            "handoff_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "commands": {
                "cg2all": list(cg_command),
                "restore_chemistry": list(chemistry_command),
            },
            "qc": {
                "receptor_ca_rmsd_A": round(metrics.receptor_ca_rmsd_A, 6),
                "peptide_ca_rmsd_A": round(metrics.peptide_ca_rmsd_A, 6),
                "fixed_histidine_pose_indices": list(fixed_histidines),
                "n_terminus": config.chemistry.n_terminus,
                "c_terminus": config.chemistry.c_terminus,
                "disulfide_bonds": [
                    [bond.first, bond.second]
                    for bond in config.chemistry.disulfide_bonds
                ],
            },
            "inputs": {
                "source_model": {
                    "path": task.pose.model_path.as_posix(),
                    "sha256": task.pose.model_sha256,
                    "model_index": task.pose.model_index,
                },
                "reference_receptor": {
                    "path": task.reference_receptor_path.as_posix(),
                    "sha256": task.reference_receptor_sha256,
                },
                "cg2all_input": file_record(cg_input_path, root=task_dir),
                "restore_chemistry_protocol": file_record(protocol_path, root=task_dir),
                "fix_disulfide": file_record(disulfide_path, root=task_dir),
            },
            "intermediates": {
                "cg2all_raw": file_record(cg_output_path, root=task_dir),
                "all_atom_reference_receptor": file_record(grafted_path, root=task_dir),
            },
            "logs": {
                "cg2all": file_record(cg_log_path, root=task_dir),
                "restore_chemistry": file_record(chemistry_log_path, root=task_dir),
            },
            "flexpepdock_input": file_record(final_path, root=task_dir),
        },
    )
    return HandoffTaskResult(task, result_path, final_path)


def run_handoff(*, config: AppConfig, run_dir: Path) -> HandoffOutcome:
    """执行或恢复全部冻结 candidate 重建; 写出唯一交接 manifest。"""
    manifest_path = run_dir / "handoff_manifest.json"
    if manifest_path.exists():
        raise ValidationError(f"handoff manifest already exists: {manifest_path}")
    _, tasks = _verify_plan(config=config, run_dir=run_dir)
    plan_path = run_dir / HANDOFF_PLAN_NAME
    plan_sha256 = sha256_file(plan_path)
    results = tuple(
        _run_task(
            config=config,
            task=task,
            run_dir=run_dir,
            plan_sha256=plan_sha256,
        )
        for task in tasks
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-handoff-manifest/3",
            "stage": "validation_candidate_handoff",
            "status": "completed",
            "completed_at": utc_now(),
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "main_discovery_handoff",
            "known_site_information_used": False,
            "handoff_plan": {
                "path": plan_path.name,
                "sha256": plan_sha256,
            },
            "tasks": [
                {
                    "task_id": result.task.task_id,
                    "candidate_id": result.task.candidate_id,
                    "receptor_site_id": result.task.receptor_site_id,
                    "pose_id": result.task.pose.pose_id,
                    "receptor_id": result.task.pose.receptor_id,
                    "target": result.task.pose.target,
                    "source_seed": result.task.pose.seed,
                    "source_model_index": result.task.pose.model_index,
                    "task_result": file_record(result.result_path, root=run_dir),
                    "flexpepdock_input": file_record(
                        result.flexpepdock_input_path, root=run_dir
                    ),
                }
                for result in results
            ],
        },
    )
    return HandoffOutcome(manifest_path, len(results))
