"""阶段三 candidate pose 全原子交接的执行与恢复。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_list, object_mapping
from vela.preparation.chemistry import ChemistryDefinition
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
    validate_flexpepdock_input,
    verify_cg2all_tool,
    write_chemistry_refine_protocol,
)
from vela.validation.refinement.topology_recovery import (
    assess_recovered_topology,
    recover_topology,
    topology_assessment_record,
)
from vela.validation.rosetta import (
    build_chemistry_refine_command,
    rosetta_crash_log_dir,
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
    execution_status: str
    reconstruction_status: str
    flexpepdock_input_path: Path | None


@dataclass(frozen=True, slots=True)
class HandoffOutcome:
    """完整 candidate 交接的最终 manifest。"""

    manifest_path: Path
    task_count: int


def handoff_result_record(
    *, result: HandoffTaskResult, run_dir: Path
) -> dict[str, JsonValue]:
    """建立正式和开发性交接清单共用的任务记录。"""
    input_record = (
        file_record(result.flexpepdock_input_path, root=run_dir)
        if result.flexpepdock_input_path is not None
        else None
    )
    return {
        "task_id": result.task.task_id,
        "candidate_id": result.task.candidate_id,
        "receptor_site_id": result.task.receptor_site_id,
        "pose_id": result.task.pose.pose_id,
        "receptor_id": result.task.pose.receptor_id,
        "target": result.task.pose.target,
        "source_seed": result.task.pose.seed,
        "source_model_index": result.task.pose.model_index,
        "execution_status": result.execution_status,
        "reconstruction_status": result.reconstruction_status,
        "task_result": file_record(result.result_path, root=run_dir),
        "flexpepdock_input": input_record,
    }


def _terminal_chemistry_records(
    *,
    task_dir: Path,
    protocol_path: Path,
    chemistry_log_path: Path,
) -> dict[str, JsonValue]:
    """记录端基恢复链已经真实生成的文件。"""

    def optional_record(path: Path) -> dict[str, JsonValue] | None:
        return file_record(path, root=task_dir) if path.is_file() else None

    return {
        "restore_chemistry_protocol": file_record(protocol_path, root=task_dir),
        "restore_chemistry_log": optional_record(chemistry_log_path),
    }


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
        plan.get("schema") != "vela.validation-handoff-plan/6"
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
        schema="vela.validation-handoff-task-result/5",
        identity={"task_id": task.task_id},
        plan_hash_key="handoff_plan_sha256",
        plan_hash=plan_sha256,
        records={},
        stale_label="handoff task result",
    )
    if result is None:
        return None
    execution_status = result.document.get("execution_status")
    status = result.document.get("reconstruction_status")
    if execution_status == "completed" and status == "passed":
        reconstruction_status = "passed"
    elif execution_status == "completed" and status == "failed":
        reconstruction_status = "failed"
    elif execution_status == "invalid" and status == "not_assessed":
        reconstruction_status = "not_assessed"
    else:
        raise ValidationError(f"invalid handoff task status: {task.task_id}")
    input_path: Path | None = None
    if reconstruction_status == "passed":
        input_path = validate_record(
            root=task_dir,
            raw=result.document.get("flexpepdock_input"),
            name="FlexPepDock input",
        )[0]
    return HandoffTaskResult(
        task,
        result.path,
        str(execution_status),
        reconstruction_status,
        input_path,
    )


def execute_handoff_task(
    *,
    config: AppConfig,
    chemistry: ChemistryDefinition,
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
    recovery = recover_topology(
        config=config,
        source_path=task.pose.model_path,
        model_index=task.pose.model_index,
        reference_receptor_path=task.reference_receptor_path,
        chemistry=chemistry,
        task_dir=task_dir,
        rosetta_seed=config.validation.handoff.chemistry_seed,
    )
    disulfide_path = task_dir / "fix_disulfide.txt"
    commands = dict(recovery.commands)
    result_path = task_dir / "task_result.json"
    common: dict[str, JsonValue] = {
        "schema": "vela.validation-handoff-task-result/5",
        "status": "completed",
        "task_id": task.task_id,
        "candidate_id": task.candidate_id,
        "receptor_site_id": task.receptor_site_id,
        "pose_id": task.pose.pose_id,
        "receptor_id": task.pose.receptor_id,
        "target": task.pose.target,
        "source_seed": task.pose.seed,
        "chemistry_id": chemistry.chemistry_id,
        "handoff_plan_sha256": plan_sha256,
        "started_at": started_at,
        "completed_at": utc_now(),
        "execution_status": recovery.execution_status,
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
        },
        "topology_recovery": {
            "execution_status": recovery.execution_status,
            "reconstruction_status": recovery.reconstruction_status,
            "failure_reasons": list(recovery.failure_reasons),
            "failure_detail": recovery.failure_detail,
            "metrics": recovery.metrics,
            "artifacts": recovery.artifacts,
        },
    }
    if (
        recovery.execution_status != "completed"
        or recovery.reconstruction_status != "passed"
        or recovery.final_path is None
    ):
        atomic_write_json(
            result_path,
            {
                **common,
                "reconstruction_status": recovery.reconstruction_status,
                "failure_reasons": list(recovery.failure_reasons),
                "commands": commands,
                "flexpepdock_input": None,
            },
        )
        return HandoffTaskResult(
            task,
            result_path,
            recovery.execution_status,
            recovery.reconstruction_status,
            None,
        )
    protocol_path = task_dir / "restore_chemistry.xml"
    write_chemistry_refine_protocol(
        destination=protocol_path,
        receptor_residue_count=recovery.receptor_residue_count,
        chemistry=chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    chemistry_dir = task_dir / "chemistry"
    chemistry_dir.mkdir()
    calibration = config.discovery.qualification.topology_calibration
    chemistry_command = build_chemistry_refine_command(
        settings=config.validation.rosetta,
        input_path=recovery.final_path,
        protocol_path=protocol_path,
        disulfide_path=disulfide_path,
        output_dir=chemistry_dir,
        seed=config.validation.handoff.chemistry_seed,
        fixed_histidine_pose_indices=recovery.fixed_histidine_pose_indices,
        site_constraint_path=task_dir / "site_coordinate_constraints.cst",
        site_constraint_weight=calibration.site_coordinate_constraint_weight,
    )
    commands["restore_and_refine_terminal_chemistry"] = list(chemistry_command)
    chemistry_log_path = task_dir / "restore_chemistry.log"
    final_path = task_dir / "flexpepdock_input.pdb"
    try:
        run_rosetta_command(
            command=chemistry_command,
            log_path=chemistry_log_path,
            crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        )
        chemistry_scores = read_rosetta_scorefile(chemistry_dir / "chemistry.sc")
        if len(chemistry_scores) != 1:
            raise ValidationError(
                "terminal chemistry restore must produce exactly one score row"
            )
        rosetta_output = single_rosetta_pdb_output(chemistry_dir)
        atomic_write_text(final_path, rosetta_output.read_text(encoding="utf-8"))
    except ValidationError as exc:
        atomic_write_json(
            result_path,
            {
                **common,
                "completed_at": utc_now(),
                "execution_status": "invalid",
                "reconstruction_status": "not_assessed",
                "failure_reasons": ["terminal_chemistry_execution_invalid"],
                "failure_detail": str(exc),
                "commands": commands,
                "terminal_chemistry": _terminal_chemistry_records(
                    task_dir=task_dir,
                    protocol_path=protocol_path,
                    chemistry_log_path=chemistry_log_path,
                ),
                "flexpepdock_input": None,
            },
        )
        return HandoffTaskResult(task, result_path, "invalid", "not_assessed", None)
    try:
        receptor_count, fixed_histidines = validate_flexpepdock_input(
            path=final_path,
            chemistry=chemistry,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        )
        final_assessment = assess_recovered_topology(
            config=config,
            input_path=recovery.cg_input_path,
            output_path=final_path,
            chemistry=chemistry,
        )
    except ValidationError as exc:
        atomic_write_json(
            result_path,
            {
                **common,
                "completed_at": utc_now(),
                "execution_status": "completed",
                "reconstruction_status": "failed",
                "failure_reasons": ["terminal_chemistry_qc_failed"],
                "failure_detail": str(exc),
                "commands": commands,
                "terminal_chemistry": _terminal_chemistry_records(
                    task_dir=task_dir,
                    protocol_path=protocol_path,
                    chemistry_log_path=chemistry_log_path,
                ),
                "flexpepdock_input": None,
            },
        )
        return HandoffTaskResult(task, result_path, "completed", "failed", None)
    if final_assessment.failures:
        atomic_write_json(
            result_path,
            {
                **common,
                "completed_at": utc_now(),
                "execution_status": "completed",
                "reconstruction_status": "failed",
                "failure_reasons": list(final_assessment.failures),
                "failure_detail": None,
                "commands": commands,
                "final_topology_qc": topology_assessment_record(final_assessment),
                "terminal_chemistry": _terminal_chemistry_records(
                    task_dir=task_dir,
                    protocol_path=protocol_path,
                    chemistry_log_path=chemistry_log_path,
                ),
                "flexpepdock_input": None,
            },
        )
        return HandoffTaskResult(task, result_path, "completed", "failed", None)
    if (
        receptor_count != recovery.receptor_residue_count
        or fixed_histidines != recovery.fixed_histidine_pose_indices
    ):
        raise ValidationError("FlexPepDock input pose numbering changed during handoff")
    atomic_write_json(
        result_path,
        {
            **common,
            "completed_at": utc_now(),
            "execution_status": "completed",
            "reconstruction_status": "passed",
            "failure_reasons": [],
            "commands": commands,
            "qc": {
                "fixed_histidine_pose_indices": list(fixed_histidines),
                "n_terminus": chemistry.n_terminus,
                "c_terminus": chemistry.c_terminus,
                "disulfide_bonds": [
                    [bond.first, bond.second] for bond in chemistry.disulfide_bonds
                ],
            },
            "final_topology_qc": topology_assessment_record(final_assessment),
            "terminal_chemistry": _terminal_chemistry_records(
                task_dir=task_dir,
                protocol_path=protocol_path,
                chemistry_log_path=chemistry_log_path,
            ),
            "flexpepdock_input": file_record(final_path, root=task_dir),
        },
    )
    return HandoffTaskResult(task, result_path, "completed", "passed", final_path)


def execute_handoff_tasks(
    *,
    config: AppConfig,
    chemistry: ChemistryDefinition,
    tasks: tuple[CandidateHandoffTask, ...],
    run_dir: Path,
    plan_sha256: str,
) -> tuple[HandoffTaskResult, ...]:
    """按配置的 worker 上限并行执行相互独立的交接任务。"""
    worker_count = min(config.validation.rosetta.parallel_tasks, len(tasks))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = tuple(
            executor.submit(
                execute_handoff_task,
                config=config,
                chemistry=chemistry,
                task=task,
                run_dir=run_dir,
                plan_sha256=plan_sha256,
            )
            for task in tasks
        )
        return tuple(future.result() for future in futures)


def run_handoff(*, config: AppConfig, run_dir: Path) -> HandoffOutcome:
    """执行或恢复全部冻结 candidate 重建; 写出唯一交接 manifest。"""
    manifest_path = run_dir / "handoff_manifest.json"
    if manifest_path.exists():
        raise ValidationError(f"handoff manifest already exists: {manifest_path}")
    _, tasks = _verify_plan(config=config, run_dir=run_dir)
    plan_path = run_dir / HANDOFF_PLAN_NAME
    plan_sha256 = sha256_file(plan_path)
    results = execute_handoff_tasks(
        config=config,
        chemistry=config.chemistry,
        tasks=tasks,
        run_dir=run_dir,
        plan_sha256=plan_sha256,
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-handoff-manifest/7",
            "stage": "validation_candidate_handoff",
            "status": (
                "invalid"
                if any(result.execution_status == "invalid" for result in results)
                else "completed"
            ),
            "completed_at": utc_now(),
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "main_discovery_handoff",
            "known_site_information_used": False,
            "parallel_tasks": min(config.validation.rosetta.parallel_tasks, len(tasks)),
            "handoff_plan": {
                "path": plan_path.name,
                "sha256": plan_sha256,
            },
            "passed_task_count": sum(
                result.reconstruction_status == "passed" for result in results
            ),
            "failed_task_count": sum(
                result.reconstruction_status == "failed" for result in results
            ),
            "invalid_task_count": sum(
                result.execution_status == "invalid" for result in results
            ),
            "tasks": [
                handoff_result_record(result=result, run_dir=run_dir)
                for result in results
            ],
        },
    )
    return HandoffOutcome(manifest_path, len(results))
