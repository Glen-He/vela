"""结合态 FlexPepDock 局部恢复资格任务的执行、恢复和报告。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    is_current_vela_software,
    is_vela_software_identity,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.typed_data import object_list, object_mapping
from vela.validation.bound_states.controls import (
    PLAN_NAME,
    ControlInput,
    ControlTask,
    build_control_tasks,
    rosetta_parameters,
)
from vela.validation.bound_states.recovery import analyze_control_recovery
from vela.validation.models import LocalRecoveryControl, ValidationError
from vela.validation.records import (
    file_record,
    read_document,
    resume_completed_result,
    validate_record,
)
from vela.validation.refinement.reconstruction import (
    write_chemistry_prepack_protocol,
    write_chemistry_production_refine_protocol,
)
from vela.validation.rosetta import (
    build_chemistry_flexpepdock_command,
    rosetta_crash_log_dir,
    run_rosetta_command,
    verify_rosetta_scripts_tool,
)
from vela.validation.scores import (
    index_rosetta_pdb_outputs,
    read_rosetta_scorefile,
    write_decoy_manifest,
)


@dataclass(frozen=True, slots=True)
class QualificationOutcome:
    """完整资格运行的最终状态和报告位置。"""

    report_path: Path
    qualified: bool


def _prepack_control(
    *,
    config: AppConfig,
    control: LocalRecoveryControl,
    control_input: ControlInput,
    run_dir: Path,
    plan_sha256: str,
) -> Path:
    prepack_dir = run_dir / "controls" / control.control_id / "prepack"
    result_path = prepack_dir / "prepack_result.json"
    resumed = resume_completed_result(
        directory=prepack_dir,
        filename="prepack_result.json",
        document_name="prepack result",
        schema="vela.validation-prepack-result/2",
        identity={"control_id": control.control_id},
        plan_hash_key="qualification_plan_sha256",
        plan_hash=plan_sha256,
        records={
            "output": "prepack",
            "scorefile": "prepack scorefile",
            "protocol": "prepack protocol",
            "log": "prepack log",
        },
        stale_label="prepack result",
    )
    if resumed is not None:
        return resumed.files["output"]
    if prepack_dir.exists():
        raise ValidationError(
            f"incomplete prepack directory requires review: {prepack_dir}"
        )
    prepack_dir.mkdir(parents=True)
    protocol_path = prepack_dir / "prepack.xml"
    write_chemistry_prepack_protocol(
        destination=protocol_path,
        receptor_residue_count=control_input.receptor_residue_count,
        chemistry=control_input.chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    command = build_chemistry_flexpepdock_command(
        settings=config.validation.rosetta,
        input_path=control_input.complex_path,
        protocol_path=protocol_path,
        disulfide_path=control_input.disulfide_path,
        output_dir=prepack_dir,
        seed=control.prepack_seed,
        fixed_histidine_pose_indices=control_input.fixed_histidine_pose_indices,
        nstruct=1,
        scorefile_name="prepack.sc",
        native_path=None,
    )
    log_path = prepack_dir / "prepack.log"
    started_at = utc_now()
    run_rosetta_command(
        command=command,
        log_path=log_path,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        thread_count=1,
    )
    outputs = tuple(prepack_dir.glob("*.pdb"))
    if len(outputs) != 1:
        raise ValidationError(
            f"prepack must produce exactly one PDB: {control.control_id}"
        )
    score_path = prepack_dir / "prepack.sc"
    read_rosetta_scorefile(score_path)
    output = outputs[0]
    atomic_write_json(
        result_path,
        {
            "schema": "vela.validation-prepack-result/2",
            "status": "completed",
            "control_id": control.control_id,
            "qualification_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "command": list(command),
            "output": file_record(output, root=prepack_dir),
            "scorefile": file_record(score_path, root=prepack_dir),
            "protocol": file_record(protocol_path, root=prepack_dir),
            "log": file_record(log_path, root=prepack_dir),
        },
    )
    return output


def _run_seed(
    *,
    config: AppConfig,
    task: ControlTask,
    prepacked_path: Path,
    run_dir: Path,
    plan_sha256: str,
) -> None:
    task_dir = run_dir / "tasks" / task.task_id
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="qualification task result",
        schema="vela.validation-control-task-result/2",
        identity={
            "task_id": task.task_id,
            "control_id": task.control.control_id,
            "batch_id": task.batch_id,
            "seed": task.seed,
        },
        plan_hash_key="qualification_plan_sha256",
        plan_hash=plan_sha256,
        records={
            "scorefile": "task scorefile",
            "decoy_manifest": "task decoy_manifest",
            "protocol": "refinement protocol",
            "log": "task log",
        },
        stale_label="qualification task",
    )
    if resumed is not None:
        return
    if task_dir.exists():
        raise ValidationError(
            f"incomplete qualification task requires review: {task_dir}"
        )
    task_dir.mkdir(parents=True)
    protocol_path = task_dir / "refine.xml"
    write_chemistry_production_refine_protocol(
        destination=protocol_path,
        receptor_residue_count=task.control_input.receptor_residue_count,
        chemistry=task.control_input.chemistry,
        score_function=config.validation.rosetta.score_function,
        random_translation_A=task.control.random_translation_A,
        random_rotation_degrees=task.control.random_rotation_degrees,
        lowres_preoptimize=config.validation.rosetta.lowres_preoptimize,
    )
    command = build_chemistry_flexpepdock_command(
        settings=config.validation.rosetta,
        input_path=prepacked_path,
        protocol_path=protocol_path,
        disulfide_path=task.control_input.disulfide_path,
        output_dir=task_dir,
        seed=task.seed,
        native_path=task.control_input.complex_path,
        fixed_histidine_pose_indices=task.control_input.fixed_histidine_pose_indices,
        nstruct=config.validation.rosetta.decoys_per_seed,
        scorefile_name="refine.sc",
    )
    log_path = task_dir / "refine.log"
    started_at = utc_now()
    run_rosetta_command(
        command=command,
        log_path=log_path,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
    )
    score_path = task_dir / "refine.sc"
    rows = read_rosetta_scorefile(score_path)
    if len(rows) != config.validation.rosetta.decoys_per_seed:
        raise ValidationError(
            f"{task.task_id} produced {len(rows)} score rows; expected "
            f"{config.validation.rosetta.decoys_per_seed}"
        )
    decoy_manifest = write_decoy_manifest(
        rows=rows,
        control=task.control,
        decoy_paths=index_rosetta_pdb_outputs(task_dir),
        task_dir=task_dir,
    )
    atomic_write_json(
        task_dir / "task_result.json",
        {
            "schema": "vela.validation-control-task-result/2",
            "status": "completed",
            "task_id": task.task_id,
            "control_id": task.control.control_id,
            "bound_state_id": task.control.bound_state_id,
            "batch_id": task.batch_id,
            "seed": task.seed,
            "qualification_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "command": list(command),
            "scorefile": file_record(score_path, root=task_dir),
            "decoy_manifest": file_record(decoy_manifest, root=task_dir),
            "protocol": file_record(protocol_path, root=task_dir),
            "log": file_record(log_path, root=task_dir),
        },
    )


def _planned_task_ids(document: dict[str, object]) -> tuple[str, ...]:
    try:
        rows = object_list(document.get("tasks"), name="qualification tasks")
    except TypeError as exc:
        raise ValidationError("qualification plan tasks are invalid") from exc
    task_ids: list[str] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="qualification task")
        except TypeError as exc:
            raise ValidationError("qualification plan task is invalid") from exc
        task_id = row.get("task_id")
        if not isinstance(task_id, str) or row.get("status") != "planned":
            raise ValidationError("qualification tasks must be planned")
        task_ids.append(task_id)
    return tuple(task_ids)


def _verify_plan(
    *,
    config: AppConfig,
    run_dir: Path,
    tasks: tuple[ControlTask, ...],
    require_current_software: bool,
) -> dict[str, object]:
    plan_path = run_dir / PLAN_NAME
    plan = read_document(plan_path, name="qualification plan")
    if (
        plan.get("schema") != "vela.validation-qualification-plan/2"
        or plan.get("stage") != "validation_qualification"
        or plan.get("status") != "planned"
        or not (
            is_current_vela_software(plan.get("software"))
            if require_current_software
            else is_vela_software_identity(plan.get("software"))
        )
        or plan.get("method_id") != config.validation.method_id
    ):
        raise ValidationError("qualification plan identity is invalid")
    inputs = object_mapping(plan.get("inputs"), name="qualification inputs")
    snapshot = object_mapping(inputs.get("config_snapshot"), name="config snapshot")
    snapshot_path, _ = validate_record(
        root=run_dir, raw=snapshot, name="config snapshot"
    )
    if config.source_snapshot_sha256 != sha256_file(snapshot_path):
        raise ValidationError("current project config differs from the frozen plan")
    if plan.get("rosetta_parameters") != rosetta_parameters(config):
        raise ValidationError("current Rosetta parameters differ from the frozen plan")
    if _planned_task_ids(plan) != tuple(task.task_id for task in tasks):
        raise ValidationError("current qualification tasks differ from the frozen plan")
    return plan


def _validate_completed_tasks(
    *, tasks: tuple[ControlTask, ...], run_dir: Path, plan_sha256: str
) -> None:
    """校验分析所依赖的全部历史任务记录和原始产物。"""
    for task in tasks:
        result = resume_completed_result(
            directory=run_dir / "tasks" / task.task_id,
            filename="task_result.json",
            document_name="qualification task result",
            schema="vela.validation-control-task-result/2",
            identity={
                "task_id": task.task_id,
                "control_id": task.control.control_id,
                "batch_id": task.batch_id,
                "seed": task.seed,
            },
            plan_hash_key="qualification_plan_sha256",
            plan_hash=plan_sha256,
            records={
                "scorefile": "task scorefile",
                "decoy_manifest": "task decoy_manifest",
                "protocol": "refinement protocol",
                "log": "task log",
            },
            stale_label="qualification task",
        )
        if result is None:
            raise ValidationError(f"qualification task is incomplete: {task.task_id}")


def _control_report_rows(
    *,
    controls: dict[str, LocalRecoveryControl],
    tasks: tuple[ControlTask, ...],
    run_dir: Path,
    output_dir: Path | None = None,
) -> tuple[bool, list[dict[str, JsonValue]]]:
    rows: list[dict[str, JsonValue]] = []
    qualified = True
    for control_id, control in controls.items():
        control_tasks = tuple(
            task for task in tasks if task.control.control_id == control_id
        )
        passed, row = analyze_control_recovery(
            control=control,
            tasks=control_tasks,
            run_dir=run_dir,
            output_dir=output_dir,
        )
        qualified = qualified and passed
        rows.append(row)
    return qualified, rows


def run_qualification(*, config: AppConfig, run_dir: Path) -> QualificationOutcome:
    """执行或恢复全部已冻结控制任务; 生成不修改项目配置的资格报告。"""
    report_path = run_dir / "qualification_report.json"
    if report_path.exists():
        raise ValidationError(f"qualification report already exists: {report_path}")
    tasks = build_control_tasks(config)
    plan = _verify_plan(
        config=config,
        run_dir=run_dir,
        tasks=tasks,
        require_current_software=True,
    )
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    software = object_mapping(plan.get("software"), name="qualification software")
    if (
        software.get("rosetta_version") != tool.version
        or software.get("rosetta_scripts_sha256") != tool.executable_sha256
    ):
        raise ValidationError("current FlexPepDock differs from the frozen plan")
    plan_path = run_dir / PLAN_NAME
    plan_sha256 = sha256_file(plan_path)
    controls = {
        control.control_id: control for control in config.validation.local_controls
    }
    control_inputs: dict[str, ControlInput] = {}
    for task in tasks:
        control_inputs.setdefault(task.control.control_id, task.control_input)
    prepacked = {
        control_id: _prepack_control(
            config=config,
            control=control,
            control_input=control_inputs[control_id],
            run_dir=run_dir,
            plan_sha256=plan_sha256,
        )
        for control_id, control in controls.items()
    }

    def execute(task: ControlTask) -> None:
        _run_seed(
            config=config,
            task=task,
            prepacked_path=prepacked[task.control.control_id],
            run_dir=run_dir,
            plan_sha256=plan_sha256,
        )

    with ThreadPoolExecutor(
        max_workers=config.validation.rosetta.parallel_tasks
    ) as executor:
        for _ in executor.map(execute, tasks):
            pass
    qualified, control_rows = _control_report_rows(
        controls=controls, tasks=tasks, run_dir=run_dir
    )
    atomic_write_json(
        report_path,
        {
            "schema": "vela.validation-qualification-report/3",
            "stage": "validation_qualification",
            "status": "qualified" if qualified else "unqualified",
            "completed_at": utc_now(),
            "method_id": config.validation.method_id,
            "qualification_plan": {
                "path": plan_path.name,
                "sha256": plan_sha256,
            },
            "software": {
                "rosetta_version": tool.version,
                "rosetta_scripts_sha256": tool.executable_sha256,
            },
            "evidence_category": "method_positive_control",
            "ligand_candidate_evidence": False,
            "controls": control_rows,
        },
    )
    return QualificationOutcome(report_path, qualified)


def analyze_qualification(*, config: AppConfig, run_dir: Path) -> QualificationOutcome:
    """使用当前分析器重析完整历史任务; 同时保留原资格报告。"""
    analysis_dir = run_dir / "qualification_analysis"
    report_path = analysis_dir / "qualification_report.json"
    if analysis_dir.exists():
        raise ValidationError(
            f"qualification analysis directory already exists: {analysis_dir}"
        )
    tasks = build_control_tasks(config)
    plan = _verify_plan(
        config=config,
        run_dir=run_dir,
        tasks=tasks,
        require_current_software=False,
    )
    plan_path = run_dir / PLAN_NAME
    plan_sha256 = sha256_file(plan_path)
    _validate_completed_tasks(tasks=tasks, run_dir=run_dir, plan_sha256=plan_sha256)
    controls = {
        control.control_id: control for control in config.validation.local_controls
    }
    qualified, control_rows = _control_report_rows(
        controls=controls,
        tasks=tasks,
        run_dir=run_dir,
        output_dir=analysis_dir / "controls",
    )
    raw_sampling_software = object_mapping(
        plan.get("software"), name="qualification sampling software"
    )
    sampling_software: dict[str, JsonValue] = {}
    for key in (
        "vela_version",
        "vela_source_sha256",
        "rosetta_version",
        "rosetta_scripts_sha256",
    ):
        value = raw_sampling_software.get(key)
        if not isinstance(value, str) or not value:
            raise ValidationError("qualification sampling software is invalid")
        sampling_software[key] = value
    atomic_write_json(
        report_path,
        {
            "schema": "vela.validation-qualification-report/3",
            "stage": "validation_qualification",
            "status": "qualified" if qualified else "unqualified",
            "completed_at": utc_now(),
            "method_id": config.validation.method_id,
            "qualification_plan": {
                "path": plan_path.relative_to(run_dir).as_posix(),
                "sha256": plan_sha256,
            },
            "sampling_software": sampling_software,
            "analysis_software": vela_software_identity(),
            "evidence_category": "method_positive_control",
            "ligand_candidate_evidence": False,
            "controls": control_rows,
        },
    )
    return QualificationOutcome(report_path, qualified)
