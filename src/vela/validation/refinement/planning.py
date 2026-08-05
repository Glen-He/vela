"""阶段三正式 candidate 局部精修来源与不可变计划。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

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
from vela.validation.models import ValidationError
from vela.validation.readiness import assess_validation_readiness
from vela.validation.records import (
    file_record,
    nonnegative_integer,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.guided import GUIDED_EVIDENCE
from vela.validation.refinement.reconstruction import validate_flexpepdock_input
from vela.validation.rosetta import (
    verify_flexpepdock_tool,
    verify_rosetta_scripts_tool,
)

REFINEMENT_PLAN_NAME = "refinement_plan.json"
BLIND_REFINEMENT_EVIDENCE = "main_discovery_local_refinement"
GUIDED_REFINEMENT_EVIDENCE = "guided_site_compatibility_refinement"


@dataclass(frozen=True, slots=True)
class RefinementSource:
    """局部精修起点清单及其不可改变的证据身份。"""

    kind: str
    category: str
    manifest_path: Path
    evidence_category: str
    known_site_information_used: bool


@dataclass(frozen=True, slots=True)
class RefinementStart:
    """一个通过阶段二来源和全原子化学核对的局部精修起点。"""

    start_id: str
    candidate_id: str
    receptor_site_id: str
    pose_id: str
    receptor_id: str
    target: str
    source_seed: int | None
    input_path: Path
    input_sha256: str
    receptor_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class RefinementTask:
    """一个全原子起点与独立局部精修 seed 的组合。"""

    task_id: str
    start: RefinementStart
    seed: int


@dataclass(frozen=True, slots=True)
class RefinementPlan:
    """冻结在独立目录中的正式 candidate 局部精修计划。"""

    run_dir: Path
    tasks: tuple[RefinementTask, ...]


def _optional_seed(value: object) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValidationError("source seed must be a non-negative integer or null")
    return value


def read_refinement_source(
    *, config: AppConfig, source_run_dir: Path
) -> tuple[RefinementSource, tuple[RefinementStart, ...]]:
    """读取 blind 或 guided 起点; 并由清单决定证据类别。"""
    handoff_path = source_run_dir / "handoff_manifest.json"
    guided_path = source_run_dir / "guided_manifest.json"
    existing = tuple(path for path in (handoff_path, guided_path) if path.is_file())
    if len(existing) != 1:
        raise ValidationError(
            "refinement source must contain exactly one supported completed manifest"
        )
    manifest_path = existing[0]
    manifest = read_document(manifest_path, name="refinement source manifest")
    if manifest_path == handoff_path:
        source = RefinementSource(
            kind="blind_handoff",
            category="handoffs",
            manifest_path=manifest_path,
            evidence_category=BLIND_REFINEMENT_EVIDENCE,
            known_site_information_used=False,
        )
        expected_schema = "vela.validation-handoff-manifest/7"
        source_plan_key = "handoff_plan"
        expected_source_evidence = "main_discovery_handoff"
    else:
        source = RefinementSource(
            kind="guided_handoff",
            category="guided",
            manifest_path=manifest_path,
            evidence_category=GUIDED_REFINEMENT_EVIDENCE,
            known_site_information_used=True,
        )
        expected_schema = "vela.validation-guided-manifest/1"
        source_plan_key = "guided_plan"
        expected_source_evidence = GUIDED_EVIDENCE
    if (
        manifest.get("schema") != expected_schema
        or manifest.get("status") != "completed"
        or manifest.get("chemistry_id") != config.chemistry.chemistry_id
        or manifest.get("known_site_information_used")
        is not source.known_site_information_used
        or manifest.get("evidence_category") != expected_source_evidence
    ):
        raise ValidationError("refinement source manifest identity is invalid")
    try:
        rows = object_list(manifest.get("tasks"), name="handoff manifest tasks")
    except TypeError as exc:
        raise ValidationError("handoff manifest structure is invalid") from exc
    validate_record(
        root=source_run_dir,
        raw=manifest.get(source_plan_key),
        name="refinement source plan",
    )
    starts: list[RefinementStart] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="handoff manifest task")
        except TypeError as exc:
            raise ValidationError("handoff manifest task is invalid") from exc
        start_id = safe_identifier(row.get("task_id"), name="handoff task ID")
        validate_record(
            root=source_run_dir,
            raw=row.get("task_result"),
            name=f"{start_id} task result",
        )
        if source.kind == "blind_handoff":
            execution_status = row.get("execution_status")
            reconstruction_status = row.get("reconstruction_status")
            if execution_status == "invalid":
                raise ValidationError("invalid handoff task cannot be refined")
            if execution_status != "completed":
                raise ValidationError("handoff execution status is invalid")
            if reconstruction_status == "failed":
                if row.get("flexpepdock_input") is not None:
                    raise ValidationError(
                        "failed handoff task must not declare a FlexPepDock input"
                    )
                continue
            if reconstruction_status != "passed":
                raise ValidationError("handoff reconstruction status is invalid")
        input_path, input_hash = validate_record(
            root=source_run_dir,
            raw=row.get("flexpepdock_input"),
            name=f"{start_id} FlexPepDock input",
        )
        receptor_count, histidines = validate_flexpepdock_input(
            path=input_path,
            chemistry=config.chemistry,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        )
        starts.append(
            RefinementStart(
                start_id=start_id,
                candidate_id=safe_identifier(
                    row.get("candidate_id"), name="candidate ID"
                ),
                receptor_site_id=safe_identifier(
                    row.get("receptor_site_id"), name="receptor site ID"
                ),
                pose_id=safe_identifier(row.get("pose_id"), name="pose ID"),
                receptor_id=safe_identifier(row.get("receptor_id"), name="receptor ID"),
                target=safe_identifier(row.get("target"), name="target ID"),
                source_seed=_optional_seed(row.get("source_seed")),
                input_path=input_path,
                input_sha256=input_hash,
                receptor_residue_count=receptor_count,
                fixed_histidine_pose_indices=histidines,
            )
        )
    if not starts or len({item.start_id for item in starts}) != len(starts):
        raise ValidationError("handoff manifest contains no unique starts")
    return source, tuple(starts)


def build_refinement_tasks(
    *, config: AppConfig, source_run_dir: Path
) -> tuple[RefinementTask, ...]:
    """将每个合格起点与全部已冻结局部精修 seed 组合。"""
    _, starts = read_refinement_source(config=config, source_run_dir=source_run_dir)
    tasks: list[RefinementTask] = []
    index = 1
    for start in starts:
        for seed in config.validation.seeds:
            tasks.append(RefinementTask(f"refine_{index:05d}", start, seed))
            index += 1
    if not tasks:
        raise ValidationError("Stage 3 production seeds are not frozen")
    return tuple(tasks)


def refinement_parameters(config: AppConfig) -> dict[str, JsonValue]:
    """返回写入计划并在运行前复核的全部局部精修参数。"""
    settings = config.validation
    return {
        "prepack_seed": settings.refinement.prepack_seed,
        "random_translation_A": settings.refinement.random_translation_A,
        "random_rotation_degrees": settings.refinement.random_rotation_degrees,
        "ranking_score": settings.refinement.ranking_score,
        "seeds": list(settings.seeds),
        "parallel_tasks": settings.rosetta.parallel_tasks,
        "decoys_per_seed": settings.rosetta.decoys_per_seed,
        "score_function": settings.rosetta.score_function,
        "lowres_preoptimize": settings.rosetta.lowres_preoptimize,
        "application": "rosetta_scripts_with_in_pose_chemistry_restoration",
    }


def refinement_identity(plan: dict[str, object]) -> tuple[str, bool]:
    """从已验证计划中收窄证据类别和已知位点标记。"""
    evidence = plan.get("evidence_category")
    known = plan.get("known_site_information_used")
    if not isinstance(evidence, str) or not isinstance(known, bool):
        raise ValidationError("refinement evidence identity is invalid")
    return evidence, known


def _relative_output(path: Path, *, config: AppConfig, name: str) -> str:
    try:
        return path.resolve().relative_to(config.paths.outputs_dir.resolve()).as_posix()
    except ValueError as exc:
        raise ValidationError(
            f"{name} is outside the configured outputs directory"
        ) from exc


def write_refinement_plan(
    *, config: AppConfig, source_run_dir: Path, run_id: str
) -> RefinementPlan:
    """冻结已放行方法、显式起点来源、seed 和局部精修参数。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    readiness = assess_validation_readiness(config)
    if not readiness.production_ready:
        raise ValidationError(
            "Stage 3 candidate refinement is not ready: "
            + "; ".join(issue.code for issue in readiness.issues)
        )
    source, _ = read_refinement_source(config=config, source_run_dir=source_run_dir)
    tasks = build_refinement_tasks(config=config, source_run_dir=source_run_dir)
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    run_dir = config.paths.outputs_dir / "validation" / "refinements" / run_id
    if run_dir.exists():
        raise ValidationError(f"refinement run directory already exists: {run_dir}")
    report = config.validation.qualification_report
    if report is None:
        raise ValidationError("qualified method lacks a qualification report")
    report_relative = _relative_output(
        report, config=config, name="qualification report"
    )
    source_relative = _relative_output(
        source.manifest_path, config=config, name="refinement source manifest"
    )
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    atomic_write_json(
        run_dir / REFINEMENT_PLAN_NAME,
        {
            "schema": "vela.validation-refinement-plan/2",
            "stage": "validation_local_refinement",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": source.evidence_category,
            "known_site_information_used": source.known_site_information_used,
            "software": {
                **vela_software_identity(),
                "rosetta_version": flexpepdock.version,
                "flexpepdock_sha256": flexpepdock.executable_sha256,
                "rosetta_scripts_sha256": scripts.executable_sha256,
            },
            "inputs": {
                "config_snapshot": file_record(snapshot, root=run_dir),
                "qualification_report": {
                    "path": report_relative,
                    "sha256": sha256_file(report),
                },
                "source_manifest": {
                    "kind": source.kind,
                    "category": source.category,
                    "path": source_relative,
                    "sha256": sha256_file(source.manifest_path),
                },
            },
            "parameters": refinement_parameters(config),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "start_id": task.start.start_id,
                    "candidate_id": task.start.candidate_id,
                    "receptor_site_id": task.start.receptor_site_id,
                    "pose_id": task.start.pose_id,
                    "receptor_id": task.start.receptor_id,
                    "target": task.start.target,
                    "source_seed": task.start.source_seed,
                    "refinement_seed": task.seed,
                    "input_sha256": task.start.input_sha256,
                    "status": "planned",
                }
                for task in tasks
            ],
        },
    )
    return RefinementPlan(run_dir, tasks)


def _source_manifest(
    *, config: AppConfig, record: dict[str, object], category: str, name: str
) -> Path:
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError(f"invalid {name} record")
    path = (config.paths.outputs_dir / relative).resolve()
    root = (config.paths.outputs_dir / "validation" / category).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValidationError(f"{name} is outside {category} runs") from exc
    if not path.is_file() or sha256_file(path) != expected_hash:
        raise ValidationError(f"{name} hash mismatch")
    return path


def verify_refinement_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], tuple[RefinementTask, ...]]:
    """复核计划快照、资格报告、起点来源、任务和工具身份。"""
    readiness = assess_validation_readiness(config)
    if not readiness.production_ready:
        raise ValidationError("Stage 3 candidate refinement is no longer ready")
    plan = read_document(run_dir / REFINEMENT_PLAN_NAME, name="refinement plan")
    if (
        plan.get("schema") != "vela.validation-refinement-plan/2"
        or plan.get("stage") != "validation_local_refinement"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("method_id") != config.validation.method_id
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
    ):
        raise ValidationError("refinement plan identity is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="refinement inputs")
        rows = object_list(plan.get("tasks"), name="refinement tasks")
    except TypeError as exc:
        raise ValidationError("refinement plan structure is invalid") from exc
    snapshot_path, _ = validate_record(
        root=run_dir, raw=inputs.get("config_snapshot"), name="config snapshot"
    )
    if sha256_file(snapshot_path) != config.source_snapshot_sha256:
        raise ValidationError("current project config differs from refinement plan")
    try:
        source_record = object_mapping(
            inputs.get("source_manifest"), name="refinement source manifest"
        )
        report_record = object_mapping(
            inputs.get("qualification_report"), name="qualification report"
        )
    except TypeError as exc:
        raise ValidationError("refinement source records are invalid") from exc
    source_kind = source_record.get("kind")
    source_category = source_record.get("category")
    resolved_category: str
    if source_kind == "blind_handoff" and source_category == "handoffs":
        resolved_category = "handoffs"
        expected_evidence = BLIND_REFINEMENT_EVIDENCE
        expected_known = False
    elif source_kind == "guided_handoff" and source_category == "guided":
        resolved_category = "guided"
        expected_evidence = GUIDED_REFINEMENT_EVIDENCE
        expected_known = True
    else:
        raise ValidationError("refinement source identity is invalid")
    if (
        plan.get("evidence_category") != expected_evidence
        or plan.get("known_site_information_used") is not expected_known
    ):
        raise ValidationError("refinement evidence identity is invalid")
    source_manifest = _source_manifest(
        config=config,
        record=source_record,
        category=resolved_category,
        name="refinement source manifest",
    )
    _source_manifest(
        config=config,
        record=report_record,
        category="controls",
        name="qualification report",
    )
    source, _ = read_refinement_source(
        config=config, source_run_dir=source_manifest.parent
    )
    if (
        source.kind != source_kind
        or source.category != resolved_category
        or source.manifest_path != source_manifest
        or source.evidence_category != expected_evidence
        or source.known_site_information_used is not expected_known
    ):
        raise ValidationError("current refinement source differs from the frozen plan")
    tasks = build_refinement_tasks(config=config, source_run_dir=source_manifest.parent)
    recorded: list[tuple[str, str, int]] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="refinement task")
        except TypeError as exc:
            raise ValidationError("refinement task is invalid") from exc
        if row.get("status") != "planned":
            raise ValidationError("refinement task status is invalid")
        recorded.append(
            (
                safe_identifier(row.get("task_id"), name="refinement task ID"),
                safe_identifier(row.get("start_id"), name="refinement start ID"),
                nonnegative_integer(row.get("refinement_seed"), name="refinement seed"),
            )
        )
    expected = tuple((task.task_id, task.start.start_id, task.seed) for task in tasks)
    if tuple(recorded) != expected or plan.get("parameters") != refinement_parameters(
        config
    ):
        raise ValidationError("current refinement tasks differ from frozen plan")
    verify_flexpepdock_tool(config.validation.rosetta)
    verify_rosetta_scripts_tool(config.validation.rosetta)
    return plan, tasks
