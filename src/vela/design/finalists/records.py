"""阶段四候选柔性复核计划的读取与完整性校验。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import is_current_vela_software, sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.design.models import (
    DesignError,
    DesignTemplate,
    FinalistStart,
    FinalistTask,
    SequenceCandidate,
)
from vela.design.screening.records import (
    SCREEN_PLAN_NAME,
    candidate_from_record,
    design_parameters,
)
from vela.validation.records import (
    nonnegative_integer,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.rosetta import verify_flexpepdock_tool

FINALIST_PLAN_NAME = "finalist_plan.json"


@dataclass(frozen=True, slots=True)
class FinalistPlan:
    """已冻结的候选、成对起点和 FlexPepDock 任务。"""

    run_dir: Path
    source_screen_dir: Path
    candidates: tuple[SequenceCandidate, ...]
    templates: tuple[DesignTemplate, ...]
    starts: tuple[FinalistStart, ...]
    tasks: tuple[FinalistTask, ...]


def _integers(value: object, *, name: str) -> tuple[int, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise DesignError(f"{name} is invalid") from exc
    result: list[int] = []
    for item in values:
        if not isinstance(item, int) or isinstance(item, bool):
            raise DesignError(f"{name} is invalid")
        result.append(item)
    return tuple(result)


def _source_screen(*, config: AppConfig, source_record: dict[str, object]) -> Path:
    source_path = source_record.get("path")
    if not isinstance(source_path, str):
        raise DesignError("finalist source path is invalid")
    source = (config.paths.outputs_dir / source_path).resolve()
    source_root = (config.paths.outputs_dir / "design" / "screens").resolve()
    try:
        source.relative_to(source_root)
    except ValueError as exc:
        raise DesignError("finalist source path escapes Stage 4 screen runs") from exc
    for filename, key in (
        (SCREEN_PLAN_NAME, "screen_plan_sha256"),
        ("screen_manifest.json", "screen_manifest_sha256"),
        ("screen_analysis/analysis_manifest.json", "analysis_manifest_sha256"),
    ):
        digest = source_record.get(key)
        if not isinstance(digest, str) or sha256_file(source / filename) != digest:
            raise DesignError(f"finalist source changed: {filename}")
    return source


def _qualification_report(*, config: AppConfig, record: dict[str, object]) -> None:
    path = record.get("path")
    digest = record.get("sha256")
    if (
        not isinstance(path, str)
        or not isinstance(digest, str)
        or digest != config.validation.qualification_report_sha256
    ):
        raise DesignError("Stage 3 qualification report record is invalid")
    recorded = (config.paths.outputs_dir / path).resolve()
    try:
        recorded.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DesignError("Stage 3 qualification report escapes outputs") from exc
    if (
        recorded != config.validation.qualification_report
        or not recorded.is_file()
        or sha256_file(recorded) != digest
    ):
        raise DesignError("Stage 3 qualification report changed")


def _templates(*, rows: list[object], run_dir: Path) -> tuple[DesignTemplate, ...]:
    templates: list[DesignTemplate] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="finalist template")
        except TypeError as exc:
            raise DesignError("finalist template is invalid") from exc
        path, digest = validate_record(
            root=run_dir, raw=row.get("structure"), name="finalist template structure"
        )
        role = row.get("evidence_role")
        target = row.get("target")
        if not isinstance(role, str) or not isinstance(target, str):
            raise DesignError("finalist template role or target is invalid")
        templates.append(
            DesignTemplate(
                safe_identifier(row.get("template_id"), name="template ID"),
                role,
                safe_identifier(row.get("cluster_id"), name="cluster ID"),
                safe_identifier(row.get("candidate_id"), name="Stage 3 candidate ID"),
                safe_identifier(row.get("receptor_id"), name="receptor ID"),
                safe_identifier(target, name="target ID"),
                path,
                digest,
                nonnegative_integer(
                    row.get("receptor_residue_count"), name="receptor residue count"
                ),
                _integers(
                    row.get("fixed_histidine_pose_indices"),
                    name="template histidine indices",
                ),
            )
        )
    result = tuple(templates)
    if not result or len({item.template_id for item in result}) != len(result):
        raise DesignError("finalist templates are empty or duplicated")
    return result


def _starts(
    *,
    rows: list[object],
    run_dir: Path,
    candidates: tuple[SequenceCandidate, ...],
    templates: tuple[DesignTemplate, ...],
) -> tuple[FinalistStart, ...]:
    candidate_by_id = {item.candidate_id: item for item in candidates}
    template_by_id = {item.template_id: item for item in templates}
    starts: list[FinalistStart] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="finalist start")
        except TypeError as exc:
            raise DesignError("finalist start is invalid") from exc
        candidate_id = safe_identifier(row.get("candidate_id"), name="candidate ID")
        template_id = safe_identifier(row.get("template_id"), name="template ID")
        state = row.get("state")
        if not isinstance(state, str):
            raise DesignError("finalist start state is invalid")
        try:
            candidate = candidate_by_id[candidate_id]
            template = template_by_id[template_id]
        except KeyError as exc:
            raise DesignError("finalist start references an unknown input") from exc
        path, digest = validate_record(
            root=run_dir, raw=row.get("structure"), name="finalist start structure"
        )
        starts.append(
            FinalistStart(
                safe_identifier(row.get("start_id"), name="start ID"),
                safe_identifier(row.get("pair_id"), name="start pair ID"),
                state,
                candidate,
                template,
                path,
                digest,
                _integers(
                    row.get("histidine_pose_indices"),
                    name="finalist histidine indices",
                ),
            )
        )
    result = tuple(starts)
    expected = len(candidates) * len(templates) * 2
    if len(result) != expected or len({item.start_id for item in result}) != len(
        result
    ):
        raise DesignError("finalist start coverage is incomplete or duplicated")
    return result


def _tasks(
    *, rows: list[object], starts: tuple[FinalistStart, ...], seed_count: int
) -> tuple[FinalistTask, ...]:
    start_by_id = {item.start_id: item for item in starts}
    tasks: list[FinalistTask] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="finalist task")
        except TypeError as exc:
            raise DesignError("finalist task is invalid") from exc
        if row.get("status") != "planned":
            raise DesignError("finalist task status is invalid")
        start_id = safe_identifier(row.get("start_id"), name="start ID")
        try:
            start = start_by_id[start_id]
        except KeyError as exc:
            raise DesignError("finalist task references an unknown start") from exc
        tasks.append(
            FinalistTask(
                safe_identifier(row.get("task_id"), name="finalist task ID"),
                safe_identifier(row.get("pair_id"), name="finalist pair ID"),
                start,
                nonnegative_integer(row.get("seed"), name="finalist seed"),
            )
        )
    result = tuple(tasks)
    expected = len(starts) * seed_count
    if len(result) != expected or len({item.task_id for item in result}) != len(result):
        raise DesignError("finalist task coverage is incomplete or duplicated")
    pairs: dict[str, list[FinalistTask]] = {}
    for task in result:
        pairs.setdefault(task.pair_id, []).append(task)
    if any(
        len(items) != 2
        or {item.start.state for item in items} != {"wt", "mutant"}
        or len({item.start.candidate.candidate_id for item in items}) != 1
        or len({item.start.template.template_id for item in items}) != 1
        or len({item.seed for item in items}) != 1
        for items in pairs.values()
    ):
        raise DesignError("finalist WT/mutant pair coverage is invalid")
    return result


def read_finalist_plan(*, config: AppConfig, run_dir: Path) -> FinalistPlan:
    """重新校验柔性复核计划、上游哈希、起点和完整成对任务矩阵。"""
    plan = read_document(run_dir / FINALIST_PLAN_NAME, name="design finalist plan")
    if (
        plan.get("schema") != "vela.design-finalist-plan/1"
        or plan.get("stage") != "design_flexible_verification"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("method_id") != config.design.method_id
        or plan.get("flexpepdock_method_id") != config.validation.method_id
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
        or plan.get("objective") != config.design.objective
        or plan.get("parameters") != design_parameters(config)
    ):
        raise DesignError("design finalist plan identity or parameters changed")
    try:
        inputs = object_mapping(plan.get("inputs"), name="finalist inputs")
        source_record = object_mapping(
            inputs.get("source_screen"), name="finalist source screen"
        )
        qualification_record = object_mapping(
            inputs.get("stage3_qualification_report"),
            name="Stage 3 qualification report",
        )
        candidate_rows = object_list(plan.get("candidates"), name="finalist candidates")
        template_rows = object_list(plan.get("templates"), name="finalist templates")
        start_rows = object_list(plan.get("starts"), name="finalist starts")
        task_rows = object_list(plan.get("tasks"), name="finalist tasks")
    except TypeError as exc:
        raise DesignError("design finalist plan structure is invalid") from exc
    for key in ("config_snapshot", "finalist_selection"):
        validate_record(root=run_dir, raw=inputs.get(key), name=f"finalist {key}")
    _qualification_report(config=config, record=qualification_record)
    source = _source_screen(config=config, source_record=source_record)
    candidates = tuple(
        candidate_from_record(raw=item, config=config) for item in candidate_rows
    )
    if not candidates or len({item.candidate_id for item in candidates}) != len(
        candidates
    ):
        raise DesignError("finalist candidates are empty or duplicated")
    templates = _templates(rows=template_rows, run_dir=run_dir)
    starts = _starts(
        rows=start_rows,
        run_dir=run_dir,
        candidates=candidates,
        templates=templates,
    )
    tasks = _tasks(
        rows=task_rows,
        starts=starts,
        seed_count=len(config.design.finalists.seeds),
    )
    verify_flexpepdock_tool(config.validation.rosetta)
    return FinalistPlan(run_dir, source, candidates, templates, starts, tasks)
