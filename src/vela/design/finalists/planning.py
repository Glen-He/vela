"""阶段四候选柔性复核的选择、起点冻结和计划。"""

from __future__ import annotations

import csv
import math
import shutil
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_list, object_mapping
from vela.design.finalists.records import FINALIST_PLAN_NAME, FinalistPlan
from vela.design.models import (
    DesignError,
    DesignTemplate,
    FinalistStart,
    FinalistTask,
    SequenceCandidate,
)
from vela.design.readiness import assess_design_readiness
from vela.design.screening.execution import (
    candidate_chemistry,
    task_histidine_pose_indices,
)
from vela.design.screening.planning import candidate_record, template_record
from vela.design.screening.records import (
    SCREEN_PLAN_NAME,
    ScreenPlan,
    design_parameters,
    read_screen_plan,
)
from vela.validation.records import (
    file_record,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.reconstruction import validate_flexpepdock_input
from vela.validation.rosetta import verify_flexpepdock_tool


def _table(path: Path, *, required: frozenset[str]) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise DesignError(f"finalist source columns are invalid: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows:
        raise DesignError(f"finalist source contains no rows: {path}")
    return rows


def _finite(value: str, *, name: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise DesignError(f"{name} is not a number") from exc
    if not math.isfinite(result):
        raise DesignError(f"{name} must be finite")
    return result


def _selected_candidates(
    *, config: AppConfig, plan: ScreenPlan, analysis_path: Path
) -> tuple[tuple[SequenceCandidate, ...], list[dict[str, str]]]:
    manifest = read_document(analysis_path, name="design screen analysis")
    if (
        manifest.get("schema") != "vela.design-screen-analysis-manifest/1"
        or manifest.get("status") != "completed"
        or manifest.get("design_round") != plan.design_round
        or manifest.get("objective") != config.design.objective
    ):
        raise DesignError("finalist source analysis identity is invalid")
    summary_path, _ = validate_record(
        root=analysis_path.parent,
        raw=manifest.get("candidate_summary"),
        name="finalist source candidate summary",
    )
    rows = _table(
        summary_path,
        required=frozenset(
            {
                "candidate_id",
                "positive_median_delta_score",
                "positive_worst_delta_score",
                "candidate_status",
            }
        ),
    )
    by_id = {item.candidate_id: item for item in plan.candidates}
    ranked: list[tuple[tuple[float, float, str], SequenceCandidate]] = []
    for row in rows:
        if row["candidate_status"] != "eligible":
            continue
        candidate = by_id.get(row["candidate_id"])
        if candidate is None:
            raise DesignError("finalist source references an unknown candidate")
        worst = _finite(row["positive_worst_delta_score"], name="positive worst score")
        median = _finite(
            row["positive_median_delta_score"], name="positive median score"
        )
        ranked.append(((worst, median, candidate.candidate_id), candidate))
    ranked.sort(key=lambda item: item[0])
    if not ranked:
        raise DesignError("screen analysis produced no eligible flexible finalists")
    limit = config.design.finalists.max_candidates
    selected = tuple(item[1] for item in ranked[:limit])
    selection_rows = [
        {
            "candidate_id": item[1].candidate_id,
            "sequence": item[1].sequence,
            "mutation_string": item[1].mutation_string,
            "generation": str(item[1].generation),
            "parent_candidate_ids": ";".join(
                parent.candidate_id for parent in item[1].parents
            ),
            "selection_status": "selected" if index < limit else "resource_deferred",
            "source_positive_worst_delta": f"{item[0][0]:.6f}",
            "source_positive_median_delta": f"{item[0][1]:.6f}",
        }
        for index, item in enumerate(ranked)
    ]
    return selected, selection_rows


def _screen_results(*, plan: ScreenPlan) -> dict[str, Path]:
    manifest_path = plan.run_dir / "screen_manifest.json"
    manifest = read_document(manifest_path, name="finalist source screen manifest")
    if (
        manifest.get("schema") != "vela.design-screen-manifest/1"
        or manifest.get("status") != "completed"
        or manifest.get("design_round") != plan.design_round
    ):
        raise DesignError("finalist source screen manifest identity is invalid")
    validate_record(
        root=plan.run_dir,
        raw=manifest.get("screen_plan"),
        name="finalist source screen plan",
    )
    try:
        rows = object_list(manifest.get("tasks"), name="source screen tasks")
    except TypeError as exc:
        raise DesignError("finalist source screen task list is invalid") from exc
    results: dict[str, Path] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="source screen task")
        except TypeError as exc:
            raise DesignError("finalist source screen task is invalid") from exc
        task_id = safe_identifier(row.get("task_id"), name="source screen task ID")
        if task_id in results:
            raise DesignError("finalist source screen task IDs are duplicated")
        path, _ = validate_record(
            root=plan.run_dir,
            raw=row.get("task_result"),
            name=f"{task_id} source result",
        )
        results[task_id] = path
    if set(results) != {item.task_id for item in plan.tasks}:
        raise DesignError("finalist source screen task coverage is incomplete")
    return results


def _copy_templates(
    *, run_dir: Path, templates: tuple[DesignTemplate, ...]
) -> tuple[DesignTemplate, ...]:
    destination_dir = run_dir / "templates"
    destination_dir.mkdir(parents=True)
    copied: list[DesignTemplate] = []
    for template in templates:
        destination = destination_dir / f"{template.template_id}.pdb"
        shutil.copyfile(template.path, destination)
        digest = sha256_file(destination)
        if digest != template.sha256:
            raise DesignError("copied finalist template hash changed")
        copied.append(
            DesignTemplate(
                template.template_id,
                template.evidence_role,
                template.cluster_id,
                template.candidate_id,
                template.receptor_id,
                template.target,
                destination,
                digest,
                template.receptor_residue_count,
                template.fixed_histidine_pose_indices,
            )
        )
    return tuple(copied)


def _copy_starts(
    *,
    config: AppConfig,
    run_dir: Path,
    source_plan: ScreenPlan,
    candidates: tuple[SequenceCandidate, ...],
    templates: tuple[DesignTemplate, ...],
    result_paths: dict[str, Path],
) -> tuple[FinalistStart, ...]:
    source_seed = min(config.design.seeds)
    source_tasks = {
        (item.candidate.candidate_id, item.template.template_id, item.state): item
        for item in source_plan.tasks
        if item.seed == source_seed
    }
    template_by_id = {item.template_id: item for item in templates}
    starts: list[FinalistStart] = []
    pair_index = 1
    start_index = 1
    for candidate in candidates:
        for source_template in source_plan.templates:
            template = template_by_id[source_template.template_id]
            pair_id = f"start_pair_{pair_index:06d}"
            pair_index += 1
            for state in ("wt", "mutant"):
                key = (candidate.candidate_id, template.template_id, state)
                screen_task = source_tasks.get(key)
                if screen_task is None:
                    raise DesignError(
                        "finalist source lacks its deterministic start task"
                    )
                result_path = result_paths[screen_task.task_id]
                result = read_document(result_path, name="source screen task result")
                if (
                    result.get("schema") != "vela.design-screen-task-result/1"
                    or result.get("status") != "completed"
                    or result.get("task_id") != screen_task.task_id
                    or result.get("state") != state
                    or result.get("seed") != source_seed
                ):
                    raise DesignError("finalist source task result identity is invalid")
                source_path, _ = validate_record(
                    root=result_path.parent,
                    raw=result.get("output"),
                    name="finalist source structure",
                )
                chemistry = candidate_chemistry(config=config, task=screen_task)
                receptor_count, _ = validate_flexpepdock_input(
                    path=source_path,
                    chemistry=chemistry,
                    min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
                    max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
                )
                if receptor_count != template.receptor_residue_count:
                    raise DesignError("finalist source receptor length changed")
                start_id = f"start_{start_index:06d}"
                start_index += 1
                destination = run_dir / "starts" / f"{start_id}.pdb"
                destination.parent.mkdir(parents=True, exist_ok=True)
                shutil.copyfile(source_path, destination)
                starts.append(
                    FinalistStart(
                        start_id,
                        pair_id,
                        state,
                        candidate,
                        template,
                        destination,
                        sha256_file(destination),
                        task_histidine_pose_indices(config=config, task=screen_task),
                    )
                )
    return tuple(starts)


def _tasks(
    *, config: AppConfig, starts: tuple[FinalistStart, ...]
) -> tuple[FinalistTask, ...]:
    tasks: list[FinalistTask] = []
    index = 1
    pair_index = 1
    starts_by_pair: dict[str, list[FinalistStart]] = {}
    for start in starts:
        starts_by_pair.setdefault(start.pair_id, []).append(start)
    for start_pair, paired in sorted(starts_by_pair.items()):
        if len(paired) != 2 or {item.state for item in paired} != {"wt", "mutant"}:
            raise DesignError(f"finalist start pair is incomplete: {start_pair}")
        for seed in config.design.finalists.seeds:
            pair_id = f"final_pair_{pair_index:07d}"
            pair_index += 1
            for start in sorted(paired, key=lambda item: item.state):
                tasks.append(FinalistTask(f"final_{index:07d}", pair_id, start, seed))
                index += 1
    if not tasks:
        raise DesignError("finalist task matrix is empty")
    return tuple(tasks)


def _selection_table(rows: list[dict[str, str]]) -> str:
    fields = (
        "candidate_id",
        "sequence",
        "mutation_string",
        "generation",
        "parent_candidate_ids",
        "selection_status",
        "source_positive_worst_delta",
        "source_positive_median_delta",
    )
    output = ["\t".join(fields)]
    output.extend("\t".join(row[field] for field in fields) for row in rows)
    return "\n".join(output) + "\n"


def write_finalist_plan(
    *, config: AppConfig, screen_run_dir: Path, run_id: str
) -> FinalistPlan:
    """从已分析界面筛查中冻结有限候选和独立柔性复核任务。"""
    validate_run_id(run_id)
    readiness = assess_design_readiness(config)
    if not readiness.finalist_ready:
        raise DesignError(
            "Stage 4 flexible verification is not ready: "
            + "; ".join(issue.code for issue in readiness.issues)
        )
    source = screen_run_dir.expanduser().resolve()
    source_root = (config.paths.outputs_dir / "design" / "screens").resolve()
    try:
        source_relative = source.relative_to(config.paths.outputs_dir.resolve())
        source.relative_to(source_root)
    except ValueError as exc:
        raise DesignError("finalist source is outside Stage 4 screen runs") from exc
    source_plan = read_screen_plan(config=config, run_dir=source)
    analysis_path = source / "screen_analysis" / "analysis_manifest.json"
    candidates, selection_rows = _selected_candidates(
        config=config, plan=source_plan, analysis_path=analysis_path
    )
    result_paths = _screen_results(plan=source_plan)
    stage3_report = config.validation.qualification_report
    if stage3_report is None:
        raise DesignError("qualified FlexPepDock method lacks its report")
    try:
        stage3_report_relative = stage3_report.resolve().relative_to(
            config.paths.outputs_dir.resolve()
        )
    except ValueError as exc:
        raise DesignError(
            "FlexPepDock qualification report must be inside outputs"
        ) from exc
    run_dir = config.paths.outputs_dir / "design" / "finalists" / run_id
    if run_dir.exists():
        raise DesignError(f"finalist run directory already exists: {run_dir}")
    run_dir.mkdir(parents=True)
    templates = _copy_templates(run_dir=run_dir, templates=source_plan.templates)
    starts = _copy_starts(
        config=config,
        run_dir=run_dir,
        source_plan=source_plan,
        candidates=candidates,
        templates=templates,
        result_paths=result_paths,
    )
    tasks = _tasks(config=config, starts=starts)
    snapshot = run_dir / "config.snapshot.txt"
    selection_path = run_dir / "finalist_selection.tsv"
    atomic_write_text(snapshot, config.source_snapshot_text)
    atomic_write_text(selection_path, _selection_table(selection_rows))
    tool = verify_flexpepdock_tool(config.validation.rosetta)
    plan_path = run_dir / FINALIST_PLAN_NAME
    atomic_write_json(
        plan_path,
        {
            "schema": "vela.design-finalist-plan/1",
            "stage": "design_flexible_verification",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "method_id": config.design.method_id,
            "flexpepdock_method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "source_screen_round": source_plan.design_round,
            "source_screen_seed_rule": "minimum_frozen_screen_seed",
            "source_screen_seed": min(config.design.seeds),
            "software": {
                **vela_software_identity(),
                "rosetta_version": tool.version,
                "flexpepdock_sha256": tool.executable_sha256,
            },
            "inputs": {
                "config_snapshot": file_record(snapshot, root=run_dir),
                "source_screen": {
                    "path": source_relative.as_posix(),
                    "screen_plan_sha256": sha256_file(source / SCREEN_PLAN_NAME),
                    "screen_manifest_sha256": sha256_file(
                        source / "screen_manifest.json"
                    ),
                    "analysis_manifest_sha256": sha256_file(analysis_path),
                },
                "finalist_selection": file_record(selection_path, root=run_dir),
                "stage3_qualification_report": {
                    "path": stage3_report_relative.as_posix(),
                    "sha256": config.validation.qualification_report_sha256,
                },
            },
            "parameters": design_parameters(config),
            "candidates": [candidate_record(item) for item in candidates],
            "templates": [template_record(item, run_dir=run_dir) for item in templates],
            "starts": [
                {
                    "start_id": item.start_id,
                    "pair_id": item.pair_id,
                    "state": item.state,
                    "candidate_id": item.candidate.candidate_id,
                    "template_id": item.template.template_id,
                    "structure": file_record(item.path, root=run_dir),
                    "histidine_pose_indices": list(item.histidine_pose_indices),
                }
                for item in starts
            ],
            "tasks": [
                {
                    "task_id": item.task_id,
                    "pair_id": item.pair_id,
                    "start_id": item.start.start_id,
                    "seed": item.seed,
                    "status": "planned",
                }
                for item in tasks
            ],
        },
    )
    return FinalistPlan(run_dir, source, candidates, templates, starts, tasks)
