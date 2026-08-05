"""阶段四迭代父序列选择、来源核验和邻域初筛计划。"""

from __future__ import annotations

import csv
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import JsonValue, sha256_file
from vela.core.run_identity import validate_run_id
from vela.design.finalists.records import read_finalist_plan
from vela.design.models import DesignError, SequenceCandidate
from vela.design.readiness import assess_design_readiness
from vela.design.screening.planning import write_screen_plan
from vela.design.screening.records import ScreenPlan
from vela.design.sequence.neighborhood import (
    iteration_candidate_table,
    iteration_library,
)
from vela.validation.records import read_document, validate_record


def _candidate_statuses(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"candidate_id", "candidate_status"}
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise DesignError("iteration parent summary columns are invalid")
        rows = tuple(dict(row) for row in reader)
    statuses: dict[str, str] = {}
    for row in rows:
        candidate_id = row["candidate_id"]
        if not candidate_id or candidate_id in statuses:
            raise DesignError("iteration parent summary identities are invalid")
        statuses[candidate_id] = row["candidate_status"]
    if not statuses:
        raise DesignError("iteration parent summary is empty")
    return statuses


def _selected_parents(
    *,
    candidates: tuple[SequenceCandidate, ...],
    statuses: dict[str, str],
    parent_candidate_ids: tuple[str, ...],
    max_parents: int,
) -> tuple[SequenceCandidate, ...]:
    if (
        not parent_candidate_ids
        or len(parent_candidate_ids) > max_parents
        or len(parent_candidate_ids) != len(set(parent_candidate_ids))
    ):
        raise DesignError(
            "iteration parent selection is empty, duplicated, or too large"
        )
    by_id = {item.candidate_id: item for item in candidates}
    selected: list[SequenceCandidate] = []
    for candidate_id in parent_candidate_ids:
        candidate = by_id.get(candidate_id)
        status = statuses.get(candidate_id)
        if candidate is None or status is None:
            raise DesignError(
                f"iteration parent is absent from its finalist run: {candidate_id}"
            )
        if status != "eligible":
            raise DesignError(f"iteration parent must be eligible: {candidate_id}")
        selected.append(candidate)
    if len({item.generation for item in selected}) != 1:
        raise DesignError("iteration parents must belong to the same generation")
    return tuple(sorted(selected, key=lambda item: item.candidate_id))


def write_iteration_screen_plan(
    *,
    config: AppConfig,
    finalist_run_dir: Path,
    run_id: str,
    parent_candidate_ids: tuple[str, ...],
) -> ScreenPlan:
    """从显式选择的柔性复核候选冻结下一代一步邻域筛查。"""
    validate_run_id(run_id)
    readiness = assess_design_readiness(config)
    if not readiness.finalist_ready:
        raise DesignError(
            "Stage 4 iteration is not ready: "
            + "; ".join(issue.code for issue in readiness.issues)
        )
    source = finalist_run_dir.expanduser().resolve()
    source_root = (config.paths.outputs_dir / "design" / "finalists").resolve()
    try:
        source.relative_to(source_root)
        source_path = source.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DesignError("iteration source is outside Stage 4 finalist runs") from exc
    parent_plan = read_finalist_plan(config=config, run_dir=source)
    analysis_path = source / "finalist_analysis" / "analysis_manifest.json"
    analysis = read_document(analysis_path, name="iteration parent analysis")
    if (
        analysis.get("schema") != "vela.design-finalist-analysis-manifest/1"
        or analysis.get("stage") != "design_flexible_verification_analysis"
        or analysis.get("status") != "completed"
        or analysis.get("objective") != config.design.objective
    ):
        raise DesignError("iteration parent analysis identity is invalid")
    summary_path, _ = validate_record(
        root=analysis_path.parent,
        raw=analysis.get("candidate_summary"),
        name="iteration parent candidate summary",
    )
    parents = _selected_parents(
        candidates=parent_plan.candidates,
        statuses=_candidate_statuses(summary_path),
        parent_candidate_ids=parent_candidate_ids,
        max_parents=config.design.iteration.max_parents,
    )
    library = iteration_library(
        chemistry=config.chemistry,
        settings=config.design,
        parents=parents,
    )
    source_inputs: dict[str, JsonValue] = {
        "parent_finalist": {
            "path": source_path.as_posix(),
            "selection_basis": "explicit_iteration_parent_candidate_ids",
            "selected_parent_candidate_ids": [item.candidate_id for item in parents],
            "finalist_plan_sha256": sha256_file(source / "finalist_plan.json"),
            "finalist_manifest_sha256": sha256_file(source / "finalist_manifest.json"),
            "analysis_manifest_sha256": sha256_file(analysis_path),
        }
    }
    return write_screen_plan(
        config=config,
        run_id=run_id,
        design_round="iteration",
        candidates=library.selected,
        templates=parent_plan.templates,
        source_inputs=source_inputs,
        candidate_library_text=iteration_candidate_table(library),
    )
