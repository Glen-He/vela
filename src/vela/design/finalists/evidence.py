"""阶段四候选柔性复核的结构 QC 与成对证据。"""

from __future__ import annotations

import csv
import math
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.design.finalists.execution import finalist_chemistry
from vela.design.finalists.records import FinalistPlan
from vela.design.models import DesignError, FinalistTask, SequenceCandidate
from vela.validation.records import read_document, validate_record
from vela.validation.refinement.geometry import (
    GeometryAssessment,
    assess_complex_geometry,
    read_complex_geometry,
    resolve_analysis_settings,
)
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class FinalistDecoy:
    """一份柔性复核 decoy 的分数、几何和完整来源。"""

    decoy_id: str
    task: FinalistTask
    path: Path
    sha256: str
    ranking_score: float
    interface_score: float
    geometry: GeometryAssessment


@dataclass(frozen=True, slots=True)
class StateSummary:
    """一个序列状态在单模板、单 seed 下的通过率和稳健分数。"""

    task: FinalistTask
    passed_fraction: float
    successful: bool
    ranking_median: float | None
    interface_median: float | None


@dataclass(frozen=True, slots=True)
class PairedSummary:
    """同模板、同 seed 的候选减 WT 柔性复核统计。"""

    pair_id: str
    candidate: SequenceCandidate
    template_id: str
    target: str
    evidence_role: str
    seed: int
    wt_passed_fraction: float
    mutant_passed_fraction: float
    successful: bool
    ranking_delta: float | None
    interface_delta: float | None


@dataclass(frozen=True, slots=True)
class CandidateEvidence:
    """候选跨模板不可补偿门槛及资源状态。"""

    candidate: SequenceCandidate
    status: str
    failed_gates: tuple[str, ...]
    positive_ranking_median: float | None
    positive_ranking_worst: float | None
    positive_interface_median: float | None
    positive_interface_worst: float | None


def manifest_results(*, config: AppConfig, plan: FinalistPlan) -> dict[str, Path]:
    manifest = read_document(
        plan.run_dir / "finalist_manifest.json", name="finalist manifest"
    )
    if (
        manifest.get("schema") != "vela.design-finalist-manifest/1"
        or manifest.get("stage") != "design_flexible_verification"
        or manifest.get("status") != "completed"
        or manifest.get("method_id") != config.design.method_id
        or manifest.get("flexpepdock_method_id") != config.validation.method_id
        or manifest.get("chemistry_id") != config.chemistry.chemistry_id
        or manifest.get("objective") != config.design.objective
    ):
        raise DesignError("finalist manifest identity is invalid")
    validate_record(
        root=plan.run_dir,
        raw=manifest.get("finalist_plan"),
        name="finalist plan",
    )
    try:
        rows = object_list(manifest.get("tasks"), name="finalist manifest tasks")
    except TypeError as exc:
        raise DesignError("finalist manifest tasks are invalid") from exc
    results: dict[str, Path] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="finalist manifest task")
        except TypeError as exc:
            raise DesignError("finalist manifest task is invalid") from exc
        task_id = row.get("task_id")
        if not isinstance(task_id, str) or task_id in results:
            raise DesignError(
                "finalist manifest task identity is invalid or duplicated"
            )
        result, _ = validate_record(
            root=plan.run_dir,
            raw=row.get("task_result"),
            name=f"{task_id} finalist result",
        )
        results[task_id] = result
    if set(results) != {item.task_id for item in plan.tasks}:
        raise DesignError("finalist manifest task coverage is incomplete")
    return results


def _decoy_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != ["description", "ranking_score", "path", "sha256"]:
            raise DesignError(f"invalid finalist decoy columns: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows:
        raise DesignError(f"finalist decoy manifest contains no rows: {path}")
    return rows


def task_decoys(
    *, config: AppConfig, task: FinalistTask, result_path: Path
) -> tuple[FinalistDecoy, ...]:
    result = read_document(result_path, name="finalist task result")
    if (
        result.get("schema") != "vela.design-finalist-task-result/1"
        or result.get("status") != "completed"
        or result.get("task_id") != task.task_id
        or result.get("pair_id") != task.pair_id
        or result.get("state") != task.start.state
        or result.get("candidate_id") != task.start.candidate.candidate_id
        or result.get("template_id") != task.start.template.template_id
        or result.get("seed") != task.seed
    ):
        raise DesignError(f"finalist task result identity is invalid: {task.task_id}")
    task_dir = result_path.parent
    score_path, _ = validate_record(
        root=task_dir, raw=result.get("scorefile"), name="finalist scorefile"
    )
    decoy_path, _ = validate_record(
        root=task_dir, raw=result.get("decoy_manifest"), name="finalist decoy manifest"
    )
    for key in ("fix_disulfide", "log"):
        validate_record(root=task_dir, raw=result.get(key), name=f"finalist {key}")
    ranking_name = config.design.finalists.ranking_score
    interface_name = config.design.finalists.interface_score
    scores = {row.description: row for row in read_rosetta_scorefile(score_path)}
    start_geometry = read_complex_geometry(
        path=task.start.path,
        interface_contact_A=config.validation.interface_contact_A,
    )
    geometry_settings = resolve_analysis_settings(config.validation.analysis)
    chemistry = finalist_chemistry(config=config, start=task.start)
    decoys: list[FinalistDecoy] = []
    for row in _decoy_rows(decoy_path):
        description = row["description"]
        score = scores.get(description)
        path = (decoy_path.parent / row["path"]).resolve()
        try:
            path.relative_to(decoy_path.parent.resolve())
            recorded_ranking = float(row["ranking_score"])
        except (ValueError, OverflowError) as exc:
            raise DesignError(f"invalid finalist decoy record: {description}") from exc
        if (
            score is None
            or not math.isfinite(recorded_ranking)
            or row["ranking_score"] != f"{score.score(ranking_name):.6f}"
            or not path.is_file()
            or sha256_file(path) != row["sha256"]
        ):
            raise DesignError(f"finalist decoy record mismatch: {description}")
        decoys.append(
            FinalistDecoy(
                f"{task.task_id}__{description}",
                task,
                path,
                row["sha256"],
                score.score(ranking_name),
                score.score(interface_name),
                assess_complex_geometry(
                    path=path,
                    chemistry=chemistry,
                    start=start_geometry,
                    cluster_reference=start_geometry,
                    config=config,
                    settings=geometry_settings,
                ),
            )
        )
    if len(decoys) != config.validation.rosetta.decoys_per_seed:
        raise DesignError("finalist task decoy count differs from its frozen budget")
    return tuple(decoys)


def state_summary(
    *, config: AppConfig, task: FinalistTask, decoys: tuple[FinalistDecoy, ...]
) -> StateSummary:
    passed = [item for item in decoys if item.geometry.passed]
    fraction = len(passed) / len(decoys)
    threshold = config.design.finalists.min_passed_decoy_fraction
    if threshold is None:
        raise DesignError("finalist passed-decoy threshold is unresolved")
    successful = fraction >= threshold and bool(passed)
    return StateSummary(
        task,
        fraction,
        successful,
        statistics.median(item.ranking_score for item in passed)
        if successful
        else None,
        statistics.median(item.interface_score for item in passed)
        if successful
        else None,
    )


def paired_summaries(states: tuple[StateSummary, ...]) -> tuple[PairedSummary, ...]:
    grouped: dict[str, list[StateSummary]] = defaultdict(list)
    for state in states:
        grouped[state.task.pair_id].append(state)
    pairs: list[PairedSummary] = []
    for pair_id, values in sorted(grouped.items()):
        by_state = {item.task.start.state: item for item in values}
        if set(by_state) != {"wt", "mutant"}:
            raise DesignError(f"finalist pair is incomplete: {pair_id}")
        wt = by_state["wt"]
        mutant = by_state["mutant"]
        if (
            wt.task.start.candidate != mutant.task.start.candidate
            or wt.task.start.template != mutant.task.start.template
            or wt.task.seed != mutant.task.seed
        ):
            raise DesignError(f"finalist pair inputs differ: {pair_id}")
        successful = wt.successful and mutant.successful
        ranking_delta = (
            mutant.ranking_median - wt.ranking_median
            if successful
            and mutant.ranking_median is not None
            and wt.ranking_median is not None
            else None
        )
        interface_delta = (
            mutant.interface_median - wt.interface_median
            if successful
            and mutant.interface_median is not None
            and wt.interface_median is not None
            else None
        )
        template = mutant.task.start.template
        pairs.append(
            PairedSummary(
                pair_id,
                mutant.task.start.candidate,
                template.template_id,
                template.target,
                template.evidence_role,
                mutant.task.seed,
                wt.passed_fraction,
                mutant.passed_fraction,
                successful,
                ranking_delta,
                interface_delta,
            )
        )
    return tuple(pairs)


def candidate_evidence(
    *, config: AppConfig, pairs: tuple[PairedSummary, ...]
) -> tuple[CandidateEvidence, ...]:
    settings = config.design.finalists
    required = settings.min_successful_seeds
    ranking_median_gate = settings.max_positive_median_ranking_delta
    ranking_worst_gate = settings.max_positive_worst_ranking_delta
    interface_median_gate = settings.max_positive_median_interface_delta
    interface_worst_gate = settings.max_positive_worst_interface_delta
    if (
        required is None
        or ranking_median_gate is None
        or ranking_worst_gate is None
        or interface_median_gate is None
        or interface_worst_gate is None
    ):
        raise DesignError("finalist analysis thresholds are unresolved")
    grouped: dict[str, list[PairedSummary]] = defaultdict(list)
    for pair in pairs:
        grouped[pair.candidate.candidate_id].append(pair)
    results: list[CandidateEvidence] = []
    for _, values in sorted(grouped.items()):
        candidate = values[0].candidate
        failed: list[str] = []
        by_template: dict[str, list[PairedSummary]] = defaultdict(list)
        for pair in values:
            by_template[pair.template_id].append(pair)
        if any(
            sum(item.successful for item in template_pairs) < required
            for template_pairs in by_template.values()
        ):
            failed.append("template_seed_support")
        if any(item.evidence_role != "positive" for item in values):
            raise DesignError("single-target finalist contains a non-target template")
        positive = [item for item in values if item.successful]
        positive_ranking = [
            item.ranking_delta for item in positive if item.ranking_delta is not None
        ]
        positive_interface = [
            item.interface_delta
            for item in positive
            if item.interface_delta is not None
        ]
        if not positive_ranking or not positive_interface:
            failed.append("positive_states_missing")
            positive_ranking_median = None
            positive_ranking_worst = None
            positive_interface_median = None
            positive_interface_worst = None
        else:
            positive_ranking_median = float(statistics.median(positive_ranking))
            positive_ranking_worst = max(positive_ranking)
            positive_interface_median = float(statistics.median(positive_interface))
            positive_interface_worst = max(positive_interface)
            if positive_ranking_median > ranking_median_gate:
                failed.append("positive_ranking_median")
            if positive_ranking_worst > ranking_worst_gate:
                failed.append("positive_ranking_worst")
            if positive_interface_median > interface_median_gate:
                failed.append("positive_interface_median")
            if positive_interface_worst > interface_worst_gate:
                failed.append("positive_interface_worst")
        results.append(
            CandidateEvidence(
                candidate,
                "eligible" if not failed else "rejected",
                tuple(failed),
                positive_ranking_median,
                positive_ranking_worst,
                positive_interface_median,
                positive_interface_worst,
            )
        )
    return tuple(results)
