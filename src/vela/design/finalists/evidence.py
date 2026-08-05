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
from vela.design.scores import (
    FINALIST_SCORE_COLUMNS,
    FinalistMetrics,
    finalist_metrics,
)
from vela.discovery.analysis.cluster_engine import complete_linkage
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
    metrics: FinalistMetrics
    geometry: GeometryAssessment


@dataclass(frozen=True, slots=True)
class PoseCluster:
    """单任务内不混合其他 pose basin 的通过构象簇。"""

    cluster_id: str
    members: tuple[FinalistDecoy, ...]
    medoid: FinalistDecoy
    ranking_median: float
    interface_median: float
    peptide_median: float


@dataclass(frozen=True, slots=True)
class StateSummary:
    """一个序列状态在单模板、单 seed 下的通过率和稳健分数。"""

    task: FinalistTask
    passed_fraction: float
    successful: bool
    pose_clusters: tuple[PoseCluster, ...]
    primary_pose_cluster_id: str | None
    ranking_median: float | None
    interface_median: float | None
    peptide_median: float | None


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
    pose_relation: str
    wt_pose_cluster_id: str | None
    mutant_pose_cluster_id: str | None
    paired_reweighted_sc_delta: float | None
    paired_I_sc_delta: float | None
    paired_pep_sc_delta: float | None


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
    positive_peptide_median: float | None
    successful_pair_fraction: float


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
        expected = ["description", *FINALIST_SCORE_COLUMNS, "path", "sha256"]
        if reader.fieldnames != expected:
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
        result.get("schema") != "vela.design-finalist-task-result/2"
        or result.get("status") != "completed"
        or result.get("task_id") != task.task_id
        or result.get("pair_id") != task.pair_id
        or result.get("state") != task.start.state
        or result.get("candidate_id") != task.start.candidate.candidate_id
        or result.get("template_id") != task.start.template.template_id
        or result.get("seed") != task.seed
        or result.get("score_columns") != FINALIST_SCORE_COLUMNS
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
    scores = {row.description: row for row in read_rosetta_scorefile(score_path)}
    start_geometry = read_complex_geometry(
        path=task.start.path,
        interface_contact_A=config.validation.interface_contact_A,
    )
    cluster_reference = read_complex_geometry(
        path=task.start.template.path,
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
            recorded_metrics = {
                name: float(row[name]) for name in FINALIST_SCORE_COLUMNS
            }
        except (ValueError, OverflowError) as exc:
            raise DesignError(f"invalid finalist decoy record: {description}") from exc
        if (
            score is None
            or any(not math.isfinite(value) for value in recorded_metrics.values())
            or any(
                row[name] != f"{value:.6f}"
                for name, value in finalist_metrics(score).as_dict().items()
            )
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
                finalist_metrics(score),
                assess_complex_geometry(
                    path=path,
                    chemistry=chemistry,
                    start=start_geometry,
                    cluster_reference=cluster_reference,
                    config=config,
                    settings=geometry_settings,
                ),
            )
        )
    if len(decoys) != config.validation.rosetta.decoys_per_seed:
        raise DesignError("finalist task decoy count differs from its frozen budget")
    return tuple(decoys)


def _backbone_rmsd(first: FinalistDecoy, second: FinalistDecoy) -> float:
    left = first.geometry.cluster_backbone
    right = second.geometry.cluster_backbone
    if not left or len(left) != len(right):
        raise DesignError("finalist peptide backbone correspondence is invalid")
    return math.sqrt(
        sum(
            first_position.dist(second_position) ** 2
            for first_position, second_position in zip(left, right, strict=True)
        )
        / len(left)
    )


def _pose_clusters(
    *, config: AppConfig, task: FinalistTask, decoys: tuple[FinalistDecoy, ...]
) -> tuple[PoseCluster, ...]:
    """按共同受体坐标系执行确定性 complete-linkage。"""
    if not decoys:
        return ()
    threshold = resolve_analysis_settings(
        config.validation.analysis
    ).max_cluster_backbone_rmsd_A
    groups = complete_linkage(
        decoys,
        distance=lambda first, second: _backbone_rmsd(first, second) / threshold,
        identity=lambda item: item.decoy_id,
    )
    clusters: list[PoseCluster] = []
    for index, members in enumerate(groups, 1):
        medoid = min(
            members,
            key=lambda candidate: (
                sum(_backbone_rmsd(candidate, member) for member in members),
                candidate.decoy_id,
            ),
        )
        clusters.append(
            PoseCluster(
                cluster_id=f"{task.task_id}__pose_{index:04d}",
                members=members,
                medoid=medoid,
                ranking_median=float(
                    statistics.median(item.metrics.reweighted_sc for item in members)
                ),
                interface_median=float(
                    statistics.median(item.metrics.I_sc for item in members)
                ),
                peptide_median=float(
                    statistics.median(item.metrics.pep_sc for item in members)
                ),
            )
        )
    return tuple(
        sorted(
            clusters,
            key=lambda cluster: (
                -len(cluster.members),
                cluster.ranking_median,
                cluster.cluster_id,
            ),
        )
    )


def state_summary(
    *, config: AppConfig, task: FinalistTask, decoys: tuple[FinalistDecoy, ...]
) -> StateSummary:
    passed = [item for item in decoys if item.geometry.passed]
    fraction = len(passed) / len(decoys)
    threshold = config.design.finalists.min_passed_decoy_fraction
    if threshold is None:
        raise DesignError("finalist passed-decoy threshold is unresolved")
    successful = fraction >= threshold and bool(passed)
    pose_clusters = _pose_clusters(config=config, task=task, decoys=tuple(passed))
    primary = pose_clusters[0] if successful and pose_clusters else None
    return StateSummary(
        task,
        fraction,
        successful,
        pose_clusters,
        primary.cluster_id if primary is not None else None,
        primary.ranking_median if primary is not None else None,
        primary.interface_median if primary is not None else None,
        primary.peptide_median if primary is not None else None,
    )


def _contact_similarity(first: PoseCluster, second: PoseCluster) -> float:
    left = first.medoid.geometry.receptor_contacts
    right = second.medoid.geometry.receptor_contacts
    union = left | right
    return len(left & right) / len(union) if union else 0.0


def _matched_pose_pair(
    *, config: AppConfig, wt: StateSummary, mutant: StateSummary
) -> tuple[PoseCluster, PoseCluster] | None:
    settings = resolve_analysis_settings(config.validation.analysis)
    matches = [
        (
            -min(len(wt_cluster.members), len(mutant_cluster.members)),
            _backbone_rmsd(wt_cluster.medoid, mutant_cluster.medoid),
            -_contact_similarity(wt_cluster, mutant_cluster),
            wt_cluster.cluster_id,
            mutant_cluster.cluster_id,
            wt_cluster,
            mutant_cluster,
        )
        for wt_cluster in wt.pose_clusters
        for mutant_cluster in mutant.pose_clusters
        if _backbone_rmsd(wt_cluster.medoid, mutant_cluster.medoid)
        <= settings.max_cluster_backbone_rmsd_A
        and _contact_similarity(wt_cluster, mutant_cluster)
        >= settings.min_start_contact_overlap
    ]
    if not matches:
        return None
    selected = min(matches, key=lambda item: item[:5])
    return selected[5], selected[6]


def paired_summaries(
    *, config: AppConfig, states: tuple[StateSummary, ...]
) -> tuple[PairedSummary, ...]:
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
        matched = (
            _matched_pose_pair(config=config, wt=wt, mutant=mutant)
            if wt.successful and mutant.successful
            else None
        )
        successful = matched is not None
        wt_cluster = matched[0] if matched is not None else None
        mutant_cluster = matched[1] if matched is not None else None
        pose_relation = (
            "matched_pose"
            if matched is not None
            else "alternative_pose"
            if wt.pose_clusters and mutant.pose_clusters
            else "unmatched_wt_pose"
            if wt.pose_clusters
            else "unmatched_mutant_pose"
            if mutant.pose_clusters
            else "no_valid_pose"
        )
        ranking_delta = (
            mutant_cluster.ranking_median - wt_cluster.ranking_median
            if mutant_cluster is not None and wt_cluster is not None
            else None
        )
        interface_delta = (
            mutant_cluster.interface_median - wt_cluster.interface_median
            if mutant_cluster is not None and wt_cluster is not None
            else None
        )
        peptide_delta = (
            mutant_cluster.peptide_median - wt_cluster.peptide_median
            if mutant_cluster is not None and wt_cluster is not None
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
                pose_relation,
                wt_cluster.cluster_id if wt_cluster is not None else None,
                mutant_cluster.cluster_id if mutant_cluster is not None else None,
                ranking_delta,
                interface_delta,
                peptide_delta,
            )
        )
    return tuple(pairs)


def candidate_evidence(
    *, config: AppConfig, pairs: tuple[PairedSummary, ...]
) -> tuple[CandidateEvidence, ...]:
    settings = config.design.finalists
    required = settings.min_successful_seeds
    ranking_median_gate = settings.max_positive_median_ranking_delta
    interface_median_gate = settings.max_positive_median_interface_delta
    if required is None or ranking_median_gate is None or interface_median_gate is None:
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
        for template_id, template_pairs in sorted(by_template.items()):
            successful_pairs = [item for item in template_pairs if item.successful]
            if len(successful_pairs) < required:
                failed.append(f"template_seed_support:{template_id}")
                continue
            template_ranking = [
                item.paired_reweighted_sc_delta
                for item in successful_pairs
                if item.paired_reweighted_sc_delta is not None
            ]
            template_interface = [
                item.paired_I_sc_delta
                for item in successful_pairs
                if item.paired_I_sc_delta is not None
            ]
            if (
                not template_ranking
                or statistics.median(template_ranking) > ranking_median_gate
            ):
                failed.append(f"template_ranking_median:{template_id}")
            if (
                not template_interface
                or statistics.median(template_interface) > interface_median_gate
            ):
                failed.append(f"template_interface_median:{template_id}")
        if any(item.evidence_role != "positive" for item in values):
            raise DesignError("single-target finalist contains a non-target template")
        positive = [item for item in values if item.successful]
        positive_ranking = [
            item.paired_reweighted_sc_delta
            for item in positive
            if item.paired_reweighted_sc_delta is not None
        ]
        positive_interface = [
            item.paired_I_sc_delta
            for item in positive
            if item.paired_I_sc_delta is not None
        ]
        positive_peptide = [
            item.paired_pep_sc_delta
            for item in positive
            if item.paired_pep_sc_delta is not None
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
        positive_peptide_median = (
            float(statistics.median(positive_peptide)) if positive_peptide else None
        )
        results.append(
            CandidateEvidence(
                candidate,
                "eligible" if not failed else "rejected",
                tuple(failed),
                positive_ranking_median,
                positive_ranking_worst,
                positive_interface_median,
                positive_interface_worst,
                positive_peptide_median,
                len(positive) / len(values),
            )
        )
    return tuple(results)
