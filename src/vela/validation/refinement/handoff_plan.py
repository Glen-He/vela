"""阶段二受支持 site 到局部精修交接任务的选择与冻结。"""

from __future__ import annotations

import csv
import math
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path
from statistics import median

import gemmi

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.discovery.analysis.cluster_engine import bounded_leader_clusters
from vela.discovery.analysis.evidence import PoseEvidence, ReceptorSite
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.analysis.reports import (
    ReportedCandidateSite,
    ReportedReceptorSite,
    read_site_analysis_report,
)
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.evidence import (
    align_receptor,
    read_reference_chains,
    read_structure,
    required_atom,
    split_model,
)
from vela.discovery.sampling.planning import (
    EXPLORATORY_DISCOVERY_EVIDENCE,
    MAIN_DISCOVERY_EVIDENCE,
)
from vela.validation.models import ValidationError
from vela.validation.records import read_document, safe_identifier, validate_record
from vela.validation.refinement.reconstruction import verify_cg2all_tool
from vela.validation.rosetta import verify_rosetta_scripts_tool

HANDOFF_PLAN_NAME = "handoff_plan.json"
HANDOFF_PLAN_SCHEMA = "vela.validation-handoff-plan/11"
MAIN_HANDOFF_EVIDENCE = "main_discovery_handoff"
EXPLORATORY_HANDOFF_EVIDENCE = "exploratory_discovery_handoff"
SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE = (
    "exploratory_source_seed_confirmation_handoff"
)
FUNNEL_SCREENING_HANDOFF_EVIDENCE = "exploratory_funnel_screening_handoff"
POSE_SELECTION_METHOD = (
    "ranked_pose_clusters; top_ceil_half_budget_clusters; "
    "cluster_medoid_then_deterministic_farthest_point"
)
POSE_SELECTION_SEED_POLICY = "diversity_tiebreak_not_hard_requirement"
SOURCE_SEED_CONFIRMATION_SELECTION_METHOD = (
    "ranked_pose_clusters; one_seed_specific_medoid_per_source_seed; "
    "best_four_distinct_source_seeds"
)
SOURCE_SEED_CONFIRMATION_POLICY = "distinct_source_seed_required_per_start"


@dataclass(frozen=True, slots=True)
class CandidateHandoffTask:
    """一个保留 blind 证据身份的全原子重建任务。"""

    task_id: str
    candidate_id: str
    receptor_site_id: str
    pose: PoseEvidence
    reference_receptor_path: Path
    reference_receptor_sha256: str


@dataclass(frozen=True, slots=True)
class CandidateHandoffPlan:
    """已经冻结到独立运行目录的候选交接计划。"""

    run_id: str
    run_dir: Path
    discovery_run_dir: Path
    tasks: tuple[CandidateHandoffTask, ...]


@dataclass(frozen=True, slots=True)
class RankedPoseFamilies:
    """一个受体site内已对齐、按跨source证据排序的CABS姿态家族。"""

    clusters: tuple[tuple[PoseEvidence, ...], ...]
    coordinates: dict[str, tuple[gemmi.Position, ...]]

    def distance(self, first: PoseEvidence, second: PoseEvidence) -> float:
        """计算受体对齐后、保持残基顺序的配体CA RMSD。"""
        left = self.coordinates[first.pose_id]
        right = self.coordinates[second.pose_id]
        if len(left) != len(right) or not left:
            raise ValidationError("handoff peptide CA identities differ")
        return math.sqrt(
            sum(
                first_position.dist(second_position) ** 2
                for first_position, second_position in zip(left, right, strict=True)
            )
            / len(left)
        )


@dataclass(frozen=True, slots=True)
class FunnelScreeningSelection:
    """Stage 3A-0冻结的候选、起点和完整审计记录。"""

    candidate_ids: tuple[str, ...]
    candidate_arms: dict[str, tuple[str, ...]]
    tasks: tuple[CandidateHandoffTask, ...]
    audit: dict[str, JsonValue]


def handoff_budget_record(
    *,
    config: AppConfig,
    tasks: tuple[CandidateHandoffTask, ...],
    candidate_count: int,
    refinement_seed_count: int | None = None,
    refinement_mode: str = "qualified_full_protocol_only",
) -> dict[str, JsonValue]:
    """返回候选交接到正式局部精修的完整确定性计算预算。"""
    if candidate_count < 1:
        raise ValidationError("handoff candidate count must be positive")
    receptor_site_count = len({task.receptor_site_id for task in tasks})
    start_count = len(tasks)
    seed_count = (
        len(config.validation.seeds)
        if refinement_seed_count is None
        else refinement_seed_count
    )
    if seed_count < 1:
        raise ValidationError("handoff refinement seed count must be positive")
    refinement_task_count = start_count * seed_count
    return {
        "candidate_site_count": candidate_count,
        "receptor_site_count": receptor_site_count,
        "all_atom_start_count": start_count,
        "refinement_seeds_per_start": seed_count,
        "refinement_decoys_per_seed": config.validation.rosetta.decoys_per_seed,
        "projected_refinement_task_count": refinement_task_count,
        "projected_refinement_decoy_count": (
            refinement_task_count * config.validation.rosetta.decoys_per_seed
        ),
        "refinement_mode": refinement_mode,
    }


def _identifier(value: str, *, name: str) -> str:
    return safe_identifier(value, name=name)


def _selected_candidates(
    *, candidates: dict[str, ReportedCandidateSite], requested_ids: tuple[str, ...]
) -> tuple[ReportedCandidateSite, ...]:
    requested_ids = tuple(
        _identifier(value, name="requested candidate ID") for value in requested_ids
    )
    if len(requested_ids) != len(set(requested_ids)):
        raise ValidationError("requested candidate IDs must be unique")
    if not requested_ids:
        raise ValidationError(
            "at least one explicit candidate ID is required for all-atom handoff"
        )
    unknown = sorted(set(requested_ids) - set(candidates))
    if unknown:
        raise ValidationError("unknown candidate IDs: " + ", ".join(unknown))
    selected = tuple(candidates[candidate_id] for candidate_id in requested_ids)
    if not selected:
        raise ValidationError("no handoff-eligible candidate sites are available")
    ineligible = [item.candidate_id for item in selected if not item.handoff_eligible]
    if ineligible:
        raise ValidationError(
            "candidates outside their frozen evidence-tier budget cannot enter "
            "blind handoff: " + ", ".join(ineligible)
        )
    return selected


def candidate_evidence_records_for_evidence(
    *,
    discovery_run_dir: Path,
    candidate_ids: tuple[str, ...],
    expected_evidence_category: str,
) -> list[dict[str, JsonValue]]:
    """读取冻结报告并返回显式候选的证据等级与预算内排名。"""
    try:
        report = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=expected_evidence_category,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    selected = _selected_candidates(
        candidates=report.candidate_sites,
        requested_ids=candidate_ids,
    )
    return [
        {
            "candidate_id": candidate.candidate_id,
            "evidence_tier": candidate.evidence_tier,
            "rank_within_tier": candidate.rank_within_tier,
        }
        for candidate in selected
    ]


def candidate_evidence_records(
    *, discovery_run_dir: Path, candidate_ids: tuple[str, ...]
) -> list[dict[str, JsonValue]]:
    """返回正式主发现候选的冻结证据等级和排名。"""
    return candidate_evidence_records_for_evidence(
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
        expected_evidence_category=MAIN_DISCOVERY_EVIDENCE,
    )


def select_ranked_pose_cluster_representatives(
    *,
    ranked_clusters: tuple[tuple[PoseEvidence, ...], ...],
    count: int,
    distance: Callable[[PoseEvidence, PoseEvidence], float],
) -> tuple[PoseEvidence, ...]:
    """在排名靠前的姿态簇间均衡保留中心与几何多样性。"""
    if count < 1:
        raise ValidationError("handoff pose count must be positive")
    available = sum(len(cluster) for cluster in ranked_clusters)
    if available < count:
        raise ValidationError(
            f"handoff has only {available} passed poses; requested {count}"
        )
    if any(not cluster for cluster in ranked_clusters):
        raise ValidationError("ranked handoff pose clusters must not be empty")

    active_count = min(len(ranked_clusters), math.ceil(count / 2))
    while sum(len(cluster) for cluster in ranked_clusters[:active_count]) < count:
        active_count += 1
    active_clusters = ranked_clusters[:active_count]
    selected: list[PoseEvidence] = []
    selected_by_cluster: list[list[PoseEvidence]] = [[] for _ in range(active_count)]
    selected_ids: set[str] = set()
    selected_seeds: set[int] = set()

    for cluster_index, cluster in enumerate(active_clusters):
        medoid = min(
            cluster,
            key=lambda pose: (
                sum(distance(pose, other) for other in cluster),
                pose.seed in selected_seeds,
                pose.ranking_score,
                pose.pose_id,
            ),
        )
        selected.append(medoid)
        selected_by_cluster[cluster_index].append(medoid)
        selected_ids.add(medoid.pose_id)
        selected_seeds.add(medoid.seed)

    while len(selected) < count:
        made_progress = False
        for cluster_index, cluster in enumerate(active_clusters):
            previous = selected_by_cluster[cluster_index]
            eligible = tuple(
                pose for pose in cluster if pose.pose_id not in selected_ids
            )
            if not eligible:
                continue
            representative = min(
                eligible,
                key=lambda pose: (
                    -min(distance(pose, other) for other in previous),
                    pose.seed in selected_seeds,
                    pose.ranking_score,
                    pose.pose_id,
                ),
            )
            selected.append(representative)
            previous.append(representative)
            selected_ids.add(representative.pose_id)
            selected_seeds.add(representative.seed)
            made_progress = True
            if len(selected) == count:
                break
        if not made_progress:
            raise ValidationError("handoff pose selection could not satisfy its budget")
    return tuple(selected)


def select_source_seed_representatives(
    *,
    ranked_clusters: tuple[tuple[PoseEvidence, ...], ...],
    count: int,
    distance: Callable[[PoseEvidence, PoseEvidence], float],
) -> tuple[PoseEvidence, ...]:
    """为确认试验选择来自不同 CABS seed 的确定性代表姿态。"""
    if count < 1:
        raise ValidationError("handoff pose count must be positive")
    if any(not cluster for cluster in ranked_clusters):
        raise ValidationError("ranked handoff pose clusters must not be empty")
    seed_poses: dict[int, list[tuple[int, tuple[PoseEvidence, ...], PoseEvidence]]] = {}
    for cluster_index, cluster in enumerate(ranked_clusters):
        for pose in cluster:
            seed_poses.setdefault(pose.seed, []).append((cluster_index, cluster, pose))
    if len(seed_poses) < count:
        raise ValidationError(
            f"handoff has only {len(seed_poses)} distinct source seeds; "
            f"requested {count}"
        )

    representatives: list[tuple[int, float, float, int, str, PoseEvidence]] = []
    for seed, candidates in seed_poses.items():
        cluster_index, cluster, representative = min(
            candidates,
            key=lambda item: (
                item[0],
                sum(distance(item[2], other) for other in item[1]),
                item[2].ranking_score,
                item[2].pose_id,
            ),
        )
        representatives.append(
            (
                cluster_index,
                sum(distance(representative, other) for other in cluster),
                representative.ranking_score,
                seed,
                representative.pose_id,
                representative,
            )
        )
    return tuple(
        item[5]
        for item in sorted(
            representatives,
            key=lambda item: (item[0], item[1], item[2], item[3], item[4]),
        )[:count]
    )


def rank_site_pose_families(
    *,
    site: ReceptorSite | ReportedReceptorSite,
    pose_by_id: dict[str, PoseEvidence],
    peptide_sequence: str,
    reference_receptor_path: Path,
    pose_clustering_rmsd_A: float,
) -> RankedPoseFamilies:
    """用阶段二同一4.0 Å合同重建受体site内的跨source姿态家族。"""
    if not site.supported:
        raise ValidationError(f"unsupported receptor site in candidate: {site.site_id}")
    if pose_clustering_rmsd_A <= 0:
        raise ValidationError("handoff pose selection settings must be positive")
    missing = sorted(set(site.pose_ids) - set(pose_by_id))
    if missing:
        raise ValidationError(f"receptor site refers to unknown pose: {missing[0]}")
    poses = tuple(
        pose_by_id[pose_id]
        for pose_id in site.pose_ids
        if pose_by_id[pose_id].qc_status == "passed"
    )
    reference = read_reference_chains(reference_receptor_path)
    structures: dict[Path, gemmi.Structure] = {}
    coordinates: dict[str, tuple[gemmi.Position, ...]] = {}
    for pose in poses:
        structure = structures.get(pose.model_path)
        if structure is None:
            structure = read_structure(pose.model_path)
            structures[pose.model_path] = structure
        if pose.model_index > len(structure):
            raise ValidationError(
                f"handoff pose model index is outside its structure: {pose.pose_id}"
            )
        receptor, peptide = split_model(
            structure[pose.model_index - 1], peptide_sequence=peptide_sequence
        )
        alignment = align_receptor(receptor=receptor, reference_chains=reference)
        coordinates[pose.pose_id] = tuple(
            gemmi.Position(alignment.transform.apply(required_atom(residue, "CA").pos))
            for residue in peptide
        )

    def distance(first: PoseEvidence, second: PoseEvidence) -> float:
        left = coordinates[first.pose_id]
        right = coordinates[second.pose_id]
        if len(left) != len(right) or not left:
            raise ValidationError("handoff peptide CA identities differ")
        return math.sqrt(
            sum(
                first_position.dist(second_position) ** 2
                for first_position, second_position in zip(left, right, strict=True)
            )
            / len(left)
        )

    clusters = bounded_leader_clusters(
        poses,
        distance=distance,
        identity=lambda pose: pose.pose_id,
        maximum_distance=pose_clustering_rmsd_A,
    )
    ranked_clusters = tuple(
        sorted(
            (tuple(cluster) for cluster in clusters),
            key=lambda cluster: (
                -len({pose.seed for pose in cluster}),
                -len(cluster),
                min(pose.ranking_score for pose in cluster),
                min(pose.pose_id for pose in cluster),
            ),
        )
    )
    return RankedPoseFamilies(ranked_clusters, coordinates)


def select_site_poses(
    *,
    site: ReceptorSite | ReportedReceptorSite,
    pose_by_id: dict[str, PoseEvidence],
    count: int,
    peptide_sequence: str,
    reference_receptor_path: Path,
    pose_clustering_rmsd_A: float,
    require_distinct_source_seeds: bool = False,
) -> tuple[PoseEvidence, ...]:
    """从受支持site的主要姿态家族选择中心和几何多样性起点。"""
    if count < 1:
        raise ValidationError("handoff pose count must be positive")
    families = rank_site_pose_families(
        site=site,
        pose_by_id=pose_by_id,
        peptide_sequence=peptide_sequence,
        reference_receptor_path=reference_receptor_path,
        pose_clustering_rmsd_A=pose_clustering_rmsd_A,
    )
    try:
        if require_distinct_source_seeds:
            return select_source_seed_representatives(
                ranked_clusters=families.clusters,
                count=count,
                distance=families.distance,
            )
        return select_ranked_pose_cluster_representatives(
            ranked_clusters=families.clusters,
            count=count,
            distance=families.distance,
        )
    except ValidationError as exc:
        raise ValidationError(f"{site.site_id}: {exc}") from exc


def _negative_refinement_evidence(
    *, config: AppConfig, run_dirs: tuple[Path, ...]
) -> tuple[dict[str, str], list[dict[str, JsonValue]]]:
    """从既有全原子分析中提取没有支持簇的候选及可审计来源。"""
    exclusions: dict[str, str] = {}
    records: list[dict[str, JsonValue]] = []
    if len(run_dirs) != len(set(run_dirs)):
        raise ValidationError("negative refinement runs must be unique")
    root = (config.paths.outputs_dir / "validation" / "refinements").resolve()
    for run_dir in run_dirs:
        resolved = run_dir.expanduser().resolve()
        if not resolved.is_relative_to(root):
            raise ValidationError("negative refinement run is outside refinements")
        manifest_path = resolved / "refinement_analysis" / "analysis_manifest.json"
        manifest = read_document(manifest_path, name="negative refinement analysis")
        if (
            manifest.get("stage") != "validation_candidate_refinement_analysis"
            or manifest.get("status") != "completed"
            or manifest.get("known_site_information_used") is not False
            or manifest.get("production_qualified") is not False
        ):
            raise ValidationError("negative refinement analysis identity is invalid")
        cluster_path, cluster_hash = validate_record(
            root=manifest_path.parent,
            raw=manifest.get("refined_clusters"),
            name="negative refinement clusters",
        )
        with cluster_path.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"candidate_id", "supported"}
            if reader.fieldnames is None or not required.issubset(reader.fieldnames):
                raise ValidationError("negative refinement cluster table is invalid")
            support: dict[str, bool] = {}
            for row in reader:
                candidate_id = safe_identifier(
                    row["candidate_id"], name="negative candidate ID"
                )
                if row["supported"] not in {"true", "false"}:
                    raise ValidationError(
                        "negative refinement support value is invalid"
                    )
                support[candidate_id] = support.get(candidate_id, False) or (
                    row["supported"] == "true"
                )
        category = manifest.get("evidence_category")
        if not isinstance(category, str):
            raise ValidationError("negative refinement evidence category is invalid")
        reason = (
            "source_seed_confirmation_failed"
            if category == "exploratory_source_seed_confirmation_refinement"
            else "atomic_convergence_not_observed"
        )
        for candidate_id, supported in support.items():
            if not supported:
                exclusions[candidate_id] = reason
        records.append(
            {
                "run_dir": resolved.relative_to(config.paths.outputs_dir).as_posix(),
                "analysis_manifest_sha256": sha256_file(manifest_path),
                "refined_clusters_sha256": cluster_hash,
                "evidence_category": category,
            }
        )
    return exclusions, records


def _family_record(
    families: RankedPoseFamilies,
) -> tuple[tuple[PoseEvidence, ...] | None, dict[str, JsonValue]]:
    """返回首个跨source家族及该site的native-free审计摘要。"""
    records: list[dict[str, JsonValue]] = []
    selected: tuple[PoseEvidence, ...] | None = None
    for index, family in enumerate(families.clusters, 1):
        seeds = tuple(sorted({pose.seed for pose in family}))
        distances = [
            families.distance(first, second)
            for left_index, first in enumerate(family)
            for second in family[left_index + 1 :]
        ]
        records.append(
            {
                "family_rank": index,
                "pose_count": len(family),
                "source_seeds": list(seeds),
                "source_seed_support": len(seeds),
                "median_pairwise_backbone_rmsd_A": (
                    median(distances) if distances else 0.0
                ),
                "best_cabs_score": min(pose.ranking_score for pose in family),
                "representative_pose_id": min(pose.pose_id for pose in family),
            }
        )
        if selected is None and len(seeds) >= 2:
            selected = family
    return selected, {"pose_families": records}


def build_funnel_screening_selection(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    negative_refinement_runs: tuple[Path, ...],
) -> FunnelScreeningSelection:
    """执行零Rosetta的Stage 3A-0审计并选择跨source筛选起点。"""
    try:
        report = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=EXPLORATORY_DISCOVERY_EVIDENCE,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    poses = read_pose_evidence(path=report.pose_path, run_dir=discovery_run_dir)
    pose_by_id = {pose.pose_id: pose for pose in poses}
    if len(pose_by_id) != len(poses):
        raise ValidationError("discovery pose IDs are not unique")
    exclusions, negative_records = _negative_refinement_evidence(
        config=config, run_dirs=negative_refinement_runs
    )
    candidates = tuple(
        candidate
        for candidate in report.candidate_sites.values()
        if candidate.handoff_eligible
    )
    audit_rows: list[dict[str, JsonValue]] = []
    task_pool: dict[str, tuple[CandidateHandoffTask, ...]] = {}
    sortable: dict[str, tuple[int, int, int, str]] = {}
    for candidate in candidates:
        site_records: list[dict[str, JsonValue]] = []
        candidate_tasks: list[CandidateHandoffTask] = []
        supports: list[int] = []
        eligible = candidate.candidate_id not in exclusions
        for site_id in candidate.receptor_site_ids:
            site = report.receptor_sites.get(site_id)
            if site is None:
                raise ValidationError(f"candidate refers to unknown site: {site_id}")
            receptor_path = (
                config.paths.data_dir
                / "receptors"
                / "prepared"
                / f"{site.receptor_id}.cif"
            )
            if not receptor_path.is_file():
                raise ValidationError(
                    f"prepared reference receptor is missing: {site.receptor_id}"
                )
            families = rank_site_pose_families(
                site=site,
                pose_by_id=pose_by_id,
                peptide_sequence=config.chemistry.sequence,
                reference_receptor_path=receptor_path,
                pose_clustering_rmsd_A=config.discovery.cabsdock.pose_clustering_rmsd_A,
            )
            family, family_record = _family_record(families)
            support = 0 if family is None else len({pose.seed for pose in family})
            supports.append(support)
            if family is None:
                eligible = False
                selected_poses: tuple[PoseEvidence, ...] = ()
            else:
                selected_poses = select_source_seed_representatives(
                    ranked_clusters=(family,),
                    count=config.validation.funnel.screening_starts_per_receptor_site,
                    distance=families.distance,
                )
            for pose in selected_poses:
                task_id = safe_identifier(
                    f"{candidate.candidate_id}__{pose.pose_id}", name="handoff task ID"
                )
                candidate_tasks.append(
                    CandidateHandoffTask(
                        task_id,
                        candidate.candidate_id,
                        site.site_id,
                        pose,
                        receptor_path,
                        sha256_file(receptor_path),
                    )
                )
            site_records.append(
                {
                    "receptor_site_id": site.site_id,
                    "receptor_id": site.receptor_id,
                    **family_record,
                    "selected_family_source_support": support,
                    "screening_pose_ids": [pose.pose_id for pose in selected_poses],
                }
            )
        minimum_support = min(supports)
        total_support = sum(supports)
        sortable[candidate.candidate_id] = (
            -minimum_support,
            -total_support,
            candidate.rank_within_tier,
            candidate.candidate_id,
        )
        task_pool[candidate.candidate_id] = tuple(candidate_tasks)
        audit_rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "evidence_tier": candidate.evidence_tier,
                "stage2_rank_within_tier": candidate.rank_within_tier,
                "excluded_by_prior_evidence": candidate.candidate_id in exclusions,
                "exclusion_reason": exclusions.get(candidate.candidate_id),
                "minimum_selected_family_source_support": minimum_support,
                "total_selected_family_source_support": total_support,
                "stage3a0_eligible": eligible,
                "receptor_sites": site_records,
            }
        )

    by_tier: dict[str, list[str]] = {
        "ensemble_consensus": [],
        "conformation_specific": [],
    }
    eligibility = {row["candidate_id"]: row["stage3a0_eligible"] for row in audit_rows}
    for candidate in candidates:
        if (
            eligibility[candidate.candidate_id] is True
            and candidate.evidence_tier in by_tier
        ):
            by_tier[candidate.evidence_tier].append(candidate.candidate_id)
    for values in by_tier.values():
        values.sort(key=sortable.__getitem__)
    ensemble = tuple(
        by_tier["ensemble_consensus"][
            : config.validation.funnel.ensemble_screening_budget
        ]
    )
    specific = tuple(
        by_tier["conformation_specific"][
            : config.validation.funnel.conformation_specific_screening_budget
        ]
    )
    selected_ids = ensemble + specific
    if not selected_ids:
        raise ValidationError("Stage 3A-0 produced no screening candidates")
    tasks = tuple(
        task for candidate_id in selected_ids for task in task_pool[candidate_id]
    )
    audit_by_id = {
        row["candidate_id"]: row
        for row in audit_rows
        if isinstance(row["candidate_id"], str)
    }
    for rank, candidate_id in enumerate(ensemble, 1):
        audit_by_id[candidate_id]["stage3a0_selected"] = True
        audit_by_id[candidate_id]["stage3a0_rank_within_tier"] = rank
    for rank, candidate_id in enumerate(specific, 1):
        audit_by_id[candidate_id]["stage3a0_selected"] = True
        audit_by_id[candidate_id]["stage3a0_rank_within_tier"] = rank
    for row in audit_rows:
        row.setdefault("stage3a0_selected", False)
        row.setdefault("stage3a0_rank_within_tier", None)
    return FunnelScreeningSelection(
        selected_ids,
        {
            "ensemble_consensus_arm": ensemble,
            "conformation_specific_arm": specific,
        },
        tasks,
        {
            "schema": "vela.validation-refinement-funnel-audit/1",
            "phase": "stage3a0_cross_source_pose_family_audit",
            "pose_family_clustering": {
                "algorithm": "bounded_leader_complete_diameter",
                "receptor_alignment": "prepared_receptor_CA",
                "ordered_peptide_atoms": "CA",
                "maximum_rmsd_A": config.discovery.cabsdock.pose_clustering_rmsd_A,
            },
            "ranking": [
                "minimum_selected_family_source_support_desc",
                "total_selected_family_source_support_desc",
                "stage2_rank_within_tier_asc",
                "candidate_id_asc",
            ],
            "budgets": {
                "ensemble_consensus": config.validation.funnel.ensemble_screening_budget,
                "conformation_specific": config.validation.funnel.conformation_specific_screening_budget,
            },
            "negative_refinement_analyses": negative_records,
            "excluded_candidates": exclusions,
            "selected_candidate_ids": list(selected_ids),
            "candidate_audit": audit_rows,
        },
    )


def build_handoff_tasks_for_evidence(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    candidate_ids: tuple[str, ...],
    expected_evidence_category: str,
    require_distinct_source_seeds: bool = False,
) -> tuple[CandidateHandoffTask, ...]:
    """只从受支持 blind site 选择主要簇内非冗余起点。"""
    try:
        report = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=expected_evidence_category,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    poses = read_pose_evidence(path=report.pose_path, run_dir=discovery_run_dir)
    pose_by_id = {pose.pose_id: pose for pose in poses}
    if len(pose_by_id) != len(poses):
        raise ValidationError("discovery pose IDs are not unique")
    candidates = _selected_candidates(
        candidates=report.candidate_sites, requested_ids=candidate_ids
    )
    tasks: list[CandidateHandoffTask] = []
    for candidate in candidates:
        for site_id in candidate.receptor_site_ids:
            site = report.receptor_sites.get(site_id)
            if site is None:
                raise ValidationError(f"candidate refers to unknown site: {site_id}")
            reference_receptor_path = (
                config.paths.data_dir
                / "receptors"
                / "prepared"
                / f"{site.receptor_id}.cif"
            )
            if not reference_receptor_path.is_file():
                raise ValidationError(
                    f"prepared reference receptor is missing: {site.receptor_id}"
                )
            for pose in select_site_poses(
                site=site,
                pose_by_id=pose_by_id,
                count=config.validation.handoff.poses_per_receptor_site,
                peptide_sequence=config.chemistry.sequence,
                reference_receptor_path=reference_receptor_path,
                pose_clustering_rmsd_A=(
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
                require_distinct_source_seeds=require_distinct_source_seeds,
            ):
                _identifier(pose.pose_id, name="pose ID")
                task_id = _identifier(
                    f"{candidate.candidate_id}__{pose.pose_id}", name="handoff task ID"
                )
                tasks.append(
                    CandidateHandoffTask(
                        task_id=task_id,
                        candidate_id=candidate.candidate_id,
                        receptor_site_id=site.site_id,
                        pose=pose,
                        reference_receptor_path=reference_receptor_path,
                        reference_receptor_sha256=sha256_file(reference_receptor_path),
                    )
                )
    return tuple(tasks)


def build_handoff_tasks(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    candidate_ids: tuple[str, ...],
) -> tuple[CandidateHandoffTask, ...]:
    """从正式主发现候选选择全原子交接起点。"""
    return build_handoff_tasks_for_evidence(
        config=config,
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
        expected_evidence_category=MAIN_DISCOVERY_EVIDENCE,
    )


def _cg2all_parameters(config: AppConfig) -> dict[str, JsonValue]:
    settings = config.validation.cg2all
    return {
        "representation": settings.representation,
        "receptor_histidine_state": settings.receptor_histidine_state,
        "device": settings.device,
        "processes": settings.processes,
        "batch_size": settings.batch_size,
        "chain_break_cutoff_A": settings.chain_break_cutoff_A,
        "max_ca_rmsd_A": settings.max_ca_rmsd_A,
    }


def handoff_task_records(
    *,
    tasks: tuple[CandidateHandoffTask, ...],
    discovery_run_dir: Path,
    data_dir: Path,
) -> list[dict[str, JsonValue]]:
    """生成计划写入和执行复核共用的完整任务合同。"""
    discovery_root = discovery_run_dir.resolve()
    records: list[dict[str, JsonValue]] = []
    for task in tasks:
        model_path = task.pose.model_path.resolve()
        if not model_path.is_relative_to(discovery_root):
            raise ValidationError(
                f"handoff source model is outside discovery run: {model_path}"
            )
        records.append(
            {
                "task_id": task.task_id,
                "candidate_id": task.candidate_id,
                "receptor_site_id": task.receptor_site_id,
                "pose_id": task.pose.pose_id,
                "receptor_id": task.pose.receptor_id,
                "target": task.pose.target,
                "seed": task.pose.seed,
                "source_model": {
                    "path": model_path.relative_to(discovery_root).as_posix(),
                    "sha256": task.pose.model_sha256,
                    "model_index": task.pose.model_index,
                },
                "reference_receptor": {
                    "path": task.reference_receptor_path.relative_to(
                        data_dir
                    ).as_posix(),
                    "sha256": task.reference_receptor_sha256,
                },
                "status": "planned",
            }
        )
    return records


def _write_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    run_id: str,
    candidate_ids: tuple[str, ...],
    expected_discovery_evidence: str,
    handoff_evidence: str,
    production_qualified: bool,
    candidate_arms: dict[str, tuple[str, ...]],
    promotion_contract: dict[str, JsonValue],
    require_distinct_source_seeds: bool = False,
    selected_tasks: tuple[CandidateHandoffTask, ...] | None = None,
    selection_details: dict[str, JsonValue] | None = None,
) -> CandidateHandoffPlan:
    """冻结阶段二来源、候选选择、输入哈希和重建工具身份。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    run_dir = config.paths.outputs_dir / "validation" / "handoffs" / run_id
    if run_dir.exists():
        raise ValidationError(f"handoff run directory already exists: {run_dir}")
    flattened_arms = tuple(
        candidate_id
        for arm_candidates in candidate_arms.values()
        for candidate_id in arm_candidates
    )
    if (
        set(flattened_arms) != set(candidate_ids)
        or len(flattened_arms) != len(set(flattened_arms))
        or len(candidate_ids) != len(set(candidate_ids))
    ):
        raise ValidationError("candidate arms must partition the requested candidates")
    tasks = selected_tasks or build_handoff_tasks_for_evidence(
        config=config,
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
        expected_evidence_category=expected_discovery_evidence,
        require_distinct_source_seeds=require_distinct_source_seeds,
    )
    if {task.candidate_id for task in tasks} != set(candidate_ids):
        raise ValidationError("selected handoff tasks do not cover frozen candidates")
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    rosetta = verify_rosetta_scripts_tool(config.validation.rosetta)
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    discovery_files = (
        *(
            discovery_run_dir / name
            for name in ("run_manifest.json", "sampling_manifest.json")
        ),
        discovery_run_dir / "site_analysis" / "analysis_manifest.json",
    )
    if any(not path.is_file() for path in discovery_files):
        raise ValidationError("discovery run manifests are incomplete")
    candidate_evidence = candidate_evidence_records_for_evidence(
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
        expected_evidence_category=expected_discovery_evidence,
    )
    atomic_write_json(
        run_dir / HANDOFF_PLAN_NAME,
        {
            "schema": HANDOFF_PLAN_SCHEMA,
            "stage": "validation_candidate_handoff",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "evidence_category": handoff_evidence,
            "source_evidence_category": expected_discovery_evidence,
            "production_qualified": production_qualified,
            "known_site_information_used": False,
            "chemistry_id": config.chemistry.chemistry_id,
            "software": {
                **vela_software_identity(),
                "cg2all_version": cg2all.version,
                "cg2all_executable_sha256": cg2all.executable_sha256,
                "cg2all_checkpoint_sha256": cg2all.checkpoint_sha256,
                "rosetta_version": rosetta.version,
                "rosetta_scripts_sha256": rosetta.executable_sha256,
            },
            "inputs": {
                "config_snapshot": {
                    "path": snapshot.name,
                    "sha256": sha256_file(snapshot),
                },
                "discovery_run": {
                    "path": discovery_run_dir.resolve()
                    .relative_to(config.paths.outputs_dir.resolve())
                    .as_posix(),
                    "manifests": [
                        {
                            "path": path.relative_to(discovery_run_dir).as_posix(),
                            "sha256": sha256_file(path),
                        }
                        for path in discovery_files
                    ],
                },
            },
            "selection": {
                "requested_candidate_ids": list(candidate_ids),
                "candidate_arms": {
                    name: list(values) for name, values in candidate_arms.items()
                },
                "candidate_evidence": candidate_evidence,
                "poses_per_receptor_site": (
                    config.validation.funnel.screening_starts_per_receptor_site
                    if handoff_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE
                    else config.validation.handoff.poses_per_receptor_site
                ),
                "pose_clustering_rmsd_A": (
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
                "pose_selection": (
                    "top_cross_source_pose_family_seed_medoids"
                    if handoff_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE
                    else (
                        SOURCE_SEED_CONFIRMATION_SELECTION_METHOD
                        if require_distinct_source_seeds
                        else POSE_SELECTION_METHOD
                    )
                ),
                "source_seed_policy": (
                    "two_distinct_source_seeds_within_one_pose_family"
                    if handoff_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE
                    else (
                        SOURCE_SEED_CONFIRMATION_POLICY
                        if require_distinct_source_seeds
                        else POSE_SELECTION_SEED_POLICY
                    )
                ),
                "frozen_evidence_tier_budgets_enforced": True,
                "promotion_contract": promotion_contract,
                "funnel_audit": selection_details,
                "budget": handoff_budget_record(
                    config=config,
                    tasks=tasks,
                    candidate_count=len(candidate_evidence),
                    refinement_seed_count=(
                        config.validation.refinement.seed_batch_sizes[0]
                        if handoff_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE
                        else None
                    ),
                    refinement_mode=(
                        "stage3a_incremental_screening"
                        if handoff_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE
                        else "qualified_full_protocol_only"
                    ),
                ),
            },
            "cg2all_parameters": _cg2all_parameters(config),
            "tasks": handoff_task_records(
                tasks=tasks,
                discovery_run_dir=discovery_run_dir,
                data_dir=config.paths.data_dir,
            ),
        },
    )
    return CandidateHandoffPlan(run_id, run_dir, discovery_run_dir, tasks)


def write_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    run_id: str,
    candidate_ids: tuple[str, ...],
) -> CandidateHandoffPlan:
    """冻结正式主发现候选的全原子交接。"""
    return _write_handoff_plan(
        config=config,
        discovery_run_dir=discovery_run_dir,
        run_id=run_id,
        candidate_ids=candidate_ids,
        expected_discovery_evidence=MAIN_DISCOVERY_EVIDENCE,
        handoff_evidence=MAIN_HANDOFF_EVIDENCE,
        production_qualified=True,
        candidate_arms={"blind_discovery_arm": candidate_ids},
        promotion_contract={
            "mode": "explicit_formal_candidate_selection",
            "qc_metrics_used_for_candidate_reranking": False,
        },
    )


def exploration_promotion_contract(config: AppConfig) -> dict[str, JsonValue]:
    """冻结探索交接后只按 pass/fail 晋级的确定性规则。"""
    minimum_starts = config.validation.analysis.min_refinement_start_support
    minimum_seeds = config.validation.analysis.min_refinement_seed_support
    minimum_source_seeds = config.validation.analysis.min_refinement_source_seed_support
    if minimum_starts is None or minimum_seeds is None or minimum_source_seeds is None:
        raise ValidationError("Stage 3 support thresholds are unresolved")
    return {
        "mode": "qc_only_deterministic_arm_selection",
        "candidate_eligibility": {
            "required_receptor_sites": "all_source_receptor_sites",
            "minimum_passed_starts_per_receptor_site": minimum_starts,
            "allowed_task_status": "execution_completed_and_reconstruction_passed",
            "qc_metrics_used_for_candidate_reranking": False,
        },
        "deep_refinement_selection": {
            "blind_discovery_arm_budget": 2,
            "functional_annotation_arm_budget": 1,
            "order": "frozen_candidate_arm_order",
            "failover": "next_eligible_candidate_within_same_arm_only",
            "insufficient_arm_support": "stop_arm_without_cross_arm_substitution",
        },
        "deep_refinement_success": {
            "minimum_independent_cabs_starts": minimum_starts,
            "minimum_independent_cabs_source_seeds": minimum_source_seeds,
            "minimum_independent_rosetta_seeds": minimum_seeds,
            "supported_cluster_required": True,
            "single_low_score_outlier_is_insufficient": True,
        },
    }


def write_exploration_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    run_id: str,
    blind_candidate_ids: tuple[str, ...],
    functional_candidate_ids: tuple[str, ...],
) -> CandidateHandoffPlan:
    """冻结永久分支的开发候选和48起点QC漏斗。"""
    if len(blind_candidate_ids) < 2 or not functional_candidate_ids:
        raise ValidationError(
            "exploration handoff requires at least two blind and one functional candidate"
        )
    candidate_ids = (*blind_candidate_ids, *functional_candidate_ids)
    evidence = candidate_evidence_records_for_evidence(
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
        expected_evidence_category=EXPLORATORY_DISCOVERY_EVIDENCE,
    )
    by_id: dict[str, dict[str, JsonValue]] = {}
    for record in evidence:
        candidate_id = record.get("candidate_id")
        if not isinstance(candidate_id, str):
            raise ValidationError("exploration candidate evidence is invalid")
        by_id[candidate_id] = record
    if any(
        by_id[candidate_id]["evidence_tier"] != "ensemble_consensus"
        for candidate_id in candidate_ids
    ):
        raise ValidationError(
            "exploration handoff arms currently require ensemble-consensus candidates"
        )

    def rank(candidate_id: str) -> int:
        value = by_id[candidate_id].get("rank_within_tier")
        if not isinstance(value, int) or isinstance(value, bool):
            raise ValidationError("exploration candidate rank is invalid")
        return value

    for arm in (blind_candidate_ids, functional_candidate_ids):
        ranks = tuple(rank(candidate_id) for candidate_id in arm)
        if ranks != tuple(sorted(ranks)):
            raise ValidationError(
                "candidate arm order must follow the frozen blind rank"
            )
    return _write_handoff_plan(
        config=config,
        discovery_run_dir=discovery_run_dir,
        run_id=run_id,
        candidate_ids=candidate_ids,
        expected_discovery_evidence=EXPLORATORY_DISCOVERY_EVIDENCE,
        handoff_evidence=EXPLORATORY_HANDOFF_EVIDENCE,
        production_qualified=False,
        candidate_arms={
            "blind_discovery_arm": blind_candidate_ids,
            "functional_annotation_arm": functional_candidate_ids,
        },
        promotion_contract=exploration_promotion_contract(config),
    )


def funnel_screening_contract(config: AppConfig) -> dict[str, JsonValue]:
    """冻结P15 Stage 3A-0及后续增量精修的完整预算。"""
    funnel = config.validation.funnel
    first, second, deep = config.validation.refinement.seed_batch_sizes
    return {
        "mode": "cross_source_incremental_refinement_funnel",
        "stage3a0": {
            "pose_family_rmsd_A": config.discovery.cabsdock.pose_clustering_rmsd_A,
            "minimum_source_seed_support": 2,
            "ensemble_candidate_budget": funnel.ensemble_screening_budget,
            "conformation_specific_candidate_budget": (
                funnel.conformation_specific_screening_budget
            ),
        },
        "stage3a_screening": {
            "source_starts_per_receptor_site": funnel.screening_starts_per_receptor_site,
            "rosetta_seed_count": first,
            "promotion_status": "cross_source_screening_hit",
            "promotion_budget": funnel.screening_promotion_budget,
        },
        "stage3b_confirmation": {
            "additional_rosetta_seed_count": second,
            "minimum_source_starts": 2,
            "minimum_rosetta_seeds": 2,
            "minimum_source_seed_task_cells": funnel.confirmation_min_task_cells,
            "promotion_budget": funnel.confirmation_promotion_budget,
        },
        "stage3c_deep_confirmation": {
            "additional_rosetta_seed_count": deep,
            "minimum_source_starts": 2,
            "minimum_rosetta_seeds_per_source": 2,
            "source_pool_per_receptor_site": funnel.source_pool_per_receptor_site,
            "final_hypothesis_budget": funnel.final_hypothesis_budget,
        },
        "cross_receptor_atomic_compatibility": {
            "receptor_alignment": "common_receptor_CA",
            "maximum_peptide_backbone_rmsd_A": (
                funnel.cross_receptor_max_backbone_rmsd_A
            ),
            "minimum_receptor_contact_jaccard": (
                funnel.cross_receptor_min_contact_jaccard
            ),
            "both_metrics_required": True,
        },
        "task_reuse": "all_existing_start_seed_tasks_are_reused",
        "production_qualified": False,
    }


def write_funnel_screening_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    negative_refinement_runs: tuple[Path, ...],
    run_id: str,
) -> CandidateHandoffPlan:
    """自动完成Stage 3A-0并冻结Stage 3A全原子筛选起点。"""
    selection = build_funnel_screening_selection(
        config=config,
        discovery_run_dir=discovery_run_dir,
        negative_refinement_runs=negative_refinement_runs,
    )
    return _write_handoff_plan(
        config=config,
        discovery_run_dir=discovery_run_dir,
        run_id=run_id,
        candidate_ids=selection.candidate_ids,
        expected_discovery_evidence=EXPLORATORY_DISCOVERY_EVIDENCE,
        handoff_evidence=FUNNEL_SCREENING_HANDOFF_EVIDENCE,
        production_qualified=False,
        candidate_arms=selection.candidate_arms,
        promotion_contract=funnel_screening_contract(config),
        selected_tasks=selection.tasks,
        selection_details=selection.audit,
    )


def source_seed_confirmation_contract(config: AppConfig) -> dict[str, JsonValue]:
    """冻结单候选跨 CABS source seed 重复性确认规则。"""
    return {
        "mode": "independent_source_seed_confirmation",
        "candidate_selection": "prespecified_development_candidate_only",
        "distinct_source_seeds_per_receptor_site": (
            config.validation.handoff.poses_per_receptor_site
        ),
        "qc_metrics_used_for_candidate_reranking": False,
        "production_qualified": False,
    }


def write_source_seed_confirmation_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    run_id: str,
    candidate_id: str,
) -> CandidateHandoffPlan:
    """冻结一个探索候选的跨 CABS source seed 确认起点。"""
    return _write_handoff_plan(
        config=config,
        discovery_run_dir=discovery_run_dir,
        run_id=run_id,
        candidate_ids=(candidate_id,),
        expected_discovery_evidence=EXPLORATORY_DISCOVERY_EVIDENCE,
        handoff_evidence=SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE,
        production_qualified=False,
        candidate_arms={"source_seed_confirmation_arm": (candidate_id,)},
        promotion_contract=source_seed_confirmation_contract(config),
        require_distinct_source_seeds=True,
    )
