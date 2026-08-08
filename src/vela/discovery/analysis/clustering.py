"""阶段二不依赖具体对接引擎的单受体与跨构象 site 分析。"""

from __future__ import annotations

from bisect import bisect_left, bisect_right
from collections import defaultdict
from dataclasses import dataclass
from statistics import median

from vela.core.provenance import JsonValue
from vela.discovery.analysis.cluster_engine import (
    bounded_leader_clusters,
    complete_linkage,
    normalized_site_distance,
)
from vela.discovery.analysis.evidence import (
    CandidateSite,
    PoseEvidence,
    ReceptorSite,
    SiteAnalysisResult,
)
from vela.discovery.models import DiscoveryError, SiteAnalysisSettings

CANDIDATE_MATCHING_METHOD = "deterministic_complete_linkage_distinct_receptors"
CANDIDATE_RANKING_METHOD = "lexicographic_within_evidence_tier"
CANDIDATE_SCORE_NORMALIZATION = "receptor_supported_site_median_midrank_quantile"


@dataclass(frozen=True, slots=True)
class _CandidateDraft:
    candidate_id: str
    target: str
    coordinate_frame_id: str
    receptor_site_ids: tuple[str, ...]
    receptor_ids: tuple[str, ...]
    seed_support_by_receptor: tuple[str, ...]
    representative_site_id: str
    evidence_tier: str
    minimum_seed_support: int
    total_seed_support: int
    maximum_normalized_site_distance: float
    minimum_selected_pose_fraction: float
    total_selected_pose_fraction: float
    median_receptor_score_quantile: float


def candidate_analysis_contract(
    settings: SiteAnalysisSettings,
) -> dict[str, JsonValue]:
    """返回跨受体匹配、分级、排序和预算的唯一冻结合同。"""
    if not settings.complete:
        raise DiscoveryError("candidate analysis settings are unresolved")
    return {
        "parameters": {
            "contact_jaccard_distance": settings.contact_jaccard_distance,
            "position_distance_A": settings.position_distance_A,
            "min_seed_support": settings.min_seed_support,
            "min_receptor_support": settings.min_receptor_support,
            "min_conformation_specific_seed_support": (
                settings.min_conformation_specific_seed_support
            ),
            "ensemble_candidate_budget": settings.ensemble_candidate_budget,
            "conformation_specific_candidate_budget": (
                settings.conformation_specific_candidate_budget
            ),
        },
        "matching": {
            "method": CANDIDATE_MATCHING_METHOD,
            "input": "supported_receptor_sites",
            "coordinate_scope": "same_target_and_aligned_coordinate_frame",
            "distance": "max_normalized_contact_jaccard_and_centroid_distance",
            "maximum_normalized_distance": 1.0,
            "linkage": "complete",
            "receptor_uniqueness": "at_most_one_site_per_receptor",
            "merge_order": [
                "maximum_intercluster_distance_asc",
                "first_minimum_site_id_asc",
                "second_minimum_site_id_asc",
            ],
            "representative": "minimum_total_normalized_distance_medoid",
        },
        "evidence_tiers": {
            "ensemble_consensus": "minimum_receptor_support",
            "conformation_specific": "minimum_conformation_specific_seed_support",
            "insufficient_evidence": "not_handoff_eligible",
        },
        "ranking": {
            "method": CANDIDATE_RANKING_METHOD,
            "scope": "within_target_coordinate_frame_and_evidence_tier",
            "keys": [
                "minimum_seed_support_desc",
                "total_seed_support_desc",
                "maximum_normalized_site_distance_asc",
                "minimum_selected_pose_fraction_desc",
                "total_selected_pose_fraction_desc",
                "median_receptor_score_quantile_asc",
                "candidate_id_asc",
            ],
            "score_normalization": {
                "method": CANDIDATE_SCORE_NORMALIZATION,
                "scope": "within_receptor_across_supported_sites",
                "raw_score_direction": "lower_is_better",
                "tie_method": "midrank",
                "range": [0.0, 1.0],
                "single_site_value": 0.5,
            },
        },
        "budgets": {
            "unit": "candidate_site",
            "ensemble_consensus": settings.ensemble_candidate_budget,
            "conformation_specific": settings.conformation_specific_candidate_budget,
            "insufficient_evidence": 0,
            "automatic_handoff": False,
        },
        "known_site_information_used": False,
    }


def candidate_ranking_key(
    *,
    minimum_seed_support: int,
    total_seed_support: int,
    maximum_normalized_site_distance: float,
    minimum_selected_pose_fraction: float,
    total_selected_pose_fraction: float,
    median_receptor_score_quantile: float,
    candidate_id: str,
) -> tuple[int, int, float, float, float, float, str]:
    """返回正式候选在证据等级内部使用的确定性字典序键。"""
    return (
        -minimum_seed_support,
        -total_seed_support,
        maximum_normalized_site_distance,
        -minimum_selected_pose_fraction,
        -total_selected_pose_fraction,
        median_receptor_score_quantile,
        candidate_id,
    )


def _midrank_quantiles(scores: dict[str, float]) -> dict[str, float]:
    """把同一受体内的 site 分数转换为保留并列关系的经验分位。"""
    if not scores:
        return {}
    ordered = sorted(scores.values())
    if len(ordered) == 1:
        return {site_id: 0.5 for site_id in scores}
    denominator = len(ordered) - 1
    return {
        site_id: (
            (bisect_left(ordered, score) + bisect_right(ordered, score) - 1)
            / 2
            / denominator
        )
        for site_id, score in scores.items()
    }


def _pose_medoid(
    members: tuple[PoseEvidence, ...], *, settings: SiteAnalysisSettings
) -> PoseEvidence:
    contact_limit = settings.contact_jaccard_distance
    position_limit = settings.position_distance_A
    if contact_limit is None or position_limit is None:
        raise DiscoveryError("site distance thresholds are unresolved")

    def total_distance(candidate: PoseEvidence) -> float:
        return sum(
            normalized_site_distance(
                first_contacts=candidate.contact_residues,
                first_position=candidate.local_position,
                second_contacts=other.contact_residues,
                second_position=other.local_position,
                contact_limit=contact_limit,
                position_limit=position_limit,
            )
            for other in members
        )

    return min(members, key=lambda item: (total_distance(item), item.pose_id))


def _receptor_sites(
    *, poses: tuple[PoseEvidence, ...], settings: SiteAnalysisSettings
) -> list[ReceptorSite]:
    contact_limit = settings.contact_jaccard_distance
    position_limit = settings.position_distance_A
    min_seed_support = settings.min_seed_support
    if contact_limit is None or position_limit is None or min_seed_support is None:
        raise DiscoveryError("site analysis settings are unresolved")
    grouped: dict[tuple[str, str, str], list[PoseEvidence]] = defaultdict(list)
    for pose in poses:
        if pose.qc_status == "passed":
            grouped[(pose.target, pose.receptor_id, pose.coordinate_frame_id)].append(
                pose
            )
    sites: list[ReceptorSite] = []
    for (target, receptor_id, frame_id), members in sorted(grouped.items()):
        clusters = bounded_leader_clusters(
            members,
            distance=lambda first, second: normalized_site_distance(
                first_contacts=first.contact_residues,
                first_position=first.local_position,
                second_contacts=second.contact_residues,
                second_position=second.local_position,
                contact_limit=contact_limit,
                position_limit=position_limit,
            ),
            identity=lambda item: item.pose_id,
            maximum_distance=1.0,
        )
        for index, cluster in enumerate(clusters, 1):
            medoid = _pose_medoid(cluster, settings=settings)
            seeds = tuple(sorted({item.seed for item in cluster}))
            sites.append(
                ReceptorSite(
                    site_id=f"{receptor_id}_S{index:03d}",
                    receptor_id=receptor_id,
                    target=target,
                    coordinate_frame_id=frame_id,
                    pose_ids=tuple(item.pose_id for item in cluster),
                    supporting_seeds=seeds,
                    representative_pose_id=medoid.pose_id,
                    representative_contacts=medoid.contact_residues,
                    representative_position=medoid.local_position,
                    pose_count=len(cluster),
                    supported=len(seeds) >= min_seed_support,
                )
            )
    return sites


def _site_medoid(
    members: tuple[ReceptorSite, ...], *, settings: SiteAnalysisSettings
) -> ReceptorSite:
    contact_limit = settings.contact_jaccard_distance
    position_limit = settings.position_distance_A
    if contact_limit is None or position_limit is None:
        raise DiscoveryError("site distance thresholds are unresolved")

    def total_distance(candidate: ReceptorSite) -> float:
        return sum(
            normalized_site_distance(
                first_contacts=candidate.representative_contacts,
                first_position=candidate.representative_position,
                second_contacts=other.representative_contacts,
                second_position=other.representative_position,
                contact_limit=contact_limit,
                position_limit=position_limit,
            )
            for other in members
        )

    return min(members, key=lambda item: (total_distance(item), item.site_id))


def _candidate_sites(
    *,
    receptor_sites: list[ReceptorSite],
    poses: tuple[PoseEvidence, ...],
    settings: SiteAnalysisSettings,
) -> list[CandidateSite]:
    contact_limit = settings.contact_jaccard_distance
    position_limit = settings.position_distance_A
    min_receptor_support = settings.min_receptor_support
    min_specific_seed_support = settings.min_conformation_specific_seed_support
    ensemble_budget = settings.ensemble_candidate_budget
    specific_budget = settings.conformation_specific_candidate_budget
    if (
        contact_limit is None
        or position_limit is None
        or min_receptor_support is None
        or min_specific_seed_support is None
        or ensemble_budget is None
        or specific_budget is None
    ):
        raise DiscoveryError("candidate site settings are unresolved")
    pose_by_id = {pose.pose_id: pose for pose in poses}
    if len(pose_by_id) != len(poses):
        raise DiscoveryError("pose_id values must be globally unique")
    selected_pose_count_by_receptor: dict[str, int] = defaultdict(int)
    median_score_by_site: dict[str, float] = {}
    score_names_by_receptor: dict[str, set[str]] = defaultdict(set)
    for site in receptor_sites:
        if site.supported:
            selected_pose_count_by_receptor[site.receptor_id] += site.pose_count
            site_poses = tuple(pose_by_id[pose_id] for pose_id in site.pose_ids)
            score_names_by_receptor[site.receptor_id].update(
                pose.score_name for pose in site_poses
            )
            median_score_by_site[site.site_id] = float(
                median(pose.ranking_score for pose in site_poses)
            )
    inconsistent_score_names = {
        receptor_id: names
        for receptor_id, names in score_names_by_receptor.items()
        if len(names) != 1
    }
    if inconsistent_score_names:
        raise DiscoveryError(
            "supported sites for one receptor must share one ranking score: "
            + "; ".join(
                f"{receptor_id}={sorted(names)}"
                for receptor_id, names in sorted(inconsistent_score_names.items())
            )
        )
    score_quantile_by_site: dict[str, float] = {}
    for receptor_id in sorted(selected_pose_count_by_receptor):
        score_quantile_by_site.update(
            _midrank_quantiles(
                {
                    site.site_id: median_score_by_site[site.site_id]
                    for site in receptor_sites
                    if site.supported and site.receptor_id == receptor_id
                }
            )
        )
    grouped: dict[tuple[str, str], list[ReceptorSite]] = defaultdict(list)
    for site in receptor_sites:
        if site.supported:
            grouped[(site.target, site.coordinate_frame_id)].append(site)
    drafts: list[_CandidateDraft] = []
    for (target, frame_id), members in sorted(grouped.items()):
        clusters = complete_linkage(
            members,
            distance=lambda first, second: normalized_site_distance(
                first_contacts=first.representative_contacts,
                first_position=first.representative_position,
                second_contacts=second.representative_contacts,
                second_position=second.representative_position,
                contact_limit=contact_limit,
                position_limit=position_limit,
            ),
            identity=lambda item: item.site_id,
            can_merge=lambda first, second: (
                not (
                    {item.receptor_id for item in first}
                    & {item.receptor_id for item in second}
                )
            ),
        )
        prefix = target.removeprefix("ck2_").upper()
        for index, cluster in enumerate(clusters, 1):
            medoid = _site_medoid(cluster, settings=settings)
            receptor_ids = tuple(sorted(item.receptor_id for item in cluster))
            seed_support = tuple(len(item.supporting_seeds) for item in cluster)
            evidence_tier = (
                "ensemble_consensus"
                if len(receptor_ids) >= min_receptor_support
                else (
                    "conformation_specific"
                    if max(seed_support) >= min_specific_seed_support
                    else "insufficient_evidence"
                )
            )
            pairwise_distances = tuple(
                normalized_site_distance(
                    first_contacts=first.representative_contacts,
                    first_position=first.representative_position,
                    second_contacts=second.representative_contacts,
                    second_position=second.representative_position,
                    contact_limit=contact_limit,
                    position_limit=position_limit,
                )
                for first_index, first in enumerate(cluster)
                for second in cluster[first_index + 1 :]
            )
            selected_pose_fractions = tuple(
                item.pose_count / selected_pose_count_by_receptor[item.receptor_id]
                for item in cluster
            )
            receptor_score_quantiles = tuple(
                score_quantile_by_site[item.site_id] for item in cluster
            )
            drafts.append(
                _CandidateDraft(
                    candidate_id=f"{prefix}_C{index:03d}",
                    target=target,
                    coordinate_frame_id=frame_id,
                    receptor_site_ids=tuple(item.site_id for item in cluster),
                    receptor_ids=receptor_ids,
                    seed_support_by_receptor=tuple(
                        f"{item.receptor_id}:{','.join(map(str, item.supporting_seeds))}"
                        for item in cluster
                    ),
                    representative_site_id=medoid.site_id,
                    evidence_tier=evidence_tier,
                    minimum_seed_support=min(seed_support),
                    total_seed_support=sum(seed_support),
                    maximum_normalized_site_distance=max(
                        pairwise_distances, default=0.0
                    ),
                    minimum_selected_pose_fraction=min(selected_pose_fractions),
                    total_selected_pose_fraction=sum(selected_pose_fractions),
                    median_receptor_score_quantile=float(
                        median(receptor_score_quantiles)
                    ),
                )
            )
    candidates: list[CandidateSite] = []
    ranking_groups = sorted(
        {(item.target, item.coordinate_frame_id, item.evidence_tier) for item in drafts}
    )
    for target, frame_id, evidence_tier in ranking_groups:
        ranked = sorted(
            (
                item
                for item in drafts
                if item.target == target
                and item.coordinate_frame_id == frame_id
                and item.evidence_tier == evidence_tier
            ),
            key=lambda item: candidate_ranking_key(
                minimum_seed_support=item.minimum_seed_support,
                total_seed_support=item.total_seed_support,
                maximum_normalized_site_distance=(
                    item.maximum_normalized_site_distance
                ),
                minimum_selected_pose_fraction=(item.minimum_selected_pose_fraction),
                total_selected_pose_fraction=item.total_selected_pose_fraction,
                median_receptor_score_quantile=(item.median_receptor_score_quantile),
                candidate_id=item.candidate_id,
            ),
        )
        budget = (
            ensemble_budget
            if evidence_tier == "ensemble_consensus"
            else specific_budget
        )
        for rank, item in enumerate(ranked, 1):
            candidates.append(
                CandidateSite(
                    candidate_id=item.candidate_id,
                    target=item.target,
                    coordinate_frame_id=item.coordinate_frame_id,
                    receptor_site_ids=item.receptor_site_ids,
                    receptor_ids=item.receptor_ids,
                    seed_support_by_receptor=item.seed_support_by_receptor,
                    representative_site_id=item.representative_site_id,
                    receptor_support=len(item.receptor_ids),
                    evidence_tier=item.evidence_tier,
                    rank_within_tier=rank,
                    minimum_seed_support=item.minimum_seed_support,
                    total_seed_support=item.total_seed_support,
                    maximum_normalized_site_distance=(
                        item.maximum_normalized_site_distance
                    ),
                    minimum_selected_pose_fraction=(
                        item.minimum_selected_pose_fraction
                    ),
                    total_selected_pose_fraction=item.total_selected_pose_fraction,
                    median_receptor_score_quantile=(
                        item.median_receptor_score_quantile
                    ),
                    handoff_eligible=(
                        evidence_tier != "insufficient_evidence" and rank <= budget
                    ),
                )
            )
    return candidates


def analyze_sites(
    *, poses: tuple[PoseEvidence, ...], settings: SiteAnalysisSettings
) -> SiteAnalysisResult:
    """先在单受体内聚类; 再在同一亚型的受体构象间建立对应 site。"""
    if not settings.complete:
        raise DiscoveryError("site analysis settings are unresolved")
    pose_ids = [pose.pose_id for pose in poses]
    if len(pose_ids) != len(set(pose_ids)):
        raise DiscoveryError("pose_id values must be globally unique")
    frames_by_target: dict[str, set[str]] = defaultdict(set)
    for pose in poses:
        if pose.qc_status == "passed":
            frames_by_target[pose.target].add(pose.coordinate_frame_id)
    inconsistent = {
        target: frames for target, frames in frames_by_target.items() if len(frames) > 1
    }
    if inconsistent:
        raise DiscoveryError(
            "passed poses for one target must use one aligned coordinate frame: "
            + "; ".join(
                f"{target}={sorted(frames)}" for target, frames in inconsistent.items()
            )
        )
    receptor_sites = _receptor_sites(poses=poses, settings=settings)
    candidates = _candidate_sites(
        receptor_sites=receptor_sites,
        poses=poses,
        settings=settings,
    )
    return SiteAnalysisResult(tuple(receptor_sites), tuple(candidates))
