"""阶段二不依赖具体对接引擎的单受体与跨构象 site 分析。"""

from __future__ import annotations

from collections import defaultdict

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
    *, receptor_sites: list[ReceptorSite], settings: SiteAnalysisSettings
) -> list[CandidateSite]:
    contact_limit = settings.contact_jaccard_distance
    position_limit = settings.position_distance_A
    min_receptor_support = settings.min_receptor_support
    if contact_limit is None or position_limit is None or min_receptor_support is None:
        raise DiscoveryError("candidate site settings are unresolved")
    grouped: dict[tuple[str, str], list[ReceptorSite]] = defaultdict(list)
    for site in receptor_sites:
        if site.supported:
            grouped[(site.target, site.coordinate_frame_id)].append(site)
    candidates: list[CandidateSite] = []
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
            candidates.append(
                CandidateSite(
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
                    receptor_support=len(receptor_ids),
                    supported=len(receptor_ids) >= min_receptor_support,
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
    candidates = _candidate_sites(receptor_sites=receptor_sites, settings=settings)
    return SiteAnalysisResult(tuple(receptor_sites), tuple(candidates))
