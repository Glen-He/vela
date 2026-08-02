"""apo 主发现与结合态复现之间的跨运行 site 对应。"""

from __future__ import annotations

import csv
import io
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.discovery.analysis.cluster_engine import normalized_site_distance
from vela.discovery.analysis.reports import (
    ReportedReceptorSite,
    SiteAnalysisReport,
    read_site_analysis_report,
)
from vela.discovery.analysis.workflow import discovery_run_target
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.planning import MAIN_DISCOVERY_EVIDENCE
from vela.validation.bound_states.replication import REPLICATION_EVIDENCE
from vela.validation.models import ValidationError


@dataclass(frozen=True, slots=True)
class ReplicationComparisonOutcome:
    """跨状态比较报告的位置和匹配规模。"""

    manifest_path: Path
    candidate_count: int
    matched_candidate_count: int
    replication_only_site_count: int


@dataclass(frozen=True, slots=True)
class CandidateReplication:
    """一个主发现 candidate 与结合态受体 site 的对应证据。"""

    candidate_id: str
    target: str
    main_supported: bool
    main_receptor_support: int
    replication_site_ids: tuple[str, ...]
    replication_state_ids: tuple[str, ...]
    minimum_normalized_distance: float | None

    @property
    def matched(self) -> bool:
        return bool(self.replication_site_ids)


def _resolved_distances(config: AppConfig, *, target_id: str) -> tuple[float, float]:
    target_settings = config.discovery.target(target_id)
    contact_limit = target_settings.analysis.contact_jaccard_distance
    position_limit = target_settings.analysis.position_distance_A
    if contact_limit is None or position_limit is None:
        raise ValidationError("discovery site distance thresholds are unresolved")
    return contact_limit, position_limit


def _distance(
    first: ReportedReceptorSite,
    second: ReportedReceptorSite,
    *,
    contact_limit: float,
    position_limit: float,
) -> float:
    return normalized_site_distance(
        first_contacts=first.representative_contacts,
        first_position=first.representative_position,
        second_contacts=second.representative_contacts,
        second_position=second.representative_position,
        contact_limit=contact_limit,
        position_limit=position_limit,
    )


def compare_sites(
    *,
    main: SiteAnalysisReport,
    replication: SiteAnalysisReport,
    contact_limit: float,
    position_limit: float,
) -> tuple[tuple[CandidateReplication, ...], tuple[ReportedReceptorSite, ...]]:
    """使用冻结的 site 距离定义建立跨受体状态对应, 不改变原始聚类。"""
    supported_replication = tuple(
        site for site in replication.receptor_sites.values() if site.supported
    )
    matches: list[CandidateReplication] = []
    matched_replication_ids: set[str] = set()
    for candidate in sorted(
        main.candidate_sites.values(), key=lambda item: item.candidate_id
    ):
        representative = main.receptor_sites[candidate.representative_site_id]
        distances = sorted(
            [
                (
                    _distance(
                        representative,
                        site,
                        contact_limit=contact_limit,
                        position_limit=position_limit,
                    ),
                    site,
                )
                for site in supported_replication
                if site.target == candidate.target
                and site.coordinate_frame_id == candidate.coordinate_frame_id
            ],
            key=lambda item: (item[0], item[1].site_id),
        )
        corresponding = tuple(site for distance, site in distances if distance <= 1.0)
        matched_replication_ids.update(site.site_id for site in corresponding)
        matches.append(
            CandidateReplication(
                candidate_id=candidate.candidate_id,
                target=candidate.target,
                main_supported=candidate.supported,
                main_receptor_support=candidate.receptor_support,
                replication_site_ids=tuple(site.site_id for site in corresponding),
                replication_state_ids=tuple(
                    sorted({site.receptor_id for site in corresponding})
                ),
                minimum_normalized_distance=(
                    distances[0][0] if corresponding else None
                ),
            )
        )
    replication_only = tuple(
        sorted(
            (
                site
                for site in supported_replication
                if site.site_id not in matched_replication_ids
            ),
            key=lambda item: item.site_id,
        )
    )
    return tuple(matches), replication_only


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _output_record(path: Path) -> dict[str, str]:
    return {"path": path.name, "sha256": sha256_file(path)}


def compare_replication_run(
    *, config: AppConfig, discovery_run_dir: Path, replication_run_dir: Path
) -> ReplicationComparisonOutcome:
    """核对两类分析并生成不可覆盖的跨状态对应报告。"""
    output_dir = replication_run_dir / "comparison"
    if output_dir.exists():
        raise ValidationError(f"replication comparison already exists: {output_dir}")
    target_id = discovery_run_target(discovery_run_dir)
    contact_limit, position_limit = _resolved_distances(config, target_id=target_id)
    try:
        main = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=MAIN_DISCOVERY_EVIDENCE,
        )
        replication = read_site_analysis_report(
            run_dir=replication_run_dir,
            expected_evidence_category=REPLICATION_EVIDENCE,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    matches, replication_only = compare_sites(
        main=main,
        replication=replication,
        contact_limit=contact_limit,
        position_limit=position_limit,
    )
    candidate_path = output_dir / "candidate_replication.tsv"
    replication_only_path = output_dir / "replication_only_sites.tsv"
    manifest_path = output_dir / "comparison_manifest.json"
    atomic_write_text(
        candidate_path,
        _tsv(
            (
                "candidate_id",
                "target",
                "main_supported",
                "main_receptor_support",
                "replication_site_ids",
                "replication_state_ids",
                "replication_state_support",
                "minimum_normalized_distance",
                "matched",
            ),
            [
                {
                    "candidate_id": item.candidate_id,
                    "target": item.target,
                    "main_supported": str(item.main_supported).lower(),
                    "main_receptor_support": str(item.main_receptor_support),
                    "replication_site_ids": ";".join(item.replication_site_ids),
                    "replication_state_ids": ";".join(item.replication_state_ids),
                    "replication_state_support": str(len(item.replication_state_ids)),
                    "minimum_normalized_distance": (
                        ""
                        if item.minimum_normalized_distance is None
                        else f"{item.minimum_normalized_distance:.6f}"
                    ),
                    "matched": str(item.matched).lower(),
                }
                for item in matches
            ],
        ),
    )
    atomic_write_text(
        replication_only_path,
        _tsv(
            (
                "site_id",
                "target",
                "state_id",
                "supporting_seeds",
                "representative_contacts",
                "representative_position_A",
            ),
            [
                {
                    "site_id": site.site_id,
                    "target": site.target,
                    "state_id": site.receptor_id,
                    "supporting_seeds": ";".join(map(str, site.supporting_seeds)),
                    "representative_contacts": ";".join(
                        sorted(site.representative_contacts)
                    ),
                    "representative_position_A": ";".join(
                        f"{value:.6f}" for value in site.representative_position
                    ),
                }
                for site in replication_only
            ],
        ),
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-replication-comparison-manifest/1",
            "stage": "validation_replication_comparison",
            "status": "completed",
            "generated_at": utc_now(),
            "evidence_category": "cross_state_replication_comparison",
            "known_site_information_used": False,
            "classification_applied": False,
            "parameters": {
                "contact_jaccard_distance": contact_limit,
                "position_distance_A": position_limit,
                "normalized_match_limit": 1.0,
            },
            "inputs": {
                "main_discovery_analysis": {
                    "path": main.manifest_path.resolve()
                    .relative_to(config.paths.outputs_dir.resolve())
                    .as_posix(),
                    "sha256": sha256_file(main.manifest_path),
                },
                "bound_state_replication_analysis": {
                    "path": replication.manifest_path.resolve()
                    .relative_to(config.paths.outputs_dir.resolve())
                    .as_posix(),
                    "sha256": sha256_file(replication.manifest_path),
                },
            },
            "candidate_replication": {
                **_output_record(candidate_path),
                "count": len(matches),
                "matched_count": sum(item.matched for item in matches),
            },
            "replication_only_sites": {
                **_output_record(replication_only_path),
                "count": len(replication_only),
            },
        },
    )
    return ReplicationComparisonOutcome(
        manifest_path=manifest_path,
        candidate_count=len(matches),
        matched_candidate_count=sum(item.matched for item in matches),
        replication_only_site_count=len(replication_only),
    )
