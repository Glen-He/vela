"""阶段二分析产生的单受体 site 与跨构象 candidate 报告。"""

import csv
import io
import json
import math
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import (
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.core.run_identity import RUN_ID_PATTERN
from vela.core.typed_data import object_mapping
from vela.discovery.analysis.clustering import (
    candidate_analysis_contract,
    candidate_ranking_key,
)
from vela.discovery.analysis.evidence import SiteAnalysisResult
from vela.discovery.models import (
    SITE_EVIDENCE_TIERS,
    DiscoveryError,
    SiteAnalysisSettings,
)


@dataclass(frozen=True, slots=True)
class ReportedReceptorSite:
    """从冻结报告读取的单受体 site。"""

    site_id: str
    target: str
    receptor_id: str
    coordinate_frame_id: str
    pose_count: int
    pose_ids: tuple[str, ...]
    supporting_seeds: tuple[int, ...]
    representative_pose_id: str
    representative_contacts: frozenset[str]
    representative_position: tuple[float, float, float]
    supported: bool


@dataclass(frozen=True, slots=True)
class ReportedCandidateSite:
    """从冻结报告读取的跨受体 candidate site。"""

    candidate_id: str
    target: str
    coordinate_frame_id: str
    receptor_ids: tuple[str, ...]
    receptor_site_ids: tuple[str, ...]
    representative_site_id: str
    receptor_support: int
    evidence_tier: str
    rank_within_tier: int
    minimum_seed_support: int
    total_seed_support: int
    maximum_normalized_site_distance: float
    minimum_selected_pose_fraction: float
    total_selected_pose_fraction: float
    median_receptor_score_quantile: float
    handoff_eligible: bool


@dataclass(frozen=True, slots=True)
class SiteAnalysisReport:
    """经过 manifest、路径和哈希校验的完整 site 报告。"""

    evidence_category: str
    pose_path: Path
    receptor_sites: dict[str, ReportedReceptorSite]
    candidate_sites: dict[str, ReportedCandidateSite]
    ensemble_candidate_budget: int
    conformation_specific_candidate_budget: int
    manifest_path: Path


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _document(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"site analysis manifest does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name="site analysis manifest")
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid site analysis manifest: {path}") from exc


def _read_table(
    path: Path, *, required: frozenset[str], allow_empty: bool = False
) -> tuple[dict[str, str], ...]:
    if not path.is_file():
        raise DiscoveryError(f"site analysis table does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if (
            fields is None
            or len(fields) != len(set(fields))
            or not required.issubset(fields)
        ):
            raise DiscoveryError(f"site analysis table has invalid columns: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows and not allow_empty:
        raise DiscoveryError(f"site analysis table contains no rows: {path}")
    return rows


def _identifier(value: str, *, name: str) -> str:
    if not RUN_ID_PATTERN.fullmatch(value):
        raise DiscoveryError(f"{name} is not a safe identifier: {value}")
    return value


def _items(value: str, *, name: str) -> tuple[str, ...]:
    items = tuple(item for item in value.split(";") if item)
    if not items or len(items) != len(set(items)):
        raise DiscoveryError(f"{name} must contain unique non-empty values")
    return items


def _boolean(value: str, *, name: str) -> bool:
    if value == "true":
        return True
    if value == "false":
        return False
    raise DiscoveryError(f"{name} must be true or false")


def _integer(value: str, *, name: str) -> int:
    try:
        result = int(value)
    except ValueError as exc:
        raise DiscoveryError(f"{name} must be an integer") from exc
    if result < 0:
        raise DiscoveryError(f"{name} must not be negative")
    return result


def _float(value: str, *, name: str) -> float:
    try:
        result = float(value)
    except ValueError as exc:
        raise DiscoveryError(f"{name} must be a number") from exc
    if not math.isfinite(result):
        raise DiscoveryError(f"{name} must be finite")
    return result


def _contract_float(value: object, *, name: str) -> float:
    if (
        not isinstance(value, int | float)
        or isinstance(value, bool)
        or not math.isfinite(value)
    ):
        raise DiscoveryError(f"{name} must be a finite number")
    return float(value)


def _contract_integer(value: object, *, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 1:
        raise DiscoveryError(f"{name} must be a positive integer")
    return value


def _settings_from_contract(value: object) -> SiteAnalysisSettings:
    try:
        contract = object_mapping(value, name="candidate analysis contract")
        parameters = object_mapping(
            contract.get("parameters"), name="candidate analysis parameters"
        )
    except TypeError as exc:
        raise DiscoveryError("candidate analysis contract is invalid") from exc
    settings = SiteAnalysisSettings(
        contact_jaccard_distance=_contract_float(
            parameters.get("contact_jaccard_distance"),
            name="contact_jaccard_distance",
        ),
        position_distance_A=_contract_float(
            parameters.get("position_distance_A"), name="position_distance_A"
        ),
        min_seed_support=_contract_integer(
            parameters.get("min_seed_support"), name="min_seed_support"
        ),
        min_receptor_support=_contract_integer(
            parameters.get("min_receptor_support"), name="min_receptor_support"
        ),
        min_conformation_specific_seed_support=_contract_integer(
            parameters.get("min_conformation_specific_seed_support"),
            name="min_conformation_specific_seed_support",
        ),
        ensemble_candidate_budget=_contract_integer(
            parameters.get("ensemble_candidate_budget"),
            name="ensemble_candidate_budget",
        ),
        conformation_specific_candidate_budget=_contract_integer(
            parameters.get("conformation_specific_candidate_budget"),
            name="conformation_specific_candidate_budget",
        ),
    )
    if contract != candidate_analysis_contract(settings):
        raise DiscoveryError("candidate analysis contract is not canonical")
    return settings


def _position(value: str, *, name: str) -> tuple[float, float, float]:
    fields = value.split(";")
    if len(fields) != 3:
        raise DiscoveryError(f"{name} must contain three coordinates")
    try:
        return float(fields[0]), float(fields[1]), float(fields[2])
    except ValueError as exc:
        raise DiscoveryError(f"{name} contains an invalid coordinate") from exc


def _report_path(*, run_dir: Path, root: Path, raw: object, expected_name: str) -> Path:
    try:
        record = object_mapping(raw, name=f"site analysis {expected_name}")
    except TypeError as exc:
        raise DiscoveryError(f"invalid site analysis record: {expected_name}") from exc
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise DiscoveryError(f"invalid site analysis record: {expected_name}")
    path = (root / relative).resolve()
    try:
        path.relative_to(run_dir.resolve())
    except ValueError as exc:
        raise DiscoveryError(f"site analysis path escapes its run: {path}") from exc
    if (
        path.name != expected_name
        or not path.is_file()
        or sha256_file(path) != expected_hash
    ):
        raise DiscoveryError(f"site analysis file hash mismatch: {path}")
    return path


def _reported_receptor_sites(path: Path) -> dict[str, ReportedReceptorSite]:
    rows = _read_table(
        path,
        required=frozenset(
            {
                "site_id",
                "target",
                "receptor_id",
                "coordinate_frame_id",
                "pose_count",
                "pose_ids",
                "supporting_seeds",
                "representative_pose_id",
                "representative_contacts",
                "representative_position_A",
                "supported",
            }
        ),
    )
    sites: dict[str, ReportedReceptorSite] = {}
    for row in rows:
        site_id = _identifier(row["site_id"], name="receptor site ID")
        pose_ids = tuple(
            _identifier(item, name=f"{site_id} pose ID")
            for item in _items(row["pose_ids"], name=f"{site_id}.pose_ids")
        )
        representative = _identifier(
            row["representative_pose_id"],
            name=f"{site_id} representative pose ID",
        )
        seeds = tuple(
            _integer(item, name=f"{site_id} supporting seed")
            for item in _items(
                row["supporting_seeds"], name=f"{site_id}.supporting_seeds"
            )
        )
        if representative not in pose_ids or site_id in sites:
            raise DiscoveryError(f"invalid or duplicate receptor site: {site_id}")
        pose_count = _integer(row["pose_count"], name=f"{site_id}.pose_count")
        if pose_count != len(pose_ids):
            raise DiscoveryError(f"{site_id}.pose_count does not match pose_ids")
        sites[site_id] = ReportedReceptorSite(
            site_id=site_id,
            target=row["target"],
            receptor_id=_identifier(row["receptor_id"], name=f"{site_id} receptor ID"),
            coordinate_frame_id=row["coordinate_frame_id"],
            pose_count=pose_count,
            pose_ids=pose_ids,
            supporting_seeds=seeds,
            representative_pose_id=representative,
            representative_contacts=frozenset(
                _items(
                    row["representative_contacts"],
                    name=f"{site_id}.representative_contacts",
                )
            ),
            representative_position=_position(
                row["representative_position_A"],
                name=f"{site_id}.representative_position_A",
            ),
            supported=_boolean(row["supported"], name=f"{site_id}.supported"),
        )
    return sites


def _reported_candidates(path: Path) -> dict[str, ReportedCandidateSite]:
    rows = _read_table(
        path,
        allow_empty=True,
        required=frozenset(
            {
                "candidate_id",
                "target",
                "coordinate_frame_id",
                "receptor_support",
                "receptor_ids",
                "receptor_site_ids",
                "representative_site_id",
                "evidence_tier",
                "rank_within_tier",
                "minimum_seed_support",
                "total_seed_support",
                "maximum_normalized_site_distance",
                "minimum_selected_pose_fraction",
                "total_selected_pose_fraction",
                "median_receptor_score_quantile",
                "handoff_eligible",
            }
        ),
    )
    candidates: dict[str, ReportedCandidateSite] = {}
    for row in rows:
        candidate_id = _identifier(row["candidate_id"], name="candidate ID")
        site_ids = tuple(
            _identifier(item, name=f"{candidate_id} receptor site ID")
            for item in _items(
                row["receptor_site_ids"],
                name=f"{candidate_id}.receptor_site_ids",
            )
        )
        representative = _identifier(
            row["representative_site_id"],
            name=f"{candidate_id} representative site ID",
        )
        receptor_ids = tuple(
            _identifier(item, name=f"{candidate_id} receptor ID")
            for item in _items(row["receptor_ids"], name=f"{candidate_id}.receptor_ids")
        )
        receptor_support = _integer(
            row["receptor_support"], name=f"{candidate_id}.receptor_support"
        )
        evidence_tier = row["evidence_tier"]
        rank_within_tier = _integer(
            row["rank_within_tier"], name=f"{candidate_id}.rank_within_tier"
        )
        if (
            representative not in site_ids
            or candidate_id in candidates
            or receptor_support != len(receptor_ids)
            or evidence_tier not in SITE_EVIDENCE_TIERS
            or rank_within_tier < 1
        ):
            raise DiscoveryError(f"invalid or duplicate candidate site: {candidate_id}")
        candidate = ReportedCandidateSite(
            candidate_id=candidate_id,
            target=row["target"],
            coordinate_frame_id=row["coordinate_frame_id"],
            receptor_ids=receptor_ids,
            receptor_site_ids=site_ids,
            representative_site_id=representative,
            receptor_support=receptor_support,
            evidence_tier=evidence_tier,
            rank_within_tier=rank_within_tier,
            minimum_seed_support=_integer(
                row["minimum_seed_support"],
                name=f"{candidate_id}.minimum_seed_support",
            ),
            total_seed_support=_integer(
                row["total_seed_support"],
                name=f"{candidate_id}.total_seed_support",
            ),
            maximum_normalized_site_distance=_float(
                row["maximum_normalized_site_distance"],
                name=f"{candidate_id}.maximum_normalized_site_distance",
            ),
            minimum_selected_pose_fraction=_float(
                row["minimum_selected_pose_fraction"],
                name=f"{candidate_id}.minimum_selected_pose_fraction",
            ),
            total_selected_pose_fraction=_float(
                row["total_selected_pose_fraction"],
                name=f"{candidate_id}.total_selected_pose_fraction",
            ),
            median_receptor_score_quantile=_float(
                row["median_receptor_score_quantile"],
                name=f"{candidate_id}.median_receptor_score_quantile",
            ),
            handoff_eligible=_boolean(
                row["handoff_eligible"],
                name=f"{candidate_id}.handoff_eligible",
            ),
        )
        if (
            candidate.minimum_seed_support < 1
            or candidate.total_seed_support
            < candidate.minimum_seed_support * candidate.receptor_support
            or not 0.0 <= candidate.maximum_normalized_site_distance <= 1.0
            or not 0.0 < candidate.minimum_selected_pose_fraction <= 1.0
            or not 0.0
            < candidate.total_selected_pose_fraction
            <= float(candidate.receptor_support)
            or not 0.0 <= candidate.median_receptor_score_quantile <= 1.0
        ):
            raise DiscoveryError(f"candidate metrics are invalid: {candidate_id}")
        candidates[candidate_id] = candidate
    return candidates


def read_site_analysis_report(
    *, run_dir: Path, expected_evidence_category: str
) -> SiteAnalysisReport:
    """读取一个类别明确且完整的 site 分析报告。"""
    analysis_dir = run_dir / "site_analysis"
    manifest_path = analysis_dir / "analysis_manifest.json"
    manifest = _document(manifest_path)
    try:
        object_mapping(
            manifest.get("analysis_contract"), name="candidate analysis contract"
        )
    except TypeError as exc:
        raise DiscoveryError("site analysis contract is invalid") from exc
    if (
        manifest.get("schema") != "vela.discovery-site-analysis-manifest/3"
        or manifest.get("evidence_category") != expected_evidence_category
        or manifest.get("known_site_information_used") is not False
    ):
        raise DiscoveryError("site analysis manifest identity is invalid")
    settings = _settings_from_contract(manifest.get("analysis_contract"))
    ensemble_budget = settings.ensemble_candidate_budget
    specific_budget = settings.conformation_specific_candidate_budget
    min_receptor_support = settings.min_receptor_support
    min_specific_seed_support = settings.min_conformation_specific_seed_support
    if (
        ensemble_budget is None
        or specific_budget is None
        or min_receptor_support is None
        or min_specific_seed_support is None
    ):
        raise DiscoveryError("site analysis contract is unresolved")
    pose_path = _report_path(
        run_dir=run_dir,
        root=run_dir,
        raw=manifest.get("pose_evidence"),
        expected_name="pose_evidence.tsv",
    )
    receptor_path = _report_path(
        run_dir=run_dir,
        root=analysis_dir,
        raw=manifest.get("receptor_sites"),
        expected_name="receptor_sites.tsv",
    )
    candidate_path = _report_path(
        run_dir=run_dir,
        root=analysis_dir,
        raw=manifest.get("candidate_sites"),
        expected_name="candidate_sites.tsv",
    )
    receptor_sites = _reported_receptor_sites(receptor_path)
    candidates = _reported_candidates(candidate_path)
    missing = {
        site_id
        for candidate in candidates.values()
        for site_id in candidate.receptor_site_ids
        if site_id not in receptor_sites
    }
    if missing:
        raise DiscoveryError(
            "candidate report refers to unknown receptor sites: "
            + ", ".join(sorted(missing))
        )
    for candidate in candidates.values():
        sites = tuple(
            receptor_sites[site_id] for site_id in candidate.receptor_site_ids
        )
        seed_counts = tuple(len(site.supporting_seeds) for site in sites)
        expected_tier = (
            "ensemble_consensus"
            if candidate.receptor_support >= min_receptor_support
            else (
                "conformation_specific"
                if max(seed_counts) >= min_specific_seed_support
                else "insufficient_evidence"
            )
        )
        if (
            candidate.minimum_seed_support != min(seed_counts)
            or candidate.total_seed_support != sum(seed_counts)
            or candidate.evidence_tier != expected_tier
        ):
            raise DiscoveryError(
                f"candidate evidence metrics are inconsistent: {candidate.candidate_id}"
            )
        selected_pose_totals = {
            receptor_id: sum(
                site.pose_count
                for site in receptor_sites.values()
                if site.supported and site.receptor_id == receptor_id
            )
            for receptor_id in candidate.receptor_ids
        }
        selected_pose_fractions = tuple(
            receptor_sites[site_id].pose_count
            / selected_pose_totals[receptor_sites[site_id].receptor_id]
            for site_id in candidate.receptor_site_ids
        )
        if not math.isclose(
            candidate.minimum_selected_pose_fraction,
            min(selected_pose_fractions),
            abs_tol=5e-8,
        ) or not math.isclose(
            candidate.total_selected_pose_fraction,
            sum(selected_pose_fractions),
            abs_tol=5e-8,
        ):
            raise DiscoveryError(
                f"candidate selected-pose fractions are inconsistent: "
                f"{candidate.candidate_id}"
            )
    ranking_groups = sorted(
        {
            (item.target, item.coordinate_frame_id, item.evidence_tier)
            for item in candidates.values()
        }
    )
    for target, frame_id, tier in ranking_groups:
        tier_candidates = sorted(
            (
                item
                for item in candidates.values()
                if item.target == target
                and item.coordinate_frame_id == frame_id
                and item.evidence_tier == tier
            ),
            key=lambda item: item.rank_within_tier,
        )
        if tuple(item.rank_within_tier for item in tier_candidates) != tuple(
            range(1, len(tier_candidates) + 1)
        ):
            raise DiscoveryError(f"candidate ranks are not contiguous for {tier}")
        expected_order = sorted(
            tier_candidates,
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
        if [item.candidate_id for item in tier_candidates] != [
            item.candidate_id for item in expected_order
        ]:
            raise DiscoveryError(f"candidate ranks disagree with the contract: {tier}")
        budget = ensemble_budget if tier == "ensemble_consensus" else specific_budget
        for candidate in tier_candidates:
            expected_eligible = (
                tier != "insufficient_evidence" and candidate.rank_within_tier <= budget
            )
            if candidate.handoff_eligible != expected_eligible:
                raise DiscoveryError(
                    f"candidate handoff eligibility is inconsistent: "
                    f"{candidate.candidate_id}"
                )
    for candidate in candidates.values():
        representative = receptor_sites[candidate.representative_site_id]
        if (
            representative.target != candidate.target
            or representative.coordinate_frame_id != candidate.coordinate_frame_id
            or tuple(
                sorted(
                    receptor_sites[site_id].receptor_id
                    for site_id in candidate.receptor_site_ids
                )
            )
            != tuple(sorted(candidate.receptor_ids))
        ):
            raise DiscoveryError(
                f"candidate site metadata is inconsistent: {candidate.candidate_id}"
            )
    return SiteAnalysisReport(
        evidence_category=expected_evidence_category,
        pose_path=pose_path,
        receptor_sites=receptor_sites,
        candidate_sites=candidates,
        ensemble_candidate_budget=ensemble_budget,
        conformation_specific_candidate_budget=specific_budget,
        manifest_path=manifest_path,
    )


def write_site_reports(
    *,
    result: SiteAnalysisResult,
    pose_table: Path,
    output_dir: Path,
    evidence_category: str,
    settings: SiteAnalysisSettings,
) -> None:
    """写出不覆盖已有结果的 site 分析报告。"""
    try:
        pose_reference = pose_table.resolve().relative_to(output_dir.parent.resolve())
    except ValueError as exc:
        raise DiscoveryError(
            f"pose evidence is outside its discovery run: {pose_table}"
        ) from exc
    output_paths = (
        output_dir / "receptor_sites.tsv",
        output_dir / "candidate_sites.tsv",
        output_dir / "analysis_manifest.json",
    )
    existing = [path for path in output_paths if path.exists()]
    if existing:
        raise DiscoveryError(
            "site analysis outputs already exist: "
            + ", ".join(str(path) for path in existing)
        )
    if (
        settings.ensemble_candidate_budget is None
        or settings.conformation_specific_candidate_budget is None
    ):
        raise DiscoveryError("candidate delivery budgets are unresolved")
    receptor_fields = (
        "site_id",
        "target",
        "receptor_id",
        "coordinate_frame_id",
        "pose_count",
        "pose_ids",
        "supporting_seeds",
        "representative_pose_id",
        "representative_contacts",
        "representative_position_A",
        "supported",
    )
    receptor_rows = [
        {
            "site_id": site.site_id,
            "target": site.target,
            "receptor_id": site.receptor_id,
            "coordinate_frame_id": site.coordinate_frame_id,
            "pose_count": str(site.pose_count),
            "pose_ids": ";".join(site.pose_ids),
            "supporting_seeds": ";".join(map(str, site.supporting_seeds)),
            "representative_pose_id": site.representative_pose_id,
            "representative_contacts": ";".join(sorted(site.representative_contacts)),
            "representative_position_A": ";".join(
                f"{value:.6f}" for value in site.representative_position
            ),
            "supported": str(site.supported).lower(),
        }
        for site in result.receptor_sites
    ]
    candidate_fields = (
        "candidate_id",
        "target",
        "coordinate_frame_id",
        "receptor_support",
        "receptor_ids",
        "receptor_site_ids",
        "seed_support_by_receptor",
        "representative_site_id",
        "evidence_tier",
        "rank_within_tier",
        "minimum_seed_support",
        "total_seed_support",
        "maximum_normalized_site_distance",
        "minimum_selected_pose_fraction",
        "total_selected_pose_fraction",
        "median_receptor_score_quantile",
        "handoff_eligible",
    )
    candidate_rows = [
        {
            "candidate_id": site.candidate_id,
            "target": site.target,
            "coordinate_frame_id": site.coordinate_frame_id,
            "receptor_support": str(site.receptor_support),
            "receptor_ids": ";".join(site.receptor_ids),
            "receptor_site_ids": ";".join(site.receptor_site_ids),
            "seed_support_by_receptor": ";".join(site.seed_support_by_receptor),
            "representative_site_id": site.representative_site_id,
            "evidence_tier": site.evidence_tier,
            "rank_within_tier": str(site.rank_within_tier),
            "minimum_seed_support": str(site.minimum_seed_support),
            "total_seed_support": str(site.total_seed_support),
            "maximum_normalized_site_distance": (
                f"{site.maximum_normalized_site_distance:.6f}"
            ),
            "minimum_selected_pose_fraction": (
                f"{site.minimum_selected_pose_fraction:.8f}"
            ),
            "total_selected_pose_fraction": (
                f"{site.total_selected_pose_fraction:.8f}"
            ),
            "median_receptor_score_quantile": (
                f"{site.median_receptor_score_quantile:.8f}"
            ),
            "handoff_eligible": str(site.handoff_eligible).lower(),
        }
        for site in result.candidate_sites
    ]
    receptor_path, candidate_path, manifest_path = output_paths
    atomic_write_text(receptor_path, _tsv(receptor_fields, receptor_rows))
    atomic_write_text(candidate_path, _tsv(candidate_fields, candidate_rows))
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.discovery-site-analysis-manifest/3",
            "generated_at": utc_now(),
            "evidence_category": evidence_category,
            "known_site_information_used": False,
            "analysis_contract": candidate_analysis_contract(settings),
            "pose_evidence": {
                "path": pose_reference.as_posix(),
                "sha256": sha256_file(pose_table),
            },
            "receptor_sites": {
                "path": receptor_path.name,
                "sha256": sha256_file(receptor_path),
                "count": len(result.receptor_sites),
                "supported_count": sum(
                    site.supported for site in result.receptor_sites
                ),
            },
            "candidate_sites": {
                "path": candidate_path.name,
                "sha256": sha256_file(candidate_path),
                "count": len(result.candidate_sites),
                "handoff_eligible_count": sum(
                    site.handoff_eligible for site in result.candidate_sites
                ),
                "evidence_tier_counts": {
                    tier: sum(
                        site.evidence_tier == tier for site in result.candidate_sites
                    )
                    for tier in SITE_EVIDENCE_TIERS
                },
            },
        },
    )
