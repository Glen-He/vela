"""阶段二分析产生的单受体 site 与跨构象 candidate 报告。"""

import csv
import io
import json
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
from vela.discovery.analysis.evidence import SiteAnalysisResult
from vela.discovery.models import DiscoveryError


@dataclass(frozen=True, slots=True)
class ReportedReceptorSite:
    """从冻结报告读取的单受体 site。"""

    site_id: str
    target: str
    receptor_id: str
    coordinate_frame_id: str
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
    supported: bool


@dataclass(frozen=True, slots=True)
class SiteAnalysisReport:
    """经过 manifest、路径和哈希校验的完整 site 报告。"""

    evidence_category: str
    pose_path: Path
    receptor_sites: dict[str, ReportedReceptorSite]
    candidate_sites: dict[str, ReportedCandidateSite]
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
        sites[site_id] = ReportedReceptorSite(
            site_id=site_id,
            target=row["target"],
            receptor_id=_identifier(row["receptor_id"], name=f"{site_id} receptor ID"),
            coordinate_frame_id=row["coordinate_frame_id"],
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
                "supported",
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
        if (
            representative not in site_ids
            or candidate_id in candidates
            or receptor_support != len(receptor_ids)
        ):
            raise DiscoveryError(f"invalid or duplicate candidate site: {candidate_id}")
        candidates[candidate_id] = ReportedCandidateSite(
            candidate_id=candidate_id,
            target=row["target"],
            coordinate_frame_id=row["coordinate_frame_id"],
            receptor_ids=receptor_ids,
            receptor_site_ids=site_ids,
            representative_site_id=representative,
            receptor_support=receptor_support,
            supported=_boolean(row["supported"], name=f"{candidate_id}.supported"),
        )
    return candidates


def read_site_analysis_report(
    *, run_dir: Path, expected_evidence_category: str
) -> SiteAnalysisReport:
    """读取一个类别明确且完整的 site 分析报告。"""
    analysis_dir = run_dir / "site_analysis"
    manifest_path = analysis_dir / "analysis_manifest.json"
    manifest = _document(manifest_path)
    if (
        manifest.get("schema") != "vela.discovery-site-analysis-manifest/1"
        or manifest.get("evidence_category") != expected_evidence_category
        or manifest.get("known_site_information_used") is not False
    ):
        raise DiscoveryError("site analysis manifest identity is invalid")
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
        manifest_path=manifest_path,
    )


def write_site_reports(
    *,
    result: SiteAnalysisResult,
    pose_table: Path,
    output_dir: Path,
    evidence_category: str,
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
        "supported",
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
            "supported": str(site.supported).lower(),
        }
        for site in result.candidate_sites
    ]
    receptor_path, candidate_path, manifest_path = output_paths
    atomic_write_text(receptor_path, _tsv(receptor_fields, receptor_rows))
    atomic_write_text(candidate_path, _tsv(candidate_fields, candidate_rows))
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.discovery-site-analysis-manifest/1",
            "generated_at": utc_now(),
            "evidence_category": evidence_category,
            "known_site_information_used": False,
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
                "supported_count": sum(
                    site.supported for site in result.candidate_sites
                ),
            },
        },
    )
