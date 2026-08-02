"""阶段一产物完整性检查。"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.preparation.chemistry import (
    ChemistryDefinition,
    ChemistryError,
    assess_chemistry,
    chemistry_record,
    chemistry_record_relative_path,
)
from vela.preparation.receptors.models import (
    ReceptorAuditConfig,
    ReceptorDefinition,
    ReceptorPreparationConfig,
)

type Problem = tuple[str, str]
type JsonObject = dict[str, object]


@dataclass(frozen=True, slots=True)
class PreparationIssue:
    """阶段一阻塞项。"""

    code: str
    message: str


@dataclass(frozen=True, slots=True)
class PreparationReadiness:
    """阶段一准备状态。"""

    issues: tuple[PreparationIssue, ...]
    audited_receptor_ids: tuple[str, ...]
    prepared_receptor_ids: tuple[str, ...]

    @property
    def ready(self) -> bool:
        """是否可以把阶段一产物交给阶段二。"""

        return not self.issues


def _json_object(path: Path) -> JsonObject | None:
    if not path.is_file():
        return None
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, UnicodeError, json.JSONDecodeError):
        return None
    return _object_mapping(raw)


def _object_mapping(value: object) -> JsonObject | None:
    try:
        return object_mapping(value, name="manifest object")
    except TypeError:
        return None


def _object_list(value: object) -> list[object] | None:
    try:
        return object_list(value, name="manifest list")
    except TypeError:
        return None


def _chemistry_problems(
    chemistry: ChemistryDefinition,
    data_dir: Path,
) -> list[Problem]:
    problems: list[Problem] = []
    assessment = assess_chemistry(chemistry)
    if assessment.errors:
        problems.append(
            (
                "chemistry_invalid",
                "Ligand chemistry definition is structurally invalid.",
            )
        )
    if assessment.unresolved_fields:
        problems.append(
            (
                "chemistry_unresolved",
                "Ligand chemistry decisions are unresolved: "
                + ", ".join(assessment.unresolved_fields),
            )
        )

    try:
        record_path = data_dir / chemistry_record_relative_path(chemistry)
    except ChemistryError:
        return problems
    record = _json_object(record_path)
    if record is None:
        problems.append(
            (
                "chemistry_record_missing",
                f"Ligand chemistry record is missing or invalid: {record_path}",
            )
        )
        return problems

    expected = chemistry_record(chemistry, assessment)
    if any(
        record.get(key) != expected[key]
        for key in ("schema", "chemistry", "assessment")
    ):
        problems.append(
            (
                "chemistry_record_stale",
                "Ligand chemistry record does not match the active config.",
            )
        )
    return problems


def _record_problem(
    raw: object,
    data_dir: Path,
    expected_relative_path: Path,
    code_prefix: str,
) -> Problem | None:
    record = _object_mapping(raw)
    if record is None:
        return (
            f"{code_prefix}_record_invalid",
            f"File record is invalid: {expected_relative_path}",
        )
    raw_path = record.get("path")
    digest = record.get("sha256")
    if raw_path != expected_relative_path.as_posix() or not isinstance(digest, str):
        return (
            f"{code_prefix}_record_invalid",
            f"File record does not match the expected path: {expected_relative_path}",
        )
    candidate = (data_dir / expected_relative_path).resolve()
    try:
        candidate.relative_to(data_dir.resolve())
    except ValueError:
        return (
            f"{code_prefix}_path_invalid",
            f"File record escapes the data directory: {expected_relative_path}",
        )
    if not candidate.is_file():
        return (
            f"{code_prefix}_file_missing",
            f"Recorded file is missing: {candidate}",
        )
    if sha256_file(candidate) != digest:
        return (
            f"{code_prefix}_hash_mismatch",
            f"Recorded file hash does not match: {candidate}",
        )
    return None


def _manifest(
    path: Path,
    expected_schema: str,
    code_prefix: str,
) -> tuple[JsonObject | None, list[Problem]]:
    payload = _json_object(path)
    if payload is None:
        return None, [
            (
                f"{code_prefix}_manifest_missing",
                f"Manifest is missing or invalid: {path}",
            )
        ]
    if payload.get("schema") != expected_schema:
        return payload, [
            (
                f"{code_prefix}_manifest_schema_invalid",
                f"Manifest schema is invalid: {path}",
            )
        ]
    return payload, []


def _entries_by_receptor(
    raw: object,
    code_prefix: str,
) -> tuple[dict[str, JsonObject], list[Problem]]:
    rows = _object_list(raw)
    if rows is None:
        return {}, [
            (
                f"{code_prefix}_entries_invalid",
                "Manifest receptor entries are missing or invalid.",
            )
        ]
    entries: dict[str, JsonObject] = {}
    problems: list[Problem] = []
    for raw_row in rows:
        row = _object_mapping(raw_row)
        if row is None or not isinstance(row.get("receptor_id"), str):
            problems.append(
                (
                    f"{code_prefix}_entry_invalid",
                    "Manifest contains an invalid receptor entry.",
                )
            )
            continue
        receptor_id = row["receptor_id"]
        if not isinstance(receptor_id, str):
            continue
        if receptor_id in entries:
            problems.append(
                (
                    f"{code_prefix}_entry_duplicate",
                    f"Manifest contains a duplicate receptor entry: {receptor_id}",
                )
            )
            continue
        entries[receptor_id] = row
    return entries, problems


def _entry_set_problems(
    entries: dict[str, JsonObject],
    expected_ids: set[str],
    code_prefix: str,
) -> list[Problem]:
    problems: list[Problem] = []
    missing = sorted(expected_ids - entries.keys())
    unexpected = sorted(entries.keys() - expected_ids)
    if missing:
        problems.append(
            (
                f"{code_prefix}_entries_missing",
                f"Manifest is missing receptor entries: {', '.join(missing)}",
            )
        )
    if unexpected:
        problems.append(
            (
                f"{code_prefix}_entries_unexpected",
                f"Manifest contains unexpected receptor entries: {', '.join(unexpected)}",
            )
        )
    return problems


def _download_problems(
    receptors: tuple[ReceptorDefinition, ...],
    data_dir: Path,
) -> list[Problem]:
    manifest_path = data_dir / "receptors" / "raw" / "download_manifest.json"
    manifest, problems = _manifest(
        manifest_path,
        "vela.receptor-download-manifest/1",
        "download",
    )
    if manifest is None:
        return problems
    entries, entry_problems = _entries_by_receptor(manifest.get("entries"), "download")
    problems.extend(entry_problems)
    expected_ids = {receptor.receptor_id for receptor in receptors}
    problems.extend(_entry_set_problems(entries, expected_ids, "download"))
    for receptor in receptors:
        entry = entries.get(receptor.receptor_id)
        if entry is None:
            continue
        if entry.get("pdb_id") != receptor.pdb_id:
            problems.append(
                (
                    "download_pdb_id_mismatch",
                    f"Downloaded PDB ID does not match {receptor.receptor_id}.",
                )
            )
        raw_files = _object_list(entry.get("files"))
        if raw_files is None:
            problems.append(
                (
                    "download_files_invalid",
                    f"Download file records are invalid for {receptor.receptor_id}.",
                )
            )
            continue
        records: dict[str, object] = {}
        for raw_file in raw_files:
            file_record = _object_mapping(raw_file)
            if file_record is None or not isinstance(file_record.get("path"), str):
                continue
            path = file_record["path"]
            if isinstance(path, str):
                records[path] = file_record
        expected_paths = (
            Path("receptors") / "raw" / f"{receptor.pdb_id}.cif",
            Path("receptors") / "raw" / f"{receptor.pdb_id}.entry.json",
        )
        if len(raw_files) != len(expected_paths) or set(records) != {
            path.as_posix() for path in expected_paths
        }:
            problems.append(
                (
                    "download_files_incomplete",
                    f"Download file records are incomplete for {receptor.receptor_id}.",
                )
            )
        for expected_path in expected_paths:
            raw_record = records.get(expected_path.as_posix())
            problem = _record_problem(
                raw_record,
                data_dir,
                expected_path,
                "download",
            )
            if problem is not None:
                problems.append(problem)
    return problems


def _audit_problems(
    receptors: tuple[ReceptorDefinition, ...],
    audit: ReceptorAuditConfig,
    data_dir: Path,
) -> tuple[list[Problem], tuple[str, ...]]:
    audit_dir = data_dir / "receptors" / "audit"
    manifest_path = audit_dir / "audit_manifest.json"
    manifest, problems = _manifest(
        manifest_path,
        "vela.receptor-audit-manifest/1",
        "audit",
    )
    if manifest is None:
        return problems, ()
    if manifest.get("parameters") != {
        "crystal_contacts": {
            "distance_A": audit.crystal_contacts.distance_A,
            "min_occupancy": audit.crystal_contacts.min_occupancy,
            "include_hydrogens": audit.crystal_contacts.include_hydrogens,
        },
    }:
        problems.append(
            (
                "audit_parameters_stale",
                "Receptor audit manifest does not match the active config.",
            )
        )
    failures = _object_list(manifest.get("identity_failures"))
    if failures is None:
        problems.append(
            (
                "audit_identity_failures_invalid",
                "Receptor identity failure records are invalid.",
            )
        )
    elif failures:
        problems.append(
            (
                "audit_identity_failed",
                "One or more receptor identity checks failed.",
            )
        )

    output_paths = (
        Path("receptors") / "audit" / "structure_summary.tsv",
        Path("receptors") / "audit" / "chain_summary.tsv",
        Path("receptors") / "audit" / "component_summary.tsv",
        Path("receptors") / "audit" / "missing_summary.tsv",
        Path("receptors") / "audit" / "sequence_difference_summary.tsv",
    )
    raw_outputs = _object_list(manifest.get("outputs"))
    output_records: dict[str, object] = {}
    if raw_outputs is None:
        problems.append(("audit_outputs_invalid", "Audit output records are invalid."))
    else:
        for raw_output in raw_outputs:
            record = _object_mapping(raw_output)
            if record is None or not isinstance(record.get("path"), str):
                continue
            path = record["path"]
            if isinstance(path, str):
                output_records[path] = record
    if raw_outputs is not None and (
        len(raw_outputs) != len(output_paths)
        or set(output_records) != {path.as_posix() for path in output_paths}
    ):
        problems.append(
            ("audit_outputs_incomplete", "Audit output records are incomplete.")
        )
    for expected_path in output_paths:
        problem = _record_problem(
            output_records.get(expected_path.as_posix()),
            data_dir,
            expected_path,
            "audit",
        )
        if problem is not None:
            problems.append(problem)

    audited_ids: list[str] = []
    summary_path = audit_dir / "structure_summary.tsv"
    if summary_path.is_file():
        try:
            with summary_path.open(encoding="utf-8", newline="") as handle:
                rows = csv.DictReader(handle, delimiter="\t")
                for row in rows:
                    receptor_id = row.get("receptor_id")
                    if receptor_id:
                        audited_ids.append(receptor_id)
        except (OSError, UnicodeError, csv.Error):
            problems.append(
                ("audit_summary_invalid", f"Audit summary is invalid: {summary_path}")
            )
    expected_ids = {receptor.receptor_id for receptor in receptors}
    audited_set = set(audited_ids)
    if audited_set != expected_ids or len(audited_ids) != len(audited_set):
        problems.append(
            (
                "audit_receptor_set_mismatch",
                "Audit summary does not contain every configured receptor exactly once.",
            )
        )
    return problems, tuple(sorted(audited_set))


def _preparation_problems(
    receptors: tuple[ReceptorDefinition, ...],
    preparation: ReceptorPreparationConfig,
    data_dir: Path,
) -> tuple[list[Problem], tuple[str, ...]]:
    prepared_dir = data_dir / "receptors" / "prepared"
    manifest_path = prepared_dir / "preparation_manifest.json"
    manifest, problems = _manifest(
        manifest_path,
        "vela.receptor-preparation-manifest/1",
        "preparation",
    )
    if manifest is None:
        return problems, ()
    if manifest.get("parameters") != {
        "altloc": {"preferred_label": preparation.altloc.preferred_label}
    }:
        problems.append(
            (
                "preparation_parameters_stale",
                "Receptor preparation manifest does not match the active config.",
            )
        )
    entries, entry_problems = _entries_by_receptor(
        manifest.get("entries"),
        "preparation",
    )
    problems.extend(entry_problems)
    selected = tuple(receptor for receptor in receptors if receptor.prepare)
    expected_ids = {receptor.receptor_id for receptor in selected}
    problems.extend(_entry_set_problems(entries, expected_ids, "preparation"))
    for receptor in selected:
        entry = entries.get(receptor.receptor_id)
        if entry is None:
            continue
        if entry.get("pdb_id") != receptor.pdb_id:
            problems.append(
                (
                    "preparation_pdb_id_mismatch",
                    f"Prepared PDB ID does not match {receptor.receptor_id}.",
                )
            )
        if entry.get("author_chain_id") != receptor.author_chain_id:
            problems.append(
                (
                    "preparation_chain_mismatch",
                    f"Prepared chain does not match {receptor.receptor_id}.",
                )
            )
        source_path = Path("receptors") / "raw" / f"{receptor.pdb_id}.cif"
        output_path = Path("receptors") / "prepared" / f"{receptor.receptor_id}.cif"
        source_problem = _record_problem(
            entry.get("source"),
            data_dir,
            source_path,
            "preparation_source",
        )
        if source_problem is not None:
            problems.append(source_problem)
        output_problem = _record_problem(
            entry.get("output"),
            data_dir,
            output_path,
            "preparation_output",
        )
        if output_problem is not None:
            problems.append(output_problem)

    report_paths = (
        Path("receptors") / "prepared" / "preparation_summary.tsv",
        Path("receptors") / "prepared" / "altloc_decisions.tsv",
    )
    raw_reports = _object_list(manifest.get("reports"))
    report_records: dict[str, object] = {}
    if raw_reports is None:
        problems.append(
            ("preparation_reports_invalid", "Preparation report records are invalid.")
        )
    else:
        for raw_report in raw_reports:
            record = _object_mapping(raw_report)
            if record is None or not isinstance(record.get("path"), str):
                continue
            path = record["path"]
            if isinstance(path, str):
                report_records[path] = record
    if raw_reports is not None and (
        len(raw_reports) != len(report_paths)
        or set(report_records) != {path.as_posix() for path in report_paths}
    ):
        problems.append(
            (
                "preparation_reports_incomplete",
                "Preparation report records are incomplete.",
            )
        )
    for expected_path in report_paths:
        problem = _record_problem(
            report_records.get(expected_path.as_posix()),
            data_dir,
            expected_path,
            "preparation_report",
        )
        if problem is not None:
            problems.append(problem)
    return problems, tuple(sorted(entries))


def assess_preparation_readiness(
    chemistry: ChemistryDefinition,
    receptors: tuple[ReceptorDefinition, ...],
    audit: ReceptorAuditConfig,
    preparation: ReceptorPreparationConfig,
    data_dir: Path,
) -> PreparationReadiness:
    """检查阶段一声明数据和落盘产物是否一致且完整。"""

    problems = _chemistry_problems(chemistry, data_dir)
    problems.extend(_download_problems(receptors, data_dir))
    audit_problems, audited_ids = _audit_problems(receptors, audit, data_dir)
    problems.extend(audit_problems)
    preparation_problems, prepared_ids = _preparation_problems(
        receptors,
        preparation,
        data_dir,
    )
    problems.extend(preparation_problems)
    unique = {
        (code, message): PreparationIssue(code=code, message=message)
        for code, message in problems
    }
    return PreparationReadiness(
        issues=tuple(unique.values()),
        audited_receptor_ids=audited_ids,
        prepared_receptor_ids=prepared_ids,
    )
