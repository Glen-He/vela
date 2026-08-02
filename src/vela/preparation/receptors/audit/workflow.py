"""多条目结构审计与报告 manifest 编排。"""

from pathlib import Path

from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.preparation.receptors.audit.entry import audit_entry
from vela.preparation.receptors.audit.models import AuditResult
from vela.preparation.receptors.audit.reports import (
    CHAIN_FIELDS,
    COMPONENT_FIELDS,
    DIFFERENCE_FIELDS,
    MISSING_FIELDS,
    STRUCTURE_FIELDS,
    tsv,
)
from vela.preparation.receptors.models import ReceptorAuditConfig, ReceptorDefinition


def audit_receptors(
    *,
    definitions: tuple[ReceptorDefinition, ...],
    settings: ReceptorAuditConfig,
    data_dir: Path,
) -> tuple[AuditResult, ...]:
    """审计全部登记结构并写出稳定 TSV 和 manifest。"""
    raw_dir = data_dir / "receptors" / "raw"
    audit_dir = data_dir / "receptors" / "audit"
    entries = [
        audit_entry(
            definition=definition,
            raw_dir=raw_dir,
            contact_settings=settings.crystal_contacts,
        )
        for definition in definitions
    ]
    outputs = (
        (
            "structure_summary.tsv",
            STRUCTURE_FIELDS,
            [entry.summary for entry in entries],
        ),
        (
            "chain_summary.tsv",
            CHAIN_FIELDS,
            [row for entry in entries for row in entry.chains],
        ),
        (
            "component_summary.tsv",
            COMPONENT_FIELDS,
            [row for entry in entries for row in entry.components],
        ),
        (
            "missing_summary.tsv",
            MISSING_FIELDS,
            [row for entry in entries for row in entry.missing],
        ),
        (
            "sequence_difference_summary.tsv",
            DIFFERENCE_FIELDS,
            [row for entry in entries for row in entry.differences],
        ),
    )
    output_records: list[dict[str, JsonValue]] = []
    for name, fields, rows in outputs:
        path = audit_dir / name
        atomic_write_text(path, tsv(fields, rows))
        output_records.append(
            {
                "path": path.relative_to(data_dir).as_posix(),
                "sha256": sha256_file(path),
            }
        )

    results = tuple(entry.result for entry in entries)
    atomic_write_json(
        audit_dir / "audit_manifest.json",
        {
            "schema": "vela.receptor-audit-manifest/1",
            "generated_at": utc_now(),
            "parameters": {
                "crystal_contacts": {
                    "distance_A": settings.crystal_contacts.distance_A,
                    "min_occupancy": settings.crystal_contacts.min_occupancy,
                    "include_hydrogens": settings.crystal_contacts.include_hydrogens,
                }
            },
            "outputs": output_records,
            "identity_failures": [
                result.receptor_id
                for result in results
                if result.identity_status != "passed"
            ],
            "manual_review_required": [
                result.receptor_id
                for result in results
                if result.manual_review_required
            ],
        },
    )
    return results
