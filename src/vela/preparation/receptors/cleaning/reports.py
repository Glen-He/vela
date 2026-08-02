"""基础受体制备的 TSV 和 provenance manifest。"""

import csv
import io
from pathlib import Path

from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.preparation.receptors.cleaning.models import (
    AltlocDecision,
    PreparationResult,
)
from vela.preparation.receptors.models import (
    ReceptorDefinition,
    ReceptorPreparationConfig,
)


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def write_preparation_reports(
    *,
    results: list[PreparationResult],
    decisions: list[AltlocDecision],
    definitions: tuple[ReceptorDefinition, ...],
    settings: ReceptorPreparationConfig,
    data_dir: Path,
) -> None:
    """写出制备统计、替代构象决定和派生关系。"""
    prepared_dir = data_dir / "receptors" / "prepared"
    summary_fields = (
        "receptor_id",
        "pdb_id",
        "residue_count",
        "atom_count",
        "removed_water_residues",
        "removed_nonpolymer_residues",
        "altloc_residues_resolved",
        "altloc_atoms_removed",
        "source_sha256",
        "output_sha256",
        "output_path",
    )
    summary_rows = [
        {
            "receptor_id": result.receptor_id,
            "pdb_id": result.pdb_id,
            "residue_count": str(result.residue_count),
            "atom_count": str(result.atom_count),
            "removed_water_residues": str(result.removed_water_residues),
            "removed_nonpolymer_residues": str(result.removed_nonpolymer_residues),
            "altloc_residues_resolved": str(result.altloc_residues_resolved),
            "altloc_atoms_removed": str(result.altloc_atoms_removed),
            "source_sha256": result.source_sha256,
            "output_sha256": result.output_sha256,
            "output_path": result.output_path.relative_to(data_dir).as_posix(),
        }
        for result in results
    ]
    altloc_fields = (
        "receptor_id",
        "residue_id",
        "residue_name",
        "available_altlocs",
        "selected_altloc",
        "removed_atom_count",
    )
    altloc_rows = [
        {
            "receptor_id": decision.receptor_id,
            "residue_id": decision.residue_id,
            "residue_name": decision.residue_name,
            "available_altlocs": decision.available_altlocs,
            "selected_altloc": decision.selected_altloc,
            "removed_atom_count": str(decision.removed_atom_count),
        }
        for decision in decisions
    ]
    summary_path = prepared_dir / "preparation_summary.tsv"
    altloc_path = prepared_dir / "altloc_decisions.tsv"
    atomic_write_text(summary_path, _tsv(summary_fields, summary_rows))
    atomic_write_text(altloc_path, _tsv(altloc_fields, altloc_rows))

    by_id = {item.receptor_id: item for item in definitions if item.prepare}
    records: list[dict[str, JsonValue]] = []
    for result in results:
        definition = by_id[result.receptor_id]
        records.append(
            {
                "receptor_id": result.receptor_id,
                "pdb_id": result.pdb_id,
                "author_chain_id": definition.author_chain_id,
                "selection_reason": definition.selection_reason,
                "water_policy": definition.water_policy,
                "remove_components": list(definition.remove_components),
                "retain_components": list(definition.retain_components),
                "source": {
                    "path": f"receptors/raw/{result.pdb_id}.cif",
                    "sha256": result.source_sha256,
                },
                "output": {
                    "path": result.output_path.relative_to(data_dir).as_posix(),
                    "sha256": result.output_sha256,
                },
            }
        )
    atomic_write_json(
        prepared_dir / "preparation_manifest.json",
        {
            "schema": "vela.receptor-preparation-manifest/1",
            "generated_at": utc_now(),
            "parameters": {
                "altloc": {
                    "preferred_label": settings.altloc.preferred_label,
                }
            },
            "entries": records,
            "reports": [
                {
                    "path": summary_path.relative_to(data_dir).as_posix(),
                    "sha256": sha256_file(summary_path),
                },
                {
                    "path": altloc_path.relative_to(data_dir).as_posix(),
                    "sha256": sha256_file(altloc_path),
                },
            ],
        },
    )
