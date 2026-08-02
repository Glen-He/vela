"""组合单个 RCSB 条目的结构审计。"""

from pathlib import Path

import gemmi

from vela.core.provenance import sha256_file
from vela.preparation.receptors.audit.cif import (
    category_rows,
    difference_rows,
    first_value,
    joined,
    mapping_by_author_chain,
    missing_rows,
    read_metadata,
    revision_date,
)
from vela.preparation.receptors.audit.contacts import crystal_contacts
from vela.preparation.receptors.audit.entities import extract_entity_rows
from vela.preparation.receptors.audit.models import AuditResult, EntryAudit
from vela.preparation.receptors.models import (
    CrystalContactConfig,
    ReceptorDefinition,
    ReceptorError,
)


def audit_entry(
    *,
    definition: ReceptorDefinition,
    raw_dir: Path,
    contact_settings: CrystalContactConfig,
) -> EntryAudit:
    """读取并审计一个结构条目。"""
    cif_path = raw_dir / f"{definition.pdb_id}.cif"
    metadata_path = raw_dir / f"{definition.pdb_id}.entry.json"
    for path in (cif_path, metadata_path):
        if not path.is_file():
            raise ReceptorError(f"missing raw receptor input: {path}")

    block = gemmi.cif.read_file(str(cif_path)).sole_block()
    metadata = read_metadata(metadata_path)
    missing = missing_rows(block=block, definition=definition)
    differences = difference_rows(block=block, definition=definition)
    entities = extract_entity_rows(
        atom_rows=category_rows(block, "atom_site"),
        entity_rows={row.get("id", ""): row for row in category_rows(block, "entity")},
        polymer_rows={
            row.get("entity_id", ""): row for row in category_rows(block, "entity_poly")
        },
        mappings=mapping_by_author_chain(block),
        missing_rows=missing,
        difference_rows=differences,
        definition=definition,
    )
    contacts = crystal_contacts(
        path=cif_path,
        author_chain_id=definition.author_chain_id,
        settings=contact_settings,
    )
    identity_ok = (
        len(entities.configured_matches) == 1
        and entities.configured_matches[0]["configured_accession_match"] == "true"
    )
    identity_status = "passed" if identity_ok else "failed"
    configured_components = [
        row for row in entities.components if row["configured_chain"] == "true"
    ]
    configured_chain = (
        entities.configured_matches[0] if entities.configured_matches else {}
    )
    other_target_chains = sum(
        row["auth_chain_id"] != definition.author_chain_id
        and definition.uniprot_accession in row["uniprot_accessions"].split(";")
        for row in entities.chains
    )
    non_target_polymer_chains = sum(
        definition.uniprot_accession not in row["uniprot_accessions"].split(";")
        for row in entities.chains
    )
    summary = {
        "receptor_id": definition.receptor_id,
        "pdb_id": definition.pdb_id,
        "target": definition.target,
        "structure_state": definition.structure_state,
        "roles": ";".join(definition.roles),
        "prepare": str(definition.prepare).lower(),
        "title": first_value(block, ("_struct.title",)),
        "experimental_method": first_value(block, ("_exptl.method",)),
        "resolution_A": first_value(
            block,
            ("_refine.ls_d_res_high", "_em_3d_reconstruction.resolution"),
        ),
        "crystal_ph": first_value(block, ("_exptl_crystal_grow.pH",)),
        "revision_date": revision_date(metadata),
        "model_count": str(len({value for value in entities.model_ids if value})),
        "protein_chain_count": str(len(entities.chains)),
        "other_target_chain_count": str(other_target_chains),
        "non_target_polymer_chain_count": str(non_target_polymer_chains),
        "configured_author_chain": definition.author_chain_id,
        "configured_uniprot_accession": definition.uniprot_accession,
        "identity_status": identity_status,
        "configured_modeled_residue_count": configured_chain.get(
            "modeled_residue_count", ""
        ),
        "configured_mean_b_factor": configured_chain.get("mean_b_factor", ""),
        "configured_missing_residue_records": configured_chain.get(
            "missing_residue_records", ""
        ),
        "configured_sequence_difference_records": configured_chain.get(
            "sequence_difference_records", ""
        ),
        "alternate_location_atom_records": configured_chain.get(
            "alternate_location_atom_records", ""
        ),
        "water_atom_count": str(
            sum(
                int(row["atom_count"])
                for row in configured_components
                if row["component_class"] == "water"
            )
        ),
        "nonpolymer_component_ids": joined(
            {
                row["component_id"]
                for row in configured_components
                if row["component_class"] != "water"
            }
        ),
        "asu_other_polymer_contact_atom_pairs": str(
            contacts.asu_other_polymer_atom_pairs
        ),
        "asu_contacted_target_residues": str(contacts.asu_target_residues),
        "symmetry_contact_atom_pairs": str(contacts.symmetry_polymer_atom_pairs),
        "symmetry_contacted_target_residues": str(contacts.symmetry_target_residues),
        "minimum_symmetry_contact_distance_A": (
            f"{contacts.minimum_symmetry_distance:.3f}"
            if contacts.minimum_symmetry_distance is not None
            else ""
        ),
        "crystal_contact_review": "quantified_manual_interpretation_required",
        "component_disposition_review": (
            "pending" if configured_components else "not_applicable"
        ),
        "audit_status": "manual_review_required" if identity_ok else "failed",
        "raw_cif_sha256": sha256_file(cif_path),
        "raw_metadata_sha256": sha256_file(metadata_path),
    }
    return EntryAudit(
        summary=summary,
        chains=entities.chains,
        components=entities.components,
        missing=missing,
        differences=differences,
        result=AuditResult(
            receptor_id=definition.receptor_id,
            pdb_id=definition.pdb_id,
            identity_status=identity_status,
            manual_review_required=True,
        ),
    )
