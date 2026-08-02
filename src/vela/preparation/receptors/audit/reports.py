"""审计 TSV schema 与稳定序列化。"""

import csv
import io

STRUCTURE_FIELDS = (
    "receptor_id",
    "pdb_id",
    "target",
    "structure_state",
    "roles",
    "prepare",
    "title",
    "experimental_method",
    "resolution_A",
    "crystal_ph",
    "revision_date",
    "model_count",
    "protein_chain_count",
    "other_target_chain_count",
    "non_target_polymer_chain_count",
    "configured_author_chain",
    "configured_uniprot_accession",
    "identity_status",
    "configured_modeled_residue_count",
    "configured_mean_b_factor",
    "configured_missing_residue_records",
    "configured_sequence_difference_records",
    "alternate_location_atom_records",
    "water_atom_count",
    "nonpolymer_component_ids",
    "asu_other_polymer_contact_atom_pairs",
    "asu_contacted_target_residues",
    "symmetry_contact_atom_pairs",
    "symmetry_contacted_target_residues",
    "minimum_symmetry_contact_distance_A",
    "crystal_contact_review",
    "component_disposition_review",
    "audit_status",
    "raw_cif_sha256",
    "raw_metadata_sha256",
)
CHAIN_FIELDS = (
    "receptor_id",
    "pdb_id",
    "auth_chain_id",
    "label_chain_id",
    "entity_id",
    "polymer_type",
    "modeled_residue_count",
    "modeled_atom_count",
    "mean_b_factor",
    "minimum_occupancy",
    "alternate_location_atom_records",
    "missing_residue_records",
    "sequence_difference_records",
    "uniprot_accessions",
    "configured_chain",
    "configured_accession_match",
)
COMPONENT_FIELDS = (
    "receptor_id",
    "pdb_id",
    "auth_chain_id",
    "label_chain_id",
    "entity_id",
    "component_id",
    "component_name",
    "component_class",
    "atom_count",
    "configured_chain",
    "disposition",
)
MISSING_FIELDS = (
    "receptor_id",
    "pdb_id",
    "auth_chain_id",
    "label_chain_id",
    "residue_name",
    "auth_seq_id",
    "label_seq_id",
    "polymer_flag",
    "occupancy_flag",
    "configured_chain",
)
DIFFERENCE_FIELDS = (
    "receptor_id",
    "pdb_id",
    "auth_chain_id",
    "auth_seq_id",
    "pdb_residue",
    "reference_residue",
    "reference_seq_id",
    "details",
    "difference_class",
    "configured_chain",
)


def tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    """按固定字段顺序序列化 TSV。"""
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()
