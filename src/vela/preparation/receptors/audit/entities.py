"""从 atom_site 组装聚合物链和非聚合物组分审计行。"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

from vela.preparation.receptors.audit.cif import joined, optional_float
from vela.preparation.receptors.models import ReceptorDefinition


@dataclass(frozen=True, slots=True)
class EntityRows:
    """atom_site 聚合后的链、组分和模型事实。"""

    chains: list[dict[str, str]]
    configured_matches: list[dict[str, str]]
    components: list[dict[str, str]]
    model_ids: set[str]


def extract_entity_rows(
    *,
    atom_rows: list[dict[str, str]],
    entity_rows: dict[str, dict[str, str]],
    polymer_rows: dict[str, dict[str, str]],
    mappings: dict[str, list[dict[str, str]]],
    missing_rows: list[dict[str, str]],
    difference_rows: list[dict[str, str]],
    definition: ReceptorDefinition,
) -> EntityRows:
    """按 author/label/entity 身份聚合坐标记录。"""
    chain_atoms: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    component_counts: dict[tuple[str, str, str, str], int] = defaultdict(int)
    model_ids: set[str] = set()
    for row in atom_rows:
        model_ids.add(row.get("pdbx_PDB_model_num", ""))
        key = (
            row.get("auth_asym_id", ""),
            row.get("label_asym_id", ""),
            row.get("label_entity_id", ""),
        )
        if row.get("group_PDB", "") == "ATOM":
            chain_atoms[key].append(row)
        else:
            component_counts[(*key, row.get("label_comp_id", ""))] += 1

    chains: list[dict[str, str]] = []
    configured_matches: list[dict[str, str]] = []
    for (auth_chain, label_chain, entity_id), atoms in sorted(chain_atoms.items()):
        polymer = polymer_rows.get(entity_id, {})
        accessions = {
            item.get("accession", "") for item in mappings.get(auth_chain, [])
        }
        residues = {
            (
                atom.get("auth_seq_id", ""),
                atom.get("pdbx_PDB_ins_code", ""),
                atom.get("label_comp_id", ""),
            )
            for atom in atoms
        }
        b_factors = [
            value
            for atom in atoms
            if (value := optional_float(atom.get("B_iso_or_equiv", ""))) is not None
        ]
        occupancies = [
            value
            for atom in atoms
            if (value := optional_float(atom.get("occupancy", ""))) is not None
        ]
        row = {
            "receptor_id": definition.receptor_id,
            "pdb_id": definition.pdb_id,
            "auth_chain_id": auth_chain,
            "label_chain_id": label_chain,
            "entity_id": entity_id,
            "polymer_type": polymer.get("type", ""),
            "modeled_residue_count": str(len(residues)),
            "modeled_atom_count": str(len(atoms)),
            "mean_b_factor": (
                f"{sum(b_factors) / len(b_factors):.3f}" if b_factors else ""
            ),
            "minimum_occupancy": f"{min(occupancies):.3f}" if occupancies else "",
            "alternate_location_atom_records": str(
                sum(bool(atom.get("label_alt_id", "")) for atom in atoms)
            ),
            "missing_residue_records": str(
                sum(item["auth_chain_id"] == auth_chain for item in missing_rows)
            ),
            "sequence_difference_records": str(
                sum(item["auth_chain_id"] == auth_chain for item in difference_rows)
            ),
            "uniprot_accessions": joined(accessions),
            "configured_chain": str(auth_chain == definition.author_chain_id).lower(),
            "configured_accession_match": str(
                definition.uniprot_accession in accessions
            ).lower(),
        }
        chains.append(row)
        if auth_chain == definition.author_chain_id:
            configured_matches.append(row)

    components: list[dict[str, str]] = []
    for (auth_chain, label_chain, entity_id, component), count in sorted(
        component_counts.items()
    ):
        entity = entity_rows.get(entity_id, {})
        components.append(
            {
                "receptor_id": definition.receptor_id,
                "pdb_id": definition.pdb_id,
                "auth_chain_id": auth_chain,
                "label_chain_id": label_chain,
                "entity_id": entity_id,
                "component_id": component,
                "component_name": entity.get("pdbx_description", ""),
                "component_class": (
                    "water" if component in {"HOH", "DOD"} else "nonpolymer"
                ),
                "atom_count": str(count),
                "configured_chain": str(
                    auth_chain == definition.author_chain_id
                ).lower(),
                "disposition": "pending_manual_review",
            }
        )
    return EntityRows(chains, configured_matches, components, model_ids)
