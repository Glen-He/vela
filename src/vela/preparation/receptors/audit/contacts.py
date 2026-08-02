"""目标作者链参与的 ASU 和晶体对称接触量化。"""

from pathlib import Path

import gemmi

from vela.preparation.receptors.audit.models import CrystalContacts
from vela.preparation.receptors.models import CrystalContactConfig, ReceptorError


def _residue_id(residue: gemmi.Residue) -> str:
    return f"{residue.seqid.num}{residue.seqid.icode.strip()}"


def crystal_contacts(
    *, path: Path, author_chain_id: str, settings: CrystalContactConfig
) -> CrystalContacts:
    """量化作者链与其他聚合物链及晶体对称像的接触。"""
    structure = gemmi.read_structure(str(path))
    if len(structure) != 1:
        raise ReceptorError(f"crystal-contact analysis requires one model: {path}")
    neighbor_search = gemmi.NeighborSearch(
        structure[0], structure.cell, settings.distance_A
    ).populate(include_h=settings.include_hydrogens)
    contact_search = gemmi.ContactSearch(settings.distance_A)
    contact_search.ignore = gemmi.ContactSearch.Ignore.SameChain
    contact_search.min_occupancy = settings.min_occupancy
    contacts = contact_search.find_contacts(neighbor_search)

    asu_pairs = 0
    symmetry_pairs = 0
    asu_residues: set[str] = set()
    symmetry_residues: set[str] = set()
    symmetry_distances: list[float] = []
    for contact in contacts:
        first = contact.partner1
        second = contact.partner2
        if (
            first.residue.entity_type != gemmi.EntityType.Polymer
            or second.residue.entity_type != gemmi.EntityType.Polymer
        ):
            continue
        first_target = first.chain.name == author_chain_id
        second_target = second.chain.name == author_chain_id
        if not first_target and not second_target:
            continue
        target_residue = first.residue if first_target else second.residue
        if contact.image_idx == 0:
            asu_pairs += 1
            asu_residues.add(_residue_id(target_residue))
        else:
            symmetry_pairs += 1
            symmetry_residues.add(_residue_id(target_residue))
            symmetry_distances.append(contact.dist)
    return CrystalContacts(
        asu_other_polymer_atom_pairs=asu_pairs,
        asu_target_residues=len(asu_residues),
        symmetry_polymer_atom_pairs=symmetry_pairs,
        symmetry_target_residues=len(symmetry_residues),
        minimum_symmetry_distance=(
            min(symmetry_distances) if symmetry_distances else None
        ),
    )
