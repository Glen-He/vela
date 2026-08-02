"""替代构象选择和原子规范化。"""

import gemmi

from vela.preparation.receptors.cleaning.models import AltlocDecision
from vela.preparation.receptors.models import AltlocConfig


def has_altloc(atom: gemmi.Atom) -> bool:
    return atom.altloc not in {"\x00", " ", "."}


def _residue_id(residue: gemmi.Residue) -> str:
    return f"{residue.seqid.num}{residue.seqid.icode.strip()}"


def resolve_alternate_conformation(
    *, receptor_id: str, residue: gemmi.Residue, settings: AltlocConfig
) -> AltlocDecision | None:
    """按残基内各标签的平均占有率和声明的平局标签选择构象。"""
    occupancy_by_label: dict[str, list[float]] = {}
    for atom in residue:
        if has_altloc(atom):
            occupancy_by_label.setdefault(atom.altloc, []).append(atom.occ)
    if not occupancy_by_label:
        return None

    def rank(label: str) -> tuple[float, bool, str]:
        occupancies = occupancy_by_label[label]
        return (
            -(sum(occupancies) / len(occupancies)),
            label != settings.preferred_label,
            label,
        )

    selected = min(occupancy_by_label, key=rank)
    remove_indices = [
        index
        for index, atom in enumerate(residue)
        if has_altloc(atom) and atom.altloc != selected
    ]
    for index in reversed(remove_indices):
        del residue[index]
    for atom in residue:
        if has_altloc(atom):
            atom.altloc = "\x00"
    return AltlocDecision(
        receptor_id=receptor_id,
        residue_id=_residue_id(residue),
        residue_name=residue.name,
        available_altlocs=";".join(sorted(occupancy_by_label)),
        selected_altloc=selected,
        removed_atom_count=len(remove_indices),
    )
