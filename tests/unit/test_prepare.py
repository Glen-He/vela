import gemmi

from vela.preparation.receptors.cleaning import resolve_alternate_conformation
from vela.preparation.receptors.models import AltlocConfig

ALTLOC_SETTINGS = AltlocConfig(preferred_label="A")


def _atom(*, name: str, altloc: str, occupancy: float) -> gemmi.Atom:
    atom = gemmi.Atom()
    atom.name = name
    atom.element = gemmi.Element("C")
    atom.altloc = altloc
    atom.occ = occupancy
    return atom


def test_altloc_resolution_selects_highest_mean_occupancy() -> None:
    residue = gemmi.Residue()
    residue.name = "VAL"
    residue.seqid = gemmi.SeqId(16, " ")
    residue.add_atom(_atom(name="CB", altloc="A", occupancy=0.4))
    residue.add_atom(_atom(name="CB", altloc="B", occupancy=0.6))
    residue.add_atom(_atom(name="CG1", altloc="A", occupancy=0.4))
    residue.add_atom(_atom(name="CG1", altloc="B", occupancy=0.6))

    decision = resolve_alternate_conformation(
        receptor_id="test", residue=residue, settings=ALTLOC_SETTINGS
    )

    assert decision is not None
    assert decision.selected_altloc == "B"
    assert decision.removed_atom_count == 2
    assert len(residue) == 2
    assert {atom.altloc for atom in residue} == {"\x00"}


def test_altloc_resolution_prefers_a_when_occupancies_tie() -> None:
    residue = gemmi.Residue()
    residue.name = "SER"
    residue.seqid = gemmi.SeqId(3, " ")
    residue.add_atom(_atom(name="OG", altloc="B", occupancy=0.5))
    residue.add_atom(_atom(name="OG", altloc="A", occupancy=0.5))

    decision = resolve_alternate_conformation(
        receptor_id="test", residue=residue, settings=ALTLOC_SETTINGS
    )

    assert decision is not None
    assert decision.selected_altloc == "A"
    assert len(residue) == 1


def test_altloc_resolution_uses_configured_tie_label() -> None:
    residue = gemmi.Residue()
    residue.name = "SER"
    residue.seqid = gemmi.SeqId(3, " ")
    residue.add_atom(_atom(name="OG", altloc="B", occupancy=0.5))
    residue.add_atom(_atom(name="OG", altloc="A", occupancy=0.5))
    settings = AltlocConfig(preferred_label="B")

    decision = resolve_alternate_conformation(
        receptor_id="test", residue=residue, settings=settings
    )

    assert decision is not None
    assert decision.selected_altloc == "B"
