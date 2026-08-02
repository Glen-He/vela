from vela.preparation.chemistry import (
    ChemistryDefinition,
    DisulfideBond,
    HistidineState,
    assess_chemistry,
    chemistry_record_relative_path,
    validate_chemistry,
)


def _definition(
    *,
    ligand_id: str = "test-peptide",
    disulfide_bonds: tuple[DisulfideBond, ...] = (DisulfideBond(1, 11),),
    n_terminus: str = "NH3+",
    c_terminus: str = "COO-",
    target_ph: float | None = 7.4,
    net_charge: int | None = 1,
    histidines: tuple[HistidineState, ...] = (HistidineState(7, "HIE"),),
    other_modifications_status: str = "none",
    decision_sources: tuple[str, ...] = ("test-source",),
) -> ChemistryDefinition:
    return ChemistryDefinition(
        ligand_id=ligand_id,
        chemistry_id="p15-test",
        sequence="CWMSPRHLGTC",
        chirality="L",
        disulfide_bonds=disulfide_bonds,
        n_terminus=n_terminus,
        c_terminus=c_terminus,
        target_ph=target_ph,
        net_charge=net_charge,
        histidines=histidines,
        other_modifications_status=other_modifications_status,
        other_modifications=(),
        decision_sources=decision_sources,
    )


def test_complete_ligand_definition_is_production_ready() -> None:
    definition = _definition()
    assessment = validate_chemistry(definition)
    assert assessment.schema_valid
    assert assessment.production_ready
    assert assessment.unresolved_fields == ()
    assert chemistry_record_relative_path(definition).as_posix() == (
        "chemistry/test-peptide/chemistry_record.json"
    )


def test_unresolved_scientific_fields_are_not_schema_errors() -> None:
    definition = _definition(
        n_terminus="unresolved",
        c_terminus="unresolved",
        target_ph=None,
        net_charge=None,
        histidines=(HistidineState(7, "unresolved"),),
        other_modifications_status="unresolved",
        decision_sources=(),
    )
    assessment = validate_chemistry(definition)
    assert assessment.schema_valid
    assert not assessment.production_ready
    assert set(assessment.unresolved_fields) == {
        "c_terminus",
        "decision_sources",
        "histidines.7",
        "n_terminus",
        "net_charge",
        "other_modifications_status",
        "target_ph",
    }


def test_disulfide_must_connect_cysteines() -> None:
    assessment = assess_chemistry(_definition(disulfide_bonds=(DisulfideBond(2, 11),)))
    assert "disulfide position 2 is not cysteine" in assessment.errors


def test_ligand_id_must_be_a_safe_stable_identifier() -> None:
    assessment = assess_chemistry(_definition(ligand_id="../P15"))
    assert any(error.startswith("ligand_id must start") for error in assessment.errors)
