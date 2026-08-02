"""从原始 mmCIF 建立单链 receptor-only 基础结构。"""

from pathlib import Path

import gemmi

from vela.core.provenance import atomic_write_text, sha256_file
from vela.preparation.receptors.cleaning.conformers import (
    has_altloc,
    resolve_alternate_conformation,
)
from vela.preparation.receptors.cleaning.models import (
    AltlocDecision,
    PreparationResult,
)
from vela.preparation.receptors.models import (
    AltlocConfig,
    ReceptorDefinition,
    ReceptorError,
)


def _source_chain(
    *, definition: ReceptorDefinition, raw_dir: Path
) -> tuple[gemmi.Chain, Path]:
    source_path = raw_dir / f"{definition.pdb_id}.cif"
    if not source_path.is_file():
        raise ReceptorError(f"missing raw receptor input: {source_path}")
    structure = gemmi.read_structure(str(source_path))
    if len(structure) != 1:
        raise ReceptorError(
            f"{definition.receptor_id} must contain exactly one coordinate model"
        )
    matches = [
        chain for chain in structure[0] if chain.name == definition.author_chain_id
    ]
    if len(matches) != 1:
        raise ReceptorError(
            f"{definition.receptor_id} author chain is not unique: "
            f"{definition.author_chain_id}"
        )
    return matches[0], source_path


def _validate_component_policy(
    *, definition: ReceptorDefinition, source_chain: gemmi.Chain
) -> None:
    observed = {
        residue.name
        for residue in source_chain
        if residue.entity_type == gemmi.EntityType.NonPolymer
    }
    declared = set(definition.remove_components) | set(definition.retain_components)
    if observed == declared:
        return
    details: list[str] = []
    if undeclared := sorted(observed - declared):
        details.append("undeclared=" + ",".join(undeclared))
    if absent := sorted(declared - observed):
        details.append("not_observed=" + ",".join(absent))
    raise ReceptorError(
        f"{definition.receptor_id} component policy does not match raw structure: "
        + "; ".join(details)
    )


def prepare_structure(
    *,
    definition: ReceptorDefinition,
    altloc_settings: AltlocConfig,
    raw_dir: Path,
    prepared_dir: Path,
) -> tuple[PreparationResult, list[AltlocDecision]]:
    """制备一个已声明 prepare=true 的 receptor-only 基础结构。"""
    if not definition.prepare:
        raise ReceptorError(f"{definition.receptor_id} is not a preparation target")
    source_chain, source_path = _source_chain(definition=definition, raw_dir=raw_dir)
    _validate_component_policy(definition=definition, source_chain=source_chain)

    cleaned = gemmi.Structure()
    cleaned.name = definition.receptor_id
    cleaned.add_model(gemmi.Model(1))
    cleaned_chain = gemmi.Chain(definition.author_chain_id)
    decisions: list[AltlocDecision] = []
    removed_water = 0
    removed_nonpolymer = 0
    for source_residue in source_chain:
        if source_residue.entity_type == gemmi.EntityType.Water:
            if definition.water_policy == "remove_all":
                removed_water += 1
                continue
        elif source_residue.entity_type == gemmi.EntityType.NonPolymer:
            if source_residue.name in definition.remove_components:
                removed_nonpolymer += 1
                continue
            if source_residue.name not in definition.retain_components:
                raise ReceptorError(
                    f"{definition.receptor_id} has unresolved component: "
                    f"{source_residue.name}"
                )
        elif source_residue.entity_type != gemmi.EntityType.Polymer:
            raise ReceptorError(
                f"{definition.receptor_id} contains unsupported residue entity type"
            )

        residue = source_residue.clone()
        decision = resolve_alternate_conformation(
            receptor_id=definition.receptor_id,
            residue=residue,
            settings=altloc_settings,
        )
        if decision is not None:
            decisions.append(decision)
        cleaned_chain.add_residue(residue)

    cleaned[0].add_chain(cleaned_chain)
    cleaned.setup_entities()
    destination = prepared_dir / f"{definition.receptor_id}.cif"
    atomic_write_text(destination, cleaned.make_mmcif_document().as_string())
    output_chain = _validated_output(destination=destination, definition=definition)
    result = PreparationResult(
        receptor_id=definition.receptor_id,
        pdb_id=definition.pdb_id,
        output_path=destination,
        source_sha256=sha256_file(source_path),
        output_sha256=sha256_file(destination),
        residue_count=sum(
            residue.entity_type == gemmi.EntityType.Polymer for residue in output_chain
        ),
        atom_count=sum(
            len(residue)
            for residue in output_chain
            if residue.entity_type == gemmi.EntityType.Polymer
        ),
        removed_water_residues=removed_water,
        removed_nonpolymer_residues=removed_nonpolymer,
        altloc_residues_resolved=len(decisions),
        altloc_atoms_removed=sum(item.removed_atom_count for item in decisions),
    )
    return result, decisions


def _validated_output(
    *, destination: Path, definition: ReceptorDefinition
) -> gemmi.Chain:
    reloaded = gemmi.read_structure(str(destination))
    if len(reloaded) != 1 or len(reloaded[0]) != 1:
        raise ReceptorError(
            f"prepared receptor is not single-model and single-chain: {destination}"
        )
    output_chain = reloaded[0][0]
    if any(has_altloc(atom) for residue in output_chain for atom in residue):
        raise ReceptorError(
            f"prepared receptor still contains alternate locations: {destination}"
        )
    unexpected_components = {
        residue.name
        for residue in output_chain
        if residue.entity_type == gemmi.EntityType.NonPolymer
        and residue.name not in definition.retain_components
    }
    unexpected_water = (
        any(residue.entity_type == gemmi.EntityType.Water for residue in output_chain)
        and definition.water_policy != "retain_all"
    )
    if unexpected_components or unexpected_water:
        raise ReceptorError(
            f"prepared receptor contains components outside policy: {destination}"
        )
    return output_chain
