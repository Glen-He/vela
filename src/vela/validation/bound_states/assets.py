"""实验结合态的成对参考、去配体受体和局部对照资产制备。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.preparation.receptors.cleaning.conformers import (
    resolve_alternate_conformation,
)
from vela.preparation.receptors.models import AltlocConfig, ReceptorDefinition
from vela.validation.models import (
    BoundStateDefinition,
    ValidationError,
    ValidationSettings,
)

RECEPTOR_CHAIN = "A"
PEPTIDE_CHAIN = "P"


@dataclass(frozen=True, slots=True)
class BoundStateAsset:
    """一个实验结合态经验证后的阶段三派生资产。"""

    state_id: str
    receptor_id: str
    receptor_path: Path
    pair_reference_path: Path
    local_control_path: Path | None
    disulfide_path: Path | None
    fixed_histidine_pose_indices: tuple[int, ...]
    receptor_residue_count: int
    ligand_residue_count: int
    interface_atom_pairs: int
    interface_receptor_residues: int
    minimum_interface_distance_A: float


def _raw_structure(*, pdb_id: str, data_dir: Path) -> tuple[gemmi.Structure, Path]:
    path = data_dir / "receptors" / "raw" / f"{pdb_id}.cif"
    if not path.is_file():
        raise ValidationError(f"missing raw bound-state structure: {path}")
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid raw bound-state structure: {path}") from exc
    if len(structure) != 1:
        raise ValidationError(f"bound-state structure must contain one model: {path}")
    return structure, path


def _source_chain(
    *, structure: gemmi.Structure, chain_id: str, state_id: str
) -> gemmi.Chain:
    matches = [chain for chain in structure[0] if chain.name == chain_id]
    if len(matches) != 1:
        raise ValidationError(
            f"{state_id} author chain is not unique or missing: {chain_id}"
        )
    return matches[0]


def _copy_chain(
    *,
    source: gemmi.Chain,
    output_name: str,
    state_id: str,
    altloc: AltlocConfig,
    include_non_water: bool,
) -> tuple[gemmi.Chain, int]:
    output = gemmi.Chain(output_name)
    decisions = 0
    for source_residue in source:
        if include_non_water:
            if source_residue.entity_type == gemmi.EntityType.Water:
                continue
        elif source_residue.entity_type != gemmi.EntityType.Polymer:
            continue
        residue = source_residue.clone()
        if (
            resolve_alternate_conformation(
                receptor_id=state_id,
                residue=residue,
                settings=altloc,
            )
            is not None
        ):
            decisions += 1
        output.add_residue(residue)
    if len(output) == 0:
        raise ValidationError(f"{state_id} selected chain contains no usable residues")
    return output, decisions


def _structure(
    *, source: gemmi.Structure, name: str, chains: tuple[gemmi.Chain, ...]
) -> gemmi.Structure:
    output = gemmi.Structure()
    output.name = name
    output.cell = source.cell
    output.spacegroup_hm = source.spacegroup_hm
    output.add_model(gemmi.Model(1))
    for chain in chains:
        output[0].add_chain(chain)
    output.setup_entities()
    return output


def _polymer_residues(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(
        residue for residue in chain if residue.entity_type == gemmi.EntityType.Polymer
    )


def _standard_sequence(chain: gemmi.Chain) -> str:
    sequence: list[str] = []
    for residue in _polymer_residues(chain):
        info = gemmi.find_tabulated_residue(residue.name)
        if not info.is_amino_acid() or len(info.one_letter_code) != 1:
            raise ValidationError(
                f"local peptide control contains non-standard residue: {residue.name}"
            )
        sequence.append(info.one_letter_code)
    return "".join(sequence)


def _heavy_atoms(
    residues: tuple[gemmi.Residue, ...],
) -> tuple[tuple[gemmi.Residue, gemmi.Atom], ...]:
    return tuple(
        (residue, atom)
        for residue in residues
        for atom in residue
        if atom.element.name != "H"
    )


def _interface(
    *, receptor: gemmi.Chain, ligand: gemmi.Chain, threshold_A: float
) -> tuple[int, int, float]:
    receptor_atoms = _heavy_atoms(_polymer_residues(receptor))
    ligand_residues = tuple(
        residue for residue in ligand if residue.entity_type != gemmi.EntityType.Water
    )
    ligand_atoms = _heavy_atoms(ligand_residues)
    if not receptor_atoms or not ligand_atoms:
        raise ValidationError("bound-state interface has no heavy atoms")
    pairs = 0
    receptor_contacts: set[tuple[int | None, str]] = set()
    minimum = float("inf")
    for residue, atom in receptor_atoms:
        for _, ligand_atom in ligand_atoms:
            distance = atom.pos.dist(ligand_atom.pos)
            minimum = min(minimum, distance)
            if distance <= threshold_A:
                pairs += 1
                receptor_contacts.add((residue.seqid.num, residue.seqid.icode.strip()))
    if pairs == 0:
        raise ValidationError("configured receptor and ligand chains do not contact")
    return pairs, len(receptor_contacts), minimum


def _renumbered_standard_peptide(chain: gemmi.Chain) -> gemmi.Chain:
    output = gemmi.Chain(PEPTIDE_CHAIN)
    for index, source_residue in enumerate(_polymer_residues(chain), 1):
        residue = source_residue.clone()
        residue.seqid = gemmi.SeqId(index, " ")
        output.add_residue(residue)
    return output


def _atom(residue: gemmi.Residue, name: str, *, state_id: str) -> gemmi.Atom:
    for atom in residue:
        if atom.name == name:
            return atom
    raise ValidationError(f"{state_id} disulfide endpoint lacks required {name} atom")


def _validate_control_chemistry(
    *,
    definition: BoundStateDefinition,
    peptide: gemmi.Chain,
    settings: ValidationSettings,
) -> None:
    sequence = _standard_sequence(peptide)
    if sequence != definition.ligand_sequence:
        raise ValidationError(
            f"{definition.state_id} ligand sequence mismatch: "
            f"expected {definition.ligand_sequence}, got {sequence}"
        )
    residues = _polymer_residues(peptide)
    for bond in definition.disulfide_bonds:
        if not 1 <= bond.first < bond.second <= len(residues):
            raise ValidationError(
                f"{definition.state_id} disulfide positions are outside ligand"
            )
        first = residues[bond.first - 1]
        second = residues[bond.second - 1]
        if first.name != "CYS" or second.name != "CYS":
            raise ValidationError(
                f"{definition.state_id} disulfide endpoint is not cysteine"
            )
        distance = _atom(first, "SG", state_id=definition.state_id).pos.dist(
            _atom(second, "SG", state_id=definition.state_id).pos
        )
        if not settings.min_disulfide_sg_A <= distance <= settings.max_disulfide_sg_A:
            raise ValidationError(
                f"{definition.state_id} SG distance is outside configured range: "
                f"{distance:.3f} A"
            )


def _write_local_control(
    *,
    definition: BoundStateDefinition,
    source: gemmi.Structure,
    receptor: gemmi.Chain,
    ligand: gemmi.Chain,
    output_dir: Path,
    settings: ValidationSettings,
) -> tuple[Path, Path]:
    peptide = _renumbered_standard_peptide(ligand)
    _validate_control_chemistry(
        definition=definition,
        peptide=peptide,
        settings=settings,
    )
    control = _structure(
        source=source,
        name=definition.state_id,
        chains=(receptor.clone(), peptide),
    )
    control_path = output_dir / "flexpepdock_control.pdb"
    atomic_write_text(control_path, control.make_pdb_string())
    receptor_count = len(_polymer_residues(receptor))
    disulfide_path = output_dir / "fix_disulfide.txt"
    atomic_write_text(
        disulfide_path,
        "".join(
            f"{receptor_count + bond.first} {receptor_count + bond.second}\n"
            for bond in definition.disulfide_bonds
        ),
    )
    return control_path, disulfide_path


def _prepare_one(
    *,
    definition: BoundStateDefinition,
    receptor: ReceptorDefinition,
    settings: ValidationSettings,
    altloc: AltlocConfig,
    data_dir: Path,
) -> tuple[BoundStateAsset, dict[str, JsonValue]]:
    source, source_path = _raw_structure(pdb_id=receptor.pdb_id, data_dir=data_dir)
    receptor_source = _source_chain(
        structure=source,
        chain_id=receptor.author_chain_id,
        state_id=definition.state_id,
    )
    ligand_source = _source_chain(
        structure=source,
        chain_id=definition.ligand_author_chain_id,
        state_id=definition.state_id,
    )
    receptor_chain, receptor_altlocs = _copy_chain(
        source=receptor_source,
        output_name=RECEPTOR_CHAIN,
        state_id=definition.state_id,
        altloc=altloc,
        include_non_water=False,
    )
    ligand_chain, ligand_altlocs = _copy_chain(
        source=ligand_source,
        output_name=PEPTIDE_CHAIN,
        state_id=definition.state_id,
        altloc=altloc,
        include_non_water=True,
    )
    interface_pairs, interface_residues, minimum_distance = _interface(
        receptor=receptor_chain,
        ligand=ligand_chain,
        threshold_A=settings.interface_contact_A,
    )
    output_dir = data_dir / "validation" / "bound_states" / definition.state_id
    receptor_path = output_dir / "receptor_only.cif"
    receptor_structure = _structure(
        source=source,
        name=f"{definition.state_id}_receptor",
        chains=(receptor_chain.clone(),),
    )
    atomic_write_text(
        receptor_path, receptor_structure.make_mmcif_document().as_string()
    )
    pair_path = output_dir / "pair_reference.cif"
    pair_structure = _structure(
        source=source,
        name=definition.state_id,
        chains=(receptor_chain.clone(), ligand_chain.clone()),
    )
    atomic_write_text(pair_path, pair_structure.make_mmcif_document().as_string())
    control_path: Path | None = None
    disulfide_path: Path | None = None
    if definition.local_control_kind == "standard_cyclic_peptide":
        control_path, disulfide_path = _write_local_control(
            definition=definition,
            source=source,
            receptor=receptor_chain,
            ligand=ligand_chain,
            output_dir=output_dir,
            settings=settings,
        )
    asset = BoundStateAsset(
        state_id=definition.state_id,
        receptor_id=definition.receptor_id,
        receptor_path=receptor_path,
        pair_reference_path=pair_path,
        local_control_path=control_path,
        disulfide_path=disulfide_path,
        fixed_histidine_pose_indices=tuple(
            len(_polymer_residues(receptor_chain)) + item.position
            for item in definition.histidines
        ),
        receptor_residue_count=len(_polymer_residues(receptor_chain)),
        ligand_residue_count=len(ligand_chain),
        interface_atom_pairs=interface_pairs,
        interface_receptor_residues=interface_residues,
        minimum_interface_distance_A=minimum_distance,
    )
    outputs = {
        "receptor_only": _file_record(receptor_path, data_dir=data_dir),
        "pair_reference": _file_record(pair_path, data_dir=data_dir),
    }
    if control_path is not None and disulfide_path is not None:
        outputs["flexpepdock_control"] = _file_record(control_path, data_dir=data_dir)
        outputs["fix_disulfide"] = _file_record(disulfide_path, data_dir=data_dir)
    record: dict[str, JsonValue] = {
        "state_id": definition.state_id,
        "receptor_id": definition.receptor_id,
        "pdb_id": receptor.pdb_id,
        "source_receptor_chain": receptor.author_chain_id,
        "source_ligand_chain": definition.ligand_author_chain_id,
        "standard_receptor_chain": RECEPTOR_CHAIN,
        "standard_peptide_chain": PEPTIDE_CHAIN,
        "ligand_id": definition.ligand_id,
        "local_control_kind": definition.local_control_kind,
        "control_chemistry": (
            {
                "sequence": definition.ligand_sequence,
                "disulfide_bonds": [
                    [bond.first, bond.second] for bond in definition.disulfide_bonds
                ],
                "histidines": {
                    str(item.position): item.state for item in definition.histidines
                },
                "fixed_histidine_pose_indices": [
                    len(_polymer_residues(receptor_chain)) + item.position
                    for item in definition.histidines
                ],
            }
            if definition.local_control_kind == "standard_cyclic_peptide"
            else None
        ),
        "selection_reason": definition.selection_reason,
        "source": {
            "path": source_path.relative_to(data_dir).as_posix(),
            "sha256": sha256_file(source_path),
        },
        "processing": {
            "receptor_policy": "polymer_only",
            "ligand_policy": "configured_author_chain_without_water",
            "receptor_altloc_residues_resolved": receptor_altlocs,
            "ligand_altloc_residues_resolved": ligand_altlocs,
        },
        "interface": {
            "contact_distance_A": settings.interface_contact_A,
            "atom_pairs": interface_pairs,
            "receptor_residues": interface_residues,
            "minimum_distance_A": round(minimum_distance, 6),
        },
        "outputs": outputs,
    }
    return asset, record


def _file_record(path: Path, *, data_dir: Path) -> dict[str, JsonValue]:
    return {
        "path": path.relative_to(data_dir).as_posix(),
        "sha256": sha256_file(path),
    }


def prepare_bound_state_assets(
    *,
    settings: ValidationSettings,
    receptors: tuple[ReceptorDefinition, ...],
    altloc: AltlocConfig,
    data_dir: Path,
) -> tuple[BoundStateAsset, ...]:
    """制备已声明结合态, 并把 blind、control 和 site-reference 用途分开。"""
    by_id = {item.receptor_id: item for item in receptors}
    assets: list[BoundStateAsset] = []
    records: list[dict[str, JsonValue]] = []
    for definition in settings.bound_states:
        receptor = by_id.get(definition.receptor_id)
        if receptor is None:
            raise ValidationError(
                f"unknown bound-state receptor: {definition.receptor_id}"
            )
        required_roles = {"bound_state_review", "bound_state_blind_replication"}
        if not required_roles.issubset(receptor.roles):
            raise ValidationError(
                f"{definition.receptor_id} lacks required Stage 3 roles"
            )
        asset, record = _prepare_one(
            definition=definition,
            receptor=receptor,
            settings=settings,
            altloc=altloc,
            data_dir=data_dir,
        )
        assets.append(asset)
        records.append(record)
    root = data_dir / "validation" / "bound_states"
    atomic_write_json(
        root / "preparation_manifest.json",
        {
            "schema": "vela.bound-state-preparation-manifest/1",
            "generated_at": utc_now(),
            "evidence_categories": {
                "receptor_only": "bound_state_blind_replication",
                "standard_cyclic_peptide_control": "method_positive_control",
                "pair_reference": "guided_site_compatibility_reference",
            },
            "parameters": {
                "altloc_preferred_label": altloc.preferred_label,
                "interface_contact_A": settings.interface_contact_A,
                "min_disulfide_sg_A": settings.min_disulfide_sg_A,
                "max_disulfide_sg_A": settings.max_disulfide_sg_A,
            },
            "entries": records,
        },
    )
    return tuple(assets)
