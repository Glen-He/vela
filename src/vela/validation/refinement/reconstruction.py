"""供阶段三局部精修使用的 CABS pose 全原子重建。"""

from __future__ import annotations

import math
import os
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import atomic_write_text, sha256_file
from vela.preparation.chemistry import ChemistryDefinition
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import Cg2AllSettings, ValidationError

PACKAGE_NAME_PATTERN = re.compile(r"^Name:\s*cg2all\s*$", re.MULTILINE)
PACKAGE_VERSION_PATTERN = re.compile(r"^Version:\s*(\S+)\s*$", re.MULTILINE)
HISTIDINE_NAMES = frozenset({"HIS", "HSD", "HSE"})


@dataclass(frozen=True, slots=True)
class Cg2AllToolInfo:
    """通过来源、版本、模型和命令能力检查的 cg2all 身份。"""

    version: str
    executable_sha256: str
    checkpoint_sha256: str


@dataclass(frozen=True, slots=True)
class CgInputAsset:
    """一个标准化为 A/P 链和显式二硫键的 CA/SC 输入。"""

    path: Path
    receptor_residue_count: int
    peptide_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class ReconstructionMetrics:
    """重建前后 CA 坐标的直接保真度。"""

    receptor_ca_rmsd_A: float
    peptide_ca_rmsd_A: float


@dataclass(frozen=True, slots=True)
class TopologyReconstructionAssessment:
    """粗粒化 pose 经全原子重建和 Rosetta prepack 后的组合判据。"""

    receptor_ca_rmsd_A: float
    peptide_pose_ca_rmsd_A: float
    peptide_internal_ca_rmsd_A: float
    ligand_centroid_displacement_A: float
    receptor_contact_retention_fraction: float
    disulfide_sg_distances_A: tuple[float, ...]
    min_interchain_heavy_atom_distance_A: float
    min_nonlocal_peptide_heavy_atom_distance_A: float
    failures: tuple[str, ...]

    @property
    def successful(self) -> bool:
        """全部预注册化学与几何判据是否通过。"""
        return not self.failures


def _named_atom(residue: gemmi.Residue, atom_name: str) -> gemmi.Atom | None:
    """按名称查找原子并规避 Gemmi find_atom 的空值声明偏差。"""
    for atom in residue:
        if atom.name == atom_name:
            return atom
    return None


def verify_cg2all_tool(settings: Cg2AllSettings) -> Cg2AllToolInfo:
    """核对同一环境中的入口、包元数据和显式 checkpoint。"""
    if not settings.executable.is_file() or not os.access(settings.executable, os.X_OK):
        raise ValidationError(
            f"cg2all executable is not runnable: {settings.executable}"
        )
    for label, path in (
        ("package metadata", settings.package_metadata),
        ("checkpoint", settings.checkpoint),
    ):
        if not path.is_file():
            raise ValidationError(f"cg2all {label} is missing: {path}")
    environment_root = settings.executable.resolve().parent.parent
    for label, path in (
        ("package metadata", settings.package_metadata),
        ("checkpoint", settings.checkpoint),
    ):
        try:
            path.resolve().relative_to(environment_root)
        except ValueError as exc:
            raise ValidationError(
                f"cg2all {label} is outside the executable environment"
            ) from exc
    metadata = settings.package_metadata.read_text(encoding="utf-8")
    version_match = PACKAGE_VERSION_PATTERN.search(metadata)
    if PACKAGE_NAME_PATTERN.search(metadata) is None or version_match is None:
        raise ValidationError("cg2all package metadata is not recognized")
    version = version_match.group(1)
    if version != settings.expected_version:
        raise ValidationError(
            f"cg2all version mismatch: expected {settings.expected_version}, got {version}"
        )
    checkpoint_hash = sha256_file(settings.checkpoint)
    if checkpoint_hash != settings.checkpoint_sha256:
        raise ValidationError("cg2all checkpoint hash does not match config")
    result = subprocess.run(
        [str(settings.executable), "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    help_text = result.stdout + result.stderr
    required = ("--cg", "--ckpt", "--device", "--batch", "--proc")
    if result.returncode != 0 or any(option not in help_text for option in required):
        raise ValidationError("cg2all help check lacks required command options")
    return Cg2AllToolInfo(
        version=version,
        executable_sha256=sha256_file(settings.executable),
        checkpoint_sha256=checkpoint_hash,
    )


def build_cg2all_command(
    *, settings: Cg2AllSettings, input_path: Path, output_path: Path
) -> tuple[str, ...]:
    """建立不下载模型且完整暴露运行参数的单 pose 重建命令。"""
    return (
        str(settings.executable),
        "-p",
        str(input_path),
        "-o",
        str(output_path),
        "--cg",
        settings.representation,
        "--ckpt",
        str(settings.checkpoint),
        "--chain-break-cutoff",
        str(settings.chain_break_cutoff_A),
        "--device",
        settings.device,
        "--batch",
        str(settings.batch_size),
        "--proc",
        str(settings.processes),
        "--standard-name",
    )


def _amino_acid_residues(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(
        residue
        for residue in chain
        if _named_atom(residue, "CA") is not None
        and (
            gemmi.find_tabulated_residue(residue.name).is_amino_acid()
            or residue.name in HISTIDINE_NAMES
        )
    )


def _sequence(residues: tuple[gemmi.Residue, ...]) -> str:
    sequence: list[str] = []
    for residue in residues:
        if residue.name in HISTIDINE_NAMES:
            sequence.append("H")
            continue
        code = gemmi.find_tabulated_residue(residue.name).one_letter_code
        if len(code) != 1:
            raise ValidationError(
                f"reconstruction contains unsupported residue: {residue.name}"
            )
        sequence.append(code)
    return "".join(sequence)


def _source_chains(
    *, path: Path, model_index: int, peptide_sequence: str
) -> tuple[tuple[gemmi.Residue, ...], tuple[gemmi.Residue, ...]]:
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid CABS pose: {path}") from exc
    if model_index < 1 or model_index > len(structure):
        raise ValidationError("CABS handoff model index is outside the structure")
    receptor: list[tuple[gemmi.Residue, ...]] = []
    peptide: list[tuple[gemmi.Residue, ...]] = []
    for chain in structure[model_index - 1]:
        residues = _amino_acid_residues(chain)
        if not residues:
            continue
        if _sequence(residues) == peptide_sequence:
            peptide.append(residues)
        else:
            receptor.append(residues)
    if len(receptor) != 1 or len(peptide) != 1:
        raise ValidationError(
            "CABS handoff pose must contain one receptor and one configured peptide"
        )
    return receptor[0], peptide[0]


def _cg_chain(
    *,
    residues: tuple[gemmi.Residue, ...],
    chain_id: str,
    renumber: bool,
    histidine_states: dict[int, str],
) -> gemmi.Chain:
    chain = gemmi.Chain(chain_id)
    for index, source in enumerate(residues, 1):
        residue = gemmi.Residue()
        residue.name = source.name
        residue.entity_type = gemmi.EntityType.Polymer
        residue.seqid = gemmi.SeqId(index, " ") if renumber else source.seqid
        if source.name == "HIS":
            state = histidine_states[index]
            if state not in {"HID", "HIE"}:
                raise ValidationError(
                    "cg2all CA/SC reconstruction supports configured HID or HIE histidine"
                )
            residue.name = "HSD" if state == "HID" else "HSE"
        for atom_name in ("CA", "SC"):
            atom = _named_atom(source, atom_name)
            if atom is not None:
                residue.add_atom(atom.clone())
            elif atom_name == "SC" and source.name != "GLY":
                raise ValidationError(
                    f"CABS CA/SC pose lacks SC atom at {chain_id}:{index}"
                )
        chain.add_residue(residue)
    return chain


def _ssbond_record(*, first: int, second: int) -> str:
    return f"SSBOND  {1:2d} CYS {PEPTIDE_CHAIN} {first:5d}   CYS {PEPTIDE_CHAIN} {second:5d}\n"


def write_cg2all_input(
    *,
    source_path: Path,
    model_index: int,
    destination: Path,
    chemistry: ChemistryDefinition,
    settings: Cg2AllSettings,
) -> CgInputAsset:
    """把任意合格 CABS 链名标准化为 A/P; 冻结配体微状态。"""
    receptor, peptide = _source_chains(
        path=source_path,
        model_index=model_index,
        peptide_sequence=chemistry.sequence,
    )
    histidines = {item.position: item.state for item in chemistry.histidines}
    receptor_chain = _cg_chain(
        residues=receptor,
        chain_id=RECEPTOR_CHAIN,
        renumber=False,
        histidine_states={
            index: settings.receptor_histidine_state
            for index, residue in enumerate(receptor, 1)
            if residue.name == "HIS"
        },
    )
    peptide_chain = _cg_chain(
        residues=peptide,
        chain_id=PEPTIDE_CHAIN,
        renumber=True,
        histidine_states=histidines,
    )
    structure = gemmi.Structure()
    structure.name = "cabs_handoff"
    structure.add_model(gemmi.Model(1))
    structure[0].add_chain(receptor_chain)
    structure[0].add_chain(peptide_chain)
    structure.setup_entities()
    records = "".join(
        _ssbond_record(first=bond.first, second=bond.second)
        for bond in chemistry.disulfide_bonds
    )
    atomic_write_text(destination, records + structure.make_pdb_string())
    receptor_count = len(receptor)
    return CgInputAsset(
        path=destination,
        receptor_residue_count=receptor_count,
        peptide_residue_count=len(peptide),
        fixed_histidine_pose_indices=tuple(
            receptor_count + item.position for item in chemistry.histidines
        ),
    )


def write_peptide_site_coordinate_constraints(
    *,
    source_path: Path,
    destination: Path,
    chemistry: ChemistryDefinition,
    flat_width_A: float,
    standard_deviation_A: float,
) -> int:
    """约束肽 C-alpha 留在重建起点附近。平底区仍允许环骨架调整。"""
    if flat_width_A < 0 or standard_deviation_A <= 0:
        raise ValidationError("site coordinate constraint parameters are invalid")
    receptor, peptide = _source_chains(
        path=source_path,
        model_index=1,
        peptide_sequence=chemistry.sequence,
    )
    if _named_atom(receptor[0], "CA") is None:
        raise ValidationError("site coordinate constraint reference lacks CA")
    receptor_count = len(receptor)
    lines: list[str] = []
    for peptide_index, residue in enumerate(peptide, 1):
        atom = _named_atom(residue, "CA")
        if atom is None:
            raise ValidationError("site coordinate constraint peptide lacks CA")
        lines.append(
            "CoordinateConstraint "
            f"CA {receptor_count + peptide_index} CA 1 "
            f"{atom.pos.x:.6f} {atom.pos.y:.6f} {atom.pos.z:.6f} "
            "FLAT_HARMONIC 0.0 "
            f"{standard_deviation_A:.6f} {flat_width_A:.6f}\n"
        )
    atomic_write_text(destination, "".join(lines))
    return len(lines)


def _ca_positions(
    *, path: Path, peptide_sequence: str
) -> tuple[tuple[gemmi.Position, ...], tuple[gemmi.Position, ...]]:
    receptor, peptide = _source_chains(
        path=path, model_index=1, peptide_sequence=peptide_sequence
    )

    return (
        _residue_ca_positions(receptor),
        _residue_ca_positions(peptide),
    )


def _residue_ca_positions(
    residues: tuple[gemmi.Residue, ...],
) -> tuple[gemmi.Position, ...]:
    result: list[gemmi.Position] = []
    for residue in residues:
        atom = _named_atom(residue, "CA")
        if atom is None:
            raise ValidationError("reconstruction residue lacks CA atom")
        result.append(atom.pos)
    return tuple(result)


def _rmsd(
    first: tuple[gemmi.Position, ...], second: tuple[gemmi.Position, ...]
) -> float:
    if len(first) != len(second) or not first:
        raise ValidationError("reconstruction CA identity differs from its CABS input")
    return math.sqrt(
        sum(a.dist(b) ** 2 for a, b in zip(first, second, strict=True)) / len(first)
    )


def _superpose(
    *, fixed: tuple[gemmi.Position, ...], moving: tuple[gemmi.Position, ...]
) -> gemmi.SupResult:
    if len(fixed) != len(moving) or len(fixed) < 3:
        raise ValidationError("structure identities differ during superposition")
    result = gemmi.superpose_positions(list(fixed), list(moving))
    if not math.isfinite(result.rmsd):
        raise ValidationError("structure superposition is geometrically degenerate")
    return result


def _transformed_positions(
    positions: tuple[gemmi.Position, ...], transform: gemmi.Transform
) -> tuple[gemmi.Position, ...]:
    return tuple(
        gemmi.Position(transformed.x, transformed.y, transformed.z)
        for position in positions
        for transformed in (transform.apply(position),)
    )


def _centroid(positions: tuple[gemmi.Position, ...]) -> tuple[float, float, float]:
    count = len(positions)
    return (
        sum(position.x for position in positions) / count,
        sum(position.y for position in positions) / count,
        sum(position.z for position in positions) / count,
    )


def _receptor_contacts(
    *,
    receptor: tuple[gemmi.Residue, ...],
    peptide: tuple[gemmi.Residue, ...],
    threshold_A: float,
) -> frozenset[tuple[int, str]]:
    peptide_positions = tuple(
        required.pos
        for residue in peptide
        for required in (_named_atom(residue, "CA"),)
        if required is not None
    )
    if len(peptide_positions) != len(peptide):
        raise ValidationError("contact assessment peptide lacks CA atoms")
    contacts: set[tuple[int, str]] = set()
    for residue in receptor:
        receptor_ca = _named_atom(residue, "CA")
        if receptor_ca is None:
            raise ValidationError("contact assessment receptor lacks a CA atom")
        if not any(
            receptor_ca.pos.dist(position) <= threshold_A
            for position in peptide_positions
        ):
            continue
        number = residue.seqid.num
        if number is None:
            raise ValidationError("contact assessment residue has no sequence number")
        contacts.add((number, residue.seqid.icode.strip()))
    return frozenset(contacts)


def measure_reconstruction(
    *,
    input_path: Path,
    output_path: Path,
    chemistry: ChemistryDefinition,
) -> ReconstructionMetrics:
    """测量重建前后的 CA 坐标差异, 不在测量层应用阈值。"""
    input_receptor, input_peptide = _ca_positions(
        path=input_path, peptide_sequence=chemistry.sequence
    )
    output_receptor, output_peptide = _ca_positions(
        path=output_path, peptide_sequence=chemistry.sequence
    )
    return ReconstructionMetrics(
        receptor_ca_rmsd_A=_rmsd(input_receptor, output_receptor),
        peptide_ca_rmsd_A=_rmsd(input_peptide, output_peptide),
    )


def assess_reconstruction(
    *,
    input_path: Path,
    output_path: Path,
    chemistry: ChemistryDefinition,
    settings: Cg2AllSettings,
) -> ReconstructionMetrics:
    """确认重建没有改变受体或环肽的 CA pose。"""
    metrics = measure_reconstruction(
        input_path=input_path,
        output_path=output_path,
        chemistry=chemistry,
    )
    if max(metrics.receptor_ca_rmsd_A, metrics.peptide_ca_rmsd_A) > (
        settings.max_ca_rmsd_A
    ):
        raise ValidationError(
            "cg2all reconstruction exceeds configured CA RMSD: "
            f"receptor={metrics.receptor_ca_rmsd_A:.3f} A, "
            f"peptide={metrics.peptide_ca_rmsd_A:.3f} A"
        )
    return metrics


def assess_topology_reconstruction(
    *,
    input_path: Path,
    output_path: Path,
    chemistry: ChemistryDefinition,
    settings: Cg2AllSettings,
    min_disulfide_sg_A: float,
    max_disulfide_sg_A: float,
    min_interchain_heavy_atom_distance_A: float,
    min_nonlocal_peptide_heavy_atom_distance_A: float,
    contact_ca_threshold_A: float,
    max_peptide_internal_ca_rmsd_A: float,
    max_ligand_centroid_displacement_A: float,
    min_receptor_contact_retention_fraction: float,
) -> TopologyReconstructionAssessment:
    """评价二硫键、CA 保真、主链完整性和严重界面穿插。"""
    input_receptor_positions, input_peptide_positions = _ca_positions(
        path=input_path, peptide_sequence=chemistry.sequence
    )
    output_receptor_positions, output_peptide_positions = _ca_positions(
        path=output_path, peptide_sequence=chemistry.sequence
    )
    receptor_alignment = _superpose(
        fixed=input_receptor_positions, moving=output_receptor_positions
    )
    aligned_output_peptide = _transformed_positions(
        output_peptide_positions, receptor_alignment.transform
    )
    peptide_pose_rmsd = _rmsd(input_peptide_positions, aligned_output_peptide)
    peptide_internal_rmsd = _superpose(
        fixed=input_peptide_positions, moving=output_peptide_positions
    ).rmsd
    centroid_displacement = math.dist(
        _centroid(input_peptide_positions), _centroid(aligned_output_peptide)
    )
    input_receptor, input_peptide = _source_chains(
        path=input_path, model_index=1, peptide_sequence=chemistry.sequence
    )
    receptor, peptide = _source_chains(
        path=output_path,
        model_index=1,
        peptide_sequence=chemistry.sequence,
    )
    input_contacts = _receptor_contacts(
        receptor=input_receptor,
        peptide=input_peptide,
        threshold_A=contact_ca_threshold_A,
    )
    if not input_contacts:
        raise ValidationError("topology calibration input contains no receptor contact")
    output_contacts = _receptor_contacts(
        receptor=receptor,
        peptide=peptide,
        threshold_A=contact_ca_threshold_A,
    )
    contact_retention = len(input_contacts & output_contacts) / len(input_contacts)
    failures: list[str] = []
    if receptor_alignment.rmsd > settings.max_ca_rmsd_A:
        failures.append("receptor_ca_pose_not_preserved")
    if peptide_internal_rmsd > max_peptide_internal_ca_rmsd_A:
        failures.append("peptide_internal_conformation_not_preserved")
    if centroid_displacement > max_ligand_centroid_displacement_A:
        failures.append("ligand_site_centroid_not_preserved")
    if contact_retention < min_receptor_contact_retention_fraction:
        failures.append("receptor_contact_fingerprint_not_preserved")
    for chain_id, residues in (
        (RECEPTOR_CHAIN, receptor),
        (PEPTIDE_CHAIN, peptide),
    ):
        if any(
            _named_atom(residue, atom_name) is None
            for residue in residues
            for atom_name in ("N", "CA", "C", "O")
        ):
            failures.append(f"incomplete_backbone_{chain_id}")
    sg_distances: list[float] = []
    for bond in chemistry.disulfide_bonds:
        try:
            first_sg = _named_atom(peptide[bond.first - 1], "SG")
            second_sg = _named_atom(peptide[bond.second - 1], "SG")
        except IndexError as exc:
            raise ValidationError(
                "topology reconstruction disulfide endpoint is outside peptide"
            ) from exc
        if first_sg is None or second_sg is None:
            failures.append("disulfide_endpoint_lacks_sg")
            continue
        distance = first_sg.pos.dist(second_sg.pos)
        sg_distances.append(distance)
        if not min_disulfide_sg_A <= distance <= max_disulfide_sg_A:
            failures.append("disulfide_sg_geometry_outside_limits")
    receptor_heavy = tuple(
        atom.pos for residue in receptor for atom in residue if atom.element.name != "H"
    )
    peptide_heavy = tuple(
        atom.pos for residue in peptide for atom in residue if atom.element.name != "H"
    )
    if not receptor_heavy or not peptide_heavy:
        raise ValidationError("topology reconstruction contains no heavy atoms")
    minimum_distance = min(
        receptor_atom.dist(peptide_atom)
        for receptor_atom in receptor_heavy
        for peptide_atom in peptide_heavy
    )
    if minimum_distance < min_interchain_heavy_atom_distance_A:
        failures.append("severe_interchain_heavy_atom_overlap")
    disulfide_endpoint_pairs = {
        frozenset((bond.first - 1, bond.second - 1))
        for bond in chemistry.disulfide_bonds
    }
    peptide_heavy_atoms = tuple(
        (residue_index, atom)
        for residue_index, residue in enumerate(peptide)
        for atom in residue
        if atom.element.name != "H"
    )
    nonlocal_distances = tuple(
        first_atom.pos.dist(second_atom.pos)
        for first_index, first_atom in peptide_heavy_atoms
        for second_index, second_atom in peptide_heavy_atoms
        if first_index < second_index
        and second_index - first_index > 1
        and not (
            first_atom.name == "SG"
            and second_atom.name == "SG"
            and frozenset((first_index, second_index)) in disulfide_endpoint_pairs
        )
    )
    if not nonlocal_distances:
        raise ValidationError(
            "topology reconstruction has no nonlocal peptide heavy-atom pairs"
        )
    minimum_nonlocal_peptide_distance = min(nonlocal_distances)
    if minimum_nonlocal_peptide_distance < min_nonlocal_peptide_heavy_atom_distance_A:
        failures.append("severe_nonlocal_peptide_heavy_atom_overlap")
    return TopologyReconstructionAssessment(
        receptor_ca_rmsd_A=receptor_alignment.rmsd,
        peptide_pose_ca_rmsd_A=peptide_pose_rmsd,
        peptide_internal_ca_rmsd_A=peptide_internal_rmsd,
        ligand_centroid_displacement_A=centroid_displacement,
        receptor_contact_retention_fraction=contact_retention,
        disulfide_sg_distances_A=tuple(sg_distances),
        min_interchain_heavy_atom_distance_A=minimum_distance,
        min_nonlocal_peptide_heavy_atom_distance_A=(minimum_nonlocal_peptide_distance),
        failures=tuple(dict.fromkeys(failures)),
    )


def _reference_receptor(path: Path) -> tuple[gemmi.Residue, ...]:
    """读取单模型、单蛋白链的全原子受体参考。"""
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid reference receptor: {path}") from exc
    if len(structure) != 1:
        raise ValidationError("reference receptor must contain exactly one model")
    chains = tuple(
        residues
        for chain in structure[0]
        for residues in (_amino_acid_residues(chain),)
        if residues
    )
    if len(chains) != 1:
        raise ValidationError(
            "reference receptor must contain exactly one protein chain"
        )
    return chains[0]


def _residue_identities(
    residues: tuple[gemmi.Residue, ...],
) -> tuple[tuple[int, str, str], ...]:
    identities: list[tuple[int, str, str]] = []
    for residue in residues:
        number = residue.seqid.num
        if number is None:
            raise ValidationError("reference receptor residue lacks a sequence number")
        name = "HIS" if residue.name in HISTIDINE_NAMES else residue.name
        identities.append((number, residue.seqid.icode.strip(), name))
    return tuple(identities)


def write_reference_receptor_complex(
    *,
    coarse_pose_path: Path,
    reconstructed_path: Path,
    reference_receptor_path: Path,
    destination: Path,
    chemistry: ChemistryDefinition,
) -> None:
    """用对齐后的实验受体替换神经网络受体, 仅保留重建肽链。"""
    coarse_receptor, _ = _source_chains(
        path=coarse_pose_path,
        model_index=1,
        peptide_sequence=chemistry.sequence,
    )
    reference_receptor = _reference_receptor(reference_receptor_path)
    if _residue_identities(coarse_receptor) != _residue_identities(reference_receptor):
        raise ValidationError(
            "reference receptor residues differ from the coarse docking receptor"
        )
    alignment = _superpose(
        fixed=_residue_ca_positions(coarse_receptor),
        moving=_residue_ca_positions(reference_receptor),
    )
    _, peptide = _source_chains(
        path=reconstructed_path,
        model_index=1,
        peptide_sequence=chemistry.sequence,
    )
    output = gemmi.Structure()
    output.name = "all_atom_handoff"
    output.add_model(gemmi.Model(1))
    for chain_id, residues, renumber in (
        (RECEPTOR_CHAIN, reference_receptor, False),
        (PEPTIDE_CHAIN, peptide, True),
    ):
        chain = gemmi.Chain(chain_id)
        for index, source in enumerate(residues, 1):
            residue = source.clone()
            if chain_id == RECEPTOR_CHAIN:
                for atom in residue:
                    transformed = alignment.transform.apply(atom.pos)
                    atom.pos = gemmi.Position(
                        transformed.x, transformed.y, transformed.z
                    )
            if residue.name in HISTIDINE_NAMES:
                residue.name = "HIS"
            if renumber:
                residue.seqid = gemmi.SeqId(index, " ")
            chain.add_residue(residue)
        output[0].add_chain(chain)
    output.setup_entities()
    records = "".join(
        _ssbond_record(first=bond.first, second=bond.second)
        for bond in chemistry.disulfide_bonds
    )
    atomic_write_text(destination, records + output.make_pdb_string())


def _write_chemistry_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
    flexpepdock_attributes: str | None,
) -> None:
    """生成端基、微状态、二硫键及可选同 Pose FlexPepDock 协议。"""
    if chemistry.n_terminus not in {"NH3+", "Acetyl"}:
        raise ValidationError(
            "all-atom handoff supports configured NH3+ or Acetyl N termini"
        )
    if chemistry.c_terminus != "CONH2":
        raise ValidationError(
            "all-atom handoff currently requires a configured CONH2 C terminus"
        )
    if chemistry.other_modifications_status != "none":
        raise ValidationError(
            "all-atom handoff does not support undeclared additional modifications"
        )
    if any(item.state != "HIE" for item in chemistry.histidines):
        raise ValidationError(
            "all-atom handoff currently supports configured HIE histidines only"
        )
    peptide_last_pose_index = receptor_residue_count + len(chemistry.sequence)
    peptide_first_pose_index = receptor_residue_count + 1
    disulfides = ",".join(
        f"{receptor_residue_count + bond.first}:{receptor_residue_count + bond.second}"
        for bond in chemistry.disulfide_bonds
    )
    histidine_movers = "\n".join(
        f'    <MutateResidue name="set_hie_{item.position}" '
        f'target="{receptor_residue_count + item.position}" new_res="HIS" '
        'preserve_atom_coords="true" break_disulfide_bonds="false" />'
        for item in chemistry.histidines
    )
    histidine_protocol = "\n".join(
        f'    <Add mover="set_hie_{item.position}" />' for item in chemistry.histidines
    )
    n_terminus_selector = (
        f'    <Index name="peptide_n_terminus" resnums="{peptide_first_pose_index}" />\n'
        if chemistry.n_terminus == "Acetyl"
        else ""
    )
    n_terminus_movers = (
        '    <ModifyVariantType name="remove_n_terminal_charge" '
        'remove_type="LOWER_TERMINUS_VARIANT" '
        'residue_selector="peptide_n_terminus" '
        'update_polymer_bond_dependent_atoms="true" />\n'
        '    <ModifyVariantType name="add_n_terminal_acetyl" '
        'add_type="N_ACETYLATION" residue_selector="peptide_n_terminus" '
        'update_polymer_bond_dependent_atoms="true" />\n'
        if chemistry.n_terminus == "Acetyl"
        else ""
    )
    n_terminus_protocol = (
        '    <Add mover="remove_n_terminal_charge" />\n'
        '    <Add mover="add_n_terminal_acetyl" />\n'
        if chemistry.n_terminus == "Acetyl"
        else ""
    )
    flexpepdock_mover = (
        f'    <FlexPepDock name="run_flexpepdock" {flexpepdock_attributes} />\n'
        if flexpepdock_attributes is not None
        else ""
    )
    flexpepdock_protocol = (
        '    <Add mover="run_flexpepdock" />\n'
        if flexpepdock_attributes is not None
        else ""
    )
    protocol = f"""<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="handoff_score" weights="{score_function}" />
  </SCOREFXNS>
  <RESIDUE_SELECTORS>
{n_terminus_selector}    <Index name="peptide_c_terminus" resnums="{peptide_last_pose_index}" />
  </RESIDUE_SELECTORS>
  <MOVERS>
{histidine_movers}
{n_terminus_movers}    <ForceDisulfides name="form_disulfides" scorefxn="handoff_score"
      disulfides="{disulfides}" remove_existing="false" repack="true" />
    <ModifyVariantType name="add_c_terminal_amide"
      add_type="CTERM_AMIDATION"
      residue_selector="peptide_c_terminus" />
{flexpepdock_mover}  </MOVERS>
  <PROTOCOLS>
{n_terminus_protocol}{histidine_protocol}
    <Add mover="add_c_terminal_amide" />
    <Add mover="form_disulfides" />
{flexpepdock_protocol}  </PROTOCOLS>
  <OUTPUT scorefxn="handoff_score" />
</ROSETTASCRIPTS>
"""
    atomic_write_text(destination, protocol)


def write_chemistry_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
) -> None:
    """生成只恢复声明化学身份、不执行局部 pose 搜索的协议。"""
    _write_chemistry_protocol(
        destination=destination,
        receptor_residue_count=receptor_residue_count,
        chemistry=chemistry,
        score_function=score_function,
        flexpepdock_attributes=None,
    )


def write_chemistry_prepack_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
) -> None:
    """在同一个 Rosetta Pose 内恢复化学身份并执行 FlexPepDock prepack。"""
    _write_chemistry_protocol(
        destination=destination,
        receptor_residue_count=receptor_residue_count,
        chemistry=chemistry,
        score_function=score_function,
        flexpepdock_attributes='ppk_only="true" recal_foldtree="true"',
    )


def write_chemistry_refine_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
) -> None:
    """在同一个 Rosetta Pose 内恢复化学身份并执行局部精修。"""
    _write_chemistry_protocol(
        destination=destination,
        receptor_residue_count=receptor_residue_count,
        chemistry=chemistry,
        score_function=score_function,
        flexpepdock_attributes='pep_refine="true" recal_foldtree="true"',
    )


def write_chemistry_production_refine_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
    random_translation_A: float,
    random_rotation_degrees: float,
    lowres_preoptimize: bool,
) -> None:
    """生成保持完整化学身份的正式 FlexPepDock 局部精修协议。"""
    if random_translation_A < 0 or random_rotation_degrees < 0:
        raise ValidationError("FlexPepDock initial perturbations must not be negative")
    attributes = (
        'pep_refine="true" recal_foldtree="true" '
        f'lowres_preoptimize="{str(lowres_preoptimize).lower()}" '
        f'rb_trans_size="{random_translation_A}" '
        f'rb_rot_size="{random_rotation_degrees}"'
    )
    _write_chemistry_protocol(
        destination=destination,
        receptor_residue_count=receptor_residue_count,
        chemistry=chemistry,
        score_function=score_function,
        flexpepdock_attributes=attributes,
    )


def write_topology_rebuild_protocol(
    *,
    destination: Path,
    receptor_residue_count: int,
    chemistry: ChemistryDefinition,
    score_function: str,
) -> None:
    """生成只建立二硫键并重排其 6 Å 邻域的 RosettaScripts 协议。"""
    disulfides = ",".join(
        f"{receptor_residue_count + bond.first}:{receptor_residue_count + bond.second}"
        for bond in chemistry.disulfide_bonds
    )
    if not disulfides:
        raise ValidationError("topology rebuild requires at least one disulfide")
    protocol = f"""<ROSETTASCRIPTS>
  <SCOREFXNS>
    <ScoreFunction name="topology_score" weights="{score_function}" />
  </SCOREFXNS>
  <MOVERS>
    <ForceDisulfides name="form_disulfides" scorefxn="topology_score"
      disulfides="{disulfides}" remove_existing="false" repack="true" />
  </MOVERS>
  <PROTOCOLS>
    <Add mover="form_disulfides" />
  </PROTOCOLS>
  <OUTPUT scorefxn="topology_score" />
</ROSETTASCRIPTS>
"""
    atomic_write_text(destination, protocol)


def write_disulfide_indices(
    *, destination: Path, receptor_residue_count: int, chemistry: ChemistryDefinition
) -> None:
    """把配置肽一基编号的二硫键转换为 Rosetta pose 编号。"""
    atomic_write_text(
        destination,
        "".join(
            f"{receptor_residue_count + bond.first} "
            f"{receptor_residue_count + bond.second}\n"
            for bond in chemistry.disulfide_bonds
        ),
    )


def validate_flexpepdock_input(
    *,
    path: Path,
    chemistry: ChemistryDefinition,
    min_disulfide_sg_A: float,
    max_disulfide_sg_A: float,
) -> tuple[int, tuple[int, ...]]:
    """核对链、全原子骨架、端基、His 和二硫键几何。"""
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid FlexPepDock handoff input: {path}") from exc
    if len(structure) != 1:
        raise ValidationError("FlexPepDock handoff input must contain one model")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError("FlexPepDock handoff input must contain only A/P chains")
    receptor = _amino_acid_residues(chains[RECEPTOR_CHAIN])
    peptide = _amino_acid_residues(chains[PEPTIDE_CHAIN])
    if _sequence(peptide) != chemistry.sequence:
        raise ValidationError(
            "FlexPepDock handoff peptide sequence differs from config"
        )
    for chain_id, residues in (
        (RECEPTOR_CHAIN, receptor),
        (PEPTIDE_CHAIN, peptide),
    ):
        for index, residue in enumerate(residues, 1):
            missing = [
                atom_name
                for atom_name in ("N", "CA", "C", "O")
                if _named_atom(residue, atom_name) is None
            ]
            if missing:
                raise ValidationError(
                    f"all-atom backbone is incomplete at {chain_id}:{index}: {missing}"
                )
    first = peptide[0]
    last = peptide[-1]
    if chemistry.n_terminus == "NH3+":
        if any(_named_atom(first, name) is None for name in ("1H", "2H", "3H")):
            raise ValidationError("configured peptide N terminus is not Rosetta NH3+")
    elif chemistry.n_terminus == "Acetyl":
        if any(_named_atom(first, name) is None for name in ("CP", "OCP")) or any(
            _named_atom(first, name) is not None for name in ("1H", "2H", "3H")
        ):
            raise ValidationError(
                "configured peptide N terminus is not Rosetta N-acetylated"
            )
    else:
        raise ValidationError("configured peptide N terminus is unsupported")
    if _named_atom(last, "NT") is None or _named_atom(last, "OXT") is not None:
        raise ValidationError("configured peptide C terminus is not Rosetta CONH2")
    histidine_by_position = {item.position: item.state for item in chemistry.histidines}
    for position, state in histidine_by_position.items():
        residue = peptide[position - 1]
        if state != "HIE" or _named_atom(residue, "HE2") is None:
            raise ValidationError(
                "unsupported or incorrect configured peptide histidine state at "
                f"position {position}"
            )
    for bond in chemistry.disulfide_bonds:
        first_sg = _named_atom(peptide[bond.first - 1], "SG")
        second_sg = _named_atom(peptide[bond.second - 1], "SG")
        if first_sg is None or second_sg is None:
            raise ValidationError("configured peptide disulfide endpoint lacks SG atom")
        distance = first_sg.pos.dist(second_sg.pos)
        if not min_disulfide_sg_A <= distance <= max_disulfide_sg_A:
            raise ValidationError(
                "configured peptide disulfide SG distance is outside config: "
                f"{distance:.3f} A"
            )
    receptor_count = len(receptor)
    return receptor_count, tuple(
        receptor_count + position for position in sorted(histidine_by_position)
    )
