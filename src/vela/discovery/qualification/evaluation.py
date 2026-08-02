"""仅在候选冻结后使用实验复合物评价 native recovery。"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.evidence import (
    align_receptor,
    ca_contact_residues,
    centroid,
    read_structure,
    required_atom,
    split_model,
)
from vela.discovery.sampling.trajectory import (
    CabsSequenceChain,
    iter_cabs_trajectory,
    read_cabs_sequence,
)
from vela.preparation.chemistry import ChemistryDefinition


@dataclass(frozen=True, slots=True)
class NativePoseMetrics:
    """冻结 pose 相对实验复合物的事后回收指标。"""

    ligand_ca_rmsd_A: float
    ligand_centroid_distance_A: float
    native_receptor_contact_fraction: float


@dataclass(frozen=True, slots=True)
class NativeTrajectoryAudit:
    """CABS 能量过滤前后以 L-RMSD 衡量的实验位点保留情况。"""

    trajectory_model_count: int
    trajectory_recovered_model_count: int
    trajectory_topology_feasible_recovered_model_count: int
    filtered_model_count: int
    filtered_recovered_model_count: int
    filtered_topology_feasible_recovered_model_count: int
    trajectory_qualified_recovered_model_count: int
    filtered_qualified_recovered_model_count: int
    trajectory_best_ligand_ca_rmsd_A: float | None
    trajectory_best_native_receptor_contact_fraction: float | None
    filtered_best_ligand_ca_rmsd_A: float | None
    filtered_best_native_receptor_contact_fraction: float | None


def _raw_receptor_alignment_indices(
    *, sequence: CabsSequenceChain, native_receptor: tuple[gemmi.Residue, ...]
) -> tuple[tuple[int, ...], tuple[gemmi.Position, ...]]:
    native_by_id = {
        (residue.seqid.num, residue.seqid.icode.strip()): residue
        for residue in native_receptor
    }
    indices: list[int] = []
    fixed: list[gemmi.Position] = []
    for index, (number, insertion, name) in enumerate(
        zip(
            sequence.residue_numbers,
            sequence.insertion_codes,
            sequence.residue_names,
            strict=True,
        )
    ):
        residue = native_by_id.get((number, insertion))
        if residue is None or residue.name != name:
            continue
        indices.append(index)
        fixed.append(required_atom(residue, "CA").pos)
    if len(indices) < 50:
        raise DiscoveryError(
            "fewer than 50 identity-matched receptor CA atoms are available for "
            "TRAF native evaluation"
        )
    return tuple(indices), tuple(fixed)


def _positions_rmsd(
    first: tuple[gemmi.Vec3, ...], second: tuple[gemmi.Vec3, ...]
) -> float:
    if len(first) != len(second) or not first:
        raise DiscoveryError("native and sampled peptide CA identities differ")
    return math.sqrt(
        sum(
            math.dist(
                (left.x, left.y, left.z),
                (right.x, right.y, right.z),
            )
            ** 2
            for left, right in zip(first, second, strict=True)
        )
        / len(first)
    )


def _raw_topology_distance(
    peptide: tuple[tuple[float, float, float], ...],
    *,
    chemistry: ChemistryDefinition,
) -> float:
    distances: list[float] = []
    for bond in chemistry.disulfide_bonds:
        try:
            distances.append(
                math.dist(peptide[bond.first - 1], peptide[bond.second - 1])
            )
        except IndexError as exc:
            raise DiscoveryError(
                "disulfide endpoint is outside the CABS TRAF peptide"
            ) from exc
    if not distances:
        raise DiscoveryError("native TRAF audit requires a disulfide definition")
    return max(distances)


def audit_native_recovery_before_filtering(
    *,
    archive_path: Path,
    filtered_path: Path,
    native_pair_path: Path,
    chemistry: ChemistryDefinition,
    max_ligand_ca_rmsd_A: float,
    max_disulfide_ca_distance_A: float,
    contact_ca_threshold_A: float,
    min_native_receptor_contact_fraction: float,
) -> NativeTrajectoryAudit:
    """事后比较完整 TRAF 与实际 Top-1000 的 native L-RMSD 召回。"""
    native = read_structure(native_pair_path)
    if len(native) != 1:
        raise DiscoveryError("native comparison structure must contain one model")
    native_receptor, native_peptide = split_model(
        native[0], peptide_sequence=chemistry.sequence
    )
    native_positions = tuple(
        required_atom(residue, "CA").pos for residue in native_peptide
    )
    native_contacts = ca_contact_residues(
        receptor=native_receptor,
        peptide=native_peptide,
        threshold_A=contact_ca_threshold_A,
    )
    if not native_contacts:
        raise DiscoveryError("native control contains no CA contacts")
    chains = read_cabs_sequence(archive_path)
    if chains[-1].sequence != chemistry.sequence:
        raise DiscoveryError("CABS TRAF peptide differs from native audit chemistry")
    receptor_indices, fixed_positions = _raw_receptor_alignment_indices(
        sequence=chains[0], native_receptor=native_receptor
    )
    trajectory_count = 0
    trajectory_recovered = 0
    trajectory_topology_recovered = 0
    trajectory_qualified_recovered = 0
    trajectory_best_rmsd: float | None = None
    trajectory_best_contact: float | None = None
    for frame in iter_cabs_trajectory(archive_path=archive_path, chains=chains):
        trajectory_count += 1
        receptor = frame.chain_ca[0]
        moving_receptor = tuple(
            gemmi.Position(*receptor[index]) for index in receptor_indices
        )
        alignment = gemmi.superpose_positions(fixed_positions, moving_receptor)
        moving_peptide = tuple(
            alignment.transform.apply(gemmi.Position(*position))
            for position in frame.chain_ca[-1]
        )
        ligand_rmsd = _positions_rmsd(moving_peptide, native_positions)
        topology_feasible = (
            _raw_topology_distance(frame.chain_ca[-1], chemistry=chemistry)
            <= max_disulfide_ca_distance_A
        )
        pose_contacts = {
            f"{number}{insertion}"
            for number, insertion, position in zip(
                chains[0].residue_numbers,
                chains[0].insertion_codes,
                frame.chain_ca[0],
                strict=True,
            )
            if any(
                math.dist(position, ligand) <= contact_ca_threshold_A
                for ligand in frame.chain_ca[-1]
            )
        }
        contact_fraction = len(native_contacts & pose_contacts) / len(native_contacts)
        if topology_feasible:
            trajectory_best_rmsd = (
                ligand_rmsd
                if trajectory_best_rmsd is None
                else min(trajectory_best_rmsd, ligand_rmsd)
            )
            trajectory_best_contact = (
                contact_fraction
                if trajectory_best_contact is None
                else max(trajectory_best_contact, contact_fraction)
            )
        recovered = ligand_rmsd <= max_ligand_ca_rmsd_A
        if recovered:
            trajectory_recovered += 1
            if topology_feasible:
                trajectory_topology_recovered += 1
                if contact_fraction >= min_native_receptor_contact_fraction:
                    trajectory_qualified_recovered += 1

    filtered = read_structure(filtered_path)
    filtered_recovered = 0
    filtered_topology_recovered = 0
    filtered_qualified_recovered = 0
    filtered_best_rmsd: float | None = None
    filtered_best_contact: float | None = None
    for model in filtered:
        receptor, peptide = split_model(model, peptide_sequence=chemistry.sequence)
        alignment = align_receptor(
            receptor=receptor,
            reference_chains=(native_receptor,),
        )
        moving_peptide = tuple(
            alignment.transform.apply(required_atom(residue, "CA").pos)
            for residue in peptide
        )
        ligand_rmsd = _positions_rmsd(moving_peptide, native_positions)
        pose_contacts = ca_contact_residues(
            receptor=receptor,
            peptide=peptide,
            threshold_A=contact_ca_threshold_A,
        )
        contact_fraction = len(native_contacts & pose_contacts) / len(native_contacts)
        distances = tuple(
            required_atom(peptide[bond.first - 1], "CA").pos.dist(
                required_atom(peptide[bond.second - 1], "CA").pos
            )
            for bond in chemistry.disulfide_bonds
        )
        topology_feasible = (
            bool(distances) and max(distances) <= max_disulfide_ca_distance_A
        )
        if topology_feasible:
            filtered_best_rmsd = (
                ligand_rmsd
                if filtered_best_rmsd is None
                else min(filtered_best_rmsd, ligand_rmsd)
            )
            filtered_best_contact = (
                contact_fraction
                if filtered_best_contact is None
                else max(filtered_best_contact, contact_fraction)
            )
        recovered = ligand_rmsd <= max_ligand_ca_rmsd_A
        if recovered:
            filtered_recovered += 1
            if topology_feasible:
                filtered_topology_recovered += 1
                if contact_fraction >= min_native_receptor_contact_fraction:
                    filtered_qualified_recovered += 1
    return NativeTrajectoryAudit(
        trajectory_model_count=trajectory_count,
        trajectory_recovered_model_count=trajectory_recovered,
        trajectory_topology_feasible_recovered_model_count=(
            trajectory_topology_recovered
        ),
        filtered_model_count=len(filtered),
        filtered_recovered_model_count=filtered_recovered,
        filtered_topology_feasible_recovered_model_count=(filtered_topology_recovered),
        trajectory_qualified_recovered_model_count=(trajectory_qualified_recovered),
        filtered_qualified_recovered_model_count=filtered_qualified_recovered,
        trajectory_best_ligand_ca_rmsd_A=trajectory_best_rmsd,
        trajectory_best_native_receptor_contact_fraction=trajectory_best_contact,
        filtered_best_ligand_ca_rmsd_A=filtered_best_rmsd,
        filtered_best_native_receptor_contact_fraction=filtered_best_contact,
    )


def compare_poses_to_native(
    *,
    poses: tuple[PoseEvidence, ...],
    native_pair_path: Path,
    peptide_sequence: str,
    contact_ca_threshold_A: float,
) -> dict[str, NativePoseMetrics]:
    """事后计算回收指标; 该函数不参与采样或候选选择。"""
    native = read_structure(native_pair_path)
    if len(native) != 1:
        raise DiscoveryError("native comparison structure must contain one model")
    native_receptor, native_peptide = split_model(
        native[0], peptide_sequence=peptide_sequence
    )
    native_contacts = ca_contact_residues(
        receptor=native_receptor,
        peptide=native_peptide,
        threshold_A=contact_ca_threshold_A,
    )
    if not native_contacts:
        raise DiscoveryError("native control contains no receptor contacts")
    native_positions = tuple(
        required_atom(residue, "CA").pos for residue in native_peptide
    )
    native_centroid = gemmi.Position(
        *centroid(
            tuple((position.x, position.y, position.z) for position in native_positions)
        )
    )
    structures: dict[Path, gemmi.Structure] = {}
    metrics: dict[str, NativePoseMetrics] = {}
    for pose in poses:
        if pose.pose_id in metrics:
            raise DiscoveryError(f"duplicate native comparison pose: {pose.pose_id}")
        structure = structures.get(pose.model_path)
        if structure is None:
            structure = read_structure(pose.model_path)
            structures[pose.model_path] = structure
        if pose.model_index > len(structure):
            raise DiscoveryError(
                f"pose model index is outside its structure: {pose.pose_id}"
            )
        receptor, peptide = split_model(
            structure[pose.model_index - 1], peptide_sequence=peptide_sequence
        )
        if len(peptide) != len(native_peptide):
            raise DiscoveryError("native and sampled ligands have different lengths")
        alignment = align_receptor(
            receptor=receptor,
            reference_chains=(native_receptor,),
        )
        moving_positions = tuple(
            alignment.transform.apply(required_atom(residue, "CA").pos)
            for residue in peptide
        )
        squared = sum(
            math.dist(
                (moving.x, moving.y, moving.z),
                (reference.x, reference.y, reference.z),
            )
            ** 2
            for moving, reference in zip(
                moving_positions, native_positions, strict=True
            )
        )
        moving_centroid = gemmi.Position(
            *centroid(
                tuple(
                    (position.x, position.y, position.z)
                    for position in moving_positions
                )
            )
        )
        pose_contacts = ca_contact_residues(
            receptor=receptor,
            peptide=peptide,
            threshold_A=contact_ca_threshold_A,
        )
        metrics[pose.pose_id] = NativePoseMetrics(
            ligand_ca_rmsd_A=math.sqrt(squared / len(native_positions)),
            ligand_centroid_distance_A=moving_centroid.dist(native_centroid),
            native_receptor_contact_fraction=len(native_contacts & pose_contacts)
            / len(native_contacts),
        )
    return metrics
