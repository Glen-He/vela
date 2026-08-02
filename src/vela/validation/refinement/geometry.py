"""阶段三局部精修 decoy 的几何判定和构象聚类。"""

from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.discovery.analysis.cluster_engine import complete_linkage
from vela.preparation.chemistry import ChemistryDefinition
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import ValidationAnalysisSettings, ValidationError
from vela.validation.refinement.planning import RefinementTask
from vela.validation.refinement.reconstruction import validate_flexpepdock_input


@dataclass(frozen=True, slots=True)
class ResolvedAnalysisSettings:
    """已经关闭未决值的候选精修分析门槛。"""

    min_interface_contact_pairs: int
    min_interface_receptor_residues: int
    max_receptor_ca_rmsd_A: float
    min_start_contact_overlap: float
    max_start_site_displacement_A: float
    max_cluster_backbone_rmsd_A: float
    min_heavy_atom_distance_A: float
    min_refinement_seed_support: int
    min_refinement_start_support: int


@dataclass(frozen=True, slots=True)
class ComplexGeometry:
    """计算局部精修 QC 和跨起点对齐所需的结构几何。"""

    receptor_ca: tuple[gemmi.Position, ...]
    peptide_backbone: tuple[gemmi.Position, ...]
    peptide_ca: tuple[gemmi.Position, ...]
    receptor_contacts: frozenset[str]
    interface_contact_pairs: int
    minimum_interface_distance_A: float


@dataclass(frozen=True, slots=True)
class GeometryAssessment:
    """不依赖任务类型的复合物几何质控结果。"""

    interface_contact_pairs: int
    interface_receptor_residues: int
    minimum_interface_distance_A: float
    receptor_ca_rmsd_A: float
    start_contact_overlap: float
    start_site_displacement_A: float
    passed: bool
    cluster_backbone: tuple[gemmi.Position, ...]


@dataclass(frozen=True, slots=True)
class RefinedDecoy:
    """一个带来源身份、配置分数、QC 和聚类坐标的 decoy。"""

    decoy_id: str
    task_id: str
    start_id: str
    candidate_id: str
    receptor_id: str
    source_seed: int | None
    refinement_seed: int
    path: Path
    sha256: str
    ranking_score: float
    interface_contact_pairs: int
    interface_receptor_residues: int
    minimum_interface_distance_A: float
    receptor_ca_rmsd_A: float
    start_contact_overlap: float
    start_site_displacement_A: float
    qc_status: str
    cluster_backbone: tuple[gemmi.Position, ...]


@dataclass(frozen=True, slots=True)
class RefinedCluster:
    """同一 candidate 和受体构象内的局部精修构象簇。"""

    cluster_id: str
    candidate_id: str
    receptor_id: str
    decoy_ids: tuple[str, ...]
    refinement_seeds: tuple[int, ...]
    start_ids: tuple[str, ...]
    representative_decoy_id: str
    supported: bool


def resolve_analysis_settings(
    settings: ValidationAnalysisSettings,
) -> ResolvedAnalysisSettings:
    """把完整性已验证的可选配置收窄为运行时分析合同。"""
    integers = (
        settings.min_interface_contact_pairs,
        settings.min_interface_receptor_residues,
        settings.min_refinement_seed_support,
        settings.min_refinement_start_support,
    )
    numbers = (
        settings.max_receptor_ca_rmsd_A,
        settings.min_start_contact_overlap,
        settings.max_start_site_displacement_A,
        settings.max_cluster_backbone_rmsd_A,
        settings.min_heavy_atom_distance_A,
    )
    if any(value is None for value in (*integers, *numbers)):
        raise ValidationError("Stage 3 analysis settings are unresolved")
    pairs, residues, seeds, starts = integers
    receptor_rmsd, overlap, displacement, cluster_rmsd, distance = numbers
    if (
        not isinstance(pairs, int)
        or not isinstance(residues, int)
        or not isinstance(seeds, int)
        or not isinstance(starts, int)
        or not isinstance(receptor_rmsd, float)
        or not isinstance(overlap, float)
        or not isinstance(displacement, float)
        or not isinstance(cluster_rmsd, float)
        or not isinstance(distance, float)
    ):
        raise ValidationError("Stage 3 analysis settings have invalid types")
    return ResolvedAnalysisSettings(
        pairs,
        residues,
        receptor_rmsd,
        overlap,
        displacement,
        cluster_rmsd,
        distance,
        seeds,
        starts,
    )


def _named_atom(residue: gemmi.Residue, name: str) -> gemmi.Atom | None:
    for atom in residue:
        if atom.name == name:
            return atom
    return None


def _amino_acids(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(residue for residue in chain if _named_atom(residue, "CA") is not None)


def _required_atom(residue: gemmi.Residue, name: str, *, path: Path) -> gemmi.Atom:
    atom = _named_atom(residue, name)
    if atom is None:
        raise ValidationError(f"required atom {name} is missing: {path}")
    return atom


def read_complex_geometry(*, path: Path, interface_contact_A: float) -> ComplexGeometry:
    """读取 A/P 复合物并计算不依赖评分函数的界面几何。"""
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid refinement structure: {path}") from exc
    if len(structure) != 1:
        raise ValidationError(f"refinement structure must contain one model: {path}")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError(f"refinement structure must contain A/P chains: {path}")
    receptor = _amino_acids(chains[RECEPTOR_CHAIN])
    peptide = _amino_acids(chains[PEPTIDE_CHAIN])
    receptor_atoms = tuple(
        (residue, atom)
        for residue in receptor
        for atom in residue
        if atom.element.name != "H"
    )
    peptide_atoms = tuple(
        atom for residue in peptide for atom in residue if atom.element.name != "H"
    )
    if not receptor_atoms or not peptide_atoms:
        raise ValidationError(f"refinement structure lacks interface atoms: {path}")
    contacts: set[str] = set()
    pairs = 0
    minimum = float("inf")
    for residue, receptor_atom in receptor_atoms:
        for peptide_atom in peptide_atoms:
            distance = receptor_atom.pos.dist(peptide_atom.pos)
            minimum = min(minimum, distance)
            if distance <= interface_contact_A:
                pairs += 1
                contacts.add(f"{residue.seqid.num}{residue.seqid.icode.strip()}")
    return ComplexGeometry(
        receptor_ca=tuple(
            _required_atom(residue, "CA", path=path).pos for residue in receptor
        ),
        peptide_backbone=tuple(
            _required_atom(residue, atom_name, path=path).pos
            for residue in peptide
            for atom_name in ("N", "CA", "C")
        ),
        peptide_ca=tuple(
            _required_atom(residue, "CA", path=path).pos for residue in peptide
        ),
        receptor_contacts=frozenset(contacts),
        interface_contact_pairs=pairs,
        minimum_interface_distance_A=minimum,
    )


def _transformed(
    positions: tuple[gemmi.Position, ...], transform: gemmi.Transform
) -> tuple[gemmi.Position, ...]:
    return tuple(gemmi.Position(transform.apply(position)) for position in positions)


def _centroid(positions: tuple[gemmi.Position, ...]) -> gemmi.Position:
    if not positions:
        raise ValidationError("cannot calculate an empty coordinate centroid")
    count = len(positions)
    return gemmi.Position(
        sum(position.x for position in positions) / count,
        sum(position.y for position in positions) / count,
        sum(position.z for position in positions) / count,
    )


def _aligned_positions(
    *,
    fixed_receptor: tuple[gemmi.Position, ...],
    movable: ComplexGeometry,
    positions: tuple[gemmi.Position, ...],
) -> tuple[float, tuple[gemmi.Position, ...]]:
    if len(fixed_receptor) != len(movable.receptor_ca) or len(fixed_receptor) < 3:
        raise ValidationError("refinement receptor CA correspondence is invalid")
    result = gemmi.superpose_positions(fixed_receptor, movable.receptor_ca)
    return result.rmsd, _transformed(positions, result.transform)


def assess_refined_decoy(
    *,
    task: RefinementTask,
    decoy_id: str,
    path: Path,
    path_sha256: str,
    ranking_score: float,
    start: ComplexGeometry,
    cluster_reference: ComplexGeometry,
    config: AppConfig,
    settings: ResolvedAnalysisSettings,
) -> RefinedDecoy:
    """按配置门槛评价单个 decoy 并生成共同受体坐标系。"""
    assessment = assess_complex_geometry(
        path=path,
        chemistry=config.chemistry,
        start=start,
        cluster_reference=cluster_reference,
        config=config,
        settings=settings,
    )
    return RefinedDecoy(
        decoy_id,
        task.task_id,
        task.start.start_id,
        task.start.candidate_id,
        task.start.receptor_id,
        task.start.source_seed,
        task.seed,
        path,
        path_sha256,
        ranking_score,
        assessment.interface_contact_pairs,
        assessment.interface_receptor_residues,
        assessment.minimum_interface_distance_A,
        assessment.receptor_ca_rmsd_A,
        assessment.start_contact_overlap,
        assessment.start_site_displacement_A,
        "passed" if assessment.passed else "failed",
        assessment.cluster_backbone,
    )


def assess_complex_geometry(
    *,
    path: Path,
    chemistry: ChemistryDefinition,
    start: ComplexGeometry,
    cluster_reference: ComplexGeometry,
    config: AppConfig,
    settings: ResolvedAnalysisSettings,
) -> GeometryAssessment:
    """用共同门槛评价任意配体序列复合物; 供精修和序列设计复用。"""
    validate_flexpepdock_input(
        path=path,
        chemistry=chemistry,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    decoy = read_complex_geometry(
        path=path, interface_contact_A=config.validation.interface_contact_A
    )
    receptor_rmsd, aligned_peptide_ca = _aligned_positions(
        fixed_receptor=start.receptor_ca,
        movable=decoy,
        positions=decoy.peptide_ca,
    )
    _, cluster_backbone = _aligned_positions(
        fixed_receptor=cluster_reference.receptor_ca,
        movable=decoy,
        positions=decoy.peptide_backbone,
    )
    overlap = (
        len(start.receptor_contacts & decoy.receptor_contacts)
        / len(start.receptor_contacts)
        if start.receptor_contacts
        else 0.0
    )
    displacement = _centroid(start.peptide_ca).dist(_centroid(aligned_peptide_ca))
    passed = (
        decoy.interface_contact_pairs >= settings.min_interface_contact_pairs
        and len(decoy.receptor_contacts) >= settings.min_interface_receptor_residues
        and decoy.minimum_interface_distance_A >= settings.min_heavy_atom_distance_A
        and receptor_rmsd <= settings.max_receptor_ca_rmsd_A
        and overlap >= settings.min_start_contact_overlap
        and displacement <= settings.max_start_site_displacement_A
    )
    return GeometryAssessment(
        decoy.interface_contact_pairs,
        len(decoy.receptor_contacts),
        decoy.minimum_interface_distance_A,
        receptor_rmsd,
        overlap,
        displacement,
        passed,
        cluster_backbone,
    )


def _backbone_rmsd(first: RefinedDecoy, second: RefinedDecoy) -> float:
    if len(first.cluster_backbone) != len(second.cluster_backbone):
        raise ValidationError("refined peptide backbone lengths differ")
    return math.sqrt(
        sum(
            left.dist(right) ** 2
            for left, right in zip(
                first.cluster_backbone, second.cluster_backbone, strict=True
            )
        )
        / len(first.cluster_backbone)
    )


def cluster_refined_decoys(
    *, decoys: tuple[RefinedDecoy, ...], settings: ResolvedAnalysisSettings
) -> tuple[RefinedCluster, ...]:
    """按 candidate 和受体隔离聚类并以独立 seed/起点定义支持。"""
    grouped: dict[tuple[str, str], list[RefinedDecoy]] = defaultdict(list)
    for decoy in decoys:
        if decoy.qc_status == "passed":
            grouped[(decoy.candidate_id, decoy.receptor_id)].append(decoy)
    results: list[RefinedCluster] = []
    for (candidate_id, receptor_id), members in sorted(grouped.items()):
        clusters = complete_linkage(
            members,
            distance=lambda first, second: (
                _backbone_rmsd(first, second) / settings.max_cluster_backbone_rmsd_A
            ),
            identity=lambda item: item.decoy_id,
        )
        for index, cluster in enumerate(clusters, 1):
            representative = min(
                cluster,
                key=lambda item: (
                    sum(_backbone_rmsd(item, other) for other in cluster),
                    item.ranking_score,
                    item.decoy_id,
                ),
            )
            seeds = tuple(sorted({item.refinement_seed for item in cluster}))
            starts = tuple(sorted({item.start_id for item in cluster}))
            results.append(
                RefinedCluster(
                    f"{candidate_id}__{receptor_id}__R{index:03d}",
                    candidate_id,
                    receptor_id,
                    tuple(item.decoy_id for item in cluster),
                    seeds,
                    starts,
                    representative.decoy_id,
                    len(seeds) >= settings.min_refinement_seed_support
                    and len(starts) >= settings.min_refinement_start_support,
                )
            )
    return tuple(results)
