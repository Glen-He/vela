"""将 CABS-dock 粗粒化采样输出转换为规范候选与资格证据。"""

from __future__ import annotations

import math
import statistics
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import JsonValue, atomic_write_json, sha256_file
from vela.discovery.analysis.cluster_engine import (
    bounded_leader_clusters,
    normalized_site_distance,
)
from vela.discovery.analysis.evidence import Point3D, PoseEvidence
from vela.discovery.models import CabsDockSettings, DiscoveryError, DiscoveryTask
from vela.discovery.sampling.cabsdock import (
    cabsdock_archive_path,
    cabsdock_medoid_paths,
)
from vela.discovery.sampling.materialization import (
    CabsFrameIdentity,
    materialize_cabs_frames,
)
from vela.discovery.sampling.trajectory import (
    CabsSequenceChain,
    CabsTrajectoryAudit,
    CabsTrajectoryFrame,
    audit_cabs_trajectory,
    iter_cabs_trajectory,
    read_cabs_sequence,
    trajectory_disulfide_ca_distance,
)
from vela.preparation.chemistry import ChemistryDefinition

type Residues = tuple[gemmi.Residue, ...]
type PeptideCoordinates = tuple[Point3D, ...]


@dataclass(frozen=True, slots=True)
class CandidateSelectionSettings:
    """项目侧位点优先、位点内姿态去重规则。"""

    contact_jaccard_distance: float
    position_distance_A: float
    peptide_ca_rmsd_A: float
    max_sites: int
    max_pose_clusters_per_site: int

    def __post_init__(self) -> None:
        if not 0.0 < self.contact_jaccard_distance <= 1.0:
            raise DiscoveryError("selection contact Jaccard distance must be in (0, 1]")
        if self.position_distance_A <= 0 or self.peptide_ca_rmsd_A <= 0:
            raise DiscoveryError("selection distance thresholds must be positive")
        if self.max_sites < 1 or self.max_pose_clusters_per_site < 1:
            raise DiscoveryError("selection cluster budgets must be positive")


@dataclass(frozen=True, slots=True)
class CabsDockEvidence:
    """单任务的完整候选池、项目侧选择和上游 Top-10 基线。"""

    poses: tuple[PoseEvidence, ...]
    baseline_poses: tuple[PoseEvidence, ...]
    sampling_model_count: int
    filtered_model_count: int
    filtered_topology_feasible_model_count: int
    filtered_topology_feasible_fraction: float
    topology_feasible_model_count: int
    topology_feasible_fraction: float
    contacting_topology_feasible_model_count: int
    disulfide_ca_distance_median_A: float
    disulfide_ca_distance_p90_A: float
    disulfide_ca_distance_max_A: float
    candidate_frame_selection_path: Path
    candidate_models_path: Path | None
    trajectory_audit: CabsTrajectoryAudit
    selection_status: str
    site_cluster_count: int
    pose_cluster_count: int
    selected_pose_cluster_count: int
    receptor_alignment_max_rmsd_A: float | None

    @property
    def selection_failure_reasons(self) -> tuple[str, ...]:
        """返回候选选择的规范失败原因。"""
        return () if self.selection_status == "completed" else (self.selection_status,)


@dataclass(frozen=True, slots=True)
class _FrameEvidence:
    """候选选择内部使用的已对齐模型特征。"""

    task: DiscoveryTask
    pose_id: str
    identity: CabsFrameIdentity | None
    contact_residues: frozenset[str]
    local_position: Point3D
    ranking_score: float
    score_name: str
    qc_status: str
    peptide_ca: PeptideCoordinates


def candidate_selection_settings(
    settings: CabsDockSettings,
) -> CandidateSelectionSettings:
    """从唯一 CABS 配置源建立候选选择参数。"""
    return CandidateSelectionSettings(
        contact_jaccard_distance=settings.selection_contact_jaccard_distance,
        position_distance_A=settings.selection_position_distance_A,
        peptide_ca_rmsd_A=settings.pose_clustering_rmsd_A,
        max_sites=settings.max_sites_per_task,
        max_pose_clusters_per_site=settings.max_pose_clusters_per_site,
    )


def candidate_selection_contract(
    settings: CabsDockSettings,
) -> dict[str, JsonValue]:
    """返回计划、执行和报告共用的完整 native-free 选择合同。"""
    return {
        "mode": "site_first_then_pose",
        "native_information_used": False,
        "input_pool": "complete_traf_with_topology_feasible_geometry_and_receptor_contact",
        "minimum_input_models": settings.min_models_for_selection,
        "topology_feasibility": {
            "representation": "peptide_disulfide_endpoint_ca",
            "maximum_distance_A": (
                settings.max_reconstructable_disulfide_ca_distance_A
            ),
            "chemical_bond_claimed": False,
        },
        "receptor_contact": {
            "representation": "CA",
            "distance_A": settings.trajectory_contact_ca_threshold_A,
            "identity": "author_residue_number_and_insertion_code",
        },
        "coordinate_frame": {
            "alignment_atoms": "identity_matched_receptor_CA",
            "minimum_alignment_atom_count": 50,
            "site_position": "aligned_peptide_CA_centroid",
        },
        "site_clustering": {
            "algorithm": "deterministic_bounded_leader_complete_diameter",
            "distance": "max_normalized_contact_jaccard_and_centroid_distance",
            "contact_jaccard_distance": settings.selection_contact_jaccard_distance,
            "position_distance_A": settings.selection_position_distance_A,
            "maximum_site_count": settings.max_sites_per_task,
            "ranking": "population_desc_then_best_interaction_energy",
        },
        "pose_clustering": {
            "algorithm": "deterministic_bounded_leader_complete_diameter",
            "atom_set": "peptide_CA_after_receptor_alignment",
            "rmsd_A": settings.pose_clustering_rmsd_A,
            "cluster_coverage": "complete",
            "maximum_clusters_per_site": settings.max_pose_clusters_per_site,
            "ranking": "population_desc_then_best_interaction_energy",
        },
        "representatives": {
            "per_pose_cluster": [
                "geometric_medoid",
                "minimum_cabsdock_interaction_energy",
            ],
            "maximum_per_pose_cluster": 2,
            "total_candidate_budget_per_task": (
                settings.max_sites_per_task * settings.max_pose_clusters_per_site * 2
            ),
            "tie_break": "lexicographic_pose_id",
        },
    }


def read_structure(path: Path) -> gemmi.Structure:
    if not path.is_file():
        raise DiscoveryError(f"CABS-dock output does not exist: {path}")
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise DiscoveryError(f"invalid CABS-dock structure: {path}") from exc
    if len(structure) == 0:
        raise DiscoveryError(f"CABS-dock structure contains no models: {path}")
    return structure


def _named_atom(residue: gemmi.Residue, name: str) -> gemmi.Atom | None:
    """绕开 Gemmi 缺失原子代理对象在运行时不可安全解引用的问题。"""
    for atom in residue:
        if atom.name == name:
            return atom
    return None


def required_atom(residue: gemmi.Residue, name: str) -> gemmi.Atom:
    atom = _named_atom(residue, name)
    if atom is None:
        number = residue.seqid.num
        raise DiscoveryError(
            f"residue {residue.name} {number} lacks required {name} atom"
        )
    return atom


def _ca_residues(chain: gemmi.Chain) -> Residues:
    return tuple(
        residue
        for residue in chain
        if _named_atom(residue, "CA") is not None
        and gemmi.find_tabulated_residue(residue.name).is_amino_acid()
    )


def _sequence(residues: Residues) -> str:
    return "".join(
        gemmi.find_tabulated_residue(residue.name).one_letter_code
        for residue in residues
    )


def split_model(
    model: gemmi.Model, *, peptide_sequence: str
) -> tuple[Residues, Residues]:
    peptide_matches: list[Residues] = []
    receptor_chains: list[Residues] = []
    for chain in model:
        residues = _ca_residues(chain)
        if not residues:
            continue
        if _sequence(residues) == peptide_sequence:
            peptide_matches.append(residues)
        else:
            receptor_chains.append(residues)
    if len(peptide_matches) != 1:
        raise DiscoveryError(
            "CABS-dock model must contain exactly one chain matching the ligand"
        )
    if len(receptor_chains) != 1:
        raise DiscoveryError(
            "CABS-dock model must contain exactly one receptor protein chain"
        )
    return receptor_chains[0], peptide_matches[0]


def _residue_key(residue: gemmi.Residue) -> tuple[int, str]:
    number = residue.seqid.num
    if number is None:
        raise DiscoveryError("protein residue has no sequence number")
    return number, residue.seqid.icode.strip()


def read_reference_chains(path: Path) -> tuple[Residues, ...]:
    """读取只含受体的坐标参考, 不接触任何实验配体坐标。"""
    reference = read_structure(path)
    chains = tuple(
        residues for chain in reference[0] if (residues := _ca_residues(chain))
    )
    if not chains:
        raise DiscoveryError(f"reference contains no protein chain: {path}")
    return chains


def align_receptor(
    *, receptor: Residues, reference_chains: tuple[Residues, ...]
) -> gemmi.SupResult:
    moving = {_residue_key(residue): residue for residue in receptor}
    candidates = tuple(
        (sum(_residue_key(residue) in moving for residue in residues), residues)
        for residues in reference_chains
    )
    common_count, reference_residues = max(candidates, key=lambda item: item[0])
    if common_count == 0:
        raise DiscoveryError("no receptor chain can be aligned to the reference")
    fixed_positions: list[gemmi.Position] = []
    moving_positions: list[gemmi.Position] = []
    for residue in reference_residues:
        counterpart = moving.get(_residue_key(residue))
        if counterpart is None or counterpart.name != residue.name:
            continue
        fixed_positions.append(required_atom(residue, "CA").pos)
        moving_positions.append(required_atom(counterpart, "CA").pos)
    if len(fixed_positions) < 50:
        raise DiscoveryError(
            "fewer than 50 identity-matched receptor CA atoms are available for alignment"
        )
    return gemmi.superpose_positions(fixed_positions, moving_positions)


def disulfide_ca_distance(
    *, peptide: Residues, chemistry: ChemistryDefinition
) -> float:
    distances: list[float] = []
    for bond in chemistry.disulfide_bonds:
        try:
            first = required_atom(peptide[bond.first - 1], "CA").pos
            second = required_atom(peptide[bond.second - 1], "CA").pos
        except IndexError as exc:
            raise DiscoveryError(
                "disulfide endpoint is outside the CABS-dock peptide"
            ) from exc
        distances.append(first.dist(second))
    if not distances:
        raise DiscoveryError("cyclic-peptide sampling requires a disulfide definition")
    return max(distances)


def disulfide_cabs_sc_distance(
    *, peptide: Residues, chemistry: ChemistryDefinition
) -> float:
    distances: list[float] = []
    for bond in chemistry.disulfide_bonds:
        try:
            first = _side_chain_position(peptide[bond.first - 1])
            second = _side_chain_position(peptide[bond.second - 1])
        except IndexError as exc:
            raise DiscoveryError(
                "disulfide endpoint is outside the CABS-dock peptide"
            ) from exc
        distances.append(first.dist(second))
    if not distances:
        raise DiscoveryError("cyclic-peptide sampling requires a disulfide definition")
    return max(distances)


def _side_chain_position(residue: gemmi.Residue) -> gemmi.Position:
    side_chain = _named_atom(residue, "SC")
    if side_chain is not None:
        return side_chain.pos
    beta_carbon = _named_atom(residue, "CB")
    return (
        beta_carbon.pos if beta_carbon is not None else required_atom(residue, "CA").pos
    )


def contact_residues(
    *, receptor: Residues, peptide: Residues, threshold_A: float
) -> frozenset[str]:
    peptide_positions = tuple(_side_chain_position(residue) for residue in peptide)
    contacts: set[str] = set()
    for residue in receptor:
        position = _side_chain_position(residue)
        if any(
            position.dist(peptide_position) <= threshold_A
            for peptide_position in peptide_positions
        ):
            number, insertion = _residue_key(residue)
            contacts.add(f"{number}{insertion}")
    return frozenset(contacts)


def ca_contact_residues(
    *, receptor: Residues, peptide: Residues, threshold_A: float
) -> frozenset[str]:
    """以两条链的 CA 距离定义粗粒化接触。"""
    peptide_positions = tuple(required_atom(residue, "CA").pos for residue in peptide)
    return frozenset(
        f"{_residue_key(residue)[0]}{_residue_key(residue)[1]}"
        for residue in receptor
        if any(
            required_atom(residue, "CA").pos.dist(position) <= threshold_A
            for position in peptide_positions
        )
    )


def _transformed_ca(
    *, peptide: Residues, transform: gemmi.Transform
) -> PeptideCoordinates:
    return tuple(
        (
            transformed.x,
            transformed.y,
            transformed.z,
        )
        for residue in peptide
        for transformed in (transform.apply(required_atom(residue, "CA").pos),)
    )


def centroid(positions: PeptideCoordinates) -> Point3D:
    count = len(positions)
    return (
        sum(position[0] for position in positions) / count,
        sum(position[1] for position in positions) / count,
        sum(position[2] for position in positions) / count,
    )


def _peptide_rmsd(first: _FrameEvidence, second: _FrameEvidence) -> float:
    if len(first.peptide_ca) != len(second.peptide_ca):
        raise DiscoveryError("candidate peptides have different CA counts")
    return math.sqrt(
        sum(
            math.dist(left, right) ** 2
            for left, right in zip(first.peptide_ca, second.peptide_ca, strict=True)
        )
        / len(first.peptide_ca)
    )


def _geometric_medoid(members: tuple[_FrameEvidence, ...]) -> _FrameEvidence:
    return min(
        members,
        key=lambda candidate: (
            sum(_peptide_rmsd(candidate, other) for other in members),
            candidate.pose_id,
        ),
    )


def _select_candidates(
    *,
    frames: tuple[_FrameEvidence, ...],
    settings: CandidateSelectionSettings,
) -> tuple[tuple[_FrameEvidence, ...], int, int, int]:
    site_clusters = bounded_leader_clusters(
        frames,
        distance=lambda first, second: normalized_site_distance(
            first_contacts=first.contact_residues,
            first_position=first.local_position,
            second_contacts=second.contact_residues,
            second_position=second.local_position,
            contact_limit=settings.contact_jaccard_distance,
            position_limit=settings.position_distance_A,
        ),
        identity=lambda item: item.pose_id,
        maximum_distance=1.0,
    )
    ranked_sites = sorted(
        site_clusters,
        key=lambda cluster: (
            -len(cluster),
            min(item.ranking_score for item in cluster),
            min(item.pose_id for item in cluster),
        ),
    )
    selected: dict[str, _FrameEvidence] = {}
    pose_cluster_count = 0
    selected_pose_cluster_count = 0
    for site_cluster in ranked_sites[: settings.max_sites]:
        pose_clusters = bounded_leader_clusters(
            site_cluster,
            distance=_peptide_rmsd,
            identity=lambda item: item.pose_id,
            maximum_distance=settings.peptide_ca_rmsd_A,
        )
        pose_cluster_count += len(pose_clusters)
        ranked_poses = sorted(
            pose_clusters,
            key=lambda cluster: (
                -len(cluster),
                min(item.ranking_score for item in cluster),
                min(item.pose_id for item in cluster),
            ),
        )
        retained_poses = ranked_poses[: settings.max_pose_clusters_per_site]
        selected_pose_cluster_count += len(retained_poses)
        for pose_cluster in retained_poses:
            medoid = _geometric_medoid(pose_cluster)
            rank_best = min(
                pose_cluster,
                key=lambda item: (item.ranking_score, item.pose_id),
            )
            selected[medoid.pose_id] = medoid
            selected[rank_best.pose_id] = rank_best
    return (
        tuple(selected[pose_id] for pose_id in sorted(selected)),
        len(site_clusters),
        pose_cluster_count,
        selected_pose_cluster_count,
    )


def _percentile(values: tuple[float, ...], fraction: float) -> float:
    ordered = sorted(values)
    index = max(0, math.ceil(fraction * len(ordered)) - 1)
    return ordered[index]


def _frame_evidence(
    *,
    task: DiscoveryTask,
    receptor: Residues,
    peptide: Residues,
    contacts: frozenset[str],
    pose_id: str,
    score: float,
    score_name: str,
    reference_chains: tuple[Residues, ...],
    reference_receptor_id: str,
    qc_status: str,
) -> tuple[_FrameEvidence, float]:
    alignment = align_receptor(receptor=receptor, reference_chains=reference_chains)
    peptide_ca = _transformed_ca(peptide=peptide, transform=alignment.transform)
    return (
        _FrameEvidence(
            task=task,
            pose_id=pose_id,
            identity=None,
            contact_residues=contacts,
            local_position=centroid(peptide_ca),
            ranking_score=score,
            score_name=score_name,
            qc_status=qc_status,
            peptide_ca=peptide_ca,
        ),
        alignment.rmsd,
    )


def _pose_evidence(
    *,
    frame: _FrameEvidence,
    model_path: Path,
    model_sha256: str,
    model_index: int,
    reference_receptor_id: str,
) -> PoseEvidence:
    """在模型产物完成并取得哈希后建立可交付 pose。"""
    return PoseEvidence(
        task_id=frame.task.task_id,
        pose_id=frame.pose_id,
        receptor_id=frame.task.receptor_id,
        target=frame.task.target,
        seed=frame.task.seed,
        model_path=model_path,
        model_sha256=model_sha256,
        model_index=model_index,
        contact_residues=frame.contact_residues,
        local_position=frame.local_position,
        coordinate_frame_id=f"{frame.task.target}:{reference_receptor_id}",
        ranking_score=frame.ranking_score,
        score_name=frame.score_name,
        qc_status=frame.qc_status,
    )


def align_trajectory_receptor(
    *,
    sequence: CabsSequenceChain,
    positions: tuple[Point3D, ...],
    reference_chains: tuple[Residues, ...],
) -> gemmi.SupResult:
    """按作者残基编号和残基名把原始 TRAF 受体对齐到公共坐标系。"""
    if len(positions) != len(sequence.residue_names):
        raise DiscoveryError("CABS receptor coordinates differ from SEQ")
    moving = {
        (number, insertion): (name, position)
        for number, insertion, name, position in zip(
            sequence.residue_numbers,
            sequence.insertion_codes,
            sequence.residue_names,
            positions,
            strict=True,
        )
    }
    candidates = tuple(
        (
            sum(
                _residue_key(residue) in moving
                and moving[_residue_key(residue)][0] == residue.name
                for residue in residues
            ),
            residues,
        )
        for residues in reference_chains
    )
    common_count, reference = max(candidates, key=lambda item: item[0])
    if common_count < 50:
        raise DiscoveryError(
            "fewer than 50 identity-matched receptor CA atoms are available for alignment"
        )
    fixed: list[gemmi.Position] = []
    mobile: list[gemmi.Position] = []
    for residue in reference:
        counterpart = moving.get(_residue_key(residue))
        if counterpart is None or counterpart[0] != residue.name:
            continue
        fixed.append(required_atom(residue, "CA").pos)
        mobile.append(gemmi.Position(*counterpart[1]))
    return gemmi.superpose_positions(fixed, mobile)


def trajectory_ca_contact_residues(
    *,
    sequence: CabsSequenceChain,
    receptor: tuple[Point3D, ...],
    peptide: tuple[Point3D, ...],
    threshold_A: float,
) -> frozenset[str]:
    """以完整 TRAF 可用的 CA 表示定义受体接触残基。"""
    contacts: set[str] = set()
    for number, insertion, position in zip(
        sequence.residue_numbers,
        sequence.insertion_codes,
        receptor,
        strict=True,
    ):
        if any(math.dist(position, ligand) <= threshold_A for ligand in peptide):
            contacts.add(f"{number}{insertion}")
    return frozenset(contacts)


def _trajectory_frame_evidence(
    *,
    task: DiscoveryTask,
    frame: CabsTrajectoryFrame,
    receptor_sequence: CabsSequenceChain,
    reference_chains: tuple[Residues, ...],
    contact_threshold_A: float,
) -> tuple[_FrameEvidence | None, float]:
    receptor = frame.chain_ca[0]
    peptide = frame.chain_ca[-1]
    contacts = trajectory_ca_contact_residues(
        sequence=receptor_sequence,
        receptor=receptor,
        peptide=peptide,
        threshold_A=contact_threshold_A,
    )
    alignment = align_trajectory_receptor(
        sequence=receptor_sequence,
        positions=receptor,
        reference_chains=reference_chains,
    )
    if not contacts:
        return None, alignment.rmsd
    peptide_ca = tuple(
        (transformed.x, transformed.y, transformed.z)
        for position in peptide
        for transformed in (alignment.transform.apply(gemmi.Position(*position)),)
    )
    identity = CabsFrameIdentity(frame.replica, frame.model)
    return (
        _FrameEvidence(
            task=task,
            pose_id=(f"{task.task_id}__traf_r{frame.replica:02d}_m{frame.model:04d}"),
            identity=identity,
            contact_residues=contacts,
            local_position=centroid(peptide_ca),
            ranking_score=frame.interaction_energy,
            score_name="cabsdock_interaction_energy",
            qc_status="passed",
            peptide_ca=peptide_ca,
        ),
        alignment.rmsd,
    )


def collect_cabsdock_evidence(
    *,
    task: DiscoveryTask,
    task_dir: Path,
    reference_path: Path,
    reference_receptor_id: str,
    chemistry: ChemistryDefinition,
    settings: CabsDockSettings,
    selection: CandidateSelectionSettings,
) -> CabsDockEvidence:
    """从完整 Top-10k TRAF 选择候选, 并保留 Top-1000/Top-10 基线。"""
    filtered_path = task_dir / "output_pdbs" / "top1000.pdb"
    filtered = read_structure(filtered_path)
    filtered_count = len(filtered)
    if filtered_count != settings.filtering_count:
        raise DiscoveryError(
            "CABS-dock filtered model count differs from the frozen method: "
            f"expected {settings.filtering_count}, got {filtered_count}"
        )
    per_replica = settings.filtering_count // settings.replicas
    filtered_distances: list[float] = []
    topology_by_replica = [0] * settings.replicas
    for model_index, model in enumerate(filtered, 1):
        _, peptide = split_model(model, peptide_sequence=chemistry.sequence)
        distance = disulfide_ca_distance(peptide=peptide, chemistry=chemistry)
        filtered_distances.append(distance)
        if distance <= settings.max_reconstructable_disulfide_ca_distance_A:
            topology_by_replica[(model_index - 1) // per_replica] += 1

    archive_path = cabsdock_archive_path(task_dir)
    trajectory_audit = audit_cabs_trajectory(
        archive_path=archive_path,
        chemistry=chemistry,
        replicas=settings.replicas,
        filtering_count=settings.filtering_count,
        max_reconstructable_disulfide_ca_distance_A=(
            settings.max_reconstructable_disulfide_ca_distance_A
        ),
        filtered_topology_feasible_by_replica=tuple(topology_by_replica),
    )
    chains = read_cabs_sequence(archive_path)
    if chains[-1].sequence != chemistry.sequence:
        raise DiscoveryError("CABS TRAF peptide differs from the frozen chemistry")
    reference_chains = read_reference_chains(reference_path)
    eligible_frames: list[_FrameEvidence] = []
    trajectory_distances: list[float] = []
    alignment_rmsds: list[float] = []
    for frame in iter_cabs_trajectory(archive_path=archive_path, chains=chains):
        distance = trajectory_disulfide_ca_distance(frame, chemistry=chemistry)
        trajectory_distances.append(distance)
        if distance > settings.max_reconstructable_disulfide_ca_distance_A:
            continue
        evidence, alignment_rmsd = _trajectory_frame_evidence(
            task=task,
            frame=frame,
            receptor_sequence=chains[0],
            reference_chains=reference_chains,
            contact_threshold_A=settings.trajectory_contact_ca_threshold_A,
        )
        alignment_rmsds.append(alignment_rmsd)
        if evidence is not None:
            eligible_frames.append(evidence)
    if len(trajectory_distances) != trajectory_audit.trajectory_model_count:
        raise DiscoveryError("CABS TRAF count changed between audit and selection")

    selection_path = task_dir / "output_data" / "trajectory_candidate_frames.json"
    candidate_models_path: Path | None = None
    if len(eligible_frames) < settings.min_models_for_selection:
        selected_frames: tuple[_FrameEvidence, ...] = ()
        site_count = 0
        pose_count = 0
        selected_pose_count = 0
        selection_status = "skipped_insufficient_models_for_selection"
        atomic_write_json(selection_path, {"frames": []})
        selected_poses: tuple[PoseEvidence, ...] = ()
    else:
        selected_frames, site_count, pose_count, selected_pose_count = (
            _select_candidates(frames=tuple(eligible_frames), settings=selection)
        )
        candidate_budget = (
            selection.max_sites * selection.max_pose_clusters_per_site * 2
        )
        if len(selected_frames) > candidate_budget:
            raise DiscoveryError("candidate selection exceeded its frozen task budget")
        identities = tuple(frame.identity for frame in selected_frames)
        if any(identity is None for identity in identities):
            raise DiscoveryError("trajectory candidate lacks its frame identity")
        typed_identities = tuple(
            identity for identity in identities if identity is not None
        )
        selection_path, candidate_models_path = materialize_cabs_frames(
            archive_path=archive_path,
            identities=typed_identities,
            task_dir=task_dir,
            settings=settings,
        )
        candidate_structure = read_structure(candidate_models_path)
        if len(candidate_structure) != len(selected_frames):
            raise DiscoveryError(
                "materialized CABS candidate count differs from selection"
            )
        for model in candidate_structure:
            split_model(model, peptide_sequence=chemistry.sequence)
        candidate_hash = sha256_file(candidate_models_path)
        selected_poses = tuple(
            _pose_evidence(
                frame=frame,
                model_path=candidate_models_path,
                model_sha256=candidate_hash,
                model_index=index,
                reference_receptor_id=reference_receptor_id,
            )
            for index, frame in enumerate(selected_frames, 1)
        )
        selection_status = "completed"

    medoid_paths = cabsdock_medoid_paths(task_dir)
    if len(medoid_paths) != settings.clustering_medoids:
        raise DiscoveryError(
            "CABS-dock medoid count differs from the frozen method: "
            f"expected {settings.clustering_medoids}, got {len(medoid_paths)}"
        )
    baseline: list[PoseEvidence] = []
    for index, medoid_path in enumerate(medoid_paths, 1):
        structure = read_structure(medoid_path)
        if len(structure) != 1:
            raise DiscoveryError(
                f"CABS-dock medoid file must contain one model: {medoid_path}"
            )
        receptor, peptide = split_model(
            structure[0], peptide_sequence=chemistry.sequence
        )
        feasible = (
            disulfide_ca_distance(peptide=peptide, chemistry=chemistry)
            <= settings.max_reconstructable_disulfide_ca_distance_A
        )
        contacts = ca_contact_residues(
            receptor=receptor,
            peptide=peptide,
            threshold_A=settings.trajectory_contact_ca_threshold_A,
        )
        frame, alignment_rmsd = _frame_evidence(
            task=task,
            receptor=receptor,
            peptide=peptide,
            contacts=contacts,
            pose_id=f"{task.task_id}__baseline_{index:02d}",
            score=float(index),
            score_name="cabsdock_medoid_rank",
            reference_chains=reference_chains,
            reference_receptor_id=reference_receptor_id,
            qc_status="passed" if feasible and contacts else "failed",
        )
        baseline.append(
            _pose_evidence(
                frame=frame,
                model_path=medoid_path,
                model_sha256=sha256_file(medoid_path),
                model_index=1,
                reference_receptor_id=reference_receptor_id,
            )
        )
        alignment_rmsds.append(alignment_rmsd)

    distance_values = tuple(trajectory_distances)
    filtered_topology_count = sum(
        distance <= settings.max_reconstructable_disulfide_ca_distance_A
        for distance in filtered_distances
    )
    return CabsDockEvidence(
        poses=selected_poses,
        baseline_poses=tuple(baseline),
        sampling_model_count=len(distance_values),
        filtered_model_count=filtered_count,
        filtered_topology_feasible_model_count=filtered_topology_count,
        filtered_topology_feasible_fraction=filtered_topology_count / filtered_count,
        topology_feasible_model_count=(
            trajectory_audit.trajectory_topology_feasible_model_count
        ),
        topology_feasible_fraction=(
            trajectory_audit.trajectory_topology_feasible_fraction
        ),
        contacting_topology_feasible_model_count=len(eligible_frames),
        disulfide_ca_distance_median_A=statistics.median(distance_values),
        disulfide_ca_distance_p90_A=_percentile(distance_values, 0.90),
        disulfide_ca_distance_max_A=max(distance_values),
        candidate_frame_selection_path=selection_path,
        candidate_models_path=candidate_models_path,
        trajectory_audit=trajectory_audit,
        selection_status=selection_status,
        site_cluster_count=site_count,
        pose_cluster_count=pose_count,
        selected_pose_cluster_count=selected_pose_count,
        receptor_alignment_max_rmsd_A=(
            max(alignment_rmsds) if alignment_rmsds else None
        ),
    )
