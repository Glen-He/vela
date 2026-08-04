"""读取 CABS 原始 TRAF, 并审计能量过滤前后的粗粒化拓扑。"""

from __future__ import annotations

import io
import math
import statistics
import tarfile
from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import JsonValue
from vela.discovery.analysis.evidence import Point3D
from vela.discovery.models import DiscoveryError
from vela.preparation.chemistry import ChemistryDefinition

CABS_GRID_A = 0.61
# CABS PDB 写出保留三位小数; 两个端点逐轴舍入造成的最大距离误差小于 0.002 Å。
PDB_DISTANCE_ROUNDING_TOLERANCE_A = 0.002


@dataclass(frozen=True, slots=True)
class CabsSequenceChain:
    """TRAF 坐标顺序中的一条 SEQ 链。"""

    chain_id: str
    residue_names: tuple[str, ...]
    residue_numbers: tuple[int, ...]
    insertion_codes: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.residue_names or not (
            len(self.residue_names)
            == len(self.residue_numbers)
            == len(self.insertion_codes)
        ):
            raise DiscoveryError("CABS SEQ chain residue fields have different lengths")

    @property
    def sequence(self) -> str:
        """返回标准氨基酸单字母序列。"""
        letters: list[str] = []
        for name in self.residue_names:
            residue = gemmi.find_tabulated_residue(name)
            if not residue.is_amino_acid() or residue.one_letter_code == "X":
                raise DiscoveryError(f"unsupported residue in CABS SEQ: {name}")
            letters.append(residue.one_letter_code)
        return "".join(letters)


@dataclass(frozen=True, slots=True)
class CabsTrajectoryFrame:
    """一个 TRAF 帧的真实 replica、能量和逐链 CA 坐标。"""

    model: int
    replica: int
    interaction_energy: float
    interaction_energy_asymmetry: float
    chain_ca: tuple[tuple[Point3D, ...], ...]


@dataclass(frozen=True, slots=True)
class ReplicaTrajectoryAudit:
    """单个 replica 在完整轨迹和能量过滤后的审计结果。"""

    replica: int
    trajectory_model_count: int
    trajectory_topology_feasible_model_count: int
    filtered_model_count: int
    filtered_topology_feasible_model_count: int
    replay_topology_feasible_min: int
    replay_topology_feasible_max: int
    filtered_interaction_energy_min: float
    filtered_interaction_energy_median: float
    filtered_interaction_energy_max: float


@dataclass(frozen=True, slots=True)
class CabsTrajectoryAudit:
    """完整 TRAF 与 CABS each-replica 能量过滤的可重放审计。"""

    trajectory_model_count: int
    trajectory_topology_feasible_model_count: int
    trajectory_topology_feasible_fraction: float
    filtered_model_count: int
    filtered_topology_feasible_model_count: int
    filtered_topology_feasible_fraction: float
    topology_feasible_enrichment_ratio: float | None
    max_interaction_energy_asymmetry: float
    filtered_interaction_energies: tuple[float, ...]
    replicas: tuple[ReplicaTrajectoryAudit, ...]


@dataclass(frozen=True, slots=True)
class _TrafHeader:
    model: int
    raw_length: int
    energy_row: tuple[float, ...]
    temperature: float
    replica: int


def _archive_member(archive: tarfile.TarFile, name: str) -> tarfile.TarInfo:
    matches = tuple(member for member in archive.getmembers() if member.name == name)
    if len(matches) != 1 or not matches[0].isfile():
        raise DiscoveryError(f"CABS archive must contain one regular {name} member")
    return matches[0]


def read_cabs_sequence(archive_path: Path) -> tuple[CabsSequenceChain, ...]:
    """读取归档内 SEQ, 并保留原受体编号用于 native 事后评价。"""
    try:
        with tarfile.open(archive_path, mode="r:*") as archive:
            member = _archive_member(archive, "SEQ")
            binary = archive.extractfile(member)
            if binary is None:
                raise DiscoveryError("CABS archive SEQ member is unreadable")
            rows = binary.read().decode("utf-8").splitlines()
    except (OSError, tarfile.TarError, UnicodeDecodeError) as exc:
        raise DiscoveryError(f"invalid CABS archive: {archive_path}") from exc
    grouped: list[tuple[str, list[str], list[int], list[str]]] = []
    for line in rows:
        if len(line) < 13:
            raise DiscoveryError("CABS SEQ contains a truncated residue record")
        chain_id = line[12].strip()
        residue_name = line[8:11].strip()
        try:
            residue_number = int(line[1:5])
        except ValueError as exc:
            raise DiscoveryError("CABS SEQ residue number is invalid") from exc
        insertion_code = line[5].strip()
        if not chain_id or not residue_name:
            raise DiscoveryError("CABS SEQ residue identity is incomplete")
        if not grouped or grouped[-1][0] != chain_id:
            if any(existing == chain_id for existing, *_ in grouped):
                raise DiscoveryError("CABS SEQ chain records are not contiguous")
            grouped.append((chain_id, [], [], []))
        grouped[-1][1].append(residue_name)
        grouped[-1][2].append(residue_number)
        grouped[-1][3].append(insertion_code)
    if len(grouped) != 2:
        raise DiscoveryError(
            "CABS trajectory audit requires one receptor chain and one peptide chain"
        )
    return tuple(
        CabsSequenceChain(
            chain_id,
            tuple(names),
            tuple(numbers),
            tuple(insertion_codes),
        )
        for chain_id, names, numbers, insertion_codes in grouped
    )


def _parse_header(line: str, *, chain_count: int) -> _TrafHeader:
    fields = line.split()
    expected = chain_count + 4
    if len(fields) != expected:
        raise DiscoveryError(
            f"CABS TRAF header has {len(fields)} fields; expected {expected}"
        )
    try:
        model = int(fields[0])
        raw_length = int(fields[1])
        energy = tuple(float(value) for value in fields[2:-2])
        temperature = float(fields[-2])
        replica = int(fields[-1])
    except ValueError as exc:
        raise DiscoveryError("CABS TRAF header contains an invalid number") from exc
    if model < 1 or replica < 1 or raw_length < 3:
        raise DiscoveryError("CABS TRAF header contains an invalid identity or length")
    if not all(math.isfinite(value) for value in (*energy, temperature)):
        raise DiscoveryError("CABS TRAF header contains a non-finite value")
    return _TrafHeader(model, raw_length, energy, temperature, replica)


def _coordinate_block(
    values: list[int], *, header: _TrafHeader, expected_residues: int
) -> tuple[Point3D, ...]:
    if len(values) != header.raw_length * 3:
        raise DiscoveryError(
            "CABS TRAF coordinate count differs from its declared chain length"
        )
    physical = values[3:-3]
    if len(physical) != expected_residues * 3:
        raise DiscoveryError("CABS TRAF and SEQ chain lengths differ")
    return tuple(
        (
            physical[index] * CABS_GRID_A,
            physical[index + 1] * CABS_GRID_A,
            physical[index + 2] * CABS_GRID_A,
        )
        for index in range(0, len(physical), 3)
    )


def _blocks(
    handle: io.TextIOBase, *, chain_count: int
) -> Iterator[tuple[_TrafHeader, list[int]]]:
    current: _TrafHeader | None = None
    coordinates: list[int] = []
    for line in handle:
        if "." in line:
            if current is not None:
                yield current, coordinates
            current = _parse_header(line, chain_count=chain_count)
            coordinates = []
            continue
        if current is None:
            raise DiscoveryError("CABS TRAF coordinates precede the first header")
        try:
            coordinates.extend(int(value) for value in line.split())
        except ValueError as exc:
            raise DiscoveryError("CABS TRAF contains an invalid coordinate") from exc
    if current is None:
        raise DiscoveryError("CABS TRAF is empty")
    yield current, coordinates


def iter_cabs_trajectory(
    *, archive_path: Path, chains: tuple[CabsSequenceChain, ...]
) -> Iterator[CabsTrajectoryFrame]:
    """逐帧读取 TRAF, 不把完整坐标轨迹常驻内存。"""
    if not archive_path.is_file():
        raise DiscoveryError(f"CABS archive does not exist: {archive_path}")
    try:
        with tarfile.open(archive_path, mode="r:*") as archive:
            member = _archive_member(archive, "TRAF")
            binary = archive.extractfile(member)
            if binary is None:
                raise DiscoveryError("CABS archive TRAF member is unreadable")
            with io.TextIOWrapper(binary, encoding="utf-8") as handle:
                pending_headers: list[_TrafHeader] = []
                pending_coordinates: list[tuple[Point3D, ...]] = []
                for header, values in _blocks(handle, chain_count=len(chains)):
                    if pending_headers and (
                        header.model != pending_headers[0].model
                        or header.replica != pending_headers[0].replica
                    ):
                        raise DiscoveryError(
                            "CABS TRAF chain blocks do not share a frame identity"
                        )
                    chain_index = len(pending_headers)
                    if chain_index >= len(chains):
                        raise DiscoveryError("CABS TRAF frame contains too many chains")
                    pending_headers.append(header)
                    pending_coordinates.append(
                        _coordinate_block(
                            values,
                            header=header,
                            expected_residues=len(chains[chain_index].residue_names),
                        )
                    )
                    if len(pending_headers) != len(chains):
                        continue
                    temperatures = {item.temperature for item in pending_headers}
                    if len(temperatures) != 1:
                        raise DiscoveryError(
                            "CABS TRAF chain blocks have different temperatures"
                        )
                    matrix = tuple(item.energy_row for item in pending_headers)
                    if any(len(row) != len(chains) for row in matrix):
                        raise DiscoveryError("CABS TRAF energy matrix is not square")
                    yield CabsTrajectoryFrame(
                        model=pending_headers[0].model,
                        replica=pending_headers[0].replica,
                        # 与上游 Header.get_energy 完全一致: 取受体行、肽列。
                        interaction_energy=matrix[0][-1],
                        interaction_energy_asymmetry=abs(matrix[0][-1] - matrix[-1][0]),
                        chain_ca=tuple(pending_coordinates),
                    )
                    pending_headers = []
                    pending_coordinates = []
                if pending_headers:
                    raise DiscoveryError("CABS TRAF ends inside an incomplete frame")
    except (OSError, tarfile.TarError, UnicodeDecodeError) as exc:
        raise DiscoveryError(f"invalid CABS archive: {archive_path}") from exc


def trajectory_disulfide_ca_distance(
    frame: CabsTrajectoryFrame, *, chemistry: ChemistryDefinition
) -> float:
    peptide = frame.chain_ca[-1]
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
        raise DiscoveryError("CABS cyclic-peptide audit requires a disulfide")
    return max(distances)


def audit_cabs_trajectory(
    *,
    archive_path: Path,
    chemistry: ChemistryDefinition,
    replicas: int,
    filtering_count: int,
    max_reconstructable_disulfide_ca_distance_A: float,
    filtered_topology_feasible_by_replica: tuple[int, ...],
) -> CabsTrajectoryAudit:
    """重放 CABS each-replica 过滤, 并返回真实能量和拓扑统计。"""
    chains = read_cabs_sequence(archive_path)
    if chains[-1].sequence != chemistry.sequence:
        raise DiscoveryError("CABS SEQ peptide differs from the frozen chemistry")
    if filtering_count % replicas != 0:
        raise DiscoveryError("CABS filtering count must be divisible by replicas")
    if len(filtered_topology_feasible_by_replica) != replicas or any(
        value < 0 for value in filtered_topology_feasible_by_replica
    ):
        raise DiscoveryError(
            "filtered topology counts must contain one non-negative value per replica"
        )
    by_replica: dict[int, list[tuple[int, float, float]]] = defaultdict(list)
    identities: set[tuple[int, int]] = set()
    max_energy_asymmetry = 0.0
    for frame in iter_cabs_trajectory(archive_path=archive_path, chains=chains):
        identity = (frame.replica, frame.model)
        if identity in identities:
            raise DiscoveryError(
                f"CABS TRAF contains a duplicate frame: {identity[0]}/{identity[1]}"
            )
        identities.add(identity)
        max_energy_asymmetry = max(
            max_energy_asymmetry, frame.interaction_energy_asymmetry
        )
        by_replica[frame.replica].append(
            (
                frame.model,
                frame.interaction_energy,
                trajectory_disulfide_ca_distance(frame, chemistry=chemistry),
            )
        )
    expected_replicas = set(range(1, replicas + 1))
    if set(by_replica) != expected_replicas:
        raise DiscoveryError("CABS TRAF replica identities differ from the method")
    model_counts = {len(frames) for frames in by_replica.values()}
    if len(model_counts) != 1:
        raise DiscoveryError("CABS TRAF replicas contain different model counts")
    per_replica_filter = filtering_count // replicas
    replica_audits: list[ReplicaTrajectoryAudit] = []
    filtered_energies: list[float] = []
    trajectory_topology_count = 0
    filtered_topology_count = 0
    for replica in sorted(by_replica):
        frames = by_replica[replica]
        if len(frames) < per_replica_filter:
            raise DiscoveryError("CABS TRAF is smaller than the frozen filter budget")
        topology_count = sum(
            distance <= max_reconstructable_disulfide_ca_distance_A
            for _, _, distance in frames
        )
        # 上游 np.argsort 的并列顺序不稳定, 但并列项能量相同; 这里只映射真实能量值。
        ordered = sorted(frames, key=lambda item: item[1])
        filtered = ordered[:per_replica_filter]
        energies = tuple(item[1] for item in filtered)
        cutoff = energies[-1]
        mandatory = tuple(item for item in frames if item[1] < cutoff)
        tied = tuple(item for item in frames if item[1] == cutoff)
        tied_needed = per_replica_filter - len(mandatory)
        lower = (
            max_reconstructable_disulfide_ca_distance_A
            - PDB_DISTANCE_ROUNDING_TOLERANCE_A
        )
        upper = (
            max_reconstructable_disulfide_ca_distance_A
            + PDB_DISTANCE_ROUNDING_TOLERANCE_A
        )
        mandatory_definite = sum(distance <= lower for _, _, distance in mandatory)
        mandatory_possible = sum(distance <= upper for _, _, distance in mandatory)
        tied_definite = sum(distance <= lower for _, _, distance in tied)
        tied_possible = sum(distance <= upper for _, _, distance in tied)
        replay_min = mandatory_definite + max(
            0, tied_needed - (len(tied) - tied_definite)
        )
        replay_max = mandatory_possible + min(tied_needed, tied_possible)
        filtered_topology = filtered_topology_feasible_by_replica[replica - 1]
        if not replay_min <= filtered_topology <= replay_max:
            raise DiscoveryError(
                "top1000 topology count is inconsistent with CABS energy filtering: "
                f"replica={replica}, actual={filtered_topology}, "
                f"expected={replay_min}..{replay_max}"
            )
        filtered_energies.extend(energies)
        trajectory_topology_count += topology_count
        filtered_topology_count += filtered_topology
        replica_audits.append(
            ReplicaTrajectoryAudit(
                replica=replica,
                trajectory_model_count=len(frames),
                trajectory_topology_feasible_model_count=topology_count,
                filtered_model_count=len(filtered),
                filtered_topology_feasible_model_count=filtered_topology,
                replay_topology_feasible_min=replay_min,
                replay_topology_feasible_max=replay_max,
                filtered_interaction_energy_min=min(energies),
                filtered_interaction_energy_median=statistics.median(energies),
                filtered_interaction_energy_max=max(energies),
            )
        )
    trajectory_count = sum(item.trajectory_model_count for item in replica_audits)
    trajectory_fraction = trajectory_topology_count / trajectory_count
    filtered_fraction = filtered_topology_count / filtering_count
    return CabsTrajectoryAudit(
        trajectory_model_count=trajectory_count,
        trajectory_topology_feasible_model_count=trajectory_topology_count,
        trajectory_topology_feasible_fraction=trajectory_fraction,
        filtered_model_count=filtering_count,
        filtered_topology_feasible_model_count=filtered_topology_count,
        filtered_topology_feasible_fraction=filtered_fraction,
        topology_feasible_enrichment_ratio=(
            filtered_fraction / trajectory_fraction
            if trajectory_fraction > 0.0
            else None
        ),
        max_interaction_energy_asymmetry=max_energy_asymmetry,
        filtered_interaction_energies=tuple(filtered_energies),
        replicas=tuple(replica_audits),
    )


def trajectory_audit_record(audit: CabsTrajectoryAudit) -> dict[str, JsonValue]:
    """将完整轨迹审计转换为稳定 JSON 合同。"""
    return {
        "schema": "vela.cabs-trajectory-audit/2",
        "energy_definition": "CABS_Header.get_energy_receptor_row_peptide_column",
        "max_reverse_matrix_cell_difference": round(
            audit.max_interaction_energy_asymmetry, 8
        ),
        "filtering_mode": "lowest_interaction_energy_per_replica",
        "trajectory": {
            "model_count": audit.trajectory_model_count,
            "topology_feasible_model_count": (
                audit.trajectory_topology_feasible_model_count
            ),
            "topology_feasible_fraction": round(
                audit.trajectory_topology_feasible_fraction, 8
            ),
        },
        "filtered": {
            "model_count": audit.filtered_model_count,
            "topology_feasible_model_count": (
                audit.filtered_topology_feasible_model_count
            ),
            "topology_feasible_fraction": round(
                audit.filtered_topology_feasible_fraction, 8
            ),
            "topology_feasible_enrichment_ratio": (
                round(audit.topology_feasible_enrichment_ratio, 8)
                if audit.topology_feasible_enrichment_ratio is not None
                else None
            ),
        },
        "replicas": [
            {
                "replica": item.replica,
                "trajectory_model_count": item.trajectory_model_count,
                "trajectory_topology_feasible_model_count": (
                    item.trajectory_topology_feasible_model_count
                ),
                "filtered_model_count": item.filtered_model_count,
                "filtered_topology_feasible_model_count": (
                    item.filtered_topology_feasible_model_count
                ),
                "filter_replay_topology_feasible_range": [
                    item.replay_topology_feasible_min,
                    item.replay_topology_feasible_max,
                ],
                "filtered_interaction_energy_min": round(
                    item.filtered_interaction_energy_min, 6
                ),
                "filtered_interaction_energy_median": round(
                    item.filtered_interaction_energy_median, 6
                ),
                "filtered_interaction_energy_max": round(
                    item.filtered_interaction_energy_max, 6
                ),
            }
            for item in audit.replicas
        ],
    }
