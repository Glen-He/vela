"""局部恢复控制的批内汇总、构象聚类和跨批次判定。"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import JsonValue, atomic_write_text, sha256_file
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.bound_states.controls import ControlTask
from vela.validation.models import LocalRecoveryControl, ValidationError
from vela.validation.scores import (
    RosettaScoreRow,
    index_rosetta_pdb_outputs,
    read_rosetta_scorefile,
)


@dataclass(frozen=True, slots=True)
class ControlDecoy:
    """一个已转换到实验受体坐标系的控制 decoy。"""

    decoy_id: str
    task_id: str
    seed: int
    description: str
    path: Path
    sha256: str
    ranking_score: float
    recovery_rmsd_A: float
    peptide_backbone: tuple[gemmi.Position, ...]


@dataclass(frozen=True, slots=True)
class RecoveryCluster:
    """批内按肽骨架 RMSD 聚合的一组 decoy。"""

    rank: int
    members: tuple[ControlDecoy, ...]

    @property
    def energy_representative(self) -> ControlDecoy:
        """返回用于分数排序的最低能量成员。"""
        return self.members[0]

    @property
    def geometric_medoid(self) -> ControlDecoy:
        """返回到簇内其他成员总骨架 RMSD 最小的几何代表。"""
        return min(
            self.members,
            key=lambda candidate: (
                sum(
                    _backbone_rmsd(candidate.peptide_backbone, member.peptide_backbone)
                    for member in self.members
                ),
                candidate.decoy_id,
            ),
        )

    @property
    def seeds(self) -> tuple[int, ...]:
        return tuple(sorted({member.seed for member in self.members}))

    @property
    def minimum_recovery_rmsd_A(self) -> float:
        return min(member.recovery_rmsd_A for member in self.members)

    def recovered_seeds(self, *, max_recovery_rmsd_A: float) -> tuple[int, ...]:
        return tuple(
            sorted(
                {
                    member.seed
                    for member in self.members
                    if member.recovery_rmsd_A <= max_recovery_rmsd_A
                }
            )
        )

    def recovered_energy_representative(
        self, *, max_recovery_rmsd_A: float
    ) -> ControlDecoy:
        """返回近天然成员中用于证明排序能力的最低能量模型。"""
        recovered = tuple(
            member
            for member in self.members
            if member.recovery_rmsd_A <= max_recovery_rmsd_A
        )
        if not recovered:
            raise ValidationError("near-native cluster contains no recovered decoy")
        return min(recovered, key=lambda item: (item.ranking_score, item.decoy_id))


@dataclass(frozen=True, slots=True)
class BatchRecoveryAssessment:
    """一个完整正式预算批次的采样、排序和聚类证据。"""

    batch_id: str
    seeds: tuple[int, ...]
    decoy_count: int
    cluster_count: int
    sampling_success: bool
    ranking_success: bool
    selected_cluster_rank: int | None
    selected_cluster_seed_support: int
    selected_energy_representative: ControlDecoy | None
    selected_geometric_medoid: ControlDecoy | None
    decoy_table_path: Path
    cluster_table_path: Path


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


def _alignment_coordinates(
    path: Path,
) -> tuple[tuple[gemmi.Position, ...], tuple[gemmi.Position, ...]]:
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(f"invalid control structure: {path}") from exc
    if len(structure) != 1:
        raise ValidationError(f"control structure must contain one model: {path}")
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError(f"control structure must contain A/P chains: {path}")
    receptor = _amino_acids(chains[RECEPTOR_CHAIN])
    peptide = _amino_acids(chains[PEPTIDE_CHAIN])
    receptor_ca = tuple(
        _required_atom(residue, "CA", path=path).pos for residue in receptor
    )
    peptide_backbone = tuple(
        _required_atom(residue, atom_name, path=path).pos
        for residue in peptide
        for atom_name in ("N", "CA", "C")
    )
    if len(receptor_ca) < 3 or not peptide_backbone:
        raise ValidationError(f"control alignment coordinates are incomplete: {path}")
    return receptor_ca, peptide_backbone


def _aligned_backbone(
    *, path: Path, native_receptor_ca: tuple[gemmi.Position, ...]
) -> tuple[gemmi.Position, ...]:
    receptor_ca, peptide_backbone = _alignment_coordinates(path)
    if len(receptor_ca) != len(native_receptor_ca):
        raise ValidationError(f"control receptor CA correspondence is invalid: {path}")
    alignment = gemmi.superpose_positions(native_receptor_ca, receptor_ca)
    return tuple(
        gemmi.Position(alignment.transform.apply(position))
        for position in peptide_backbone
    )


def _backbone_rmsd(
    first: tuple[gemmi.Position, ...], second: tuple[gemmi.Position, ...]
) -> float:
    if len(first) != len(second) or not first:
        raise ValidationError("control peptide backbone correspondence is invalid")
    return math.sqrt(
        sum(left.dist(right) ** 2 for left, right in zip(first, second, strict=True))
        / len(first)
    )


def _decoys_for_task(
    *,
    task: ControlTask,
    task_dir: Path,
    native_receptor_ca: tuple[gemmi.Position, ...],
) -> tuple[ControlDecoy, ...]:
    rows = read_rosetta_scorefile(task_dir / "refine.sc")
    paths = index_rosetta_pdb_outputs(task_dir)
    if {row.description for row in rows} != set(paths):
        raise ValidationError(f"control score/PDB identities differ: {task.task_id}")
    return tuple(
        _control_decoy(
            task=task,
            row=row,
            path=paths[row.description],
            native_receptor_ca=native_receptor_ca,
        )
        for row in rows
    )


def _control_decoy(
    *,
    task: ControlTask,
    row: RosettaScoreRow,
    path: Path,
    native_receptor_ca: tuple[gemmi.Position, ...],
) -> ControlDecoy:
    control = task.control
    return ControlDecoy(
        decoy_id=f"{task.task_id}__{row.description}",
        task_id=task.task_id,
        seed=task.seed,
        description=row.description,
        path=path,
        sha256=sha256_file(path),
        ranking_score=row.score(control.ranking_score),
        recovery_rmsd_A=row.score(control.recovery_rmsd_score),
        peptide_backbone=_aligned_backbone(
            path=path, native_receptor_ca=native_receptor_ca
        ),
    )


def _cluster_decoys(
    *, decoys: tuple[ControlDecoy, ...], threshold_A: float
) -> tuple[RecoveryCluster, ...]:
    """按分数顺序执行确定性的 complete-link 阈值聚类。"""
    groups: list[list[ControlDecoy]] = []
    ranked = sorted(decoys, key=lambda item: (item.ranking_score, item.decoy_id))
    for decoy in ranked:
        for group in groups:
            if all(
                _backbone_rmsd(decoy.peptide_backbone, member.peptide_backbone)
                <= threshold_A
                for member in group
            ):
                group.append(decoy)
                break
        else:
            groups.append([decoy])
    return tuple(
        RecoveryCluster(index, tuple(group)) for index, group in enumerate(groups, 1)
    )


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    lines = ["\t".join(fields) + "\n"]
    lines.extend("\t".join(row[field] for field in fields) + "\n" for row in rows)
    return "".join(lines)


def _write_batch_tables(
    *,
    control: LocalRecoveryControl,
    batch_id: str,
    clusters: tuple[RecoveryCluster, ...],
    output_dir: Path,
    run_dir: Path,
) -> tuple[Path, Path]:
    cluster_by_decoy = {
        member.decoy_id: cluster.rank
        for cluster in clusters
        for member in cluster.members
    }
    decoy_fields = (
        "decoy_id",
        "task_id",
        "seed",
        "ranking_score",
        "recovery_rmsd_A",
        "recovered",
        "cluster_rank",
        "path",
        "sha256",
    )
    decoy_rows = [
        {
            "decoy_id": member.decoy_id,
            "task_id": member.task_id,
            "seed": str(member.seed),
            "ranking_score": f"{member.ranking_score:.6f}",
            "recovery_rmsd_A": f"{member.recovery_rmsd_A:.6f}",
            "recovered": str(
                member.recovery_rmsd_A <= control.max_recovery_rmsd_A
            ).lower(),
            "cluster_rank": str(cluster_by_decoy[member.decoy_id]),
            "path": member.path.relative_to(run_dir).as_posix(),
            "sha256": member.sha256,
        }
        for cluster in clusters
        for member in cluster.members
    ]
    decoy_path = output_dir / f"{batch_id}_decoys.tsv"
    atomic_write_text(decoy_path, _tsv(decoy_fields, decoy_rows))

    cluster_fields = (
        "cluster_rank",
        "decoy_count",
        "all_seeds",
        "recovered_seeds",
        "recovered_seed_support",
        "best_ranking_score",
        "minimum_recovery_rmsd_A",
        "near_native",
        "supported",
        "in_ranking_window",
        "energy_representative_decoy_id",
        "geometric_medoid_decoy_id",
    )
    cluster_rows = [
        {
            "cluster_rank": str(cluster.rank),
            "decoy_count": str(len(cluster.members)),
            "all_seeds": ",".join(str(seed) for seed in cluster.seeds),
            "recovered_seeds": ",".join(
                str(seed)
                for seed in cluster.recovered_seeds(
                    max_recovery_rmsd_A=control.max_recovery_rmsd_A
                )
            ),
            "recovered_seed_support": str(
                len(
                    cluster.recovered_seeds(
                        max_recovery_rmsd_A=control.max_recovery_rmsd_A
                    )
                )
            ),
            "best_ranking_score": (
                f"{cluster.energy_representative.ranking_score:.6f}"
            ),
            "minimum_recovery_rmsd_A": f"{cluster.minimum_recovery_rmsd_A:.6f}",
            "near_native": str(
                cluster.minimum_recovery_rmsd_A <= control.max_recovery_rmsd_A
            ).lower(),
            "supported": str(
                len(
                    cluster.recovered_seeds(
                        max_recovery_rmsd_A=control.max_recovery_rmsd_A
                    )
                )
                >= control.min_cluster_seed_support
            ).lower(),
            "in_ranking_window": str(cluster.rank <= control.top_clusters).lower(),
            "energy_representative_decoy_id": (cluster.energy_representative.decoy_id),
            "geometric_medoid_decoy_id": cluster.geometric_medoid.decoy_id,
        }
        for cluster in clusters
    ]
    cluster_path = output_dir / f"{batch_id}_clusters.tsv"
    atomic_write_text(cluster_path, _tsv(cluster_fields, cluster_rows))
    return decoy_path, cluster_path


def _assess_batch(
    *,
    control: LocalRecoveryControl,
    batch_id: str,
    tasks: tuple[ControlTask, ...],
    run_dir: Path,
    output_dir: Path,
    native_receptor_ca: tuple[gemmi.Position, ...],
) -> BatchRecoveryAssessment:
    decoys = tuple(
        decoy
        for task in tasks
        for decoy in _decoys_for_task(
            task=task,
            task_dir=run_dir / "tasks" / task.task_id,
            native_receptor_ca=native_receptor_ca,
        )
    )
    clusters = _cluster_decoys(
        decoys=decoys, threshold_A=control.max_cluster_backbone_rmsd_A
    )
    selected = next(
        (
            cluster
            for cluster in clusters[: control.top_clusters]
            if cluster.minimum_recovery_rmsd_A <= control.max_recovery_rmsd_A
            and len(
                cluster.recovered_seeds(max_recovery_rmsd_A=control.max_recovery_rmsd_A)
            )
            >= control.min_cluster_seed_support
        ),
        None,
    )
    decoy_path, cluster_path = _write_batch_tables(
        control=control,
        batch_id=batch_id,
        clusters=clusters,
        output_dir=output_dir,
        run_dir=run_dir,
    )
    return BatchRecoveryAssessment(
        batch_id=batch_id,
        seeds=tuple(task.seed for task in tasks),
        decoy_count=len(decoys),
        cluster_count=len(clusters),
        sampling_success=any(
            decoy.recovery_rmsd_A <= control.max_recovery_rmsd_A for decoy in decoys
        ),
        ranking_success=selected is not None,
        selected_cluster_rank=selected.rank if selected is not None else None,
        selected_cluster_seed_support=(
            len(
                selected.recovered_seeds(
                    max_recovery_rmsd_A=control.max_recovery_rmsd_A
                )
            )
            if selected is not None
            else 0
        ),
        selected_energy_representative=(
            selected.recovered_energy_representative(
                max_recovery_rmsd_A=control.max_recovery_rmsd_A
            )
            if selected is not None
            else None
        ),
        selected_geometric_medoid=(
            selected.geometric_medoid if selected is not None else None
        ),
        decoy_table_path=decoy_path,
        cluster_table_path=cluster_path,
    )


def analyze_control_recovery(
    *,
    control: LocalRecoveryControl,
    tasks: tuple[ControlTask, ...],
    run_dir: Path,
    output_dir: Path | None = None,
) -> tuple[bool, dict[str, JsonValue]]:
    """合并每批全部 seed; 要求两个完整预算批次重复恢复同一姿态。"""
    expected_tasks = {
        (f"batch_{batch_index:02d}", seed)
        for batch_index, batch in enumerate(control.seed_batches, 1)
        for seed in batch
    }
    actual_tasks = {(task.batch_id, task.seed) for task in tasks}
    if not tasks or actual_tasks != expected_tasks or len(tasks) != len(actual_tasks):
        raise ValidationError(
            f"control task matrix is incomplete or duplicated: {control.control_id}"
        )
    native_receptor_ca, _ = _alignment_coordinates(tasks[0].control_input.complex_path)
    control_output_dir = (
        output_dir / control.control_id
        if output_dir is not None
        else run_dir / "analysis" / control.control_id
    )
    control_output_dir.mkdir(parents=True, exist_ok=False)
    batches = tuple(
        _assess_batch(
            control=control,
            batch_id=batch_id,
            tasks=tuple(task for task in tasks if task.batch_id == batch_id),
            run_dir=run_dir,
            output_dir=control_output_dir,
            native_receptor_ca=native_receptor_ca,
        )
        for batch_id in tuple(dict.fromkeys(task.batch_id for task in tasks))
    )
    geometric_medoids = tuple(
        batch.selected_geometric_medoid
        for batch in batches
        if batch.selected_geometric_medoid is not None
    )
    pairwise_rmsd = tuple(
        _backbone_rmsd(first.peptide_backbone, second.peptide_backbone)
        for first_index, first in enumerate(geometric_medoids)
        for second in geometric_medoids[first_index + 1 :]
    )
    pose_consistent = len(geometric_medoids) == len(batches) and all(
        value <= control.max_batch_pose_rmsd_A for value in pairwise_rmsd
    )
    passed = (
        all(batch.sampling_success for batch in batches)
        and all(batch.ranking_success for batch in batches)
        and pose_consistent
    )
    return passed, {
        "control_id": control.control_id,
        "bound_state_id": control.bound_state_id,
        "passed": passed,
        "all_batches_sampled": all(batch.sampling_success for batch in batches),
        "all_batches_ranked": all(batch.ranking_success for batch in batches),
        "batch_poses_consistent": pose_consistent,
        "maximum_observed_batch_pose_rmsd_A": (
            max(pairwise_rmsd) if pairwise_rmsd else None
        ),
        "batches": [
            {
                "batch_id": batch.batch_id,
                "seeds": list(batch.seeds),
                "decoy_count": batch.decoy_count,
                "cluster_count": batch.cluster_count,
                "sampling_success": batch.sampling_success,
                "ranking_success": batch.ranking_success,
                "selected_cluster_rank": batch.selected_cluster_rank,
                "selected_cluster_seed_support": batch.selected_cluster_seed_support,
                "selected_energy_representative_decoy_id": (
                    batch.selected_energy_representative.decoy_id
                    if batch.selected_energy_representative is not None
                    else None
                ),
                "selected_geometric_medoid_decoy_id": (
                    batch.selected_geometric_medoid.decoy_id
                    if batch.selected_geometric_medoid is not None
                    else None
                ),
                "decoy_table": {
                    "path": batch.decoy_table_path.relative_to(run_dir).as_posix(),
                    "sha256": sha256_file(batch.decoy_table_path),
                },
                "cluster_table": {
                    "path": batch.cluster_table_path.relative_to(run_dir).as_posix(),
                    "sha256": sha256_file(batch.cluster_table_path),
                },
            }
            for batch in batches
        ],
    }
