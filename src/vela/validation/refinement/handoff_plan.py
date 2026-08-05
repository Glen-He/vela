"""阶段二受支持 site 到局部精修交接任务的选择与冻结。"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.discovery.analysis.cluster_engine import bounded_leader_clusters
from vela.discovery.analysis.evidence import PoseEvidence, ReceptorSite
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.analysis.reports import (
    ReportedCandidateSite,
    ReportedReceptorSite,
    read_site_analysis_report,
)
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.evidence import (
    align_receptor,
    read_reference_chains,
    read_structure,
    required_atom,
    split_model,
)
from vela.discovery.sampling.planning import MAIN_DISCOVERY_EVIDENCE
from vela.validation.models import ValidationError
from vela.validation.records import safe_identifier
from vela.validation.refinement.reconstruction import verify_cg2all_tool
from vela.validation.rosetta import verify_rosetta_scripts_tool

HANDOFF_PLAN_NAME = "handoff_plan.json"


@dataclass(frozen=True, slots=True)
class CandidateHandoffTask:
    """一个保留 blind 证据身份的全原子重建任务。"""

    task_id: str
    candidate_id: str
    receptor_site_id: str
    pose: PoseEvidence
    reference_receptor_path: Path
    reference_receptor_sha256: str


@dataclass(frozen=True, slots=True)
class CandidateHandoffPlan:
    """已经冻结到独立运行目录的候选交接计划。"""

    run_id: str
    run_dir: Path
    discovery_run_dir: Path
    tasks: tuple[CandidateHandoffTask, ...]


def _identifier(value: str, *, name: str) -> str:
    return safe_identifier(value, name=name)


def _selected_candidates(
    *, candidates: dict[str, ReportedCandidateSite], requested_ids: tuple[str, ...]
) -> tuple[ReportedCandidateSite, ...]:
    requested_ids = tuple(
        _identifier(value, name="requested candidate ID") for value in requested_ids
    )
    if len(requested_ids) != len(set(requested_ids)):
        raise ValidationError("requested candidate IDs must be unique")
    if not requested_ids:
        raise ValidationError(
            "at least one explicit candidate ID is required for all-atom handoff"
        )
    unknown = sorted(set(requested_ids) - set(candidates))
    if unknown:
        raise ValidationError("unknown candidate IDs: " + ", ".join(unknown))
    selected = tuple(candidates[candidate_id] for candidate_id in requested_ids)
    if not selected:
        raise ValidationError("no supported candidate sites are available for handoff")
    unsupported = [item.candidate_id for item in selected if not item.supported]
    if unsupported:
        raise ValidationError(
            "unsupported candidate sites cannot enter blind handoff: "
            + ", ".join(unsupported)
        )
    return selected


def select_site_poses(
    *,
    site: ReceptorSite | ReportedReceptorSite,
    pose_by_id: dict[str, PoseEvidence],
    count: int,
    peptide_sequence: str,
    reference_receptor_path: Path,
    pose_clustering_rmsd_A: float,
) -> tuple[PoseEvidence, ...]:
    """从不同姿态簇选择跨 seed 的几何 medoid。"""
    if not site.supported:
        raise ValidationError(f"unsupported receptor site in candidate: {site.site_id}")
    if count < 1 or pose_clustering_rmsd_A <= 0:
        raise ValidationError("handoff pose selection settings must be positive")
    missing = sorted(set(site.pose_ids) - set(pose_by_id))
    if missing:
        raise ValidationError(f"receptor site refers to unknown pose: {missing[0]}")
    poses = tuple(
        pose_by_id[pose_id]
        for pose_id in site.pose_ids
        if pose_by_id[pose_id].qc_status == "passed"
    )
    reference = read_reference_chains(reference_receptor_path)
    structures: dict[Path, gemmi.Structure] = {}
    coordinates: dict[str, tuple[gemmi.Position, ...]] = {}
    for pose in poses:
        structure = structures.get(pose.model_path)
        if structure is None:
            structure = read_structure(pose.model_path)
            structures[pose.model_path] = structure
        if pose.model_index > len(structure):
            raise ValidationError(
                f"handoff pose model index is outside its structure: {pose.pose_id}"
            )
        receptor, peptide = split_model(
            structure[pose.model_index - 1], peptide_sequence=peptide_sequence
        )
        alignment = align_receptor(receptor=receptor, reference_chains=reference)
        coordinates[pose.pose_id] = tuple(
            gemmi.Position(alignment.transform.apply(required_atom(residue, "CA").pos))
            for residue in peptide
        )

    def distance(first: PoseEvidence, second: PoseEvidence) -> float:
        left = coordinates[first.pose_id]
        right = coordinates[second.pose_id]
        if len(left) != len(right) or not left:
            raise ValidationError("handoff peptide CA identities differ")
        return math.sqrt(
            sum(
                first_position.dist(second_position) ** 2
                for first_position, second_position in zip(left, right, strict=True)
            )
            / len(left)
        )

    clusters = bounded_leader_clusters(
        poses,
        distance=distance,
        identity=lambda pose: pose.pose_id,
        maximum_distance=pose_clustering_rmsd_A,
    )
    ranked_clusters = sorted(
        clusters,
        key=lambda cluster: (
            -len({pose.seed for pose in cluster}),
            -len(cluster),
            min(pose.ranking_score for pose in cluster),
            min(pose.pose_id for pose in cluster),
        ),
    )
    selected: list[PoseEvidence] = []
    seeds: set[int] = set()
    for cluster in ranked_clusters:
        eligible = tuple(pose for pose in cluster if pose.seed not in seeds)
        if eligible:
            medoid = min(
                eligible,
                key=lambda pose: (
                    sum(distance(pose, other) for other in cluster),
                    pose.ranking_score,
                    pose.pose_id,
                ),
            )
            selected.append(medoid)
            seeds.add(medoid.seed)
        if len(selected) == count:
            break
    if len(selected) < count:
        selected_ids = {pose.pose_id for pose in selected}
        for cluster in ranked_clusters:
            eligible = tuple(
                pose
                for pose in cluster
                if pose.pose_id not in selected_ids and pose.seed not in seeds
            )
            if eligible:
                medoid = min(
                    eligible,
                    key=lambda pose: (
                        sum(distance(pose, other) for other in cluster),
                        pose.ranking_score,
                        pose.pose_id,
                    ),
                )
                selected.append(medoid)
                selected_ids.add(medoid.pose_id)
                seeds.add(medoid.seed)
            if len(selected) == count:
                break
    if len(selected) != count:
        raise ValidationError(
            f"{site.site_id} lacks {count} passed poses from distinct seeds"
        )
    return tuple(selected)


def build_handoff_tasks(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    candidate_ids: tuple[str, ...],
) -> tuple[CandidateHandoffTask, ...]:
    """只从受支持 blind site 选择跨 seed 的非冗余起点。"""
    try:
        report = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=MAIN_DISCOVERY_EVIDENCE,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    poses = read_pose_evidence(path=report.pose_path, run_dir=discovery_run_dir)
    pose_by_id = {pose.pose_id: pose for pose in poses}
    if len(pose_by_id) != len(poses):
        raise ValidationError("discovery pose IDs are not unique")
    candidates = _selected_candidates(
        candidates=report.candidate_sites, requested_ids=candidate_ids
    )
    tasks: list[CandidateHandoffTask] = []
    for candidate in candidates:
        for site_id in candidate.receptor_site_ids:
            site = report.receptor_sites.get(site_id)
            if site is None:
                raise ValidationError(f"candidate refers to unknown site: {site_id}")
            reference_receptor_path = (
                config.paths.data_dir
                / "receptors"
                / "prepared"
                / f"{site.receptor_id}.cif"
            )
            if not reference_receptor_path.is_file():
                raise ValidationError(
                    f"prepared reference receptor is missing: {site.receptor_id}"
                )
            for pose in select_site_poses(
                site=site,
                pose_by_id=pose_by_id,
                count=config.validation.handoff.poses_per_receptor_site,
                peptide_sequence=config.chemistry.sequence,
                reference_receptor_path=reference_receptor_path,
                pose_clustering_rmsd_A=(
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
            ):
                _identifier(pose.pose_id, name="pose ID")
                task_id = _identifier(
                    f"{candidate.candidate_id}__{pose.pose_id}", name="handoff task ID"
                )
                tasks.append(
                    CandidateHandoffTask(
                        task_id=task_id,
                        candidate_id=candidate.candidate_id,
                        receptor_site_id=site.site_id,
                        pose=pose,
                        reference_receptor_path=reference_receptor_path,
                        reference_receptor_sha256=sha256_file(reference_receptor_path),
                    )
                )
    return tuple(tasks)


def _cg2all_parameters(config: AppConfig) -> dict[str, JsonValue]:
    settings = config.validation.cg2all
    return {
        "representation": settings.representation,
        "receptor_histidine_state": settings.receptor_histidine_state,
        "device": settings.device,
        "processes": settings.processes,
        "batch_size": settings.batch_size,
        "chain_break_cutoff_A": settings.chain_break_cutoff_A,
        "max_ca_rmsd_A": settings.max_ca_rmsd_A,
    }


def handoff_task_records(
    *,
    tasks: tuple[CandidateHandoffTask, ...],
    discovery_run_dir: Path,
    data_dir: Path,
) -> list[dict[str, JsonValue]]:
    """生成计划写入和执行复核共用的完整任务合同。"""
    discovery_root = discovery_run_dir.resolve()
    records: list[dict[str, JsonValue]] = []
    for task in tasks:
        model_path = task.pose.model_path.resolve()
        if not model_path.is_relative_to(discovery_root):
            raise ValidationError(
                f"handoff source model is outside discovery run: {model_path}"
            )
        records.append(
            {
                "task_id": task.task_id,
                "candidate_id": task.candidate_id,
                "receptor_site_id": task.receptor_site_id,
                "pose_id": task.pose.pose_id,
                "receptor_id": task.pose.receptor_id,
                "target": task.pose.target,
                "seed": task.pose.seed,
                "source_model": {
                    "path": model_path.relative_to(discovery_root).as_posix(),
                    "sha256": task.pose.model_sha256,
                    "model_index": task.pose.model_index,
                },
                "reference_receptor": {
                    "path": task.reference_receptor_path.relative_to(
                        data_dir
                    ).as_posix(),
                    "sha256": task.reference_receptor_sha256,
                },
                "status": "planned",
            }
        )
    return records


def write_handoff_plan(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    run_id: str,
    candidate_ids: tuple[str, ...],
) -> CandidateHandoffPlan:
    """冻结阶段二来源、候选选择、输入哈希和重建工具身份。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    run_dir = config.paths.outputs_dir / "validation" / "handoffs" / run_id
    if run_dir.exists():
        raise ValidationError(f"handoff run directory already exists: {run_dir}")
    tasks = build_handoff_tasks(
        config=config,
        discovery_run_dir=discovery_run_dir,
        candidate_ids=candidate_ids,
    )
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    rosetta = verify_rosetta_scripts_tool(config.validation.rosetta)
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    discovery_files = (
        *(
            discovery_run_dir / name
            for name in ("run_manifest.json", "sampling_manifest.json")
        ),
        discovery_run_dir / "site_analysis" / "analysis_manifest.json",
    )
    if any(not path.is_file() for path in discovery_files):
        raise ValidationError("discovery run manifests are incomplete")
    atomic_write_json(
        run_dir / HANDOFF_PLAN_NAME,
        {
            "schema": "vela.validation-handoff-plan/6",
            "stage": "validation_candidate_handoff",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "evidence_category": "main_discovery_handoff",
            "known_site_information_used": False,
            "chemistry_id": config.chemistry.chemistry_id,
            "software": {
                **vela_software_identity(),
                "cg2all_version": cg2all.version,
                "cg2all_executable_sha256": cg2all.executable_sha256,
                "cg2all_checkpoint_sha256": cg2all.checkpoint_sha256,
                "rosetta_version": rosetta.version,
                "rosetta_scripts_sha256": rosetta.executable_sha256,
            },
            "inputs": {
                "config_snapshot": {
                    "path": snapshot.name,
                    "sha256": sha256_file(snapshot),
                },
                "discovery_run": {
                    "path": discovery_run_dir.resolve()
                    .relative_to(config.paths.outputs_dir.resolve())
                    .as_posix(),
                    "manifests": [
                        {
                            "path": path.relative_to(discovery_run_dir).as_posix(),
                            "sha256": sha256_file(path),
                        }
                        for path in discovery_files
                    ],
                },
            },
            "selection": {
                "requested_candidate_ids": list(candidate_ids),
                "poses_per_receptor_site": (
                    config.validation.handoff.poses_per_receptor_site
                ),
                "pose_clustering_rmsd_A": (
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
                "pose_selection": (
                    "pose_cluster_seed_support_desc,population_desc,energy_asc; "
                    "cluster_medoid; distinct_seed; within_cluster_fill_if_needed"
                ),
                "distinct_seed_required": True,
                "supported_candidates_only": True,
            },
            "cg2all_parameters": _cg2all_parameters(config),
            "tasks": handoff_task_records(
                tasks=tasks,
                discovery_run_dir=discovery_run_dir,
                data_dir=config.paths.data_dir,
            ),
        },
    )
    return CandidateHandoffPlan(run_id, run_dir, discovery_run_dir, tasks)
