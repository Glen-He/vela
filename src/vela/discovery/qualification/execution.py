"""执行或恢复阶段二资格采样并冻结原始证据。"""

from __future__ import annotations

import json
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_mapping
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.analysis.pose_table import write_pose_evidence
from vela.discovery.models import DiscoveryError
from vela.discovery.qualification.evaluation import (
    audit_native_recovery_before_filtering,
    compare_poses_to_native,
)
from vela.discovery.qualification.planning import (
    CONTROL_RECOVERY,
    QualificationCase,
    build_qualification_cases,
    case_records,
)
from vela.discovery.qualification.schemas import PLAN_SCHEMA, SAMPLING_SCHEMA
from vela.discovery.sampling.cabsdock import (
    cabsdock_archive_path,
    verify_cabsdock_tool,
)
from vela.discovery.sampling.evidence import (
    CabsDockEvidence,
    CandidateSelectionSettings,
    candidate_selection_settings,
)
from vela.discovery.sampling.execution import run_cabsdock_task
from vela.discovery.sampling.planning import cabsdock_parameters


def _document(path: Path) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"qualification plan does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name="qualification plan")
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid qualification plan: {path}") from exc


def _run_seed_cases(
    *,
    cases: tuple[QualificationCase, ...],
    config: AppConfig,
    run_dir: Path,
    plan_hash: str,
    tool_version: str,
    selection: CandidateSelectionSettings,
) -> tuple[tuple[str, CabsDockEvidence], ...]:
    """一个 worker 顺序完成同一 seed 在全部控制受体上的任务。"""
    return tuple(
        (
            case.task.task_id,
            run_cabsdock_task(
                task=case.task,
                run_dir=run_dir,
                run_manifest_sha256=plan_hash,
                chemistry=case.chemistry,
                secondary_structure=case.secondary_structure,
                reference_receptor_id=case.reference_receptor_id,
                reference_path=case.reference_path,
                config=config,
                tool_version=tool_version,
                selection=selection,
            ),
        )
        for case in cases
    )


def _metrics_table(rows: list[dict[str, str]]) -> str:
    fields: tuple[str, ...] = (
        "evidence_set",
        "pose_id",
        "task_id",
        "seed",
        "model_index",
        "qc_status",
        "ligand_ca_rmsd_A",
        "ligand_centroid_distance_A",
        "native_receptor_contact_fraction",
    )
    lines = ["\t".join(fields)]
    lines.extend("\t".join(row[field] for field in fields) for row in rows)
    return "\n".join(lines) + "\n"


def run_qualification(*, config: AppConfig, run_dir: Path) -> None:
    """执行资格计划, 写出多受体位点控制的原始证据。"""
    plan_path = run_dir / "qualification_plan.json"
    plan = _document(plan_path)
    target_id = plan.get("target_id")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or not isinstance(target_id, str)
    ):
        raise DiscoveryError("qualification plan identity is invalid")
    if (run_dir / "qualification_sampling.json").exists():
        raise DiscoveryError(f"qualification sampling already exists: {run_dir}")
    snapshot = object_mapping(
        plan.get("config_snapshot"), name="qualification config snapshot"
    )
    snapshot_path = run_dir / str(snapshot.get("path"))
    if snapshot.get("sha256") != sha256_file(
        snapshot_path
    ) or config.source_snapshot_sha256 != sha256_file(snapshot_path):
        raise DiscoveryError("current config differs from the qualification plan")
    method_parameters = object_mapping(
        plan.get("method_parameters"), name="qualification method parameters"
    )
    recorded_cabsdock = object_mapping(
        method_parameters.get("cabsdock"), name="qualification CABS-dock parameters"
    )
    if recorded_cabsdock != cabsdock_parameters(config):
        raise DiscoveryError(
            "current CABS-dock parameters differ from the qualification plan"
        )
    cases = build_qualification_cases(config=config, target_id=target_id)
    records = case_records(cases=cases, data_dir=config.paths.data_dir)
    if plan.get("tasks") != records or plan.get("task_count") != len(cases):
        raise DiscoveryError("qualification tasks differ from the frozen plan")
    tool_version = verify_cabsdock_tool(config.discovery.cabsdock)
    plan_hash = sha256_file(plan_path)
    selection = candidate_selection_settings(config.discovery.cabsdock)
    grouped: dict[int, list[QualificationCase]] = {}
    for case in cases:
        grouped.setdefault(case.task.seed, []).append(case)
    batches = tuple(tuple(batch) for batch in grouped.values())
    worker_count = min(config.discovery.cabsdock.seed_workers, len(batches))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = tuple(
            executor.submit(
                _run_seed_cases,
                cases=batch,
                config=config,
                run_dir=run_dir,
                plan_hash=plan_hash,
                tool_version=tool_version,
                selection=selection,
            )
            for batch in batches
        )
        evidence_by_task = {
            task_id: evidence
            for future in futures
            for task_id, evidence in future.result()
        }
    if set(evidence_by_task) != {case.task.task_id for case in cases}:
        raise DiscoveryError("qualification workers returned incomplete evidence")
    poses: list[PoseEvidence] = []
    baseline_poses: list[PoseEvidence] = []
    metric_rows: list[dict[str, str]] = []
    task_rows: list[dict[str, JsonValue]] = []
    for case in cases:
        evidence = evidence_by_task[case.task.task_id]
        poses.extend(evidence.poses)
        baseline_poses.extend(evidence.baseline_poses)
        task_row: dict[str, JsonValue] = {
            "task_id": case.task.task_id,
            "case": case.task.evidence_category,
            "seed": case.task.seed,
            "execution_status": "completed",
            "selection_status": evidence.selection_status,
            "selection_failure_reasons": list(evidence.selection_failure_reasons),
            "sampling_model_count": evidence.sampling_model_count,
            "filtered_model_count": evidence.filtered_model_count,
            "filtered_topology_feasible_model_count": (
                evidence.filtered_topology_feasible_model_count
            ),
            "topology_feasible_model_count": (evidence.topology_feasible_model_count),
            "topology_feasible_fraction": evidence.topology_feasible_fraction,
            "contacting_topology_feasible_model_count": (
                evidence.contacting_topology_feasible_model_count
            ),
            "trajectory_model_count": (
                evidence.trajectory_audit.trajectory_model_count
            ),
            "trajectory_topology_feasible_model_count": (
                evidence.trajectory_audit.trajectory_topology_feasible_model_count
            ),
            "trajectory_topology_feasible_fraction": (
                evidence.trajectory_audit.trajectory_topology_feasible_fraction
            ),
            "topology_feasible_enrichment_ratio": (
                evidence.trajectory_audit.topology_feasible_enrichment_ratio
            ),
            "selected_model_count": len(evidence.poses),
            "selected_model_budget": (
                config.discovery.cabsdock.max_sites_per_task
                * config.discovery.cabsdock.max_pose_clusters_per_site
                * 2
            ),
            "baseline_model_count": len(evidence.baseline_poses),
        }
        if case.task.evidence_category == CONTROL_RECOVERY:
            if case.native_pair_path is None:
                raise DiscoveryError("control recovery case lacks its native pair")
            task_dir = run_dir / "tasks" / case.task.task_id
            native_filter_audit = audit_native_recovery_before_filtering(
                archive_path=cabsdock_archive_path(task_dir),
                filtered_path=task_dir / "output_pdbs" / "top1000.pdb",
                native_pair_path=case.native_pair_path,
                chemistry=case.chemistry,
                max_ligand_ca_rmsd_A=(
                    config.discovery.qualification.max_native_ligand_rmsd_A
                ),
                max_native_site_centroid_distance_A=(
                    config.discovery.qualification.max_native_site_centroid_distance_A
                ),
                max_reconstructable_disulfide_ca_distance_A=(
                    config.discovery.cabsdock.max_reconstructable_disulfide_ca_distance_A
                ),
                contact_ca_threshold_A=(
                    config.discovery.cabsdock.trajectory_contact_ca_threshold_A
                ),
                min_native_receptor_contact_fraction=(
                    config.discovery.qualification.min_native_receptor_contact_fraction
                ),
            )
            task_row["native_filter_audit"] = {
                "metric": "ligand_CA_RMSD_after_receptor_CA_alignment",
                "qualification_gate": False,
                "trajectory_model_count": (native_filter_audit.trajectory_model_count),
                "trajectory_recovered_model_count": (
                    native_filter_audit.trajectory_recovered_model_count
                ),
                "trajectory_topology_feasible_recovered_model_count": (
                    native_filter_audit.trajectory_topology_feasible_recovered_model_count
                ),
                "filtered_model_count": native_filter_audit.filtered_model_count,
                "filtered_recovered_model_count": (
                    native_filter_audit.filtered_recovered_model_count
                ),
                "filtered_topology_feasible_recovered_model_count": (
                    native_filter_audit.filtered_topology_feasible_recovered_model_count
                ),
                "trajectory_site_recovered_model_count": (
                    native_filter_audit.trajectory_site_recovered_model_count
                ),
                "trajectory_topology_feasible_site_recovered_model_count": (
                    native_filter_audit.trajectory_topology_feasible_site_recovered_model_count
                ),
                "filtered_site_recovered_model_count": (
                    native_filter_audit.filtered_site_recovered_model_count
                ),
                "filtered_topology_feasible_site_recovered_model_count": (
                    native_filter_audit.filtered_topology_feasible_site_recovered_model_count
                ),
                "trajectory_best_ligand_ca_rmsd_A": (
                    native_filter_audit.trajectory_best_ligand_ca_rmsd_A
                ),
                "trajectory_best_native_receptor_contact_fraction": (
                    native_filter_audit.trajectory_best_native_receptor_contact_fraction
                ),
                "trajectory_best_ligand_centroid_distance_A": (
                    native_filter_audit.trajectory_best_ligand_centroid_distance_A
                ),
                "filtered_best_ligand_ca_rmsd_A": (
                    native_filter_audit.filtered_best_ligand_ca_rmsd_A
                ),
                "filtered_best_native_receptor_contact_fraction": (
                    native_filter_audit.filtered_best_native_receptor_contact_fraction
                ),
                "filtered_best_ligand_centroid_distance_A": (
                    native_filter_audit.filtered_best_ligand_centroid_distance_A
                ),
            }
            comparison_sets = (
                ("selected_candidates", evidence.poses),
                ("cabsdock_top10_baseline", evidence.baseline_poses),
            )
            for evidence_set, comparison_poses in comparison_sets:
                metrics_by_pose = compare_poses_to_native(
                    poses=comparison_poses,
                    native_pair_path=case.native_pair_path,
                    peptide_sequence=case.chemistry.sequence,
                    contact_ca_threshold_A=(
                        config.discovery.cabsdock.trajectory_contact_ca_threshold_A
                    ),
                )
                for pose in comparison_poses:
                    metrics = metrics_by_pose[pose.pose_id]
                    metric_rows.append(
                        {
                            "evidence_set": evidence_set,
                            "pose_id": pose.pose_id,
                            "task_id": pose.task_id,
                            "seed": str(pose.seed),
                            "model_index": str(pose.model_index),
                            "qc_status": pose.qc_status,
                            "ligand_ca_rmsd_A": f"{metrics.ligand_ca_rmsd_A:.6f}",
                            "ligand_centroid_distance_A": (
                                f"{metrics.ligand_centroid_distance_A:.6f}"
                            ),
                            "native_receptor_contact_fraction": (
                                f"{metrics.native_receptor_contact_fraction:.6f}"
                            ),
                        }
                    )
        task_rows.append(task_row)
    pose_path = run_dir / "pose_evidence.tsv"
    baseline_path = run_dir / "baseline_pose_evidence.tsv"
    metrics_path = run_dir / "native_recovery.tsv"
    write_pose_evidence(poses=tuple(poses), path=pose_path, run_dir=run_dir)
    write_pose_evidence(
        poses=tuple(baseline_poses), path=baseline_path, run_dir=run_dir
    )
    atomic_write_text(metrics_path, _metrics_table(metric_rows))
    atomic_write_json(
        run_dir / "qualification_sampling.json",
        {
            "schema": SAMPLING_SCHEMA,
            "stage": "discovery_qualification",
            "status": "sampling_completed",
            "target_id": target_id,
            "completed_at": utc_now(),
            "qualification_plan_sha256": plan_hash,
            "software": {
                "method_version": tool_version,
                "adapter_version": version("vela"),
                "cabsdock_source_revision": (config.discovery.cabsdock.source_revision),
            },
            "execution": {
                "unit": "qualification_seed_batch",
                "seed_workers": worker_count,
                "seed_batch_count": len(batches),
                "cases_per_seed": len(batches[0]),
            },
            "pose_evidence": {
                "path": pose_path.name,
                "sha256": sha256_file(pose_path),
            },
            "baseline_pose_evidence": {
                "path": baseline_path.name,
                "sha256": sha256_file(baseline_path),
            },
            "native_recovery": {
                "path": metrics_path.name,
                "sha256": sha256_file(metrics_path),
            },
            "tasks": task_rows,
        },
    )
