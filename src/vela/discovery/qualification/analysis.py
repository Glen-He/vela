"""按受体分离评价位点访问、预算化交付与精确姿态诊断。"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    is_vela_software_identity,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.analysis.clustering import (
    analyze_sites,
    candidate_analysis_contract,
)
from vela.discovery.analysis.evidence import CandidateSite, PoseEvidence, ReceptorSite
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.models import DiscoveryError
from vela.discovery.qualification.planning import (
    CONTROL_RECOVERY,
    qualification_decision_rules,
    topology_calibration_record,
)
from vela.discovery.qualification.schemas import (
    PLAN_SCHEMA,
    REPORT_SCHEMA,
    SAMPLING_SCHEMA,
)

SELECTED_EVIDENCE = "selected_candidates"
BASELINE_EVIDENCE = "cabsdock_top10_baseline"


@dataclass(frozen=True, slots=True)
class RecoveryMetric:
    """一个候选姿态的事后 native 评价。"""

    evidence_set: str
    pose_id: str
    seed: int
    qc_status: str
    ligand_ca_rmsd_A: float
    ligand_centroid_distance_A: float
    native_receptor_contact_fraction: float


@dataclass(frozen=True, slots=True)
class TaskIdentity:
    """资格任务的冻结身份。"""

    task_id: str
    receptor_id: str
    seed: int


def _document(path: Path, *, name: str) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"{name} does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name=name)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid {name}: {path}") from exc


def _recovery_metrics(path: Path) -> tuple[RecoveryMetric, ...]:
    if not path.is_file():
        raise DiscoveryError(f"native recovery table does not exist: {path}")
    metrics: list[RecoveryMetric] = []
    identities: set[tuple[str, str]] = set()
    try:
        with path.open("r", encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                metric = RecoveryMetric(
                    evidence_set=row["evidence_set"],
                    pose_id=row["pose_id"],
                    seed=int(row["seed"]),
                    qc_status=row["qc_status"],
                    ligand_ca_rmsd_A=float(row["ligand_ca_rmsd_A"]),
                    ligand_centroid_distance_A=float(row["ligand_centroid_distance_A"]),
                    native_receptor_contact_fraction=float(
                        row["native_receptor_contact_fraction"]
                    ),
                )
                identity = (metric.evidence_set, metric.pose_id)
                if identity in identities:
                    raise DiscoveryError(
                        f"duplicate native recovery pose: {metric.evidence_set}/"
                        f"{metric.pose_id}"
                    )
                identities.add(identity)
                metrics.append(metric)
    except (KeyError, ValueError, UnicodeDecodeError) as exc:
        raise DiscoveryError(f"invalid native recovery table: {path}") from exc
    return tuple(metrics)


def _task_identities(plan: dict[str, object]) -> dict[str, TaskIdentity]:
    try:
        raw_tasks = object_list(plan.get("tasks"), name="qualification tasks")
    except TypeError as exc:
        raise DiscoveryError("qualification tasks are invalid") from exc
    identities: dict[str, TaskIdentity] = {}
    for raw in raw_tasks:
        try:
            task = object_mapping(raw, name="qualification task")
        except TypeError as exc:
            raise DiscoveryError("qualification task is invalid") from exc
        task_id = task.get("task_id")
        receptor_id = task.get("receptor_id")
        seed = task.get("seed")
        if (
            not isinstance(task_id, str)
            or not isinstance(receptor_id, str)
            or not isinstance(seed, int)
            or isinstance(seed, bool)
            or task.get("case") != CONTROL_RECOVERY
            or task_id in identities
        ):
            raise DiscoveryError("qualification task identity is invalid")
        identities[task_id] = TaskIdentity(task_id, receptor_id, seed)
    return identities


def _sampling_rows(
    *, sampling: dict[str, object], tasks: dict[str, TaskIdentity]
) -> dict[str, dict[str, object]]:
    try:
        raw_rows = object_list(sampling.get("tasks"), name="qualification sampling")
    except TypeError as exc:
        raise DiscoveryError("qualification sampling tasks are invalid") from exc
    rows: dict[str, dict[str, object]] = {}
    for raw in raw_rows:
        try:
            row = object_mapping(raw, name="qualification sampling task")
        except TypeError as exc:
            raise DiscoveryError("qualification sampling task is invalid") from exc
        task_id = row.get("task_id")
        if (
            not isinstance(task_id, str)
            or task_id not in tasks
            or task_id in rows
            or row.get("execution_status") != "completed"
        ):
            raise DiscoveryError("qualification sampling task identity is invalid")
        rows[task_id] = row
    if set(rows) != set(tasks):
        raise DiscoveryError("qualification sampling tasks do not match the plan")
    return rows


def _nonnegative_integer(value: object, *, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise DiscoveryError(f"{name} must be a non-negative integer")
    return value


def _native_filter_audit(row: dict[str, object]) -> dict[str, object]:
    try:
        audit = object_mapping(row.get("native_filter_audit"), name="native audit")
    except TypeError as exc:
        raise DiscoveryError("control task lacks its native filter audit") from exc
    if audit.get("qualification_gate") is not False:
        raise DiscoveryError("native evaluation must occur after candidate freezing")
    return audit


def _site_recovered(metric: RecoveryMetric, *, config: AppConfig) -> bool:
    rules = config.discovery.qualification
    return (
        metric.qc_status == "passed"
        and metric.ligand_centroid_distance_A
        <= rules.max_native_site_centroid_distance_A
        and metric.native_receptor_contact_fraction
        >= rules.min_native_receptor_contact_fraction
    )


def _pose_recovered(metric: RecoveryMetric, *, config: AppConfig) -> bool:
    rules = config.discovery.qualification
    return (
        metric.qc_status == "passed"
        and metric.ligand_ca_rmsd_A <= rules.max_native_ligand_rmsd_A
        and metric.native_receptor_contact_fraction
        >= rules.min_native_receptor_contact_fraction
    )


def _rank_sites(sites: tuple[ReceptorSite, ...]) -> tuple[ReceptorSite, ...]:
    return tuple(
        sorted(
            (site for site in sites if site.supported),
            key=lambda site: (
                -len(site.supporting_seeds),
                -site.pose_count,
                site.site_id,
            ),
        )
    )


def _best_values(
    metrics: tuple[RecoveryMetric, ...],
) -> tuple[float | None, float | None, float | None]:
    if not metrics:
        return None, None, None
    return (
        round(min(metric.ligand_ca_rmsd_A for metric in metrics), 6),
        round(min(metric.ligand_centroid_distance_A for metric in metrics), 6),
        round(max(metric.native_receptor_contact_fraction for metric in metrics), 6),
    )


def _receptor_evidence(
    *,
    config: AppConfig,
    receptor_id: str,
    tasks: dict[str, TaskIdentity],
    rows: dict[str, dict[str, object]],
    selected: tuple[PoseEvidence, ...],
    baseline: tuple[PoseEvidence, ...],
    metrics: tuple[RecoveryMetric, ...],
) -> dict[str, JsonValue]:
    rules = config.discovery.qualification
    receptor_task_ids = {
        task_id for task_id, task in tasks.items() if task.receptor_id == receptor_id
    }
    receptor_rows = tuple(rows[task_id] for task_id in sorted(receptor_task_ids))
    receptor_poses = tuple(
        pose for pose in selected if pose.task_id in receptor_task_ids
    )
    receptor_baseline = tuple(
        pose for pose in baseline if pose.task_id in receptor_task_ids
    )
    selected_ids = {pose.pose_id for pose in receptor_poses}
    baseline_ids = {pose.pose_id for pose in receptor_baseline}
    selected_metrics = tuple(
        metric
        for metric in metrics
        if metric.evidence_set == SELECTED_EVIDENCE and metric.pose_id in selected_ids
    )
    baseline_metrics = tuple(
        metric
        for metric in metrics
        if metric.evidence_set == BASELINE_EVIDENCE and metric.pose_id in baseline_ids
    )
    if {metric.pose_id for metric in selected_metrics} != selected_ids:
        raise DiscoveryError(f"selected metrics are incomplete for {receptor_id}")
    if {metric.pose_id for metric in baseline_metrics} != baseline_ids:
        raise DiscoveryError(f"baseline metrics are incomplete for {receptor_id}")

    trajectory_site_seeds: list[int] = []
    trajectory_site_models = 0
    for row in receptor_rows:
        audit = _native_filter_audit(row)
        count = _nonnegative_integer(
            audit.get("trajectory_topology_feasible_site_recovered_model_count"),
            name="trajectory_topology_feasible_site_recovered_model_count",
        )
        trajectory_site_models += count
        if count > 0:
            trajectory_site_seeds.append(
                _nonnegative_integer(row.get("seed"), name="seed")
            )
    trajectory_site_seeds = sorted(set(trajectory_site_seeds))
    sampling_passed = len(trajectory_site_seeds) >= rules.min_native_site_seed_support

    recovered_site_ids = {
        metric.pose_id
        for metric in selected_metrics
        if _site_recovered(metric, config=config)
    }
    exact_ids = {
        metric.pose_id
        for metric in selected_metrics
        if _pose_recovered(metric, config=config)
    }
    baseline_exact_ids = {
        metric.pose_id
        for metric in baseline_metrics
        if _pose_recovered(metric, config=config)
    }
    target = config.discovery.target(rules.control_target_id)
    site_result = analyze_sites(poses=receptor_poses, settings=target.analysis)
    ranked_sites = _rank_sites(site_result.receptor_sites)
    recovered_sites: list[dict[str, JsonValue]] = []
    qualifying_ranks: list[int] = []
    for rank, site in enumerate(ranked_sites, 1):
        recovered = recovered_site_ids & set(site.pose_ids)
        if not recovered:
            continue
        recovered_seeds = sorted(
            {pose.seed for pose in receptor_poses if pose.pose_id in recovered}
        )
        recovered_sites.append(
            {
                "rank": rank,
                "site_id": site.site_id,
                "supporting_seeds": list(site.supporting_seeds),
                "recovered_seeds": recovered_seeds,
                "recovered_seed_count": len(recovered_seeds),
                "recovered_model_count": len(recovered),
            }
        )
        if len(recovered_seeds) >= rules.min_native_site_seed_support:
            qualifying_ranks.append(rank)
    first_rank = min(qualifying_ranks, default=None)
    diagnostic_delivery_passed = (
        first_rank is not None and first_rank <= rules.receptor_site_diagnostic_budget
    )
    selected_best = _best_values(selected_metrics)
    baseline_best = _best_values(baseline_metrics)
    return {
        "receptor_id": receptor_id,
        "sampling_pool_site_access": {
            "passed": sampling_passed,
            "topology_feasible_site_recovered_model_count": trajectory_site_models,
            "successful_seeds": trajectory_site_seeds,
            "successful_seed_count": len(trajectory_site_seeds),
        },
        "receptor_site_diagnostic": {
            "passed": diagnostic_delivery_passed,
            "qualification_gate": False,
            "ranking": "supporting_seed_count_desc,pose_count_desc,site_id_asc",
            "site_budget": rules.receptor_site_diagnostic_budget,
            "eligible_supported_site_count": len(ranked_sites),
            "first_qualifying_native_site_rank": first_rank,
            "native_sites": recovered_sites,
            "maximum_handoff_start_count": (
                rules.receptor_site_diagnostic_budget
                * config.validation.handoff.poses_per_receptor_site
            ),
        },
        "exact_pose_recovery": {
            "qualification_gate": False,
            "selected_model_count": len(exact_ids),
            "selected_successful_seeds": sorted(
                {pose.seed for pose in receptor_poses if pose.pose_id in exact_ids}
            ),
            "baseline_model_count": len(baseline_exact_ids),
            "baseline_successful_seeds": sorted(
                {
                    pose.seed
                    for pose in receptor_baseline
                    if pose.pose_id in baseline_exact_ids
                }
            ),
            "selected_best_ligand_ca_rmsd_A": selected_best[0],
            "selected_best_centroid_distance_A": selected_best[1],
            "selected_best_contact_fraction": selected_best[2],
            "baseline_best_ligand_ca_rmsd_A": baseline_best[0],
            "baseline_best_centroid_distance_A": baseline_best[1],
            "baseline_best_contact_fraction": baseline_best[2],
        },
        "candidate_count": len(receptor_poses),
    }


def _candidate_delivery_evidence(
    *,
    config: AppConfig,
    selected: tuple[PoseEvidence, ...],
    metrics: tuple[RecoveryMetric, ...],
) -> dict[str, JsonValue]:
    rules = config.discovery.qualification
    target = config.discovery.target(rules.control_target_id)
    result = analyze_sites(poses=selected, settings=target.analysis)
    pose_by_id = {pose.pose_id: pose for pose in selected}
    site_by_id = {site.site_id: site for site in result.receptor_sites}
    recovered_ids = {
        metric.pose_id
        for metric in metrics
        if metric.evidence_set == SELECTED_EVIDENCE
        and _site_recovered(metric, config=config)
    }
    recovered_candidates: list[dict[str, JsonValue]] = []
    delivered_candidates: list[CandidateSite] = []
    for candidate in sorted(
        result.candidate_sites,
        key=lambda item: (
            item.evidence_tier,
            item.rank_within_tier,
            item.candidate_id,
        ),
    ):
        receptor_recovery: list[dict[str, JsonValue]] = []
        recovered_receptor_count = 0
        for site_id in candidate.receptor_site_ids:
            site = site_by_id[site_id]
            site_recovered_ids = recovered_ids & set(site.pose_ids)
            recovered_seeds = sorted(
                {pose_by_id[pose_id].seed for pose_id in site_recovered_ids}
            )
            receptor_passed = len(recovered_seeds) >= rules.min_native_site_seed_support
            recovered_receptor_count += receptor_passed
            receptor_recovery.append(
                {
                    "receptor_id": site.receptor_id,
                    "site_id": site.site_id,
                    "recovered_model_count": len(site_recovered_ids),
                    "recovered_seeds": recovered_seeds,
                    "recovered_seed_count": len(recovered_seeds),
                    "passed": receptor_passed,
                }
            )
        if recovered_receptor_count < 1:
            continue
        native_supported = recovered_receptor_count >= rules.min_native_receptor_support
        formally_delivered = native_supported and candidate.handoff_eligible
        if formally_delivered:
            delivered_candidates.append(candidate)
        recovered_candidates.append(
            {
                "candidate_id": candidate.candidate_id,
                "evidence_tier": candidate.evidence_tier,
                "rank_within_tier": candidate.rank_within_tier,
                "handoff_eligible": candidate.handoff_eligible,
                "native_supported": native_supported,
                "formally_delivered": formally_delivered,
                "recovered_receptor_count": recovered_receptor_count,
                "minimum_seed_support": candidate.minimum_seed_support,
                "total_seed_support": candidate.total_seed_support,
                "maximum_normalized_site_distance": round(
                    candidate.maximum_normalized_site_distance, 6
                ),
                "receptor_recovery": receptor_recovery,
            }
        )
    tier_counts = {
        tier: sum(
            candidate.evidence_tier == tier for candidate in result.candidate_sites
        )
        for tier in (
            "ensemble_consensus",
            "conformation_specific",
            "insufficient_evidence",
        )
    }
    eligible_candidates = tuple(
        candidate for candidate in result.candidate_sites if candidate.handoff_eligible
    )
    start_count = sum(
        candidate.receptor_support * config.validation.handoff.poses_per_receptor_site
        for candidate in eligible_candidates
    )
    return {
        "passed": bool(delivered_candidates),
        "qualification_gate": True,
        "analysis_contract": candidate_analysis_contract(target.analysis),
        "ensemble_candidate_budget": target.analysis.ensemble_candidate_budget,
        "conformation_specific_candidate_budget": (
            target.analysis.conformation_specific_candidate_budget
        ),
        "candidate_count": len(result.candidate_sites),
        "evidence_tier_counts": tier_counts,
        "handoff_eligible_candidate_count": len(eligible_candidates),
        "maximum_handoff_start_count": start_count,
        "maximum_flexpepdock_decoy_count": (
            start_count
            * len(config.validation.seeds)
            * config.validation.rosetta.decoys_per_seed
        ),
        "native_candidates": recovered_candidates,
        "formally_delivered_candidate_ids": [
            candidate.candidate_id for candidate in delivered_candidates
        ],
    }


def analyze_qualification(*, config: AppConfig, run_dir: Path) -> Path:
    """以冻结的跨受体候选分级与预算合同生成 Stage 2 资格结论。"""
    report_path = run_dir / "qualification_report.json"
    if report_path.exists():
        raise DiscoveryError(f"qualification report already exists: {report_path}")
    plan_path = run_dir / "qualification_plan.json"
    sampling_path = run_dir / "qualification_sampling.json"
    plan = _document(plan_path, name="qualification plan")
    sampling = _document(sampling_path, name="qualification sampling")
    target_id = plan.get("target_id")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or not is_vela_software_identity(plan.get("software"))
        or sampling.get("schema") != SAMPLING_SCHEMA
        or sampling.get("qualification_plan_sha256") != sha256_file(plan_path)
        or sampling.get("target_id") != target_id
        or not isinstance(target_id, str)
    ):
        raise DiscoveryError("qualification evidence identity is invalid")
    if target_id != config.discovery.qualification.control_target_id:
        raise DiscoveryError("qualification target lacks a target-matched control")
    match plan.get("software"):
        case {
            "vela_version": str(sampling_version),
            "vela_source_sha256": str(sampling_source_sha256),
        }:
            sampling_software: dict[str, JsonValue] = {
                "vela_version": sampling_version,
                "vela_source_sha256": sampling_source_sha256,
            }
        case _:
            raise DiscoveryError("qualification sampling software is invalid")
    if plan.get("decision_rules") != qualification_decision_rules(
        config=config, target_id=target_id
    ):
        raise DiscoveryError("qualification decision rules have changed")
    tasks = _task_identities(plan)
    rows = _sampling_rows(sampling=sampling, tasks=tasks)
    selected = read_pose_evidence(path=run_dir / "pose_evidence.tsv", run_dir=run_dir)
    baseline = read_pose_evidence(
        path=run_dir / "baseline_pose_evidence.tsv", run_dir=run_dir
    )
    metrics = _recovery_metrics(run_dir / "native_recovery.tsv")
    known_sets = {SELECTED_EVIDENCE, BASELINE_EVIDENCE}
    if {metric.evidence_set for metric in metrics} - known_sets:
        raise DiscoveryError("native recovery contains an unknown evidence set")
    all_selected_ids = {pose.pose_id for pose in selected}
    all_baseline_ids = {pose.pose_id for pose in baseline}
    if {
        metric.pose_id for metric in metrics if metric.evidence_set == SELECTED_EVIDENCE
    } != all_selected_ids or {
        metric.pose_id for metric in metrics if metric.evidence_set == BASELINE_EVIDENCE
    } != all_baseline_ids:
        raise DiscoveryError("native recovery metrics do not match frozen candidates")

    rules = config.discovery.qualification
    scope = object_mapping(plan.get("control_scope"), name="control scope")
    expected_scope: dict[str, JsonValue] = {
        "control_target_id": rules.control_target_id,
        "control_receptor_ids": list(rules.control_receptor_ids),
        "benchmark_receptor_id": rules.benchmark_receptor_id,
        "reference_receptor_id": config.discovery.target(target_id).reference_receptor,
        "native_bound_state_id": rules.control_bound_state_id,
        "requested_target_id": target_id,
        "protein_target_matched": True,
        "receptor_relation": "cross_receptor_target_domain",
    }
    if scope != expected_scope:
        raise DiscoveryError("qualification control scope has changed")
    receptor_evidence = tuple(
        _receptor_evidence(
            config=config,
            receptor_id=receptor_id,
            tasks=tasks,
            rows=rows,
            selected=selected,
            baseline=baseline,
            metrics=metrics,
        )
        for receptor_id in rules.control_receptor_ids
    )
    sampling_successful_receptors = tuple(
        str(evidence["receptor_id"])
        for evidence in receptor_evidence
        if object_mapping(
            evidence["sampling_pool_site_access"], name="sampling pool site access"
        ).get("passed")
        is True
    )
    candidate_delivery = _candidate_delivery_evidence(
        config=config,
        selected=selected,
        metrics=metrics,
    )
    topology = topology_calibration_record(config)
    if plan.get("topology_calibration") != topology:
        raise DiscoveryError("topology calibration differs from the frozen plan")
    failure_reasons: list[str] = []
    if topology.get("status") != "qualified":
        failure_reasons.append("all_atom_topology_recoverability_not_qualified")
    if len(sampling_successful_receptors) < rules.min_native_receptor_support:
        failure_reasons.append("insufficient_receptor_sampling_support")
    if candidate_delivery["passed"] is not True:
        failure_reasons.append("cross_receptor_candidate_delivery_failed")
    status = "unqualified" if failure_reasons else "qualified"
    atomic_write_json(
        report_path,
        {
            "schema": REPORT_SCHEMA,
            "stage": "discovery_qualification",
            "status": status,
            "evidence_level": (
                "site_hypothesis_ready" if status == "qualified" else "not_qualified"
            ),
            "atomic_pose_design_ready": False,
            "failure_reasons": failure_reasons,
            "target_id": target_id,
            "generated_at": utc_now(),
            "sampling_software": sampling_software,
            "analysis_software": vela_software_identity(),
            "qualification_plan": {
                "path": plan_path.name,
                "sha256": sha256_file(plan_path),
            },
            "qualification_sampling": {
                "path": sampling_path.name,
                "sha256": sha256_file(sampling_path),
            },
            "control_scope": expected_scope,
            "topology_calibration": topology,
            "verified_decision_rules": qualification_decision_rules(
                config=config, target_id=target_id
            ),
            "receptor_evidence": list(receptor_evidence),
            "candidate_delivery": candidate_delivery,
            "sampling_successful_receptors": list(sampling_successful_receptors),
            "sampling_successful_receptor_count": len(sampling_successful_receptors),
            "ensemble_robust": (
                len(sampling_successful_receptors) == len(rules.control_receptor_ids)
                and candidate_delivery["passed"] is True
            ),
            "interpretation": (
                "qualification supports native-free cross-receptor candidate "
                "delivery within frozen evidence-tier budgets; exact-pose recovery "
                "and per-receptor ranks remain descriptive"
            ),
            "recommended_target_config": (
                {
                    "qualification_status": "qualified",
                    "qualification_report": report_path.as_posix(),
                    "contact_jaccard_distance": (
                        config.discovery.target(
                            target_id
                        ).analysis.contact_jaccard_distance
                    ),
                    "position_distance_A": (
                        config.discovery.target(target_id).analysis.position_distance_A
                    ),
                    "min_seed_support": (
                        config.discovery.target(target_id).analysis.min_seed_support
                    ),
                    "min_receptor_support": (
                        config.discovery.target(target_id).analysis.min_receptor_support
                    ),
                    "min_conformation_specific_seed_support": (
                        config.discovery.target(
                            target_id
                        ).analysis.min_conformation_specific_seed_support
                    ),
                    "ensemble_candidate_budget": (
                        config.discovery.target(
                            target_id
                        ).analysis.ensemble_candidate_budget
                    ),
                    "conformation_specific_candidate_budget": (
                        config.discovery.target(
                            target_id
                        ).analysis.conformation_specific_candidate_budget
                    ),
                }
                if status == "qualified"
                else None
            ),
        },
    )
    return report_path
