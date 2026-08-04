"""分离采样、候选选择与位点聚类, 生成阶段二资格结论。"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    atomic_write_json,
    is_current_vela_software,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.analysis.cluster_engine import normalized_site_distance
from vela.discovery.analysis.clustering import analyze_sites
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.models import DiscoveryError
from vela.discovery.qualification.planning import (
    CONTROL_RECOVERY,
    TARGET_PILOT,
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
    """一个阳性控制模型的实验回收指标。"""

    evidence_set: str
    pose_id: str
    seed: int
    qc_status: str
    ligand_ca_rmsd_A: float
    ligand_centroid_distance_A: float
    native_receptor_contact_fraction: float


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
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
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


def _task_cases(plan: dict[str, object]) -> dict[str, str]:
    try:
        rows = object_list(plan.get("tasks"), name="qualification tasks")
    except TypeError as exc:
        raise DiscoveryError("qualification tasks are invalid") from exc
    cases: dict[str, str] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="qualification task")
        except TypeError as exc:
            raise DiscoveryError("qualification task is invalid") from exc
        task_id = row.get("task_id")
        case = row.get("case")
        if (
            not isinstance(task_id, str)
            or not isinstance(case, str)
            or case not in {CONTROL_RECOVERY, TARGET_PILOT}
        ):
            raise DiscoveryError("qualification task identity is invalid")
        cases[task_id] = case
    return cases


def _shared_control_paths(*, config: AppConfig, raw: object) -> dict[str, Path] | None:
    if raw is None:
        return None
    try:
        shared = object_mapping(raw, name="shared control")
        files = object_mapping(shared.get("files"), name="shared control files")
    except TypeError as exc:
        raise DiscoveryError("shared control record is invalid") from exc
    relative_run = shared.get("run_path")
    if not isinstance(relative_run, str) or not relative_run:
        raise DiscoveryError("shared control run path is invalid")
    root = (config.paths.outputs_dir / "discovery" / "qualifications").resolve()
    run_dir = (root / relative_run).resolve()
    if not run_dir.is_relative_to(root):
        raise DiscoveryError("shared control run escapes qualification outputs")
    paths: dict[str, Path] = {}
    required = {
        "qualification_plan": "qualification_plan.json",
        "qualification_sampling": "qualification_sampling.json",
        "pose_evidence": "pose_evidence.tsv",
        "baseline_pose_evidence": "baseline_pose_evidence.tsv",
        "native_recovery": "native_recovery.tsv",
        "qualification_report": "qualification_report.json",
    }
    for name, filename in required.items():
        try:
            record = object_mapping(files.get(name), name=f"shared control {name}")
        except TypeError as exc:
            raise DiscoveryError(f"shared control {name} record is invalid") from exc
        if record.get("path") != filename:
            raise DiscoveryError(f"shared control {name} path is invalid")
        digest = record.get("sha256")
        path = run_dir / filename
        if not isinstance(digest, str) or not path.is_file():
            raise DiscoveryError(f"shared control {name} is missing")
        if sha256_file(path) != digest:
            raise DiscoveryError(f"shared control {name} hash mismatch")
        paths[name] = path
    report = _document(
        paths["qualification_report"], name="shared qualification report"
    )
    plan = _document(paths["qualification_plan"], name="shared qualification plan")
    control = object_mapping(
        report.get("control_recovery"), name="shared control recovery"
    )
    selection = object_mapping(
        control.get("candidate_selection"), name="shared candidate selection"
    )
    if (
        report.get("schema") != REPORT_SCHEMA
        or report.get("status") != "qualified"
        or selection.get("passed") is not True
        or not is_current_vela_software(plan.get("software"))
    ):
        raise DiscoveryError("shared qualification control did not pass")
    return paths


def best_native_coherent_site_evidence(
    *,
    recovered: set[str],
    poses: dict[str, PoseEvidence],
    contact_limit: float,
    position_limit: float,
) -> tuple[int, float]:
    """求跨 seed 的两两相容 native 集合, 避免星形邻域虚增支持度。"""
    unknown = recovered - set(poses)
    if unknown:
        raise DiscoveryError("native recovery references unknown poses")
    if not recovered:
        return 0, 0.0

    compatibility: dict[tuple[str, str], bool] = {}

    def compatible(first_id: str, second_id: str) -> bool:
        if first_id == second_id:
            return True
        key = (first_id, second_id) if first_id < second_id else (second_id, first_id)
        cached = compatibility.get(key)
        if cached is not None:
            return cached
        first = poses[first_id]
        second = poses[second_id]
        result = (
            normalized_site_distance(
                first_contacts=first.contact_residues,
                first_position=first.local_position,
                second_contacts=second.contact_residues,
                second_position=second.local_position,
                contact_limit=contact_limit,
                position_limit=position_limit,
            )
            <= 1.0
        )
        compatibility[key] = result
        return result

    by_seed: dict[int, tuple[str, ...]] = {}
    for seed in sorted({poses[pose_id].seed for pose_id in recovered}):
        by_seed[seed] = tuple(
            sorted(pose_id for pose_id in recovered if poses[pose_id].seed == seed)
        )

    all_pose_ids = tuple(sorted(poses))
    seeds = tuple(sorted(by_seed))
    for support in range(len(seeds), 0, -1):
        for seed_subset in combinations(seeds, support):
            ordered_seeds = tuple(
                sorted(seed_subset, key=lambda seed: len(by_seed[seed]))
            )

            def search(
                index: int,
                chosen: tuple[str, ...],
                seed_order: tuple[int, ...] = ordered_seeds,
            ) -> tuple[str, ...] | None:
                if index == len(seed_order):
                    return chosen
                for pose_id in by_seed[seed_order[index]]:
                    extended = (*chosen, pose_id)
                    if not all(compatible(member, pose_id) for member in chosen):
                        continue
                    if any(
                        not any(
                            all(compatible(member, option) for member in extended)
                            for option in by_seed[remaining_seed]
                        )
                        for remaining_seed in seed_order[index + 1 :]
                    ):
                        continue
                    result = search(index + 1, extended, seed_order)
                    if result is not None:
                        return result
                return None

            chosen = search(0, ())
            if chosen is None:
                continue
            neighborhood = tuple(
                pose_id
                for pose_id in all_pose_ids
                if all(compatible(member, pose_id) for member in chosen)
            )
            return support, len(chosen) / len(neighborhood)
    return 0, 0.0


def _sampling_task_rows(
    *, sampling: dict[str, object], task_cases: dict[str, str]
) -> tuple[dict[str, object], ...]:
    try:
        raw_rows = object_list(sampling.get("tasks"), name="qualification sampling")
    except TypeError as exc:
        raise DiscoveryError("qualification sampling tasks are invalid") from exc
    rows: list[dict[str, object]] = []
    for raw in raw_rows:
        try:
            row = object_mapping(raw, name="qualification sampling task")
        except TypeError as exc:
            raise DiscoveryError("qualification sampling task is invalid") from exc
        task_id = row.get("task_id")
        if not isinstance(task_id, str) or task_id not in task_cases:
            raise DiscoveryError("qualification sampling task identity is invalid")
        if row.get("execution_status") != "completed":
            raise DiscoveryError("qualification sampling contains incomplete execution")
        rows.append(row)
    if {str(row["task_id"]) for row in rows} != set(task_cases):
        raise DiscoveryError("qualification sampling tasks do not match the plan")
    return tuple(rows)


def _nonnegative_integer(value: object, *, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise DiscoveryError(f"{name} must be a non-negative integer")
    return value


def _optional_nonnegative_number(value: object, *, name: str) -> float | None:
    if value is None:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool) or value < 0:
        raise DiscoveryError(f"{name} must be a non-negative number or null")
    return float(value)


def _pose_is_recovered(metric: RecoveryMetric, *, config: AppConfig) -> bool:
    rules = config.discovery.qualification
    return (
        metric.qc_status == "passed"
        and metric.ligand_ca_rmsd_A <= rules.max_native_ligand_rmsd_A
        and metric.native_receptor_contact_fraction
        >= rules.min_native_receptor_contact_fraction
    )


def _site_is_recovered(metric: RecoveryMetric, *, config: AppConfig) -> bool:
    rules = config.discovery.qualification
    return (
        metric.qc_status == "passed"
        and metric.ligand_centroid_distance_A
        <= rules.max_native_site_centroid_distance_A
        and metric.native_receptor_contact_fraction
        >= rules.min_native_receptor_contact_fraction
    )


def _best_metric_values(
    metrics: tuple[RecoveryMetric, ...],
) -> tuple[float | None, float | None]:
    if not metrics:
        return None, None
    return (
        round(min(metric.ligand_ca_rmsd_A for metric in metrics), 6),
        round(max(metric.native_receptor_contact_fraction for metric in metrics), 6),
    )


def analyze_qualification(*, config: AppConfig, run_dir: Path) -> Path:
    """生成分层资格结论和可直接冻结到目标配置的聚类参数。"""
    report_path = run_dir / "qualification_report.json"
    if report_path.exists():
        raise DiscoveryError(f"qualification report already exists: {report_path}")
    plan_path = run_dir / "qualification_plan.json"
    sampling_path = run_dir / "qualification_sampling.json"
    plan = _document(plan_path, name="qualification plan")
    sampling = _document(sampling_path, name="qualification sampling manifest")
    target_id = plan.get("target_id")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or not is_current_vela_software(plan.get("software"))
        or sampling.get("schema") != SAMPLING_SCHEMA
        or sampling.get("qualification_plan_sha256") != sha256_file(plan_path)
        or sampling.get("target_id") != target_id
        or not isinstance(target_id, str)
    ):
        raise DiscoveryError("qualification evidence identity is invalid")
    task_cases = _task_cases(plan)
    _sampling_task_rows(sampling=sampling, task_cases=task_cases)
    poses = read_pose_evidence(path=run_dir / "pose_evidence.tsv", run_dir=run_dir)
    baseline = read_pose_evidence(
        path=run_dir / "baseline_pose_evidence.tsv", run_dir=run_dir
    )
    shared_paths = _shared_control_paths(config=config, raw=plan.get("shared_control"))
    if shared_paths is None:
        control_plan = plan
        control_sampling = sampling
        control_cases = task_cases
        control_run_dir = run_dir
        control_poses_source = poses
        baseline_source = baseline
        metrics_path = run_dir / "native_recovery.tsv"
    else:
        control_plan = _document(
            shared_paths["qualification_plan"], name="shared qualification plan"
        )
        control_sampling = _document(
            shared_paths["qualification_sampling"],
            name="shared qualification sampling",
        )
        control_cases = _task_cases(control_plan)
        _sampling_task_rows(sampling=control_sampling, task_cases=control_cases)
        control_run_dir = shared_paths["qualification_plan"].parent
        control_poses_source = read_pose_evidence(
            path=shared_paths["pose_evidence"], run_dir=control_run_dir
        )
        baseline_source = read_pose_evidence(
            path=shared_paths["baseline_pose_evidence"], run_dir=control_run_dir
        )
        metrics_path = shared_paths["native_recovery"]
    control_poses = tuple(
        pose
        for pose in control_poses_source
        if control_cases.get(pose.task_id) == CONTROL_RECOVERY
    )
    control_baseline = tuple(
        pose
        for pose in baseline_source
        if control_cases.get(pose.task_id) == CONTROL_RECOVERY
    )
    pilot_poses = tuple(
        pose for pose in poses if task_cases.get(pose.task_id) == TARGET_PILOT
    )
    metrics = _recovery_metrics(metrics_path)
    selected_metrics = tuple(
        metric for metric in metrics if metric.evidence_set == SELECTED_EVIDENCE
    )
    baseline_metrics = tuple(
        metric for metric in metrics if metric.evidence_set == BASELINE_EVIDENCE
    )
    unknown_sets = {metric.evidence_set for metric in metrics} - {
        SELECTED_EVIDENCE,
        BASELINE_EVIDENCE,
    }
    if unknown_sets:
        raise DiscoveryError(
            "native recovery contains unknown evidence sets: "
            + ", ".join(sorted(unknown_sets))
        )
    selected_by_pose = {metric.pose_id: metric for metric in selected_metrics}
    baseline_by_pose = {metric.pose_id: metric for metric in baseline_metrics}
    control_pose_ids = {pose.pose_id for pose in control_poses}
    baseline_pose_ids = {pose.pose_id for pose in control_baseline}
    if control_pose_ids != set(selected_by_pose):
        raise DiscoveryError("selected recovery metrics do not match control poses")
    if baseline_pose_ids != set(baseline_by_pose):
        raise DiscoveryError("baseline recovery metrics do not match baseline poses")

    selection_pose_recovered = {
        metric.pose_id
        for metric in selected_metrics
        if _pose_is_recovered(metric, config=config)
    }
    selection_site_recovered = {
        metric.pose_id
        for metric in selected_metrics
        if _site_is_recovered(metric, config=config)
    }
    baseline_recovered = {
        metric.pose_id
        for metric in baseline_metrics
        if _pose_is_recovered(metric, config=config)
    }
    selection_seeds = sorted(
        {selected_by_pose[pose_id].seed for pose_id in selection_site_recovered}
    )
    selection_pose_seeds = sorted(
        {selected_by_pose[pose_id].seed for pose_id in selection_pose_recovered}
    )
    baseline_seeds = sorted(
        {baseline_by_pose[pose_id].seed for pose_id in baseline_recovered}
    )
    rules = config.discovery.qualification
    decision_rules = qualification_decision_rules(config=config, target_id=target_id)
    if plan.get("decision_rules") != decision_rules:
        raise DiscoveryError(
            "qualification decision rules differ from the frozen current protocol"
        )
    control_rows = tuple(
        row
        for row in _sampling_task_rows(
            sampling=control_sampling, task_cases=control_cases
        )
        if control_cases[str(row["task_id"])] == CONTROL_RECOVERY
    )
    trajectory_count = sum(
        _nonnegative_integer(
            row.get("trajectory_model_count"), name="trajectory_model_count"
        )
        for row in control_rows
    )
    trajectory_topology_count = sum(
        _nonnegative_integer(
            row.get("trajectory_topology_feasible_model_count"),
            name="trajectory_topology_feasible_model_count",
        )
        for row in control_rows
    )
    filtered_count = sum(
        _nonnegative_integer(
            row.get("filtered_model_count"), name="filtered_model_count"
        )
        for row in control_rows
    )
    sampling_count = sum(
        _nonnegative_integer(
            row.get("sampling_model_count"), name="sampling_model_count"
        )
        for row in control_rows
    )
    if sampling_count != trajectory_count:
        raise DiscoveryError("sampling pool count differs from the complete TRAF")
    filtered_topology_count = sum(
        _nonnegative_integer(
            row.get("filtered_topology_feasible_model_count"),
            name="filtered_topology_feasible_model_count",
        )
        for row in control_rows
    )
    topology_count = sum(
        _nonnegative_integer(
            row.get("topology_feasible_model_count"),
            name="topology_feasible_model_count",
        )
        for row in control_rows
    )
    if topology_count != trajectory_topology_count:
        raise DiscoveryError("sampling topology count differs from TRAF audit")
    pool_count = sum(
        _nonnegative_integer(
            row.get("contacting_topology_feasible_model_count"),
            name="contacting_topology_feasible_model_count",
        )
        for row in control_rows
    )
    native_filter_rows: list[tuple[dict[str, object], dict[str, object]]] = []
    selected_model_budget = (
        config.discovery.cabsdock.max_sites_per_task
        * config.discovery.cabsdock.max_pose_clusters_per_site
        * 2
    )
    for row in control_rows:
        try:
            audit = object_mapping(
                row.get("native_filter_audit"), name="native filter audit"
            )
        except TypeError as exc:
            raise DiscoveryError("control task lacks its native filter audit") from exc
        if audit.get("qualification_gate") is not False:
            raise DiscoveryError("native filter audit must remain descriptive")
        selected_count = _nonnegative_integer(
            row.get("selected_model_count"), name="selected_model_count"
        )
        if (
            row.get("selected_model_budget") != selected_model_budget
            or selected_count > selected_model_budget
        ):
            raise DiscoveryError(
                "control candidate selection exceeded its frozen budget"
            )
        native_filter_rows.append((row, audit))
    full_native_recovered = sum(
        _nonnegative_integer(
            audit.get("trajectory_recovered_model_count"),
            name="trajectory_recovered_model_count",
        )
        for _, audit in native_filter_rows
    )
    full_topology_native_recovered = sum(
        _nonnegative_integer(
            audit.get("trajectory_topology_feasible_recovered_model_count"),
            name="trajectory_topology_feasible_recovered_model_count",
        )
        for _, audit in native_filter_rows
    )
    filtered_native_recovered = sum(
        _nonnegative_integer(
            audit.get("filtered_recovered_model_count"),
            name="filtered_recovered_model_count",
        )
        for _, audit in native_filter_rows
    )
    filtered_topology_native_recovered = sum(
        _nonnegative_integer(
            audit.get("filtered_topology_feasible_recovered_model_count"),
            name="filtered_topology_feasible_recovered_model_count",
        )
        for _, audit in native_filter_rows
    )
    pool_recovered_count = sum(
        _nonnegative_integer(
            audit.get("trajectory_qualified_recovered_model_count"),
            name="trajectory_qualified_recovered_model_count",
        )
        for _, audit in native_filter_rows
    )
    pool_seeds = sorted(
        _nonnegative_integer(row.get("seed"), name="seed")
        for row, audit in native_filter_rows
        if _nonnegative_integer(
            audit.get("trajectory_qualified_recovered_model_count"),
            name="trajectory_qualified_recovered_model_count",
        )
        > 0
    )
    pool_rmsd_values_list: list[float] = []
    pool_contact_values_list: list[float] = []
    for _, audit in native_filter_rows:
        rmsd = _optional_nonnegative_number(
            audit.get("trajectory_best_ligand_ca_rmsd_A"),
            name="trajectory_best_ligand_ca_rmsd_A",
        )
        contact = _optional_nonnegative_number(
            audit.get("trajectory_best_native_receptor_contact_fraction"),
            name="trajectory_best_native_receptor_contact_fraction",
        )
        if rmsd is not None:
            pool_rmsd_values_list.append(rmsd)
        if contact is not None:
            pool_contact_values_list.append(contact)
    pool_rmsd_values = tuple(pool_rmsd_values_list)
    pool_contact_values = tuple(pool_contact_values_list)

    target = config.discovery.target(target_id)
    control_target = config.discovery.target(rules.control_target_id)
    if not target.analysis.complete or not control_target.analysis.complete:
        raise DiscoveryError(
            "qualification site analysis parameters must be frozen before execution"
        )
    pilot_sites = analyze_sites(poses=pilot_poses, settings=target.analysis)
    control_pose_by_id = {pose.pose_id: pose for pose in control_poses}
    contact_limit = control_target.analysis.contact_jaccard_distance
    position_limit = control_target.analysis.position_distance_A
    if contact_limit is None or position_limit is None:
        raise DiscoveryError("control site distance thresholds are unresolved")
    native_neighborhood_seed_support, native_neighborhood_precision = (
        best_native_coherent_site_evidence(
            recovered=selection_site_recovered,
            poses=control_pose_by_id,
            contact_limit=contact_limit,
            position_limit=position_limit,
        )
    )
    pilot_supported = tuple(
        site for site in pilot_sites.receptor_sites if site.supported
    )
    pilot_max_seed_support = max(
        (len(site.supporting_seeds) for site in pilot_supported), default=0
    )
    if not set(selection_seeds).issubset(pool_seeds):
        raise DiscoveryError("selected native seeds are absent from the sampling pool")
    selection_seed_recall_fraction = (
        len(selection_seeds) / len(pool_seeds) if pool_seeds else 0.0
    )
    sampling_passed = len(pool_seeds) >= rules.min_native_sampling_seed_support
    selection_passed = (
        bool(control_poses)
        and selection_seed_recall_fraction
        >= rules.min_selection_native_seed_recall_fraction
    )
    site_validation_passed = (
        native_neighborhood_seed_support >= rules.min_native_site_seed_support
    )
    topology_calibration = topology_calibration_record(config)
    if plan.get("topology_calibration") != topology_calibration:
        raise DiscoveryError(
            "topology calibration evidence differs from the frozen plan"
        )
    topology_calibrated = topology_calibration.get("status") == "qualified"
    control_scope = object_mapping(plan.get("control_scope"), name="control scope")
    expected_target_matched = target_id == rules.control_target_id
    scope_valid = (
        control_scope.get("control_target_id") == rules.control_target_id
        and control_scope.get("control_receptor_id") == rules.control_receptor_id
        and control_scope.get("native_bound_state_id") == rules.control_bound_state_id
        and control_scope.get("requested_target_id") == target_id
        and control_scope.get("target_matched") is expected_target_matched
    )
    if not scope_valid:
        raise DiscoveryError("qualification control scope differs from the frozen plan")
    failure_reasons: list[str] = []
    if not topology_calibrated:
        failure_reasons.append("all_atom_topology_recoverability_not_qualified")
    if not sampling_passed:
        failure_reasons.append("insufficient_native_recovery_in_sampling_pool")
    if not selection_passed:
        failure_reasons.append(
            "insufficient_native_site_recovery_in_selected_candidates"
        )
    if not control_poses:
        failure_reasons.append("control_candidate_selection_empty")
    if not site_validation_passed:
        failure_reasons.append("native_recovery_lacks_cross_seed_site_support")
    limitations: list[str] = []
    if failure_reasons:
        status = "unqualified"
    elif target_id != rules.control_target_id:
        status = "transferability_unresolved"
        limitations.append("no_target_matched_standard_cyclic_peptide_control")
    else:
        status = "qualified"
    pool_best_rmsd = round(min(pool_rmsd_values), 6) if pool_rmsd_values else None
    pool_best_contact = (
        round(max(pool_contact_values), 6) if pool_contact_values else None
    )
    selected_best_rmsd, selected_best_contact = _best_metric_values(selected_metrics)
    baseline_best_rmsd, baseline_best_contact = _best_metric_values(baseline_metrics)
    atomic_write_json(
        report_path,
        {
            "schema": REPORT_SCHEMA,
            "stage": "discovery_qualification",
            "status": status,
            "failure_reasons": failure_reasons,
            "limitations": limitations,
            "target_id": target_id,
            "generated_at": utc_now(),
            "qualification_plan": {
                "path": plan_path.name,
                "sha256": sha256_file(plan_path),
            },
            "qualification_sampling": {
                "path": sampling_path.name,
                "sha256": sha256_file(sampling_path),
            },
            "shared_control": (
                {
                    "used": True,
                    "qualification_report": {
                        "path": shared_paths["qualification_report"].as_posix(),
                        "sha256": sha256_file(shared_paths["qualification_report"]),
                    },
                }
                if shared_paths is not None
                else {"used": False}
            ),
            "control_scope": {
                "control_target_id": rules.control_target_id,
                "control_receptor_id": rules.control_receptor_id,
                "native_bound_state_id": rules.control_bound_state_id,
                "requested_target_id": target_id,
                "target_matched": target_id == rules.control_target_id,
            },
            "topology_calibration": topology_calibration,
            "verified_decision_rules": decision_rules,
            "control_recovery": {
                "max_native_ligand_rmsd_A": rules.max_native_ligand_rmsd_A,
                "max_native_site_centroid_distance_A": (
                    rules.max_native_site_centroid_distance_A
                ),
                "min_native_receptor_contact_fraction": (
                    rules.min_native_receptor_contact_fraction
                ),
                "min_native_sampling_seed_support": (
                    rules.min_native_sampling_seed_support
                ),
                "min_native_site_seed_support": rules.min_native_site_seed_support,
                "min_selection_native_seed_recall_fraction": (
                    rules.min_selection_native_seed_recall_fraction
                ),
                "topology_feasibility": {
                    "qualification_gate": False,
                    "task_count": len(control_rows),
                    "sampling_model_count": sampling_count,
                    "topology_feasible_model_count": topology_count,
                    "topology_feasible_fraction": (
                        round(topology_count / sampling_count, 6)
                        if sampling_count
                        else 0.0
                    ),
                    "filtered_model_count": filtered_count,
                    "filtered_topology_feasible_model_count": (filtered_topology_count),
                    "filtered_topology_feasible_fraction": (
                        round(filtered_topology_count / filtered_count, 6)
                        if filtered_count
                        else 0.0
                    ),
                    "topology_feasible_enrichment_ratio": (
                        round(
                            (filtered_topology_count / filtered_count)
                            / (topology_count / sampling_count),
                            6,
                        )
                        if filtered_count and sampling_count and topology_count
                        else None
                    ),
                    "minimum_models_for_selection": (
                        config.discovery.cabsdock.min_models_for_selection
                    ),
                    "selection_skipped_task_ids": [
                        str(row["task_id"])
                        for row in control_rows
                        if row.get("selection_status") != "completed"
                    ],
                },
                "energy_filter_native_recovery_audit": {
                    "qualification_gate": False,
                    "metric": "ligand_CA_RMSD_after_receptor_CA_alignment",
                    "trajectory_recovered_model_count": full_native_recovered,
                    "trajectory_topology_feasible_recovered_model_count": (
                        full_topology_native_recovered
                    ),
                    "filtered_recovered_model_count": filtered_native_recovered,
                    "filtered_topology_feasible_recovered_model_count": (
                        filtered_topology_native_recovered
                    ),
                    "interpretation": (
                        "describes whether each-replica CABS interaction-energy "
                        "filtering retains native-region samples"
                    ),
                },
                "sampling_pool": {
                    "passed": sampling_passed,
                    "model_count": pool_count,
                    "recovered_model_count": pool_recovered_count,
                    "successful_seeds": pool_seeds,
                    "successful_seed_count": len(pool_seeds),
                    "best_ligand_ca_rmsd_A": pool_best_rmsd,
                    "best_native_receptor_contact_fraction": pool_best_contact,
                },
                "candidate_selection": {
                    "passed": selection_passed,
                    "model_count": len(control_poses),
                    "maximum_models_per_task": selected_model_budget,
                    "site_recovered_model_count": len(selection_site_recovered),
                    "successful_seeds": selection_seeds,
                    "successful_seed_count": len(selection_seeds),
                    "sampling_pool_seed_recall_fraction": round(
                        selection_seed_recall_fraction, 6
                    ),
                    "missing_sampling_pool_seeds": sorted(
                        set(pool_seeds) - set(selection_seeds)
                    ),
                    "best_ligand_ca_rmsd_A": selected_best_rmsd,
                    "best_native_receptor_contact_fraction": selected_best_contact,
                    "exact_pose_recovery_is_qualification_gate": False,
                    "exact_pose_recovered_model_count": len(selection_pose_recovered),
                    "exact_pose_successful_seeds": selection_pose_seeds,
                },
                "frozen_site_analysis": {
                    "passed": site_validation_passed,
                    "method": "maximum_pairwise_compatible_native_seed_set",
                    "native_neighborhood_seed_support": (
                        native_neighborhood_seed_support
                    ),
                    "native_neighborhood_precision": round(
                        native_neighborhood_precision, 6
                    ),
                    "qualification_gate": True,
                },
                "cabsdock_top10_baseline": {
                    "passed": (
                        len(baseline_seeds) >= rules.min_native_sampling_seed_support
                    ),
                    "model_count": len(control_baseline),
                    "recovered_model_count": len(baseline_recovered),
                    "successful_seeds": baseline_seeds,
                    "successful_seed_count": len(baseline_seeds),
                    "best_ligand_ca_rmsd_A": baseline_best_rmsd,
                    "best_native_receptor_contact_fraction": baseline_best_contact,
                    "qualification_gate": False,
                },
            },
            "target_pilot": {
                "qualification_gate": False,
                "selected_model_count": len(pilot_poses),
                "receptor_site_count": len(pilot_sites.receptor_sites),
                "supported_site_count": len(pilot_supported),
                "maximum_seed_support": pilot_max_seed_support,
                "interpretation": (
                    "technical_and_descriptive_only; absence_of_a_supported_site_is_not_"
                    "a_method_failure"
                ),
            },
            "recommended_target_config": (
                {
                    "qualification_status": "qualified",
                    "qualification_report": report_path.as_posix(),
                    "contact_jaccard_distance": (
                        target.analysis.contact_jaccard_distance
                    ),
                    "position_distance_A": target.analysis.position_distance_A,
                    "min_seed_support": target.analysis.min_seed_support,
                    "min_receptor_support": target.analysis.min_receptor_support,
                }
                if status == "qualified"
                else None
            ),
        },
    )
    return report_path
