"""阶段二 CABS-dock 任务执行、恢复和采样清单编排。"""

from __future__ import annotations

import json
import subprocess
from concurrent.futures import ThreadPoolExecutor
from importlib.metadata import version
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    is_current_vela_software,
    sha256_file,
    utc_now,
)
from vela.core.run_identity import RUN_ID_PATTERN
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.analysis.pose_table import write_pose_evidence
from vela.discovery.models import DiscoveryError, DiscoveryTask
from vela.discovery.sampling.cabsdock import (
    build_cabsdock_command,
    cabsdock_output_records,
    verify_cabsdock_tool,
    verify_native_disulfide,
)
from vela.discovery.sampling.evidence import (
    CabsDockEvidence,
    CandidateSelectionSettings,
    candidate_selection_settings,
    collect_cabsdock_evidence,
)
from vela.discovery.sampling.planning import cabsdock_parameters
from vela.discovery.sampling.trajectory import trajectory_audit_record
from vela.preparation.chemistry import ChemistryDefinition

TASK_RESULT_NAME = "task_result.json"


def _json_document(path: Path, *, name: str) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"{name} does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name=name)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid {name}: {path}") from exc


def _reference_path(*, config: AppConfig, target: str) -> tuple[str, Path]:
    receptor_id = config.discovery.target(target).reference_receptor
    path = config.paths.data_dir / "receptors" / "prepared" / f"{receptor_id}.cif"
    if not path.is_file():
        raise DiscoveryError(f"coordinate-frame receptor does not exist: {path}")
    return receptor_id, path


def _collect_task(
    *,
    task: DiscoveryTask,
    task_dir: Path,
    chemistry: ChemistryDefinition,
    reference_receptor_id: str,
    reference_path: Path,
    config: AppConfig,
    selection: CandidateSelectionSettings,
) -> CabsDockEvidence:
    return collect_cabsdock_evidence(
        task=task,
        task_dir=task_dir,
        reference_path=reference_path,
        reference_receptor_id=reference_receptor_id,
        chemistry=chemistry,
        settings=config.discovery.cabsdock,
        selection=selection,
    )


def _validate_completed_task(
    *,
    task: DiscoveryTask,
    task_dir: Path,
    run_manifest_sha256: str,
    chemistry: ChemistryDefinition,
    reference_receptor_id: str,
    reference_path: Path,
    config: AppConfig,
    selection: CandidateSelectionSettings,
) -> CabsDockEvidence:
    document = _json_document(task_dir / TASK_RESULT_NAME, name="CABS-dock task result")
    if (
        document.get("schema") != "vela.cabsdock-task-result/4"
        or document.get("execution_status") != "completed"
        or document.get("task_id") != task.task_id
        or document.get("run_manifest_sha256") != run_manifest_sha256
    ):
        raise DiscoveryError(f"stale or invalid completed task: {task.task_id}")
    outputs = object_mapping(document.get("outputs"), name="task result outputs")
    for raw in outputs.values():
        record = object_mapping(raw, name="task output record")
        relative = record.get("path")
        expected = record.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected, str):
            raise DiscoveryError(f"invalid output record for task: {task.task_id}")
        path = (task_dir / relative).resolve()
        if not path.is_relative_to(task_dir.resolve()):
            raise DiscoveryError(f"task output escapes its task directory: {path}")
        if not path.is_file() or sha256_file(path) != expected:
            raise DiscoveryError(f"task output hash mismatch: {path}")
    verify_native_disulfide(task_dir=task_dir, chemistry=chemistry)
    evidence = _collect_task(
        task=task,
        task_dir=task_dir,
        chemistry=chemistry,
        reference_receptor_id=reference_receptor_id,
        reference_path=reference_path,
        config=config,
        selection=selection,
    )
    audit = _json_document(
        task_dir / "trajectory_audit.json", name="CABS trajectory audit"
    )
    if audit != trajectory_audit_record(evidence.trajectory_audit):
        raise DiscoveryError(f"stale CABS trajectory audit: {task.task_id}")
    return evidence


def run_cabsdock_task(
    *,
    task: DiscoveryTask,
    run_dir: Path,
    run_manifest_sha256: str,
    chemistry: ChemistryDefinition,
    secondary_structure: str,
    reference_receptor_id: str,
    reference_path: Path,
    config: AppConfig,
    tool_version: str,
    selection: CandidateSelectionSettings,
) -> CabsDockEvidence:
    task_dir = run_dir / "tasks" / task.task_id
    result_path = task_dir / TASK_RESULT_NAME
    if result_path.is_file():
        return _validate_completed_task(
            task=task,
            task_dir=task_dir,
            run_manifest_sha256=run_manifest_sha256,
            chemistry=chemistry,
            reference_receptor_id=reference_receptor_id,
            reference_path=reference_path,
            config=config,
            selection=selection,
        )
    if task_dir.exists():
        raise DiscoveryError(
            f"incomplete task directory must be reviewed before retry: {task_dir}"
        )
    task_dir.mkdir(parents=True)
    command = build_cabsdock_command(
        task=task,
        settings=config.discovery.cabsdock,
        chemistry=chemistry,
        secondary_structure=secondary_structure,
        task_dir=task_dir,
    )
    started_at = utc_now()
    log_path = task_dir / "cabsdock.log"
    with log_path.open("w", encoding="utf-8") as log:
        result = subprocess.run(
            command,
            check=False,
            stdout=log,
            stderr=subprocess.STDOUT,
            text=True,
        )
    if result.returncode != 0:
        raise DiscoveryError(
            f"CABS-dock task failed with exit code {result.returncode}: {task.task_id}"
        )
    verify_native_disulfide(task_dir=task_dir, chemistry=chemistry)
    evidence = _collect_task(
        task=task,
        task_dir=task_dir,
        chemistry=chemistry,
        reference_receptor_id=reference_receptor_id,
        reference_path=reference_path,
        config=config,
        selection=selection,
    )
    atomic_write_json(
        task_dir / "trajectory_audit.json",
        trajectory_audit_record(evidence.trajectory_audit),
    )
    outputs = cabsdock_output_records(task_dir)
    outputs["candidate_frame_selection"] = {
        "path": evidence.candidate_frame_selection_path.relative_to(
            task_dir
        ).as_posix(),
        "sha256": sha256_file(evidence.candidate_frame_selection_path),
    }
    if evidence.candidate_models_path is not None:
        outputs["candidate_models"] = {
            "path": evidence.candidate_models_path.relative_to(task_dir).as_posix(),
            "sha256": sha256_file(evidence.candidate_models_path),
        }
    atomic_write_json(
        result_path,
        {
            "schema": "vela.cabsdock-task-result/4",
            "execution_status": "completed",
            "task_id": task.task_id,
            "run_manifest_sha256": run_manifest_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "software": {
                "cabsdock_version": tool_version,
                "cabsdock_source_revision": (config.discovery.cabsdock.source_revision),
                "cabsdock_patch_sha256": sha256_file(
                    config.discovery.cabsdock.patch_file
                ),
            },
            "command": list(command),
            "sampling_evidence": {
                "sampling_model_count": evidence.sampling_model_count,
                "filtered_model_count": evidence.filtered_model_count,
                "filtered_topology_feasible_model_count": (
                    evidence.filtered_topology_feasible_model_count
                ),
                "filtered_topology_feasible_fraction": (
                    evidence.filtered_topology_feasible_fraction
                ),
                "topology_feasible_model_count": (
                    evidence.topology_feasible_model_count
                ),
                "topology_feasible_fraction": evidence.topology_feasible_fraction,
                "contacting_topology_feasible_model_count": (
                    evidence.contacting_topology_feasible_model_count
                ),
                "disulfide_ca_distance_median_A": (
                    evidence.disulfide_ca_distance_median_A
                ),
                "disulfide_ca_distance_p90_A": (evidence.disulfide_ca_distance_p90_A),
                "disulfide_ca_distance_max_A": (evidence.disulfide_ca_distance_max_A),
                "receptor_alignment_max_rmsd_A": (
                    evidence.receptor_alignment_max_rmsd_A
                ),
            },
            "selection": {
                "status": evidence.selection_status,
                "failure_reasons": list(evidence.selection_failure_reasons),
                "site_cluster_count": evidence.site_cluster_count,
                "pose_cluster_count": evidence.pose_cluster_count,
                "selected_pose_cluster_count": (evidence.selected_pose_cluster_count),
                "selected_model_count": len(evidence.poses),
            },
            "baseline": {
                "method": "cabsdock_top10_kmedoids",
                "model_count": len(evidence.baseline_poses),
            },
            "outputs": outputs,
        },
    )
    return evidence


def _required_text(value: object, *, name: str) -> str:
    if not isinstance(value, str) or not value.strip():
        raise DiscoveryError(f"{name} must be non-empty text")
    return value


def _required_seed(value: object, *, name: str) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise DiscoveryError(f"{name} must be a non-negative integer")
    return value


def _receptor_path(*, config: AppConfig, raw: object, task_id: str) -> tuple[Path, str]:
    try:
        record = object_mapping(raw, name=f"{task_id} receptor")
    except TypeError as exc:
        raise DiscoveryError(f"invalid receptor record for task: {task_id}") from exc
    relative = _required_text(record.get("path"), name=f"{task_id} receptor path")
    expected = _required_text(record.get("sha256"), name=f"{task_id} receptor SHA-256")
    path = (config.paths.data_dir / relative).resolve()
    try:
        path.relative_to(config.paths.data_dir.resolve())
    except ValueError as exc:
        raise DiscoveryError(
            f"task receptor escapes the data directory: {path}"
        ) from exc
    if not path.is_file() or sha256_file(path) != expected:
        raise DiscoveryError(f"task receptor hash mismatch: {path}")
    return path, expected


def _planned_tasks(
    *, config: AppConfig, document: dict[str, object]
) -> tuple[DiscoveryTask, ...]:
    try:
        rows = object_list(document.get("tasks"), name="run manifest tasks")
    except TypeError as exc:
        raise DiscoveryError("run manifest tasks are invalid") from exc
    category = _required_text(
        document.get("evidence_category"), name="run evidence_category"
    )
    target_id = _required_text(document.get("target_id"), name="run target_id")
    tasks: list[DiscoveryTask] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="run manifest task")
        except TypeError as exc:
            raise DiscoveryError("run manifest task is invalid") from exc
        task_id = _required_text(row.get("task_id"), name="task_id")
        if not RUN_ID_PATTERN.fullmatch(task_id) or row.get("status") != "planned":
            raise DiscoveryError("run manifest tasks must be planned and well-formed")
        receptor_path, receptor_sha256 = _receptor_path(
            config=config, raw=row.get("receptor"), task_id=task_id
        )
        task = DiscoveryTask(
            task_id=task_id,
            receptor_id=_required_text(
                row.get("receptor_id"), name=f"{task_id} receptor_id"
            ),
            target=_required_text(row.get("target"), name=f"{task_id} target"),
            receptor_path=receptor_path,
            receptor_sha256=receptor_sha256,
            chemistry_id=_required_text(
                row.get("chemistry_id"), name=f"{task_id} chemistry_id"
            ),
            method_id=_required_text(row.get("method_id"), name=f"{task_id} method_id"),
            adapter_id=_required_text(
                row.get("adapter_id"), name=f"{task_id} adapter_id"
            ),
            seed=_required_seed(row.get("seed"), name=f"{task_id} seed"),
            evidence_category=_required_text(
                row.get("evidence_category"),
                name=f"{task_id} evidence_category",
            ),
        )
        if (
            task.evidence_category != category
            or task.target != target_id
            or task.chemistry_id != config.chemistry.chemistry_id
            or task.method_id != config.discovery.method_id
            or task.adapter_id != config.discovery.adapter_id
            or task.seed not in config.discovery.seeds
        ):
            raise DiscoveryError(
                f"task identity differs from the current frozen config: {task_id}"
            )
        tasks.append(task)
    task_ids = tuple(task.task_id for task in tasks)
    if not tasks or len(task_ids) != len(set(task_ids)):
        raise DiscoveryError("run manifest must contain unique sampling tasks")
    if document.get("task_count") != len(tasks):
        raise DiscoveryError("run manifest task_count does not match its tasks")
    return tuple(tasks)


def group_tasks_by_seed(
    tasks: tuple[DiscoveryTask, ...],
) -> tuple[tuple[DiscoveryTask, ...], ...]:
    """把单目标任务按 seed 分组, 保留受体构象的计划顺序。"""
    if not tasks:
        raise DiscoveryError("sampling task list is empty")
    targets = {task.target for task in tasks}
    if len(targets) != 1:
        raise DiscoveryError("one sampling run must contain exactly one target")
    grouped: dict[int, list[DiscoveryTask]] = {}
    for task in tasks:
        grouped.setdefault(task.seed, []).append(task)
    batches = tuple(tuple(batch) for batch in grouped.values())
    for batch in batches:
        receptor_ids = tuple(task.receptor_id for task in batch)
        if len(receptor_ids) != len(set(receptor_ids)):
            raise DiscoveryError(
                f"seed batch contains duplicate receptors: {batch[0].seed}"
            )
    return batches


def _run_seed_batch(
    *,
    tasks: tuple[DiscoveryTask, ...],
    run_dir: Path,
    run_manifest_sha256: str,
    config: AppConfig,
    tool_version: str,
    selection: CandidateSelectionSettings,
) -> tuple[tuple[str, CabsDockEvidence], ...]:
    """在一个 worker 中顺序完成同一 seed 的全部受体构象。"""
    return tuple(
        (
            task.task_id,
            run_cabsdock_task(
                task=task,
                run_dir=run_dir,
                run_manifest_sha256=run_manifest_sha256,
                chemistry=config.chemistry,
                secondary_structure=(
                    config.discovery.cabsdock.peptide_secondary_structure
                ),
                reference_receptor_id=_reference_path(
                    config=config, target=task.target
                )[0],
                reference_path=_reference_path(config=config, target=task.target)[1],
                config=config,
                tool_version=tool_version,
                selection=selection,
            ),
        )
        for task in tasks
    )


def run_cabsdock_sampling(*, config: AppConfig, run_dir: Path) -> None:
    """按冻结任务执行或恢复 CABS-dock; 写出统一采样证据。"""
    plan_path = run_dir / "run_manifest.json"
    plan = _json_document(plan_path, name="discovery run manifest")
    if (
        plan.get("schema") != "vela.discovery-run-manifest/4"
        or plan.get("stage") != "discovery"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("known_site_information_used") is not False
    ):
        raise DiscoveryError("run manifest must be a planned discovery run")
    if (run_dir / "sampling_manifest.json").exists() or (
        run_dir / "pose_evidence.tsv"
    ).exists():
        raise DiscoveryError(f"discovery sampling outputs already exist: {run_dir}")
    inputs = object_mapping(plan.get("inputs"), name="run manifest inputs")
    snapshot = object_mapping(inputs.get("config_snapshot"), name="config snapshot")
    relative_snapshot = snapshot.get("path")
    if not isinstance(relative_snapshot, str):
        raise DiscoveryError("config snapshot path is invalid")
    snapshot_path = run_dir / relative_snapshot
    if snapshot.get("sha256") != sha256_file(
        snapshot_path
    ) or config.source_snapshot_sha256 != sha256_file(snapshot_path):
        raise DiscoveryError("current project config differs from the frozen run")
    parameters = object_mapping(plan.get("method_parameters"), name="method parameters")
    recorded = object_mapping(parameters.get("cabsdock"), name="CABS-dock parameters")
    if recorded != cabsdock_parameters(config):
        raise DiscoveryError("current CABS-dock parameters differ from the frozen run")

    tasks = _planned_tasks(config=config, document=plan)
    target_id = _required_text(plan.get("target_id"), name="run target_id")
    evidence_category = tasks[0].evidence_category
    tool_version = verify_cabsdock_tool(config.discovery.cabsdock)
    analysis = config.discovery.target(target_id).analysis
    if (
        analysis.contact_jaccard_distance is None
        or analysis.position_distance_A is None
    ):
        raise DiscoveryError("production candidate selection thresholds are unresolved")
    selection = candidate_selection_settings(config.discovery.cabsdock)
    plan_hash = sha256_file(plan_path)
    all_poses: list[PoseEvidence] = []
    task_rows: list[dict[str, JsonValue]] = []
    seed_batches = group_tasks_by_seed(tasks)
    worker_count = min(config.discovery.cabsdock.seed_workers, len(seed_batches))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = tuple(
            executor.submit(
                _run_seed_batch,
                tasks=batch,
                run_dir=run_dir,
                run_manifest_sha256=plan_hash,
                config=config,
                tool_version=tool_version,
                selection=selection,
            )
            for batch in seed_batches
        )
        evidence_by_task = {
            task_id: evidence
            for future in futures
            for task_id, evidence in future.result()
        }
    if set(evidence_by_task) != {task.task_id for task in tasks}:
        raise DiscoveryError("seed workers returned incomplete sampling evidence")
    for task in tasks:
        evidence = evidence_by_task[task.task_id]
        all_poses.extend(evidence.poses)
        task_rows.append(
            {
                "task_id": task.task_id,
                "receptor_id": task.receptor_id,
                "target": task.target,
                "seed": task.seed,
                "chemistry_id": task.chemistry_id,
                "method_id": task.method_id,
                "adapter_id": task.adapter_id,
                "evidence_category": task.evidence_category,
                "execution_status": "completed",
                "selection_status": evidence.selection_status,
                "selection_failure_reasons": list(evidence.selection_failure_reasons),
                "trajectory_model_count": (
                    evidence.trajectory_audit.trajectory_model_count
                ),
                "trajectory_topology_feasible_model_count": (
                    evidence.trajectory_audit.trajectory_topology_feasible_model_count
                ),
                "topology_feasible_model_count": (
                    evidence.topology_feasible_model_count
                ),
                "contacting_topology_feasible_model_count": (
                    evidence.contacting_topology_feasible_model_count
                ),
                "selected_model_count": len(evidence.poses),
            }
        )
    pose_path = run_dir / "pose_evidence.tsv"
    write_pose_evidence(poses=tuple(all_poses), path=pose_path, run_dir=run_dir)
    atomic_write_json(
        run_dir / "sampling_manifest.json",
        {
            "schema": "vela.discovery-sampling-manifest/4",
            "stage": "discovery",
            "target_id": target_id,
            "status": "sampling_completed",
            "completed_at": utc_now(),
            "run_manifest_sha256": plan_hash,
            "software": {
                "method_version": tool_version,
                "adapter_version": version("vela"),
                "cabsdock_source_revision": (config.discovery.cabsdock.source_revision),
                "cabsdock_patch_sha256": sha256_file(
                    config.discovery.cabsdock.patch_file
                ),
            },
            "representation": "CABS coarse-grained CA/SC",
            "evidence_category": evidence_category,
            "known_site_information_used": False,
            "all_atom_reconstruction_used": False,
            "execution": {
                "unit": "seed_batch",
                "seed_workers": worker_count,
                "seed_batch_count": len(seed_batches),
                "receptors_per_seed": len(seed_batches[0]),
            },
            "pose_evidence": {
                "path": pose_path.name,
                "sha256": sha256_file(pose_path),
            },
            "tasks": task_rows,
        },
    )
