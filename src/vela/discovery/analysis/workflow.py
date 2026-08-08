"""阶段二规范 pose 完整性检查与 site 分析编排。"""

from __future__ import annotations

import json
from pathlib import Path

from vela.core.provenance import is_current_vela_software, sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.analysis.clustering import (
    analyze_sites,
    candidate_analysis_contract,
)
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.analysis.reports import write_site_reports
from vela.discovery.models import DiscoveryError, SiteAnalysisSettings


def _document(path: Path, *, name: str) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"{name} does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        document = object_mapping(value, name=str(path))
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid {name}: {path}") from exc
    return document


type TaskIdentity = tuple[str, str, int, str, str, str]


def _task_identities(
    document: dict[str, object], *, status_field: str, required_status: str
) -> dict[str, TaskIdentity]:
    try:
        raw_tasks = object_list(document.get("tasks"), name="run_manifest.tasks")
    except TypeError as exc:
        raise DiscoveryError("run manifest tasks are invalid") from exc
    tasks: dict[str, TaskIdentity] = {}
    for raw in raw_tasks:
        try:
            task = object_mapping(raw, name="run_manifest task")
        except TypeError as exc:
            raise DiscoveryError("run manifest task is invalid") from exc
        task_id = task.get("task_id")
        receptor_id = task.get("receptor_id")
        target = task.get("target")
        seed = task.get("seed")
        chemistry_id = task.get("chemistry_id")
        method_id = task.get("method_id")
        adapter_id = task.get("adapter_id")
        status = task.get(status_field)
        if (
            not isinstance(task_id, str)
            or not isinstance(receptor_id, str)
            or not isinstance(target, str)
            or not isinstance(seed, int)
            or isinstance(seed, bool)
            or not isinstance(chemistry_id, str)
            or not chemistry_id
            or not isinstance(method_id, str)
            or not method_id
            or not isinstance(adapter_id, str)
            or not adapter_id
            or status != required_status
        ):
            raise DiscoveryError(f"all tasks must be {required_status} and well-formed")
        if task_id in tasks:
            raise DiscoveryError(f"duplicate task_id in run manifest: {task_id}")
        tasks[task_id] = (
            receptor_id,
            target,
            seed,
            chemistry_id,
            method_id,
            adapter_id,
        )
    if not tasks:
        raise DiscoveryError("run manifest contains no sampling tasks")
    return tasks


def _completed_tasks(run_dir: Path) -> tuple[dict[str, TaskIdentity], str, str]:
    plan_path = run_dir / "run_manifest.json"
    plan = _document(plan_path, name="discovery run manifest")
    if plan.get("stage") != "discovery" or plan.get("status") != "planned":
        raise DiscoveryError("run manifest must be a planned discovery run")
    evidence_category = plan.get("evidence_category")
    target_id = plan.get("target_id")
    if (
        not isinstance(evidence_category, str)
        or not evidence_category.strip()
        or not isinstance(target_id, str)
        or not target_id.strip()
        or plan.get("known_site_information_used") is not False
        or not is_current_vela_software(plan.get("software"))
    ):
        raise DiscoveryError("run manifest has an invalid blind evidence category")
    if plan.get("schema") != "vela.discovery-run-manifest/8":
        raise DiscoveryError("discovery run manifest schema is invalid")
    planned_tasks = _task_identities(
        plan, status_field="status", required_status="planned"
    )

    sampling = _document(run_dir / "sampling_manifest.json", name="sampling manifest")
    if (
        sampling.get("stage") != "discovery"
        or sampling.get("target_id") != target_id
        or sampling.get("status") != "sampling_completed"
    ):
        raise DiscoveryError("sampling manifest must report sampling_completed")
    if sampling.get("run_manifest_sha256") != sha256_file(plan_path):
        raise DiscoveryError("sampling manifest does not match the frozen run manifest")
    try:
        software = object_mapping(
            sampling.get("software"), name="sampling_manifest.software"
        )
    except TypeError as exc:
        raise DiscoveryError("sampling manifest software is invalid") from exc
    for field in ("method_version", "adapter_version"):
        value = software.get(field)
        if not isinstance(value, str) or not value.strip():
            raise DiscoveryError(f"sampling manifest must record {field}")
    if sampling.get("schema") != "vela.discovery-sampling-manifest/7":
        raise DiscoveryError("sampling manifest schema is invalid")
    completed_tasks = _task_identities(
        sampling, status_field="execution_status", required_status="completed"
    )
    if completed_tasks != planned_tasks:
        raise DiscoveryError("sampling task identities do not match the frozen plan")
    if {identity[1] for identity in planned_tasks.values()} != {target_id}:
        raise DiscoveryError("discovery run mixes targets or disagrees with target_id")
    return completed_tasks, evidence_category, target_id


def discovery_run_target(run_dir: Path) -> str:
    """从冻结清单读取并校验一个已计划运行的靶标身份。"""
    plan = _document(run_dir / "run_manifest.json", name="discovery run manifest")
    target_id = plan.get("target_id")
    if plan.get("stage") != "discovery" or not isinstance(target_id, str):
        raise DiscoveryError("discovery run target identity is invalid")
    return target_id


def analyze_discovery_run(*, run_dir: Path, settings: SiteAnalysisSettings) -> None:
    """确认正式任务完整后生成单受体和跨构象 site 报告。"""
    resolved_run_dir = run_dir.resolve()
    tasks, evidence_category, _ = _completed_tasks(resolved_run_dir)
    plan = _document(
        resolved_run_dir / "run_manifest.json", name="discovery run manifest"
    )
    if plan.get("analysis_contract") != candidate_analysis_contract(settings):
        raise DiscoveryError(
            "current candidate analysis contract differs from the frozen run"
        )
    pose_table = resolved_run_dir / "pose_evidence.tsv"
    poses = read_pose_evidence(path=pose_table, run_dir=resolved_run_dir)
    observed_tasks: set[str] = set()
    for pose in poses:
        expected = tasks.get(pose.task_id)
        if expected is None:
            raise DiscoveryError(f"pose refers to unknown task: {pose.task_id}")
        if expected[:3] != (pose.receptor_id, pose.target, pose.seed):
            raise DiscoveryError(f"pose identity disagrees with task: {pose.pose_id}")
        observed_tasks.add(pose.task_id)
    sampling = _document(
        resolved_run_dir / "sampling_manifest.json", name="sampling manifest"
    )
    selectable_tasks = {
        str(task.get("task_id"))
        for raw in object_list(sampling.get("tasks"), name="sampling tasks")
        for task in (object_mapping(raw, name="sampling task"),)
        if task.get("selection_status") == "completed"
    }
    missing = sorted(selectable_tasks - observed_tasks)
    if missing:
        raise DiscoveryError(
            "completed tasks have no pose evidence records: " + ", ".join(missing)
        )
    result = analyze_sites(poses=poses, settings=settings)
    write_site_reports(
        result=result,
        pose_table=pose_table,
        output_dir=resolved_run_dir / "site_analysis",
        evidence_category=evidence_category,
        settings=settings,
    )
