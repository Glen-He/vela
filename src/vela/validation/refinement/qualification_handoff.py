"""阶段二资格证据到阶段三的开发性全原子交接。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    is_vela_software_identity,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.analysis.clustering import analyze_sites
from vela.discovery.analysis.pose_table import read_pose_evidence
from vela.discovery.qualification.control import (
    control_bound_state,
    control_chemistry,
)
from vela.discovery.qualification.planning import CONTROL_RECOVERY
from vela.discovery.qualification.schemas import (
    PLAN_SCHEMA as DISCOVERY_PLAN_SCHEMA,
)
from vela.discovery.qualification.schemas import (
    REPORT_SCHEMA as DISCOVERY_REPORT_SCHEMA,
)
from vela.discovery.qualification.schemas import (
    SAMPLING_SCHEMA as DISCOVERY_SAMPLING_SCHEMA,
)
from vela.preparation.chemistry import chemistry_identity_record
from vela.validation.models import ValidationError
from vela.validation.records import (
    file_record,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.handoff_plan import (
    POSE_SELECTION_METHOD,
    POSE_SELECTION_SEED_POLICY,
    CandidateHandoffTask,
    handoff_task_records,
    select_site_poses,
)
from vela.validation.refinement.handoff_run import (
    execute_handoff_tasks,
    handoff_result_record,
)
from vela.validation.refinement.reconstruction import verify_cg2all_tool
from vela.validation.rosetta import verify_rosetta_scripts_tool

PLAN_NAME = "qualification_handoff_plan.json"
MANIFEST_NAME = "qualification_handoff_manifest.json"
PLAN_SCHEMA = "vela.validation-qualification-handoff-plan/5"
MANIFEST_SCHEMA = "vela.validation-qualification-handoff-manifest/5"
EVIDENCE_CATEGORY = "qualification_development_handoff"
SITE_RANKING = "supporting_seed_count_desc,pose_count_desc,site_id_asc"


@dataclass(frozen=True, slots=True)
class QualificationHandoffPlan:
    """一个冻结的开发性交接计划。"""

    run_dir: Path
    tasks: tuple[CandidateHandoffTask, ...]


@dataclass(frozen=True, slots=True)
class QualificationHandoffOutcome:
    """开发性交接完成清单及其规模。"""

    manifest_path: Path
    task_count: int


def _qualification_documents(
    qualification_run_dir: Path,
) -> tuple[dict[str, object], dict[str, object]]:
    plan = read_document(
        qualification_run_dir / "qualification_plan.json",
        name="discovery qualification plan",
    )
    sampling = read_document(
        qualification_run_dir / "qualification_sampling.json",
        name="discovery qualification sampling",
    )
    if (
        plan.get("schema") != DISCOVERY_PLAN_SCHEMA
        or plan.get("stage") != "discovery_qualification"
        or plan.get("status") != "planned"
    ):
        raise ValidationError("discovery qualification plan identity is invalid")
    if (
        sampling.get("schema") != DISCOVERY_SAMPLING_SCHEMA
        or sampling.get("stage") != "discovery_qualification"
        or sampling.get("status") != "sampling_completed"
    ):
        raise ValidationError("discovery qualification sampling identity is invalid")
    target_id = plan.get("target_id")
    if target_id != sampling.get("target_id"):
        raise ValidationError("discovery qualification target identity is inconsistent")
    try:
        sampling_rows = object_list(
            sampling.get("tasks"), name="qualification sampling tasks"
        )
    except TypeError as exc:
        raise ValidationError("qualification sampling tasks are invalid") from exc
    parsed_sampling_rows: list[dict[str, object]] = []
    for raw in sampling_rows:
        try:
            parsed_sampling_rows.append(
                object_mapping(raw, name="qualification sampling task")
            )
        except TypeError as exc:
            raise ValidationError("qualification sampling task is invalid") from exc
    if not parsed_sampling_rows or any(
        row.get("execution_status") != "completed"
        or row.get("selection_status") != "completed"
        for row in parsed_sampling_rows
    ):
        raise ValidationError("qualification sampling is incomplete")
    return plan, sampling


def _qualification_report(
    *, qualification_run_dir: Path, target_id: object
) -> dict[str, object]:
    """只在候选选择完成后核对分析报告身份和来源状态。"""
    report = read_document(
        qualification_run_dir / "qualification_report.json",
        name="discovery qualification report",
    )
    if (
        report.get("schema") != DISCOVERY_REPORT_SCHEMA
        or report.get("stage") != "discovery_qualification"
        or report.get("status") not in {"qualified", "unqualified"}
        or report.get("target_id") != target_id
        or not is_vela_software_identity(report.get("analysis_software"))
    ):
        raise ValidationError("discovery qualification report identity is invalid")
    return report


def _control_task_ids(plan: dict[str, object]) -> frozenset[str]:
    try:
        rows = object_list(plan.get("tasks"), name="qualification plan tasks")
    except TypeError as exc:
        raise ValidationError("qualification plan tasks are invalid") from exc
    task_ids: set[str] = set()
    for raw in rows:
        try:
            row = object_mapping(raw, name="qualification plan task")
        except TypeError as exc:
            raise ValidationError("qualification plan task is invalid") from exc
        task_id = row.get("task_id")
        case = row.get("case")
        if not isinstance(task_id, str) or not isinstance(case, str):
            raise ValidationError("qualification plan task identity is invalid")
        if case == CONTROL_RECOVERY:
            task_ids.add(task_id)
    if not task_ids:
        raise ValidationError("qualification plan contains no control tasks")
    return frozenset(task_ids)


def build_qualification_handoff_tasks(
    *, config: AppConfig, qualification_run_dir: Path, site_budget: int
) -> tuple[CandidateHandoffTask, ...]:
    """按冻结的 native-free 排名选择 Top-B site 及非冗余起点。"""
    if site_budget != config.discovery.qualification.receptor_site_diagnostic_budget:
        raise ValidationError(
            "qualification handoff must use the frozen site delivery budget"
        )
    plan, _ = _qualification_documents(qualification_run_dir)
    try:
        scope = object_mapping(plan.get("control_scope"), name="control scope")
    except TypeError as exc:
        raise ValidationError("qualification control scope is invalid") from exc
    target_id = scope.get("control_target_id")
    receptor_ids = scope.get("control_receptor_ids")
    if (
        target_id != config.discovery.qualification.control_target_id
        or receptor_ids != list(config.discovery.qualification.control_receptor_ids)
        or scope.get("protein_target_matched") is not True
    ):
        raise ValidationError(
            "qualification control does not match the current project"
        )
    control_ids = _control_task_ids(plan)
    poses = tuple(
        pose
        for pose in read_pose_evidence(
            path=qualification_run_dir / "pose_evidence.tsv",
            run_dir=qualification_run_dir,
        )
        if pose.task_id in control_ids
    )
    pose_by_id = {pose.pose_id: pose for pose in poses}
    analysis = config.discovery.target(str(target_id)).analysis
    chemistry = control_chemistry(control_bound_state(config))
    tasks: list[CandidateHandoffTask] = []
    for receptor_id in config.discovery.qualification.control_receptor_ids:
        receptor_poses = tuple(
            pose for pose in poses if pose.receptor_id == receptor_id
        )
        sites = tuple(
            sorted(
                (
                    site
                    for site in analyze_sites(
                        poses=receptor_poses, settings=analysis
                    ).receptor_sites
                    if site.supported
                ),
                key=lambda site: (
                    -len(site.supporting_seeds),
                    -site.pose_count,
                    site.site_id,
                ),
            )
        )
        if len(sites) < site_budget:
            raise ValidationError(
                f"qualification control {receptor_id} has only {len(sites)} "
                f"supported sites; requested {site_budget}"
            )
        reference_receptor_path = (
            config.paths.data_dir / "receptors" / "prepared" / f"{receptor_id}.cif"
        )
        if not reference_receptor_path.is_file():
            raise ValidationError(
                f"prepared control receptor is missing: {reference_receptor_path}"
            )
        reference_hash = sha256_file(reference_receptor_path)
        for site in sites[:site_budget]:
            candidate_id = safe_identifier(
                f"{receptor_id}__{site.site_id}", name="qualification site ID"
            )
            for pose in select_site_poses(
                site=site,
                pose_by_id=pose_by_id,
                count=config.validation.handoff.poses_per_receptor_site,
                peptide_sequence=chemistry.sequence,
                reference_receptor_path=reference_receptor_path,
                pose_clustering_rmsd_A=(
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
            ):
                task_id = safe_identifier(
                    f"{receptor_id}__{site.site_id}__{pose.pose_id}",
                    name="qualification handoff task ID",
                )
                tasks.append(
                    CandidateHandoffTask(
                        task_id=task_id,
                        candidate_id=candidate_id,
                        receptor_site_id=site.site_id,
                        pose=pose,
                        reference_receptor_path=reference_receptor_path,
                        reference_receptor_sha256=reference_hash,
                    )
                )
    return tuple(tasks)


def _source_records(qualification_run_dir: Path) -> list[dict[str, JsonValue]]:
    return [
        file_record(qualification_run_dir / name, root=qualification_run_dir)
        for name in (
            "qualification_plan.json",
            "qualification_sampling.json",
            "qualification_report.json",
            "pose_evidence.tsv",
        )
    ]


def write_qualification_handoff_plan(
    *,
    config: AppConfig,
    qualification_run_dir: Path,
    run_id: str,
    site_budget: int,
) -> QualificationHandoffPlan:
    """冻结开发证据来源、Top-B 预算、起点选择及工具身份。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    run_dir = (
        config.paths.outputs_dir / "validation" / "qualification_handoffs" / run_id
    )
    if run_dir.exists():
        raise ValidationError(
            f"qualification handoff run directory already exists: {run_dir}"
        )
    tasks = build_qualification_handoff_tasks(
        config=config,
        qualification_run_dir=qualification_run_dir,
        site_budget=site_budget,
    )
    source_plan, _ = _qualification_documents(qualification_run_dir)
    _qualification_report(
        qualification_run_dir=qualification_run_dir,
        target_id=source_plan.get("target_id"),
    )
    chemistry = control_chemistry(control_bound_state(config))
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    rosetta = verify_rosetta_scripts_tool(config.validation.rosetta)
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    source_root = qualification_run_dir.resolve()
    outputs_root = config.paths.outputs_dir.resolve()
    if not source_root.is_relative_to(
        (outputs_root / "discovery" / "qualifications").resolve()
    ):
        raise ValidationError("source is outside discovery qualification outputs")
    atomic_write_json(
        run_dir / PLAN_NAME,
        {
            "schema": PLAN_SCHEMA,
            "stage": "validation_qualification_handoff",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "evidence_category": EVIDENCE_CATEGORY,
            "development_only": True,
            "known_site_information_used_for_selection": False,
            "formal_qualification_gate": False,
            "chemistry": chemistry_identity_record(chemistry),
            "software": {
                **vela_software_identity(),
                "cg2all_version": cg2all.version,
                "cg2all_executable_sha256": cg2all.executable_sha256,
                "cg2all_checkpoint_sha256": cg2all.checkpoint_sha256,
                "rosetta_version": rosetta.version,
                "rosetta_scripts_sha256": rosetta.executable_sha256,
            },
            "inputs": {
                "config_snapshot": file_record(snapshot, root=run_dir),
                "qualification_run": {
                    "path": source_root.relative_to(outputs_root).as_posix(),
                    "files": _source_records(source_root),
                },
            },
            "selection": {
                "site_budget": site_budget,
                "site_ranking": SITE_RANKING,
                "poses_per_site": config.validation.handoff.poses_per_receptor_site,
                "pose_clustering_rmsd_A": (
                    config.discovery.cabsdock.pose_clustering_rmsd_A
                ),
                "pose_selection": POSE_SELECTION_METHOD,
                "source_seed_policy": POSE_SELECTION_SEED_POLICY,
                "native_metrics_read": False,
            },
            "execution": {
                "parallel_tasks": min(
                    config.validation.rosetta.parallel_tasks, len(tasks)
                )
            },
            "tasks": handoff_task_records(
                tasks=tasks,
                discovery_run_dir=source_root,
                data_dir=config.paths.data_dir,
            ),
        },
    )
    return QualificationHandoffPlan(run_dir, tasks)


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], tuple[CandidateHandoffTask, ...]]:
    plan = read_document(run_dir / PLAN_NAME, name="qualification handoff plan")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or plan.get("stage") != "validation_qualification_handoff"
        or plan.get("status") != "planned"
        or plan.get("evidence_category") != EVIDENCE_CATEGORY
        or plan.get("development_only") is not True
        or plan.get("known_site_information_used_for_selection") is not False
        or plan.get("formal_qualification_gate") is not False
        or not is_current_vela_software(plan.get("software"))
    ):
        raise ValidationError("qualification handoff plan identity is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="handoff inputs")
        snapshot = object_mapping(inputs.get("config_snapshot"), name="config snapshot")
        source = object_mapping(
            inputs.get("qualification_run"), name="qualification source"
        )
        source_files = object_list(source.get("files"), name="qualification files")
        selection = object_mapping(plan.get("selection"), name="handoff selection")
        execution = object_mapping(plan.get("execution"), name="handoff execution")
        recorded_tasks = object_list(plan.get("tasks"), name="handoff tasks")
        recorded_chemistry = object_mapping(
            plan.get("chemistry"), name="control chemistry"
        )
    except TypeError as exc:
        raise ValidationError(
            "qualification handoff plan structure is invalid"
        ) from exc
    snapshot_path, _ = validate_record(
        root=run_dir, raw=snapshot, name="config snapshot"
    )
    if sha256_file(snapshot_path) != config.source_snapshot_sha256:
        raise ValidationError("current project config differs from the handoff plan")
    chemistry = control_chemistry(control_bound_state(config))
    if recorded_chemistry != chemistry_identity_record(chemistry):
        raise ValidationError("qualification control chemistry changed")
    relative_source = source.get("path")
    if not isinstance(relative_source, str):
        raise ValidationError("qualification source path is invalid")
    source_root = (config.paths.outputs_dir / relative_source).resolve()
    qualification_root = (
        config.paths.outputs_dir / "discovery" / "qualifications"
    ).resolve()
    if not source_root.is_relative_to(qualification_root):
        raise ValidationError("qualification source is outside declared outputs")
    for record in source_files:
        validate_record(root=source_root, raw=record, name="qualification source file")
    site_budget = selection.get("site_budget")
    if not isinstance(site_budget, int) or isinstance(site_budget, bool):
        raise ValidationError("qualification handoff site budget is invalid")
    if (
        selection.get("site_ranking") != SITE_RANKING
        or selection.get("poses_per_site")
        != config.validation.handoff.poses_per_receptor_site
        or selection.get("pose_clustering_rmsd_A")
        != config.discovery.cabsdock.pose_clustering_rmsd_A
        or selection.get("pose_selection") != POSE_SELECTION_METHOD
        or selection.get("source_seed_policy") != POSE_SELECTION_SEED_POLICY
        or selection.get("native_metrics_read") is not False
    ):
        raise ValidationError("qualification handoff selection contract changed")
    if execution.get("parallel_tasks") != min(
        config.validation.rosetta.parallel_tasks, len(recorded_tasks)
    ):
        raise ValidationError("qualification handoff execution contract changed")
    tasks = build_qualification_handoff_tasks(
        config=config,
        qualification_run_dir=source_root,
        site_budget=site_budget,
    )
    expected = handoff_task_records(
        tasks=tasks,
        discovery_run_dir=source_root,
        data_dir=config.paths.data_dir,
    )
    if recorded_tasks != expected:
        raise ValidationError("qualification handoff tasks differ from the frozen plan")
    verify_cg2all_tool(config.validation.cg2all)
    verify_rosetta_scripts_tool(config.validation.rosetta)
    return plan, tasks


def run_qualification_handoff(
    *, config: AppConfig, run_dir: Path
) -> QualificationHandoffOutcome:
    """执行开发性全原子重建; 保持其不可晋升的证据身份。"""
    manifest_path = run_dir / MANIFEST_NAME
    if manifest_path.exists():
        raise ValidationError(
            f"qualification handoff manifest already exists: {manifest_path}"
        )
    _, tasks = _verify_plan(config=config, run_dir=run_dir)
    chemistry = control_chemistry(control_bound_state(config))
    plan_path = run_dir / PLAN_NAME
    plan_hash = sha256_file(plan_path)
    results = execute_handoff_tasks(
        config=config,
        chemistry=chemistry,
        tasks=tasks,
        run_dir=run_dir,
        plan_sha256=plan_hash,
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": MANIFEST_SCHEMA,
            "stage": "validation_qualification_handoff",
            "status": (
                "invalid"
                if any(result.execution_status == "invalid" for result in results)
                else "completed"
            ),
            "completed_at": utc_now(),
            "chemistry": chemistry_identity_record(chemistry),
            "evidence_category": EVIDENCE_CATEGORY,
            "development_only": True,
            "known_site_information_used_for_selection": False,
            "formal_qualification_gate": False,
            "parallel_tasks": min(config.validation.rosetta.parallel_tasks, len(tasks)),
            "qualification_handoff_plan": file_record(plan_path, root=run_dir),
            "passed_task_count": sum(
                result.reconstruction_status == "passed" for result in results
            ),
            "failed_task_count": sum(
                result.reconstruction_status == "failed" for result in results
            ),
            "invalid_task_count": sum(
                result.execution_status == "invalid" for result in results
            ),
            "tasks": [
                handoff_result_record(result=result, run_dir=run_dir)
                for result in results
            ],
        },
    )
    return QualificationHandoffOutcome(manifest_path, len(results))
