"""根据冻结输入建立阶段二全表面采样任务计划。"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

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
from vela.core.typed_data import object_mapping
from vela.discovery.analysis.clustering import candidate_analysis_contract
from vela.discovery.models import DiscoveryError, DiscoveryTask
from vela.discovery.readiness import (
    DiscoveryReadiness,
    assess_discovery_exploration_readiness,
    assess_discovery_readiness,
)
from vela.discovery.sampling.cabsdock import cabsdock_source_records
from vela.discovery.sampling.evidence import candidate_selection_contract
from vela.discovery.sampling.materialization import materializer_record
from vela.preparation.chemistry import chemistry_record_relative_path

MAIN_DISCOVERY_EVIDENCE = "main_discovery"
EXPLORATORY_DISCOVERY_EVIDENCE = "exploratory_discovery"


def cabsdock_parameters(config: AppConfig) -> dict[str, JsonValue]:
    """将阶段二实际使用的 CABS-dock 参数写入冻结清单。"""
    settings = config.discovery.cabsdock
    return {
        "executable": settings.executable.as_posix(),
        "executable_sha256": sha256_file(settings.executable),
        "source_dir": settings.source_dir.as_posix(),
        "source_revision": settings.source_revision,
        "critical_source_files": cabsdock_source_records(settings),
        "trajectory_materializer": materializer_record(),
        "disulfide_ca_restraint": {
            "representation": "CABS_CA_pseudoatom_distance_restraint",
            "distance_A": settings.disulfide_ca_restraint_distance_A,
            "weight": settings.disulfide_ca_restraint_weight,
        },
        "seed_workers": settings.seed_workers,
        "peptide_secondary_structure": settings.peptide_secondary_structure,
        "mc_annealing": settings.mc_annealing,
        "mc_cycles": settings.mc_cycles,
        "mc_steps": settings.mc_steps,
        "replicas": settings.replicas,
        "replicas_dtemp": settings.replicas_dtemp,
        "temperature_initial": settings.temperature_initial,
        "temperature_final": settings.temperature_final,
        "binding_interactions": settings.binding_interactions,
        "protein_restraints": {
            "mode": "rigid",
            "gap": settings.protein_restraint_gap,
            "min_A": settings.protein_restraint_min_A,
            "max_A": settings.protein_restraint_max_A,
        },
        "filtering_count": settings.filtering_count,
        "filtering_mode": "each_replica",
        "clustering_medoids": settings.clustering_medoids,
        "clustering_iterations": settings.clustering_iterations,
        "trajectory_contact_ca_threshold_A": (
            settings.trajectory_contact_ca_threshold_A
        ),
        "max_reconstructable_disulfide_ca_distance_A": (
            settings.max_reconstructable_disulfide_ca_distance_A
        ),
        "min_models_for_selection": settings.min_models_for_selection,
        "candidate_selection": candidate_selection_contract(settings),
        "all_atom_reconstruction": False,
    }


@dataclass(frozen=True, slots=True)
class DiscoveryPlan:
    """已写入磁盘的阶段二任务计划。"""

    run_id: str
    target_id: str
    run_dir: Path
    tasks: tuple[DiscoveryTask, ...]


def _required_text(value: str | None, *, name: str) -> str:
    if value is None:
        raise DiscoveryError(f"{name} is unresolved")
    return value


def _build_tasks(
    config: AppConfig,
    *,
    target_id: str,
    evidence_category: str,
    readiness: DiscoveryReadiness,
) -> tuple[DiscoveryTask, ...]:
    """按已核验的受体集合展开独立 seed 任务。"""
    if not readiness.ready:
        raise DiscoveryError(
            "discovery is not ready: "
            + "; ".join(f"{item.code}: {item.message}" for item in readiness.issues)
        )
    method_id = _required_text(config.discovery.method_id, name="method_id")
    adapter_id = _required_text(config.discovery.adapter_id, name="adapter_id")
    selected = {
        item.receptor_id: item
        for item in config.receptors
        if item.receptor_id in readiness.receptor_ids
    }
    tasks: list[DiscoveryTask] = []
    for seed in config.discovery.seeds:
        for receptor_id in readiness.receptor_ids:
            receptor = selected[receptor_id]
            receptor_path = (
                config.paths.data_dir
                / "receptors"
                / "prepared"
                / f"{receptor.receptor_id}.cif"
            )
            tasks.append(
                DiscoveryTask(
                    task_id=f"{receptor.receptor_id}__seed_{seed}",
                    receptor_id=receptor.receptor_id,
                    target=receptor.target,
                    receptor_path=receptor_path,
                    receptor_sha256=sha256_file(receptor_path),
                    chemistry_id=config.chemistry.chemistry_id,
                    method_id=method_id,
                    adapter_id=adapter_id,
                    seed=seed,
                    evidence_category=evidence_category,
                )
            )
    return tuple(tasks)


def build_tasks(config: AppConfig, *, target_id: str) -> tuple[DiscoveryTask, ...]:
    """为一个正式目标的受体构象展开独立 seed 任务。"""
    readiness = assess_discovery_readiness(
        target_id=target_id,
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )
    return _build_tasks(
        config,
        target_id=target_id,
        evidence_category=MAIN_DISCOVERY_EVIDENCE,
        readiness=readiness,
    )


def _production_authorization(
    *, config: AppConfig, target_id: str
) -> dict[str, JsonValue]:
    """冻结正式任务实际依赖的资格报告。"""
    target = config.discovery.target(target_id)
    report = target.qualification_report
    if report is None or not report.is_file():
        raise DiscoveryError("qualification_report is unresolved or missing")
    digest = sha256_file(report)
    if digest != target.qualification_report_sha256:
        raise DiscoveryError("qualification_report hash differs from the config")
    return {
        "scope": "production",
        "production_qualified": True,
        "qualification_report": {
            "path": report.as_posix(),
            "sha256": digest,
        },
    }


def _exploration_authorization(
    *, config: AppConfig, target_id: str, basis_run: Path
) -> dict[str, JsonValue]:
    """记录开发性运行的历史依据和不可用于正式结论的边界。"""
    resolved = basis_run.expanduser().resolve()
    root = (config.paths.outputs_dir / "discovery" / "qualifications").resolve()
    if not resolved.is_relative_to(root):
        raise DiscoveryError(
            "exploration basis run is outside outputs/discovery/qualifications: "
            f"{resolved}"
        )
    required_names = (
        "qualification_plan.json",
        "qualification_sampling.json",
        "qualification_report.json",
        "pose_evidence.tsv",
    )
    missing = tuple(name for name in required_names if not (resolved / name).is_file())
    if missing:
        raise DiscoveryError(
            "exploration basis run is incomplete: " + ", ".join(missing)
        )
    try:
        raw_plan: object = json.loads(
            (resolved / "qualification_plan.json").read_text(encoding="utf-8")
        )
        raw_report: object = json.loads(
            (resolved / "qualification_report.json").read_text(encoding="utf-8")
        )
        plan = object_mapping(raw_plan, name="exploration basis plan")
        report = object_mapping(raw_report, name="exploration basis report")
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError("exploration basis records are invalid") from exc
    if plan.get("target_id") != target_id or report.get("target_id") != target_id:
        raise DiscoveryError(
            "exploration basis target differs from the requested target"
        )
    recorded_status = report.get("status")
    if not isinstance(recorded_status, str) or not recorded_status.strip():
        raise DiscoveryError("exploration basis report status is invalid")
    artifacts: dict[str, JsonValue] = {
        name: {
            "path": (resolved / name).as_posix(),
            "sha256": sha256_file(resolved / name),
        }
        for name in required_names
    }
    basis_record: dict[str, JsonValue] = {
        "path": resolved.as_posix(),
        "recorded_status": recorded_status,
        "artifacts": artifacts,
    }
    authorization: dict[str, JsonValue] = {
        "scope": "development_only",
        "production_qualified": False,
        "independent_holdout_qualified": False,
        "current_target_qualification_status": (
            config.discovery.target(target_id).qualification_status
        ),
        "evidence_limit": (
            "post_hoc_candidate_contract_not_independently_qualified; "
            "results_must_not_be_reported_as_production_discovery"
        ),
        "basis_run": basis_record,
    }
    return authorization


def production_authorization(
    *, config: AppConfig, target_id: str
) -> dict[str, JsonValue]:
    """供其他正式全局采样入口复用同一资格事实。"""
    return _production_authorization(config=config, target_id=target_id)


def write_sampling_plan(
    *,
    config: AppConfig,
    run_id: str,
    target_id: str,
    run_dir: Path,
    tasks: tuple[DiscoveryTask, ...],
    evidence_category: str,
    receptor_selection: dict[str, JsonValue],
    method_authorization: dict[str, JsonValue],
    additional_inputs: dict[str, JsonValue] | None = None,
) -> DiscoveryPlan:
    """把任意已放行的全表面受体集合冻结为统一采样计划。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise DiscoveryError(str(exc)) from exc
    resolved_run_dir = run_dir.resolve()
    try:
        resolved_run_dir.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise DiscoveryError(
            f"sampling run directory is outside outputs: {resolved_run_dir}"
        ) from exc
    if resolved_run_dir.name != run_id:
        raise DiscoveryError("sampling run directory name must equal run_id")
    if run_dir.exists():
        raise DiscoveryError(f"run directory already exists: {run_dir}")
    if not tasks:
        raise DiscoveryError("sampling plan must contain at least one task")
    if not target_id.strip() or {task.target for task in tasks} != {target_id}:
        raise DiscoveryError("sampling plan must contain exactly one declared target")
    if not evidence_category.strip() or any(
        task.evidence_category != evidence_category for task in tasks
    ):
        raise DiscoveryError("sampling tasks must share the declared evidence category")
    authorization_scope = method_authorization.get("scope")
    expected_scope = (
        "development_only"
        if evidence_category == EXPLORATORY_DISCOVERY_EVIDENCE
        else "production"
    )
    if authorization_scope != expected_scope:
        raise DiscoveryError(
            "sampling evidence category and method authorization scope disagree"
        )
    expected_qualified = expected_scope == "production"
    if method_authorization.get("production_qualified") is not expected_qualified:
        raise DiscoveryError("method authorization qualification flag is invalid")
    task_ids = tuple(task.task_id for task in tasks)
    if len(task_ids) != len(set(task_ids)):
        raise DiscoveryError("sampling task IDs must be unique")
    chemistry_path = config.paths.data_dir / chemistry_record_relative_path(
        config.chemistry
    )
    if not chemistry_path.is_file():
        raise DiscoveryError(f"chemistry record does not exist: {chemistry_path}")
    target_settings = config.discovery.target(target_id)
    analysis = target_settings.analysis
    manifest_tasks: list[dict[str, JsonValue]] = []
    for task in tasks:
        try:
            receptor_relative = task.receptor_path.resolve().relative_to(
                config.paths.data_dir.resolve()
            )
        except ValueError as exc:
            raise DiscoveryError(
                f"sampling receptor is outside the data directory: {task.receptor_path}"
            ) from exc
        if (
            not task.receptor_path.is_file()
            or sha256_file(task.receptor_path) != task.receptor_sha256
            or task.chemistry_id != config.chemistry.chemistry_id
            or task.method_id != config.discovery.method_id
            or task.adapter_id != config.discovery.adapter_id
            or task.seed not in config.discovery.seeds
        ):
            raise DiscoveryError(
                f"sampling task differs from the current frozen config: {task.task_id}"
            )
        manifest_tasks.append(
            {
                "task_id": task.task_id,
                "receptor_id": task.receptor_id,
                "target": task.target,
                "receptor": {
                    "path": receptor_relative.as_posix(),
                    "sha256": task.receptor_sha256,
                },
                "chemistry_id": task.chemistry_id,
                "method_id": task.method_id,
                "adapter_id": task.adapter_id,
                "seed": task.seed,
                "evidence_category": task.evidence_category,
                "status": "planned",
            }
        )
    snapshot_path = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot_path, config.source_snapshot_text)
    input_records: dict[str, JsonValue] = {
        "config_snapshot": {
            "path": snapshot_path.name,
            "sha256": sha256_file(snapshot_path),
        },
        "chemistry_record": {
            "path": chemistry_path.relative_to(config.paths.data_dir).as_posix(),
            "sha256": sha256_file(chemistry_path),
        },
        "method_authorization": method_authorization,
    }
    if additional_inputs is not None:
        overlap = set(input_records) & set(additional_inputs)
        if overlap:
            raise DiscoveryError(
                "additional sampling inputs duplicate reserved names: "
                + ", ".join(sorted(overlap))
            )
        input_records.update(additional_inputs)
    atomic_write_json(
        run_dir / "run_manifest.json",
        {
            "schema": "vela.discovery-run-manifest/8",
            "run_id": run_id,
            "target_id": target_id,
            "stage": "discovery",
            "status": "planned",
            "planned_at": utc_now(),
            "software": vela_software_identity(),
            "evidence_category": evidence_category,
            "known_site_information_used": False,
            "site_definition_frozen": True,
            "inputs": input_records,
            "analysis_contract": candidate_analysis_contract(analysis),
            "method_parameters": {
                "cabsdock": cabsdock_parameters(config),
                "chemical_fidelity": (
                    "canonical_sequence_and_disulfide_topology; "
                    "terminal_charge_and_amide_state_deferred_to_stage_3"
                ),
            },
            "receptor_selection": receptor_selection,
            "adapter_handoff": {
                "sampling_manifest": "sampling_manifest.json",
                "pose_evidence": "pose_evidence.tsv",
            },
            "task_count": len(tasks),
            "tasks": manifest_tasks,
        },
    )
    return DiscoveryPlan(run_id, target_id, run_dir, tasks)


def write_discovery_plan(
    *, config: AppConfig, run_id: str, target_id: str
) -> DiscoveryPlan:
    """创建阶段二主发现目录并保存配置快照和不可变任务身份。"""
    tasks = build_tasks(config, target_id=target_id)
    ensemble = config.discovery.ensemble
    return write_sampling_plan(
        config=config,
        run_id=run_id,
        target_id=target_id,
        run_dir=config.paths.outputs_dir / "runs" / run_id,
        tasks=tasks,
        evidence_category=MAIN_DISCOVERY_EVIDENCE,
        receptor_selection={
            "mode": "configured_receptor_role",
            "role": "blind_discovery",
            "target_id": target_id,
            "min_receptors_per_target": ensemble.min_receptors_per_target,
            "allowed_structure_states": list(ensemble.allowed_structure_states),
            "reference_receptor": config.discovery.target(target_id).reference_receptor,
        },
        method_authorization=_production_authorization(
            config=config, target_id=target_id
        ),
    )


def write_exploration_plan(
    *, config: AppConfig, run_id: str, target_id: str, basis_run: Path
) -> DiscoveryPlan:
    """冻结一个不会被误认为正式证据的双受体开发性发现任务。"""
    readiness = assess_discovery_exploration_readiness(
        target_id=target_id,
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )
    tasks = _build_tasks(
        config,
        target_id=target_id,
        evidence_category=EXPLORATORY_DISCOVERY_EVIDENCE,
        readiness=readiness,
    )
    ensemble = config.discovery.ensemble
    return write_sampling_plan(
        config=config,
        run_id=run_id,
        target_id=target_id,
        run_dir=config.paths.outputs_dir / "runs" / run_id,
        tasks=tasks,
        evidence_category=EXPLORATORY_DISCOVERY_EVIDENCE,
        receptor_selection={
            "mode": "configured_receptor_role",
            "role": "blind_discovery",
            "target_id": target_id,
            "min_receptors_per_target": ensemble.min_receptors_per_target,
            "allowed_structure_states": list(ensemble.allowed_structure_states),
            "reference_receptor": config.discovery.target(target_id).reference_receptor,
        },
        method_authorization=_exploration_authorization(
            config=config,
            target_id=target_id,
            basis_run=basis_run,
        ),
    )
