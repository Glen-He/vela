"""冻结一个靶标的阶段二资格控制与技术先导。"""

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
    is_current_vela_software,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_mapping
from vela.discovery.models import DiscoveryError, DiscoveryTask
from vela.discovery.qualification.control import control_bound_state, control_chemistry
from vela.discovery.sampling.evidence import candidate_selection_contract
from vela.discovery.sampling.planning import cabsdock_parameters
from vela.preparation.chemistry import ChemistryDefinition

CONTROL_RECOVERY = "known_complex_recovery"
TARGET_PILOT = "target_apo_pilot"


@dataclass(frozen=True, slots=True)
class QualificationCase:
    """一个资格采样任务及其独立配体和坐标参考。"""

    task: DiscoveryTask
    chemistry: ChemistryDefinition
    secondary_structure: str
    reference_receptor_id: str
    reference_path: Path
    native_pair_path: Path | None


@dataclass(frozen=True, slots=True)
class QualificationPlan:
    """已写入磁盘的阶段二资格计划。"""

    run_id: str
    target_id: str
    run_dir: Path
    cases: tuple[QualificationCase, ...]


def _method_identity(config: AppConfig) -> tuple[str, str]:
    method_id = config.discovery.method_id
    adapter_id = config.discovery.adapter_id
    if method_id is None or adapter_id is None:
        raise DiscoveryError("qualification method_id and adapter_id must be resolved")
    return method_id, adapter_id


def build_qualification_cases(
    *, config: AppConfig, target_id: str, include_control: bool = True
) -> tuple[QualificationCase, ...]:
    """按资格 seed 展开实验回收控制和当前靶标 apo 技术先导。"""
    method_id, adapter_id = _method_identity(config)
    target = config.discovery.target(target_id)
    receptor_by_id = {item.receptor_id: item for item in config.receptors}
    pilot = receptor_by_id.get(target.pilot_receptor)
    if (
        pilot is None
        or pilot.target != target_id
        or "blind_discovery" not in pilot.roles
    ):
        raise DiscoveryError(
            f"qualification pilot is not a blind-discovery receptor: {target.pilot_receptor}"
        )
    pilot_path = (
        config.paths.data_dir / "receptors" / "prepared" / f"{pilot.receptor_id}.cif"
    )
    reference_path = (
        config.paths.data_dir
        / "receptors"
        / "prepared"
        / f"{target.reference_receptor}.cif"
    )
    state = control_bound_state(config)
    control_root = (
        config.paths.data_dir / "validation" / "bound_states" / state.state_id
    )
    control_receptor = control_root / "receptor_only.cif"
    native_pair = control_root / "pair_reference.cif"
    required = (pilot_path, reference_path, control_receptor, native_pair)
    missing = [path for path in required if not path.is_file()]
    if missing:
        raise DiscoveryError(
            "qualification input is missing: " + ", ".join(map(str, missing))
        )
    control_definition = control_chemistry(state)
    control_secondary = config.discovery.qualification.control_secondary_structure
    if len(control_secondary) != len(control_definition.sequence):
        raise DiscoveryError(
            "qualification control secondary structure length is invalid"
        )
    cases: list[QualificationCase] = []
    for seed in config.discovery.qualification.seeds:
        if include_control:
            cases.append(
                QualificationCase(
                    task=DiscoveryTask(
                        task_id=f"control_{state.state_id}__seed_{seed}",
                        receptor_id=state.state_id,
                        target=config.discovery.qualification.control_target_id,
                        receptor_path=control_receptor,
                        receptor_sha256=sha256_file(control_receptor),
                        chemistry_id=control_definition.chemistry_id,
                        method_id=method_id,
                        adapter_id=adapter_id,
                        seed=seed,
                        evidence_category=CONTROL_RECOVERY,
                    ),
                    chemistry=control_definition,
                    secondary_structure=control_secondary,
                    reference_receptor_id=state.state_id,
                    reference_path=control_receptor,
                    native_pair_path=native_pair,
                )
            )
        cases.append(
            QualificationCase(
                task=DiscoveryTask(
                    task_id=f"pilot_{pilot.receptor_id}__seed_{seed}",
                    receptor_id=pilot.receptor_id,
                    target=target_id,
                    receptor_path=pilot_path,
                    receptor_sha256=sha256_file(pilot_path),
                    chemistry_id=config.chemistry.chemistry_id,
                    method_id=method_id,
                    adapter_id=adapter_id,
                    seed=seed,
                    evidence_category=TARGET_PILOT,
                ),
                chemistry=config.chemistry,
                secondary_structure=(
                    config.discovery.cabsdock.peptide_secondary_structure
                ),
                reference_receptor_id=target.reference_receptor,
                reference_path=reference_path,
                native_pair_path=None,
            )
        )
    return tuple(cases)


def _shared_control_record(
    *, config: AppConfig, control_run: Path
) -> dict[str, JsonValue]:
    resolved = control_run.expanduser().resolve()
    root = (config.paths.outputs_dir / "discovery" / "qualifications").resolve()
    if not resolved.is_relative_to(root):
        raise DiscoveryError(
            f"shared control run is outside discovery qualifications: {resolved}"
        )
    relative = resolved.relative_to(root)
    paths = {
        "qualification_plan": resolved / "qualification_plan.json",
        "qualification_sampling": resolved / "qualification_sampling.json",
        "pose_evidence": resolved / "pose_evidence.tsv",
        "baseline_pose_evidence": resolved / "baseline_pose_evidence.tsv",
        "native_recovery": resolved / "native_recovery.tsv",
        "qualification_report": resolved / "qualification_report.json",
    }
    missing = [path for path in paths.values() if not path.is_file()]
    if missing:
        raise DiscoveryError(
            "shared control run is incomplete: " + ", ".join(map(str, missing))
        )
    try:
        plan_value: object = json.loads(
            paths["qualification_plan"].read_text(encoding="utf-8")
        )
        plan = object_mapping(plan_value, name="shared qualification plan")
        report_value: object = json.loads(
            paths["qualification_report"].read_text(encoding="utf-8")
        )
        report = object_mapping(report_value, name="shared qualification report")
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError("shared qualification report is invalid") from exc
    control = object_mapping(
        report.get("control_recovery"), name="shared control recovery"
    )
    selection = object_mapping(
        control.get("candidate_selection"), name="shared candidate selection"
    )
    if (
        report.get("schema") != "vela.discovery-qualification-report/4"
        or report.get("status") != "qualified"
        or selection.get("passed") is not True
        or not is_current_vela_software(plan.get("software"))
    ):
        raise DiscoveryError("shared qualification control did not pass")
    return {
        "run_path": relative.as_posix(),
        "files": {
            name: {
                "path": path.name,
                "sha256": sha256_file(path),
            }
            for name, path in paths.items()
        },
    }


def topology_calibration_record(config: AppConfig) -> dict[str, JsonValue]:
    """冻结粗粒化拓扑到全原子化学恢复的独立资格证据。"""
    from vela.discovery.qualification.topology import (
        MANIFEST_NAME,
        MANIFEST_SCHEMA,
        PLAN_NAME,
        PLAN_SCHEMA,
        REPORT_SCHEMA,
        topology_calibration_contract,
    )

    qualification = config.discovery.qualification
    report = qualification.topology_calibration_report
    digest = qualification.topology_calibration_report_sha256
    if report is None or digest is None:
        if report is not None or digest is not None:
            raise DiscoveryError(
                "topology calibration report path and SHA-256 must resolve together"
            )
        if qualification.topology_calibration_status == "qualified":
            raise DiscoveryError(
                "qualified topology calibration requires a report and SHA-256"
            )
        return {
            "status": qualification.topology_calibration_status,
            "report": None,
        }
    if not report.is_file() or sha256_file(report) != digest:
        raise DiscoveryError("topology calibration report is missing or has changed")
    resolved_report = report.resolve()
    calibration_root = (
        config.paths.outputs_dir / "discovery" / "topology_calibrations"
    ).resolve()
    if not resolved_report.is_relative_to(calibration_root):
        raise DiscoveryError(
            "topology calibration report is outside its declared output root"
        )
    try:
        report_document: object = json.loads(
            resolved_report.read_text(encoding="utf-8")
        )
        report_record = object_mapping(
            report_document, name="topology calibration report"
        )
        tools = object_mapping(report_record.get("tools"), name="calibration tools")
        cg2all = object_mapping(tools.get("cg2all"), name="calibration cg2all")
        flexpepdock = object_mapping(
            tools.get("flexpepdock"), name="calibration FlexPepDock"
        )
        rosetta_scripts = object_mapping(
            tools.get("rosetta_scripts"), name="calibration RosettaScripts"
        )
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError("topology calibration report is invalid") from exc
    if (
        report_record.get("schema") != REPORT_SCHEMA
        or report_record.get("stage") != "disulfide_topology_calibration"
        or report_record.get("status") != qualification.topology_calibration_status
        or report_record.get("contract") != topology_calibration_contract(config)
        or report_record.get("native_information_used") is not False
        or report_record.get("evidence_scope") != "development_calibration_only"
    ):
        raise DiscoveryError(
            "topology calibration report differs from the current method contract"
        )
    calibrated_threshold = report_record.get("calibrated_max_disulfide_ca_distance_A")
    if (
        qualification.topology_calibration_status == "qualified"
        and calibrated_threshold
        != config.discovery.cabsdock.max_disulfide_ca_distance_A
    ):
        raise DiscoveryError(
            "configured CABS topology threshold differs from the calibrated threshold"
        )
    if (
        qualification.topology_calibration_status == "unqualified"
        and calibrated_threshold is not None
    ):
        raise DiscoveryError(
            "unqualified topology calibration must not recommend a threshold"
        )
    if (
        cg2all.get("version") != config.validation.cg2all.expected_version
        or cg2all.get("executable_sha256")
        != sha256_file(config.validation.cg2all.executable)
        or cg2all.get("checkpoint_sha256") != config.validation.cg2all.checkpoint_sha256
        or flexpepdock.get("version") != config.validation.rosetta.expected_version
        or flexpepdock.get("executable_sha256")
        != sha256_file(config.validation.rosetta.executable)
        or rosetta_scripts.get("version") != config.validation.rosetta.expected_version
        or rosetta_scripts.get("executable_sha256")
        != sha256_file(config.validation.rosetta.scripts_executable)
    ):
        raise DiscoveryError("topology calibration tool identity has changed")
    for key, filename, schema in (
        ("topology_calibration_plan", PLAN_NAME, PLAN_SCHEMA),
        ("topology_calibration_manifest", MANIFEST_NAME, MANIFEST_SCHEMA),
    ):
        try:
            record = object_mapping(report_record.get(key), name=key)
        except TypeError as exc:
            raise DiscoveryError(f"topology calibration {key} is invalid") from exc
        path = resolved_report.parent / filename
        if (
            record.get("path") != filename
            or record.get("sha256") != sha256_file(path)
            or _document_record(path).get("schema") != schema
        ):
            raise DiscoveryError(f"topology calibration {key} has changed")
    plan = _document_record(resolved_report.parent / PLAN_NAME)
    if not is_current_vela_software(plan.get("software")):
        raise DiscoveryError("topology calibration used a different Vela source tree")
    return {
        "status": qualification.topology_calibration_status,
        "report": {"path": report.as_posix(), "sha256": digest},
    }


def _document_record(path: Path) -> dict[str, object]:
    """读取同一校准目录内已哈希 JSON 文档。"""
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name=str(path))
    except (OSError, UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid topology calibration artifact: {path}") from exc


def qualification_decision_rules(
    *, config: AppConfig, target_id: str
) -> dict[str, JsonValue]:
    """返回在验证 seed 运行前必须完整冻结的资格判据。"""
    rules = config.discovery.qualification
    target_analysis = config.discovery.target(target_id).analysis
    if not target_analysis.complete:
        raise DiscoveryError(
            "qualification requires site analysis parameters frozen before execution"
        )
    return {
        "candidate_selection": candidate_selection_contract(config.discovery.cabsdock),
        "topology_feasibility": {
            "max_disulfide_ca_distance_A": (
                config.discovery.cabsdock.max_disulfide_ca_distance_A
            ),
            "qualification_gate": False,
            "min_models_for_selection": (
                config.discovery.cabsdock.min_models_for_selection
            ),
        },
        "max_native_ligand_rmsd_A": rules.max_native_ligand_rmsd_A,
        "min_native_receptor_contact_fraction": (
            rules.min_native_receptor_contact_fraction
        ),
        "min_successful_control_seeds": rules.min_successful_control_seeds,
        "site_analysis": {
            "parameter_selection": "frozen_before_validation_seeds",
            "contact_jaccard_distance": target_analysis.contact_jaccard_distance,
            "position_distance_A": target_analysis.position_distance_A,
            "min_seed_support": target_analysis.min_seed_support,
            "min_receptor_support": target_analysis.min_receptor_support,
        },
        "target_pilot": {
            "purpose": "technical_execution_and_descriptive_site_evidence",
            "qualification_gate": False,
        },
    }


def case_records(
    *, cases: tuple[QualificationCase, ...], data_dir: Path
) -> list[dict[str, JsonValue]]:
    """生成可比较的资格任务冻结记录。"""
    records: list[dict[str, JsonValue]] = []
    for case in cases:
        task = case.task
        records.append(
            {
                "task_id": task.task_id,
                "case": task.evidence_category,
                "target": task.target,
                "receptor_id": task.receptor_id,
                "receptor": {
                    "path": task.receptor_path.relative_to(data_dir).as_posix(),
                    "sha256": task.receptor_sha256,
                },
                "reference_receptor_id": case.reference_receptor_id,
                "reference": {
                    "path": case.reference_path.relative_to(data_dir).as_posix(),
                    "sha256": sha256_file(case.reference_path),
                },
                "native_pair": (
                    {
                        "path": case.native_pair_path.relative_to(data_dir).as_posix(),
                        "sha256": sha256_file(case.native_pair_path),
                    }
                    if case.native_pair_path is not None
                    else None
                ),
                "chemistry_id": task.chemistry_id,
                "sequence": case.chemistry.sequence,
                "disulfide_bonds": [
                    [bond.first, bond.second] for bond in case.chemistry.disulfide_bonds
                ],
                "secondary_structure": case.secondary_structure,
                "method_id": task.method_id,
                "adapter_id": task.adapter_id,
                "seed": task.seed,
                "status": "planned",
            }
        )
    return records


def write_qualification_plan(
    *,
    config: AppConfig,
    run_id: str,
    target_id: str,
    control_run: Path | None = None,
) -> QualificationPlan:
    """写出一个靶标的不可变资格计划。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise DiscoveryError(str(exc)) from exc
    run_dir = config.paths.outputs_dir / "discovery" / "qualifications" / run_id
    if run_dir.exists():
        raise DiscoveryError(f"qualification run directory already exists: {run_dir}")
    shared_control = (
        _shared_control_record(config=config, control_run=control_run)
        if control_run is not None
        else None
    )
    cases = build_qualification_cases(
        config=config,
        target_id=target_id,
        include_control=shared_control is None,
    )
    snapshot_path = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot_path, config.source_snapshot_text)
    rules = config.discovery.qualification
    atomic_write_json(
        run_dir / "qualification_plan.json",
        {
            "schema": "vela.discovery-qualification-plan/4",
            "stage": "discovery_qualification",
            "status": "planned",
            "run_id": run_id,
            "target_id": target_id,
            "planned_at": utc_now(),
            "software": vela_software_identity(),
            "known_site_information_use": {
                CONTROL_RECOVERY: "evaluation_only",
                TARGET_PILOT: False,
            },
            "shared_control": shared_control,
            "control_scope": {
                "control_target_id": rules.control_target_id,
                "requested_target_id": target_id,
                "target_matched": target_id == rules.control_target_id,
            },
            "topology_calibration": topology_calibration_record(config),
            "config_snapshot": {
                "path": snapshot_path.name,
                "sha256": sha256_file(snapshot_path),
            },
            "method_parameters": {"cabsdock": cabsdock_parameters(config)},
            "decision_rules": qualification_decision_rules(
                config=config, target_id=target_id
            ),
            "task_count": len(cases),
            "tasks": case_records(cases=cases, data_dir=config.paths.data_dir),
        },
    )
    return QualificationPlan(run_id, target_id, run_dir, cases)
