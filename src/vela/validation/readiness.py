"""阶段三独立准备和正式局部精修放行门槛。"""

from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    is_vela_software_identity,
    sha256_file,
)
from vela.core.typed_data import object_mapping
from vela.discovery.readiness import assess_discovery_readiness
from vela.validation.bound_states.schemas import REPORT_SCHEMA


@dataclass(frozen=True, slots=True)
class ValidationIssue:
    """一个可定位的阶段三阻断条件。"""

    code: str
    message: str


@dataclass(frozen=True, slots=True)
class ValidationReadiness:
    """区分独立资产准备与正式候选精修的放行状态。"""

    setup_ready: bool
    replication_ready: bool
    production_ready: bool
    issues: tuple[ValidationIssue, ...]


def _asset_paths(config: AppConfig) -> tuple[Path, ...]:
    root = config.paths.data_dir / "validation" / "bound_states"
    paths: list[Path] = [root / "preparation_manifest.json"]
    for state in config.validation.bound_states:
        state_dir = root / state.state_id
        paths.extend(
            [state_dir / "receptor_only.cif", state_dir / "pair_reference.cif"]
        )
        if state.local_control_kind == "standard_cyclic_peptide":
            paths.extend(
                [
                    state_dir / "flexpepdock_control.pdb",
                    state_dir / "fix_disulfide.txt",
                ]
            )
    environment_root = config.paths.data_dir / "validation" / "environments"
    paths.append(environment_root / "preparation_manifest.json")
    paths.extend(
        environment_root / reference.reference_id / "assembly_reference.cif"
        for reference in config.validation.environment_references
    )
    return tuple(paths)


def _qualification_report_issue(config: AppConfig) -> ValidationIssue | None:
    """校验 matched local-control 报告的身份、范围与不可变内容。"""
    report = config.validation.qualification_report
    digest = config.validation.qualification_report_sha256
    if report is None or not report.is_file():
        return ValidationIssue(
            "qualification_report_missing",
            "Qualified Stage 3 method requires an existing report",
        )
    if digest is None or sha256_file(report) != digest:
        return ValidationIssue(
            "qualification_report_mismatch",
            "Stage 3 qualification report hash does not match config",
        )
    try:
        raw: object = json.loads(report.read_text(encoding="utf-8"))
        document = object_mapping(raw, name="validation qualification report")
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError):
        return ValidationIssue(
            "qualification_report_invalid",
            "Stage 3 qualification report is invalid",
        )
    if (
        document.get("schema") != REPORT_SCHEMA
        or document.get("stage") != "validation_qualification"
        or document.get("status") != "qualified"
        or document.get("method_id") != config.validation.method_id
        or document.get("evidence_category") != "method_positive_control"
        or document.get("ligand_candidate_evidence") is not False
        or not is_vela_software_identity(document.get("sampling_software"))
        or not is_vela_software_identity(document.get("analysis_software"))
    ):
        return ValidationIssue(
            "qualification_report_mismatch",
            "Stage 3 qualification report does not match the local-control contract",
        )
    return None


def assess_validation_readiness(config: AppConfig) -> ValidationReadiness:
    """检查当前可自行完成的阶段三资产、工具路径和科学门槛。"""
    asset_issues: list[ValidationIssue] = []
    receptor_by_id = {item.receptor_id: item for item in config.receptors}
    for state in config.validation.bound_states:
        receptor = receptor_by_id.get(state.receptor_id)
        if receptor is None:
            asset_issues.append(
                ValidationIssue(
                    "bound_state_registry_invalid",
                    f"Unknown receptor for {state.state_id}: {state.receptor_id}",
                )
            )
            continue
        raw_path = (
            config.paths.data_dir / "receptors" / "raw" / f"{receptor.pdb_id}.cif"
        )
        if not raw_path.is_file():
            asset_issues.append(
                ValidationIssue(
                    "bound_state_raw_missing",
                    f"Raw structure is missing for {state.state_id}: {raw_path}",
                )
            )
    for reference in config.validation.environment_references:
        receptor = receptor_by_id.get(reference.receptor_id)
        if receptor is None:
            asset_issues.append(
                ValidationIssue(
                    "environment_registry_invalid",
                    f"Unknown receptor for {reference.reference_id}: {reference.receptor_id}",
                )
            )
            continue
        raw_path = (
            config.paths.data_dir / "receptors" / "raw" / f"{receptor.pdb_id}.cif"
        )
        if not raw_path.is_file():
            asset_issues.append(
                ValidationIssue(
                    "environment_raw_missing",
                    f"Raw structure is missing for {reference.reference_id}: {raw_path}",
                )
            )
    tool_issues: list[ValidationIssue] = []
    rosetta = config.validation.rosetta
    for code, label, path, directory in (
        ("flexpepdock_missing", "FlexPepDock executable", rosetta.executable, False),
        (
            "rosetta_scripts_missing",
            "RosettaScripts executable",
            rosetta.scripts_executable,
            False,
        ),
        ("rosetta_database_missing", "Rosetta database", rosetta.database, True),
        (
            "rosetta_version_file_missing",
            "Rosetta version declaration",
            rosetta.version_file,
            False,
        ),
    ):
        exists = path.is_dir() if directory else path.is_file()
        if not exists:
            tool_issues.append(ValidationIssue(code, f"{label} is missing: {path}"))
    cg2all = config.validation.cg2all
    for code, label, path in (
        ("cg2all_missing", "cg2all executable", cg2all.executable),
        ("cg2all_metadata_missing", "cg2all package metadata", cg2all.package_metadata),
        ("cg2all_checkpoint_missing", "cg2all checkpoint", cg2all.checkpoint),
    ):
        if not path.is_file():
            tool_issues.append(ValidationIssue(code, f"{label} is missing: {path}"))
    missing_assets = [path for path in _asset_paths(config) if not path.is_file()]
    if missing_assets:
        asset_issues.append(
            ValidationIssue(
                "bound_state_assets_missing",
                f"Stage 3 bound-state assets are incomplete ({len(missing_assets)} missing)",
            )
        )

    method_issues: list[ValidationIssue] = []
    if config.validation.qualification_status != "qualified":
        method_issues.append(
            ValidationIssue(
                "method_not_qualified",
                "FlexPepDock local-refinement method has not passed its positive control",
            )
        )
    if not config.validation.seeds:
        method_issues.append(
            ValidationIssue(
                "seeds_unresolved", "Stage 3 production seeds are not frozen"
            )
        )
    if not config.validation.analysis.complete:
        method_issues.append(
            ValidationIssue(
                "analysis_rules_unresolved",
                "Stage 3 QC and clustering thresholds are not calibrated",
            )
        )
    if config.validation.qualification_status == "qualified":
        report_issue = _qualification_report_issue(config)
        if report_issue is not None:
            method_issues.append(report_issue)
    replication_issues = [
        ValidationIssue(f"global_{issue.code}", f"{target_id}: {issue.message}")
        for target_id in (target.target_id for target in config.discovery.targets)
        for issue in assess_discovery_readiness(
            target_id=target_id,
            chemistry=config.chemistry,
            settings=config.discovery,
            receptors=config.receptors,
            audit=config.audit,
            preparation=config.preparation,
            data_dir=config.paths.data_dir,
        ).issues
    ]
    all_issues = asset_issues + tool_issues + method_issues + replication_issues
    issues = tuple(
        {(issue.code, issue.message): issue for issue in all_issues}.values()
    )
    return ValidationReadiness(
        setup_ready=not asset_issues and not tool_issues,
        replication_ready=not asset_issues and not replication_issues,
        production_ready=not asset_issues and not tool_issues and not method_issues,
        issues=issues,
    )
