"""进入阶段二 production 前的放行规则与编排。"""

import json
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import is_current_vela_software, sha256_file
from vela.core.typed_data import object_mapping
from vela.discovery.models import (
    DiscoverySettings,
    DiscoveryTargetSettings,
    ReceptorEnsembleSettings,
)
from vela.discovery.qualification.schemas import REPORT_SCHEMA
from vela.preparation.chemistry import ChemistryDefinition
from vela.preparation.readiness import assess_preparation_readiness
from vela.preparation.receptors.models import (
    ReceptorAuditConfig,
    ReceptorDefinition,
    ReceptorPreparationConfig,
)


@dataclass(frozen=True, slots=True)
class ReadinessIssue:
    """一个阻止正式阶段二任务的独立原因。"""

    code: str
    message: str


@dataclass(frozen=True, slots=True)
class DiscoveryReadiness:
    """阶段二放行结果和正式主受体集合。"""

    target_id: str
    issues: tuple[ReadinessIssue, ...]
    receptor_ids: tuple[str, ...]

    @property
    def ready(self) -> bool:
        return not self.issues


def _method_issues(
    *,
    chemistry: ChemistryDefinition,
    settings: DiscoverySettings,
    target: DiscoveryTargetSettings,
) -> list[ReadinessIssue]:
    issues: list[ReadinessIssue] = []
    if target.qualification_status != "qualified":
        issues.append(
            ReadinessIssue(
                "method_not_qualified",
                f"{target.target_id} qualification status: "
                f"{target.qualification_status}",
            )
        )
    if settings.method_id is None:
        issues.append(
            ReadinessIssue("method_unresolved", "global method_id is unresolved")
        )
    if settings.adapter_id is None:
        issues.append(
            ReadinessIssue("adapter_unresolved", "global adapter_id is unresolved")
        )
    if not settings.seeds:
        issues.append(
            ReadinessIssue("seeds_unresolved", "production seeds are not frozen")
        )
    if not target.analysis.complete:
        issues.append(
            ReadinessIssue(
                "analysis_rules_unresolved", "site analysis thresholds are not frozen"
            )
        )
    cabsdock = settings.cabsdock
    if not cabsdock.executable.is_file():
        issues.append(
            ReadinessIssue(
                "cabsdock_executable_missing",
                f"CABS-dock executable is missing: {cabsdock.executable}",
            )
        )
    if not cabsdock.source_dir.is_dir():
        issues.append(
            ReadinessIssue(
                "cabsdock_source_missing",
                f"CABS source directory is missing: {cabsdock.source_dir}",
            )
        )
    if len(cabsdock.peptide_secondary_structure) != len(chemistry.sequence):
        issues.append(
            ReadinessIssue(
                "cabsdock_secondary_structure_mismatch",
                "CABS-dock peptide secondary structure length does not match the ligand",
            )
        )
    report = target.qualification_report
    digest = target.qualification_report_sha256
    if report is None or digest is None:
        issues.append(
            ReadinessIssue(
                "qualification_report_unresolved",
                "qualification report path or SHA-256 is unresolved",
            )
        )
    elif not report.is_file():
        issues.append(
            ReadinessIssue("qualification_report_missing", f"missing report: {report}")
        )
    elif sha256_file(report) != digest:
        issues.append(
            ReadinessIssue(
                "qualification_report_hash_mismatch",
                f"qualification report hash mismatch: {report}",
            )
        )
    else:
        try:
            raw: object = json.loads(report.read_text(encoding="utf-8"))
            document = object_mapping(raw, name="qualification report")
        except (
            UnicodeDecodeError,
            json.JSONDecodeError,
            TypeError,
        ):
            issues.append(
                ReadinessIssue(
                    "qualification_report_invalid",
                    f"invalid qualification report: {report}",
                )
            )
        else:
            identity_mismatch = (
                document.get("schema") != REPORT_SCHEMA
                or document.get("status") != target.qualification_status
                or document.get("target_id") != target.target_id
                or not is_current_vela_software(document.get("analysis_software"))
            )
            recommendation_mismatch = (
                document.get("recommended_target_config") is not None
            )
            if target.qualification_status == "qualified":
                try:
                    recommended = object_mapping(
                        document.get("recommended_target_config"),
                        name="recommended target config",
                    )
                except TypeError:
                    recommendation_mismatch = True
                else:
                    expected = {
                        "qualification_status": "qualified",
                        "contact_jaccard_distance": (
                            target.analysis.contact_jaccard_distance
                        ),
                        "position_distance_A": target.analysis.position_distance_A,
                        "min_seed_support": target.analysis.min_seed_support,
                        "min_receptor_support": target.analysis.min_receptor_support,
                    }
                    recommendation_mismatch = any(
                        recommended.get(key) != value for key, value in expected.items()
                    )
            if identity_mismatch or recommendation_mismatch:
                issues.append(
                    ReadinessIssue(
                        "qualification_report_mismatch",
                        "qualification report does not authorize the configured target rules",
                    )
                )
    return issues


def _receptor_issues(
    *,
    target_id: str,
    receptors: tuple[ReceptorDefinition, ...],
    ensemble: ReceptorEnsembleSettings,
    target_settings: DiscoveryTargetSettings | None,
    prepared_receptor_ids: tuple[str, ...],
) -> tuple[list[ReadinessIssue], tuple[str, ...]]:
    issues: list[ReadinessIssue] = []
    if target_settings is None:
        return [
            ReadinessIssue(
                "target_not_configured",
                f"discovery target is not configured: {target_id}",
            )
        ], ()
    selected = tuple(
        item
        for item in receptors
        if "blind_discovery" in item.roles and item.target == target_id
    )
    receptor_ids = tuple(item.receptor_id for item in selected)
    if not selected:
        return [
            ReadinessIssue(
                "main_receptors_missing",
                f"no blind_discovery receptors for target {target_id}",
            )
        ], ()
    if len(selected) < ensemble.min_receptors_per_target:
        issues.append(
            ReadinessIssue(
                "receptor_ensemble_incomplete",
                "blind_discovery receptor counts below "
                f"min_receptors_per_target={ensemble.min_receptors_per_target}: "
                f"{target_id}={len(selected)}",
            )
        )
    selected_by_id = {item.receptor_id: item for item in selected}
    reference_id = target_settings.reference_receptor
    reference = selected_by_id.get(reference_id)
    if reference is None:
        issues.append(
            ReadinessIssue(
                "reference_receptor_invalid",
                f"{target_id} reference receptor is not a matching blind-discovery "
                f"member: {reference_id}",
            )
        )
    prepared_ids = set(prepared_receptor_ids)
    for receptor in selected:
        if receptor.structure_state not in ensemble.allowed_structure_states:
            issues.append(
                ReadinessIssue(
                    "main_receptor_state_disallowed",
                    f"{receptor.receptor_id} state {receptor.structure_state} is not "
                    "allowed by discovery.ensemble.allowed_structure_states",
                )
            )
        if receptor.receptor_id not in prepared_ids:
            issues.append(
                ReadinessIssue(
                    "prepared_receptor_missing",
                    f"prepared receptor is absent from Stage 1: {receptor.receptor_id}",
                )
            )
    return issues, receptor_ids


def _analysis_issues(
    *,
    target_id: str,
    settings: DiscoverySettings,
    target: DiscoveryTargetSettings,
    receptors: tuple[ReceptorDefinition, ...],
) -> list[ReadinessIssue]:
    issues: list[ReadinessIssue] = []
    seed_support = target.analysis.min_seed_support
    if seed_support is not None and seed_support > len(settings.seeds):
        issues.append(
            ReadinessIssue(
                "seed_support_impossible",
                "min_seed_support exceeds the number of frozen independent seeds",
            )
        )
    receptor_support = target.analysis.min_receptor_support
    receptor_count = sum(
        "blind_discovery" in item.roles and item.target == target_id
        for item in receptors
    )
    if receptor_support is not None and receptor_support > receptor_count:
        issues.append(
            ReadinessIssue(
                "receptor_support_impossible",
                "min_receptor_support exceeds available conformations: "
                f"{target_id}={receptor_count}",
            )
        )
    return issues


def assess_discovery_readiness(
    *,
    target_id: str,
    chemistry: ChemistryDefinition,
    settings: DiscoverySettings,
    receptors: tuple[ReceptorDefinition, ...],
    audit: ReceptorAuditConfig,
    preparation: ReceptorPreparationConfig,
    data_dir: Path,
) -> DiscoveryReadiness:
    """汇总化学、阶段一产物、方法、分析规则和主受体门槛。"""
    target_settings = next(
        (target for target in settings.targets if target.target_id == target_id), None
    )
    preparation_readiness = assess_preparation_readiness(
        chemistry=chemistry,
        receptors=receptors,
        audit=audit,
        preparation=preparation,
        data_dir=data_dir,
    )
    receptor_problems, receptor_ids = _receptor_issues(
        target_id=target_id,
        receptors=receptors,
        ensemble=settings.ensemble,
        target_settings=target_settings,
        prepared_receptor_ids=preparation_readiness.prepared_receptor_ids,
    )
    issues = (
        [
            ReadinessIssue(issue.code, issue.message)
            for issue in preparation_readiness.issues
        ]
        + (
            _method_issues(
                chemistry=chemistry,
                settings=settings,
                target=target_settings,
            )
            if target_settings is not None
            else []
        )
        + receptor_problems
        + (
            _analysis_issues(
                target_id=target_id,
                settings=settings,
                target=target_settings,
                receptors=receptors,
            )
            if target_settings is not None
            else []
        )
    )
    unique = {(issue.code, issue.message): issue for issue in issues}
    return DiscoveryReadiness(target_id, tuple(unique.values()), receptor_ids)


__all__ = ["DiscoveryReadiness", "ReadinessIssue", "assess_discovery_readiness"]
