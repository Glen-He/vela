"""阶段四工具准备与正式序列优化放行门槛。"""

from __future__ import annotations

from dataclasses import dataclass

from vela.config import AppConfig
from vela.core.provenance import sha256_file


@dataclass(frozen=True, slots=True)
class DesignIssue:
    """一个可定位的阶段四阻断条件。"""

    code: str
    message: str


@dataclass(frozen=True, slots=True)
class DesignReadiness:
    """区分工具可用性与正式科学运行放行。"""

    setup_ready: bool
    screen_ready: bool
    finalist_ready: bool
    production_ready: bool
    issues: tuple[DesignIssue, ...]


def assess_design_readiness(config: AppConfig) -> DesignReadiness:
    """检查 RosettaScripts、配体序列空间和待冻结科学决定。"""
    setup: list[DesignIssue] = []
    rosetta = config.validation.rosetta
    for code, label, path, is_directory in (
        (
            "flexpepdock_missing",
            "FlexPepDock executable",
            rosetta.executable,
            False,
        ),
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
        exists = path.is_dir() if is_directory else path.is_file()
        if not exists:
            setup.append(DesignIssue(code, f"{label} is missing: {path}"))

    scientific: list[DesignIssue] = []
    settings = config.design
    if settings.qualification_status != "qualified":
        scientific.append(
            DesignIssue(
                "method_not_qualified",
                "Stage 4 paired interface screen has not passed independent controls",
            )
        )
    if settings.objective is None:
        scientific.append(
            DesignIssue(
                "objective_unresolved",
                "Stage 4 objective must follow supported Stage 3 target states",
            )
        )
    if not settings.seeds:
        scientific.append(
            DesignIssue("seeds_unresolved", "Stage 4 independent seeds are not frozen")
        )
    if not settings.analysis.complete:
        scientific.append(
            DesignIssue(
                "analysis_rules_unresolved",
                "Stage 4 paired multi-state thresholds are not calibrated",
            )
        )
    finalist: list[DesignIssue] = []
    if not settings.finalists.seeds:
        finalist.append(
            DesignIssue(
                "finalist_seeds_unresolved",
                "Stage 4 flexible-verification seeds are not frozen",
            )
        )
    if not settings.finalists.complete:
        finalist.append(
            DesignIssue(
                "finalist_rules_unresolved",
                "Stage 4 flexible-verification thresholds are not calibrated",
            )
        )
    if not config.validation.config_complete:
        finalist.append(
            DesignIssue(
                "flexpepdock_not_qualified",
                "Stage 4 finalists require the qualified Stage 3 FlexPepDock protocol",
            )
        )
    else:
        validation_report = config.validation.qualification_report
        validation_hash = config.validation.qualification_report_sha256
        if validation_report is None or not validation_report.is_file():
            finalist.append(
                DesignIssue(
                    "flexpepdock_report_missing",
                    "Qualified Stage 3 FlexPepDock report is missing",
                )
            )
        elif (
            validation_hash is None or sha256_file(validation_report) != validation_hash
        ):
            finalist.append(
                DesignIssue(
                    "flexpepdock_report_mismatch",
                    "Stage 3 FlexPepDock report hash does not match config",
                )
            )
    report = settings.qualification_report
    expected_hash = settings.qualification_report_sha256
    if settings.qualification_status == "qualified":
        if report is None or not report.is_file():
            scientific.append(
                DesignIssue(
                    "qualification_report_missing",
                    "Qualified Stage 4 method requires an existing report",
                )
            )
        elif expected_hash is None or sha256_file(report) != expected_hash:
            scientific.append(
                DesignIssue(
                    "qualification_report_mismatch",
                    "Stage 4 qualification report hash does not match config",
                )
            )
    issues = tuple(setup + scientific + finalist)
    setup_ready = not setup
    screen_ready = setup_ready and not scientific
    finalist_ready = screen_ready and not finalist
    return DesignReadiness(
        setup_ready=setup_ready,
        screen_ready=screen_ready,
        finalist_ready=finalist_ready,
        production_ready=finalist_ready,
        issues=issues,
    )
