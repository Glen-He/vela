"""阶段级一键入口与科学放行状态。"""

from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import sha256_file
from vela.design.finalists.execution import run_finalists
from vela.design.finalists.planning import write_finalist_plan
from vela.design.finalists.reporting import analyze_finalists
from vela.design.models import DesignError
from vela.design.readiness import DesignReadiness, assess_design_readiness
from vela.design.screening.analysis import analyze_screen
from vela.design.screening.execution import run_screen
from vela.design.screening.planning import (
    write_combination_screen_plan,
    write_single_screen_plan,
)
from vela.design.sequence.iteration import write_iteration_screen_plan
from vela.discovery.analysis.workflow import analyze_discovery_run, discovery_run_target
from vela.discovery.models import DiscoveryError
from vela.discovery.qualification.analysis import analyze_qualification
from vela.discovery.qualification.execution import (
    run_qualification as run_discovery_qualification,
)
from vela.discovery.qualification.planning import (
    write_qualification_plan as write_discovery_qualification_plan,
)
from vela.discovery.qualification.topology import (
    analyze_topology_calibration,
    resolve_topology_calibration_run,
    run_topology_calibration,
    write_topology_calibration_plan,
)
from vela.discovery.readiness import DiscoveryReadiness, assess_discovery_readiness
from vela.discovery.sampling.execution import run_cabsdock_sampling
from vela.discovery.sampling.planning import (
    write_discovery_plan,
    write_exploration_plan,
)
from vela.preparation.chemistry import (
    chemistry_record_relative_path,
    write_chemistry_record,
)
from vela.preparation.readiness import (
    PreparationReadiness,
    assess_preparation_readiness,
)
from vela.preparation.receptors.audit import audit_receptors
from vela.preparation.receptors.cleaning import prepare_receptors
from vela.preparation.receptors.download import download_receptors
from vela.preparation.receptors.models import ReceptorError
from vela.validation.assessment.candidates import write_candidate_review
from vela.validation.assessment.environment import (
    map_refinement_environment,
    prepare_environment_references,
)
from vela.validation.bound_states.assets import prepare_bound_state_assets
from vela.validation.bound_states.comparison import compare_replication_run
from vela.validation.bound_states.controls import write_qualification_plan
from vela.validation.bound_states.qualification import (
    analyze_qualification as analyze_validation_qualification,
)
from vela.validation.bound_states.qualification import run_qualification
from vela.validation.bound_states.replication import write_replication_plan
from vela.validation.models import ValidationError
from vela.validation.readiness import (
    ValidationReadiness,
    assess_validation_readiness,
)
from vela.validation.refinement.analysis import (
    analyze_funnel_confirmation,
    analyze_funnel_deep_confirmation,
    analyze_refinement_run,
)
from vela.validation.refinement.execution import run_refinement
from vela.validation.refinement.guided import run_guided_handoff, write_guided_plan
from vela.validation.refinement.handoff_plan import (
    write_exploration_handoff_plan,
    write_funnel_screening_handoff_plan,
    write_handoff_plan,
    write_source_seed_confirmation_handoff_plan,
)
from vela.validation.refinement.handoff_run import run_handoff
from vela.validation.refinement.planning import (
    write_funnel_confirmation_plan,
    write_funnel_deep_plan,
    write_refinement_plan,
)
from vela.validation.refinement.qualification_analysis import (
    analyze_qualification_refinement,
)
from vela.validation.refinement.qualification_diagnostic import (
    run_qualification_refinement,
    write_qualification_refinement_plan,
)
from vela.validation.refinement.qualification_handoff import (
    run_qualification_handoff,
    write_qualification_handoff_plan,
)
from vela.validation.refinement.reconstruction import verify_cg2all_tool
from vela.validation.rosetta import (
    verify_flexpepdock_tool,
    verify_rosetta_scripts_tool,
)


def _print_design_readiness(readiness: DesignReadiness) -> None:
    print(f"Stage 4 setup ready: {str(readiness.setup_ready).lower()}")
    print(f"Stage 4 interface screen ready: {str(readiness.screen_ready).lower()}")
    print(
        "Stage 4 flexible finalist verification ready: "
        f"{str(readiness.finalist_ready).lower()}"
    )
    print(f"Stage 4 production ready: {str(readiness.production_ready).lower()}")
    if readiness.issues:
        print("Blocking conditions:")
        for issue in readiness.issues:
            print(f"- {issue.code}: {issue.message}")


def design_status(config: AppConfig) -> int:
    """只读报告阶段四工具与科学门槛。"""
    _print_design_readiness(assess_design_readiness(config))
    return 0


def design_tool_check(config: AppConfig) -> int:
    """核对阶段四使用的 RosettaScripts 工具身份。"""
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    flexible_tool = verify_flexpepdock_tool(config.validation.rosetta)
    print(f"RosettaScripts version: {tool.version}")
    print(f"RosettaScripts executable SHA-256: {tool.executable_sha256}")
    print(f"Design score function: {config.design.screen.score_function}")
    print(f"Design ranking score: {config.design.screen.ranking_score}")
    print(f"FlexPepDock version: {flexible_tool.version}")
    print(f"FlexPepDock executable SHA-256: {flexible_tool.executable_sha256}")
    return 0


def design_single_plan(
    *,
    config: AppConfig,
    source_run: Path,
    run_id: str,
    target_cluster_ids: tuple[str, ...],
) -> int:
    """冻结确定性单点库及其 WT 成对多状态任务。"""
    plan = write_single_screen_plan(
        config=config,
        refinement_run_dir=source_run,
        run_id=run_id,
        target_cluster_ids=target_cluster_ids,
    )
    print(f"Design single screen plan: {plan.run_dir / 'screen_plan.json'}")
    print(f"Single-mutation candidates: {len(plan.candidates)}")
    print(f"Selected Stage 3 templates: {len(plan.templates)}")
    print(f"Paired screen tasks: {len(plan.tasks)}")
    return 0


def design_combination_plan(*, config: AppConfig, source_run: Path, run_id: str) -> int:
    """根据已分析单点证据冻结有限组合重新评价任务。"""
    plan = write_combination_screen_plan(
        config=config, single_run_dir=source_run, run_id=run_id
    )
    print(f"Design combination screen plan: {plan.run_dir / 'screen_plan.json'}")
    print(f"Combination candidates: {len(plan.candidates)}")
    print(f"Paired screen tasks: {len(plan.tasks)}")
    return 0


def design_iteration_plan(
    *,
    config: AppConfig,
    source_run: Path,
    run_id: str,
    parent_candidate_ids: tuple[str, ...],
) -> int:
    """冻结多个同代父序列的一步邻域及 WT 成对初筛。"""
    plan = write_iteration_screen_plan(
        config=config,
        finalist_run_dir=source_run,
        run_id=run_id,
        parent_candidate_ids=parent_candidate_ids,
    )
    print(f"Design iteration screen plan: {plan.run_dir / 'screen_plan.json'}")
    print(f"Iteration generation: {plan.candidates[0].generation}")
    print(f"Iteration candidates selected: {len(plan.candidates)}")
    print(f"Paired screen tasks: {len(plan.tasks)}")
    return 0


def _design_run_dir(*, config: AppConfig, run_dir: Path) -> Path:
    resolved = run_dir.expanduser().resolve()
    root = (config.paths.outputs_dir / "design" / "screens").resolve()
    if not resolved.is_relative_to(root):
        raise DesignError(
            f"design screen run is outside outputs/design/screens: {resolved}"
        )
    return resolved


def design_screen_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复固定的阶段四成对界面筛查。"""
    resolved = _design_run_dir(config=config, run_dir=run_dir)
    outcome = run_screen(config=config, run_dir=resolved)
    print(f"Design screen completed: {outcome.manifest_path}")
    print(f"Completed design tasks: {outcome.task_count}")
    return 0


def design_screen_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """汇总 WT 配对分数和多状态不可补偿门槛。"""
    resolved = _design_run_dir(config=config, run_dir=run_dir)
    outcome = analyze_screen(config=config, run_dir=resolved)
    print(f"Design screen analysis: {outcome.manifest_path}")
    print(f"Candidates analyzed: {outcome.candidate_count}")
    print(f"Candidates eligible for the next design step: {outcome.eligible_count}")
    return 0


def design_finalist_plan(*, config: AppConfig, source_run: Path, run_id: str) -> int:
    """冻结有限候选、配对起点和独立柔性复核任务。"""
    plan = write_finalist_plan(config=config, screen_run_dir=source_run, run_id=run_id)
    print(f"Design finalist plan: {plan.run_dir / 'finalist_plan.json'}")
    print(f"Flexible-verification candidates: {len(plan.candidates)}")
    print(f"Flexible-verification tasks: {len(plan.tasks)}")
    return 0


def _design_finalist_dir(*, config: AppConfig, run_dir: Path) -> Path:
    resolved = run_dir.expanduser().resolve()
    root = (config.paths.outputs_dir / "design" / "finalists").resolve()
    if not resolved.is_relative_to(root):
        raise DesignError(
            f"design finalist run is outside outputs/design/finalists: {resolved}"
        )
    return resolved


def design_finalist_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复阶段四候选柔性复核。"""
    resolved = _design_finalist_dir(config=config, run_dir=run_dir)
    outcome = run_finalists(config=config, run_dir=resolved)
    print(f"Design finalist run completed: {outcome.manifest_path}")
    print(f"Completed flexible-verification tasks: {outcome.task_count}")
    return 0


def design_finalist_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """分析柔性复核并生成有限阶段五队列。"""
    resolved = _design_finalist_dir(config=config, run_dir=run_dir)
    outcome = analyze_finalists(config=config, run_dir=resolved)
    print(f"Design finalist analysis: {outcome.manifest_path}")
    print(f"Finalist candidates analyzed: {outcome.candidate_count}")
    print(f"Finalist candidates eligible: {outcome.eligible_count}")
    print(f"Stage 5 systems queued: {outcome.md_system_count}")
    return 0


def _preparation_readiness(config: AppConfig) -> PreparationReadiness:
    return assess_preparation_readiness(
        chemistry=config.chemistry,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )


def _discovery_readiness(config: AppConfig, *, target_id: str) -> DiscoveryReadiness:
    return assess_discovery_readiness(
        target_id=target_id,
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )


def _print_preparation_readiness(readiness: PreparationReadiness) -> None:
    print(f"Stage 1 preparation ready: {str(readiness.ready).lower()}")
    print(f"Receptors audited: {len(readiness.audited_receptor_ids)}")
    print(f"Receptors prepared for Stage 2: {len(readiness.prepared_receptor_ids)}")
    if readiness.issues:
        print("Blocking conditions:")
        for issue in readiness.issues:
            print(f"- {issue.code}: {issue.message}")


def preparation_status(config: AppConfig) -> int:
    """只读报告阶段一数据准备是否完整。"""
    _print_preparation_readiness(_preparation_readiness(config))
    return 0


def preparation_run(config: AppConfig) -> int:
    """顺序执行阶段一任务, 并检查可交付数据是否完整。"""
    chemistry_path = config.paths.data_dir / chemistry_record_relative_path(
        config.chemistry
    )
    write_chemistry_record(definition=config.chemistry, destination=chemistry_path)
    files = download_receptors(
        definitions=config.receptors,
        settings=config.download,
        data_dir=config.paths.data_dir,
    )
    audit = audit_receptors(
        definitions=config.receptors,
        settings=config.audit,
        data_dir=config.paths.data_dir,
    )
    failures = [item.receptor_id for item in audit if item.identity_status != "passed"]
    if failures:
        raise ReceptorError("receptor identity audit failed: " + ", ".join(failures))
    prepared = prepare_receptors(
        definitions=config.receptors,
        settings=config.preparation,
        data_dir=config.paths.data_dir,
    )
    print(f"Receptor files verified: {len(files)}")
    print(f"Receptors audited: {len(audit)}")
    print(f"Receptors prepared: {len(prepared)}")
    readiness = _preparation_readiness(config)
    _print_preparation_readiness(readiness)
    return 0 if readiness.ready else 3


def discovery_status(*, config: AppConfig, target_id: str) -> int:
    """报告阶段二是否允许规划 production。"""
    readiness = _discovery_readiness(config, target_id=target_id)
    print(f"Discovery target: {target_id}")
    print(f"Discovery production ready: {str(readiness.ready).lower()}")
    if readiness.issues:
        print("Blocking conditions:")
        for issue in readiness.issues:
            print(f"- {issue.code}: {issue.message}")
    return 0


def _qualification_run_dir(*, config: AppConfig, run_dir: Path) -> Path:
    """将阶段二资格目录限制在声明的输出根下。"""
    resolved = run_dir.expanduser().resolve()
    root = (config.paths.outputs_dir / "discovery" / "qualifications").resolve()
    if not resolved.is_relative_to(root):
        raise DiscoveryError(
            f"qualification run is outside outputs/discovery/qualifications: {resolved}"
        )
    return resolved


def discovery_qualification_plan(
    *,
    config: AppConfig,
    run_id: str,
    target_id: str,
) -> int:
    plan = write_discovery_qualification_plan(
        config=config,
        run_id=run_id,
        target_id=target_id,
    )
    print(f"Discovery qualification plan: {plan.run_dir / 'qualification_plan.json'}")
    print(f"Qualification target: {plan.target_id}")
    print(f"Qualification tasks: {len(plan.cases)}")
    return 0


def discovery_qualification_run(*, config: AppConfig, run_dir: Path) -> int:
    resolved = _qualification_run_dir(config=config, run_dir=run_dir)
    run_discovery_qualification(config=config, run_dir=resolved)
    print(f"Discovery qualification sampling completed: {resolved}")
    return 0


def discovery_qualification_analyze(*, config: AppConfig, run_dir: Path) -> int:
    resolved = _qualification_run_dir(config=config, run_dir=run_dir)
    report = analyze_qualification(config=config, run_dir=resolved)
    print(f"Discovery qualification report: {report}")
    print(f"Qualification report SHA-256: {sha256_file(report)}")
    return 0


def discovery_topology_calibration_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> int:
    """冻结开发 seed 的分层全原子二硫拓扑校准。"""
    plan = write_topology_calibration_plan(
        config=config, source_run=source_run, run_id=run_id
    )
    print(
        f"Topology calibration plan: {plan.run_dir / 'topology_calibration_plan.json'}"
    )
    print(f"Topology calibration tasks: {plan.task_count}")
    return 0


def discovery_topology_calibration_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复冻结的全原子二硫拓扑校准。"""
    resolved = resolve_topology_calibration_run(config=config, run_dir=run_dir)
    manifest = run_topology_calibration(config=config, run_dir=resolved)
    print(f"Topology calibration completed: {manifest}")
    return 0


def discovery_topology_calibration_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """按预注册门槛汇总全原子二硫拓扑校准。"""
    resolved = resolve_topology_calibration_run(config=config, run_dir=run_dir)
    report = analyze_topology_calibration(config=config, run_dir=resolved)
    print(f"Topology calibration report: {report}")
    print(f"Topology calibration report SHA-256: {sha256_file(report)}")
    return 0


def discovery_plan(*, config: AppConfig, run_id: str, target_id: str) -> int:
    plan = write_discovery_plan(config=config, run_id=run_id, target_id=target_id)
    print(f"Discovery plan written: {plan.run_dir / 'run_manifest.json'}")
    print(f"Discovery target: {plan.target_id}")
    print(f"Planned tasks: {len(plan.tasks)}")
    return 0


def discovery_exploration_plan(
    *, config: AppConfig, run_id: str, target_id: str, basis_run: Path
) -> int:
    """冻结开发性盲发现, 并明确阻止其冒充正式证据。"""
    plan = write_exploration_plan(
        config=config,
        run_id=run_id,
        target_id=target_id,
        basis_run=basis_run,
    )
    print(f"Exploration plan written: {plan.run_dir / 'run_manifest.json'}")
    print(f"Exploration target: {plan.target_id}")
    print(f"Planned tasks: {len(plan.tasks)}")
    print("Evidence scope: development_only; production claims are prohibited")
    return 0


def _discovery_run_dir(*, config: AppConfig, run_dir: Path) -> Path:
    """将阶段二运行路径限制在配置声明的 outputs/runs。"""
    resolved = run_dir.expanduser().resolve()
    runs_root = (config.paths.outputs_dir / "runs").resolve()
    if not resolved.is_relative_to(runs_root):
        raise DiscoveryError(f"discovery run is outside outputs/runs: {resolved}")
    return resolved


def discovery_run(*, config: AppConfig, run_dir: Path) -> int:
    resolved = _discovery_run_dir(config=config, run_dir=run_dir)
    run_cabsdock_sampling(config=config, run_dir=resolved)
    print(f"Discovery sampling completed: {resolved}")
    return 0


def discovery_analyze(*, config: AppConfig, run_dir: Path) -> int:
    resolved = _discovery_run_dir(config=config, run_dir=run_dir)
    target_id = discovery_run_target(resolved)
    analyze_discovery_run(
        run_dir=resolved,
        settings=config.discovery.target(target_id).analysis,
    )
    print(f"Site analysis written: {resolved / 'site_analysis'}")
    return 0


def _print_validation_readiness(readiness: ValidationReadiness) -> None:
    print(f"Stage 3 independent setup ready: {str(readiness.setup_ready).lower()}")
    print(
        "Stage 3 bound-state blind replication ready: "
        f"{str(readiness.replication_ready).lower()}"
    )
    print(
        f"Stage 3 candidate refinement ready: {str(readiness.production_ready).lower()}"
    )
    if readiness.issues:
        print("Blocking conditions:")
        for issue in readiness.issues:
            print(f"- {issue.code}: {issue.message}")


def validation_status(config: AppConfig) -> int:
    """只读报告阶段三独立准备与正式局部精修门槛。"""
    _print_validation_readiness(assess_validation_readiness(config))
    return 0


def validation_prepare(config: AppConfig) -> int:
    """制备结合态派生资产; 不启动盲搜或局部精修。"""
    assets = prepare_bound_state_assets(
        settings=config.validation,
        receptors=config.receptors,
        altloc=config.preparation.altloc,
        data_dir=config.paths.data_dir,
    )
    environments = prepare_environment_references(config)
    print(f"Bound-state assets prepared: {len(assets)}")
    for asset in assets:
        kind = "local-control" if asset.local_control_path else "site-reference-only"
        print(
            f"Prepared {asset.state_id}: {kind}; "
            f"interface_pairs={asset.interface_atom_pairs}"
        )
    print(f"Full-enzyme environment references prepared: {len(environments)}")
    _print_validation_readiness(assess_validation_readiness(config))
    return 0


def validation_tool_check(config: AppConfig) -> int:
    """运行阶段三工具身份检查, 不执行结构采样。"""
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    print(f"FlexPepDock version: {flexpepdock.version}")
    print(f"FlexPepDock executable SHA-256: {flexpepdock.executable_sha256}")
    print(f"RosettaScripts executable SHA-256: {scripts.executable_sha256}")
    print(f"cg2all version: {cg2all.version}")
    print(f"cg2all executable SHA-256: {cg2all.executable_sha256}")
    print(f"cg2all checkpoint SHA-256: {cg2all.checkpoint_sha256}")
    return 0


def validation_control_plan(*, config: AppConfig, run_id: str) -> int:
    """冻结配置所选局部恢复控制; 不执行 Rosetta 采样。"""
    plan = write_qualification_plan(config=config, run_id=run_id)
    print(
        f"Validation control plan written: {plan.run_dir / 'qualification_plan.json'}"
    )
    print(f"Planned control tasks: {len(plan.tasks)}")
    return 0


def _validation_run_dir(*, config: AppConfig, run_dir: Path, category: str) -> Path:
    resolved = run_dir.expanduser().resolve()
    validation_root = (config.paths.outputs_dir / "validation" / category).resolve()
    if not resolved.is_relative_to(validation_root):
        raise ValidationError(
            f"validation run is outside outputs/validation/{category}: {resolved}"
        )
    return resolved


def validation_control_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复冻结的通用局部恢复控制。"""
    resolved = _validation_run_dir(config=config, run_dir=run_dir, category="controls")
    outcome = run_qualification(config=config, run_dir=resolved)
    print(f"Validation control qualified: {str(outcome.qualified).lower()}")
    print(f"Qualification report: {outcome.report_path}")
    print(f"Qualification report SHA-256: {sha256_file(outcome.report_path)}")
    return 0


def validation_control_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """使用当前分析合同重新分析完整控制任务; 不覆盖历史报告。"""
    resolved = _validation_run_dir(config=config, run_dir=run_dir, category="controls")
    outcome = analyze_validation_qualification(config=config, run_dir=resolved)
    print(f"Validation control qualified: {str(outcome.qualified).lower()}")
    print(f"Qualification report: {outcome.report_path}")
    print(f"Qualification report SHA-256: {sha256_file(outcome.report_path)}")
    return 0


def validation_replication_plan(
    *, config: AppConfig, run_id: str, target_id: str
) -> int:
    """冻结配置登记的去配体结合态全表面复现任务。"""
    plan = write_replication_plan(config=config, run_id=run_id, target_id=target_id)
    print(f"Validation replication plan written: {plan.run_dir / 'run_manifest.json'}")
    print(f"Validation replication target: {plan.target_id}")
    print(f"Planned replication tasks: {len(plan.tasks)}")
    return 0


def validation_replication_run(*, config: AppConfig, run_dir: Path) -> int:
    """用已资格化的阶段二方法执行结合态受体全表面复现。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="replications"
    )
    run_cabsdock_sampling(config=config, run_dir=resolved)
    print(f"Validation replication sampling completed: {resolved}")
    return 0


def validation_replication_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """独立分析结合态受体证据, 不并入 apo 主发现支持率。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="replications"
    )
    target_id = discovery_run_target(resolved)
    analyze_discovery_run(
        run_dir=resolved,
        settings=config.discovery.target(target_id).analysis,
    )
    print(f"Validation replication analysis written: {resolved / 'site_analysis'}")
    return 0


def validation_replication_compare(
    *, config: AppConfig, discovery_run: Path, replication_run: Path
) -> int:
    """比较主发现 candidate 与结合态复现 site, 不重新聚类或分级。"""
    main = _discovery_run_dir(config=config, run_dir=discovery_run)
    replication = _validation_run_dir(
        config=config, run_dir=replication_run, category="replications"
    )
    outcome = compare_replication_run(
        config=config,
        discovery_run_dir=main,
        replication_run_dir=replication,
    )
    print(f"Validation replication comparison: {outcome.manifest_path}")
    print(f"Main candidates compared: {outcome.candidate_count}")
    print(f"Candidates with bound-state matches: {outcome.matched_candidate_count}")
    print(f"Replication-only sites: {outcome.replication_only_site_count}")
    return 0


def validation_handoff_plan(
    *,
    config: AppConfig,
    discovery_run: Path,
    run_id: str,
    candidate_ids: tuple[str, ...],
) -> int:
    """冻结证据等级预算内 blind candidate 的多 seed 全原子交接任务。"""
    source = _discovery_run_dir(config=config, run_dir=discovery_run)
    plan = write_handoff_plan(
        config=config,
        discovery_run_dir=source,
        run_id=run_id,
        candidate_ids=candidate_ids,
    )
    print(f"Validation handoff plan written: {plan.run_dir / 'handoff_plan.json'}")
    print(f"Planned handoff tasks: {len(plan.tasks)}")
    return 0


def validation_exploration_handoff_plan(
    *,
    config: AppConfig,
    discovery_run: Path,
    run_id: str,
    blind_candidate_ids: tuple[str, ...],
    functional_candidate_ids: tuple[str, ...],
) -> int:
    """冻结探索候选的分支身份、48起点和事前晋级规则。"""
    source = _discovery_run_dir(config=config, run_dir=discovery_run)
    plan = write_exploration_handoff_plan(
        config=config,
        discovery_run_dir=source,
        run_id=run_id,
        blind_candidate_ids=blind_candidate_ids,
        functional_candidate_ids=functional_candidate_ids,
    )
    print(f"Exploration handoff plan written: {plan.run_dir / 'handoff_plan.json'}")
    print(f"Blind-arm candidates: {len(blind_candidate_ids)}")
    print(f"Functional-arm candidates: {len(functional_candidate_ids)}")
    print(f"Planned handoff tasks: {len(plan.tasks)}")
    print("Evidence scope: development_only; QC metrics cannot rerank candidates")
    return 0


def validation_confirmation_handoff_plan(
    *,
    config: AppConfig,
    discovery_run: Path,
    run_id: str,
    candidate_ids: tuple[str, ...],
) -> int:
    """冻结单候选的跨 CABS source seed 确认交接。"""
    if len(candidate_ids) != 1:
        raise ValidationError("confirmation handoff requires exactly one candidate ID")
    source = _discovery_run_dir(config=config, run_dir=discovery_run)
    plan = write_source_seed_confirmation_handoff_plan(
        config=config,
        discovery_run_dir=source,
        run_id=run_id,
        candidate_id=candidate_ids[0],
    )
    print(f"Confirmation handoff plan written: {plan.run_dir / 'handoff_plan.json'}")
    print(f"Planned distinct-source-seed starts: {len(plan.tasks)}")
    print("Evidence scope: development_only; production qualification unchanged")
    return 0


def validation_funnel_handoff_plan(
    *,
    config: AppConfig,
    discovery_run: Path,
    negative_refinement_runs: tuple[Path, ...],
    run_id: str,
) -> int:
    """执行Stage 3A-0审计并冻结跨source筛选的全原子起点。"""
    source = _discovery_run_dir(config=config, run_dir=discovery_run)
    plan = write_funnel_screening_handoff_plan(
        config=config,
        discovery_run_dir=source,
        negative_refinement_runs=negative_refinement_runs,
        run_id=run_id,
    )
    print(
        f"Funnel screening handoff plan written: {plan.run_dir / 'handoff_plan.json'}"
    )
    print(
        f"Stage 3A-0 candidates selected: {len({task.candidate_id for task in plan.tasks})}"
    )
    print(f"Stage 3A screening starts planned: {len(plan.tasks)}")
    print("Rosetta tasks started: 0")
    return 0


def validation_handoff_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复冻结的 candidate 全原子重建和化学恢复。"""
    resolved = _validation_run_dir(config=config, run_dir=run_dir, category="handoffs")
    outcome = run_handoff(config=config, run_dir=resolved)
    print(f"Validation handoff completed: {outcome.manifest_path}")
    print(f"Handoff tasks completed: {outcome.task_count}")
    return 0


def validation_qualification_handoff_plan(
    *, config: AppConfig, qualification_run: Path, run_id: str, site_budget: int
) -> int:
    """冻结资格控制的 native-free Top-B 开发交接。"""
    source = _qualification_run_dir(config=config, run_dir=qualification_run)
    plan = write_qualification_handoff_plan(
        config=config,
        qualification_run_dir=source,
        run_id=run_id,
        site_budget=site_budget,
    )
    print(
        "Validation qualification handoff plan written: "
        f"{plan.run_dir / 'qualification_handoff_plan.json'}"
    )
    print(f"Planned qualification handoff tasks: {len(plan.tasks)}")
    return 0


def validation_qualification_handoff_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复开发性交接; 不把结果晋升为正式资格。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="qualification_handoffs"
    )
    outcome = run_qualification_handoff(config=config, run_dir=resolved)
    print(f"Validation qualification handoff completed: {outcome.manifest_path}")
    print(f"Qualification handoff tasks completed: {outcome.task_count}")
    return 0


def validation_qualification_refinement_plan(
    *,
    config: AppConfig,
    source_run: Path,
    control_run: Path,
    run_id: str,
    start_ids: tuple[str, ...],
    receptor_backbone_mode: str,
) -> int:
    """冻结依赖方法阳性对照的 native-aware 阶段二到三恢复诊断。"""
    source = _validation_run_dir(
        config=config, run_dir=source_run, category="qualification_handoffs"
    )
    control = _validation_run_dir(
        config=config, run_dir=control_run, category="controls"
    )
    plan = write_qualification_refinement_plan(
        config=config,
        handoff_run_dir=source,
        control_run_dir=control,
        run_id=run_id,
        start_ids=start_ids,
        receptor_backbone_mode=receptor_backbone_mode,
    )
    print(
        "Validation qualification refinement plan written: "
        f"{plan.run_dir / 'qualification_refinement_plan.json'}"
    )
    print(f"Planned qualification refinement tasks: {len(plan.tasks)}")
    return 0


def validation_qualification_refinement_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行 native-aware 开发恢复诊断; 不改变正式资格状态。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="qualification_refinements"
    )
    outcome = run_qualification_refinement(config=config, run_dir=resolved)
    print(f"Validation qualification refinement completed: {outcome.manifest_path}")
    print(f"Qualification refinement tasks completed: {outcome.task_count}")
    return 0


def validation_qualification_refinement_analyze(
    *, config: AppConfig, run_dir: Path
) -> int:
    """分析 native-aware 开发诊断; 不改变任何正式资格状态。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="qualification_refinements"
    )
    outcome = analyze_qualification_refinement(config=config, run_dir=resolved)
    print(f"Qualification refinement analysis written: {outcome.report_path}")
    print(f"Chemistry-valid decoys analyzed: {outcome.valid_decoy_count}")
    print(f"Refinement clusters analyzed: {outcome.cluster_count}")
    print(f"Recovery supported: {str(outcome.recovery_supported).lower()}")
    return 0


def validation_guided_plan(*, config: AppConfig, run_id: str) -> int:
    """冻结配置声明的实验骨架配体线程化任务。"""
    plan = write_guided_plan(config=config, run_id=run_id)
    print(f"Validation guided plan written: {plan.run_dir / 'guided_plan.json'}")
    print(f"Planned guided tasks: {len(plan.tasks)}")
    return 0


def validation_guided_run(*, config: AppConfig, run_dir: Path) -> int:
    """生成带已知位点身份的配体全原子局部精修起点。"""
    resolved = _validation_run_dir(config=config, run_dir=run_dir, category="guided")
    outcome = run_guided_handoff(config=config, run_dir=resolved)
    print(f"Validation guided handoff completed: {outcome.manifest_path}")
    print(f"Guided starts completed: {outcome.task_count}")
    return 0


def _refinement_source_dir(*, config: AppConfig, source_run: Path) -> Path:
    """只接受 blind handoff 或 guided handoff 两类局部精修来源。"""
    resolved = source_run.expanduser().resolve()
    roots = tuple(
        (config.paths.outputs_dir / "validation" / category).resolve()
        for category in ("handoffs", "guided")
    )
    if not any(resolved.is_relative_to(root) for root in roots):
        raise ValidationError(
            "refinement source is outside outputs/validation/handoffs or guided"
        )
    return resolved


def validation_refinement_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> int:
    """冻结来源身份、候选晋级、起点、seed 和局部精修预算。"""
    source = _refinement_source_dir(config=config, source_run=source_run)
    plan = write_refinement_plan(config=config, source_run_dir=source, run_id=run_id)
    print(
        f"Validation refinement plan written: {plan.run_dir / 'refinement_plan.json'}"
    )
    print(f"Selected candidates: {', '.join(plan.selected_candidate_ids)}")
    print(f"Selected all-atom starts: {plan.start_count}")
    print(f"Planned refinement tasks: {len(plan.tasks)}")
    print(f"Planned refinement decoys: {plan.total_decoy_count}")
    return 0


def validation_funnel_confirmation_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> int:
    """冻结Stage 3A命中候选的第二随机流增量任务。"""
    source = _validation_run_dir(
        config=config,
        run_dir=source_run,
        category="refinements",
    )
    plan = write_funnel_confirmation_plan(
        config=config,
        screening_run_dir=source,
        run_id=run_id,
    )
    print(f"Stage 3B plan written: {plan.run_dir / 'refinement_plan.json'}")
    print(f"Selected candidates: {', '.join(plan.selected_candidate_ids)}")
    print(f"Reused all-atom starts: {plan.start_count}")
    print(f"Planned new refinement tasks: {len(plan.tasks)}")
    print(f"Planned new refinement decoys: {plan.total_decoy_count}")
    return 0


def validation_funnel_deep_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> int:
    """冻结Stage 3B确认候选的第3和第4条随机流增量任务。"""
    source = _validation_run_dir(
        config=config,
        run_dir=source_run,
        category="refinements",
    )
    plan = write_funnel_deep_plan(
        config=config,
        confirmation_run_dir=source,
        run_id=run_id,
    )
    print(f"Stage 3C plan written: {plan.run_dir / 'refinement_plan.json'}")
    print(f"Selected candidates: {', '.join(plan.selected_candidate_ids)}")
    print(f"Reused all-atom starts: {plan.start_count}")
    print(f"Planned new refinement tasks: {len(plan.tasks)}")
    print(f"Planned new refinement decoys: {plan.total_decoy_count}")
    return 0


def validation_refinement_run(*, config: AppConfig, run_dir: Path) -> int:
    """执行或恢复冻结的正式 candidate 局部精修任务。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="refinements"
    )
    outcome = run_refinement(config=config, run_dir=resolved)
    print(f"Validation refinement completed: {outcome.manifest_path}")
    print(f"Refinement tasks completed: {outcome.task_count}")
    return 0


def validation_refinement_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """对已完成的局部精修执行配置驱动的 QC 和构象聚类。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="refinements"
    )
    outcome = analyze_refinement_run(config=config, run_dir=resolved)
    print(f"Validation refinement analysis: {outcome.manifest_path}")
    print(f"Refined decoys analyzed: {outcome.decoy_count}")
    print(f"Refined clusters: {outcome.cluster_count}")
    return 0


def validation_funnel_confirmation_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """合并Stage 3A和3B并报告满足3/4任务单元的候选。"""
    resolved = _validation_run_dir(
        config=config,
        run_dir=run_dir,
        category="refinements",
    )
    outcome = analyze_funnel_confirmation(config=config, run_dir=resolved)
    print(f"Stage 3B analysis: {outcome.manifest_path}")
    print(f"Combined refined decoys: {outcome.decoy_count}")
    print(f"Combined refined clusters: {outcome.cluster_count}")
    print(
        "Confirmed candidates: "
        + (", ".join(outcome.confirmed_candidate_ids) or "none")
    )
    return 0


def validation_funnel_deep_analyze(*, config: AppConfig, run_dir: Path) -> int:
    """合并Stage 3A至3C并报告最终原子姿态假设。"""
    resolved = _validation_run_dir(
        config=config,
        run_dir=run_dir,
        category="refinements",
    )
    outcome = analyze_funnel_deep_confirmation(config=config, run_dir=resolved)
    print(f"Stage 3C analysis: {outcome.manifest_path}")
    print(f"Combined refined decoys: {outcome.decoy_count}")
    print(f"Combined refined clusters: {outcome.cluster_count}")
    print(
        "Final atomic hypotheses: "
        + (", ".join(outcome.final_hypothesis_cluster_ids) or "none")
    )
    return 0


def validation_environment_map(*, config: AppConfig, run_dir: Path) -> int:
    """把受支持的精修代表姿态映射到配置声明的全酶实验布局。"""
    resolved = _validation_run_dir(
        config=config, run_dir=run_dir, category="refinements"
    )
    outcome = map_refinement_environment(config=config, refinement_run_dir=resolved)
    print(f"Validation environment mapping: {outcome.manifest_path}")
    print(f"Environment mappings: {outcome.mapping_count}")
    return 0


def validation_candidate_review(
    *,
    config: AppConfig,
    discovery_run: Path,
    replication_run: Path,
    refinement_run: Path,
) -> int:
    """汇总 blind candidate 的发现、复现、精修和环境事实。"""
    discovery = _discovery_run_dir(config=config, run_dir=discovery_run)
    replication = _validation_run_dir(
        config=config, run_dir=replication_run, category="replications"
    )
    refinement = _validation_run_dir(
        config=config, run_dir=refinement_run, category="refinements"
    )
    outcome = write_candidate_review(
        config=config,
        discovery_run_dir=discovery,
        replication_run_dir=replication,
        refinement_run_dir=refinement,
    )
    print(f"Validation candidate review: {outcome.manifest_path}")
    print(f"Candidates reviewed: {outcome.candidate_count}")
    return 0
