"""Vela 命令行语法。"""

import argparse
from pathlib import Path

DEFAULT_CONFIG_PATH = Path("configs")


def _add_config_argument(parser: argparse.ArgumentParser) -> None:
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG_PATH,
        help="Project config directory (default: configs).",
    )


def build_parser() -> argparse.ArgumentParser:
    """建立唯一 CLI 语法树。"""
    parser = argparse.ArgumentParser(
        prog="vela",
        description="Reproducible CK2 cyclic-peptide structural workflow.",
    )
    groups = parser.add_subparsers(dest="group", required=True)

    config_group = groups.add_parser("config", help="Validate resolved config.")
    config_commands = config_group.add_subparsers(dest="command", required=True)
    config_check = config_commands.add_parser("check", help="Validate project inputs.")
    _add_config_argument(config_check)

    chemistry_group = groups.add_parser(
        "chemistry", help="Manage the configured ligand chemistry."
    )
    chemistry_commands = chemistry_group.add_subparsers(dest="command", required=True)
    chemistry_record = chemistry_commands.add_parser(
        "record", help="Validate and write the configured ligand chemistry record."
    )
    _add_config_argument(chemistry_record)

    preparation_group = groups.add_parser(
        "preparation", help="Run or inspect the Stage 1 release gate."
    )
    preparation_commands = preparation_group.add_subparsers(
        dest="command", required=True
    )
    preparation_status = preparation_commands.add_parser(
        "status", help="Inspect Stage 1 production readiness without writing files."
    )
    _add_config_argument(preparation_status)
    preparation_run = preparation_commands.add_parser(
        "run", help="Run automated Stage 1 work and evaluate its release gate."
    )
    _add_config_argument(preparation_run)

    receptor_group = groups.add_parser("receptors", help="Manage receptor structures.")
    receptor_commands = receptor_group.add_subparsers(dest="command", required=True)
    receptor_download = receptor_commands.add_parser(
        "download", help="Download and verify raw RCSB inputs."
    )
    _add_config_argument(receptor_download)
    receptor_audit = receptor_commands.add_parser(
        "audit", help="Audit identities, chains, components, and missing records."
    )
    _add_config_argument(receptor_audit)
    receptor_prepare = receptor_commands.add_parser(
        "prepare", help="Build audited receptor-only base structures."
    )
    _add_config_argument(receptor_prepare)

    discovery_group = groups.add_parser(
        "discovery", help="Plan and analyze Stage 2 blind discovery."
    )
    discovery_commands = discovery_group.add_subparsers(dest="command", required=True)
    qualification_plan = discovery_commands.add_parser(
        "qualification-plan",
        help="Freeze one target's positive control and apo pilot.",
    )
    _add_config_argument(qualification_plan)
    qualification_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the qualification run."
    )
    qualification_plan.add_argument(
        "--target", required=True, help="One configured receptor target to qualify."
    )
    qualification_plan.add_argument(
        "--control-run",
        type=Path,
        help="Completed qualification run whose verified positive control is reused.",
    )
    qualification_run = discovery_commands.add_parser(
        "qualification-run", help="Execute or resume frozen qualification sampling."
    )
    _add_config_argument(qualification_run)
    qualification_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned directory below outputs/discovery/qualifications.",
    )
    qualification_analyze = discovery_commands.add_parser(
        "qualification-analyze",
        help="Evaluate native recovery and select registered site thresholds.",
    )
    _add_config_argument(qualification_analyze)
    qualification_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed directory below outputs/discovery/qualifications.",
    )
    topology_plan = discovery_commands.add_parser(
        "topology-calibration-plan",
        help="Freeze stratified all-atom disulfide-reconstruction calibration.",
    )
    _add_config_argument(topology_plan)
    topology_plan.add_argument(
        "--source-run",
        dest="topology_source_run",
        type=Path,
        required=True,
        help="Development qualification directory containing completed control seeds.",
    )
    topology_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the calibration run."
    )
    topology_run = discovery_commands.add_parser(
        "topology-calibration-run",
        help="Execute or resume frozen all-atom topology calibration tasks.",
    )
    _add_config_argument(topology_run)
    topology_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned directory below outputs/discovery/topology_calibrations.",
    )
    topology_analyze = discovery_commands.add_parser(
        "topology-calibration-analyze",
        help="Apply the preregistered topology-reconstruction qualification rule.",
    )
    _add_config_argument(topology_analyze)
    topology_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed directory below outputs/discovery/topology_calibrations.",
    )
    discovery_status = discovery_commands.add_parser(
        "status", help="Inspect whether production discovery may be planned."
    )
    _add_config_argument(discovery_status)
    discovery_status.add_argument(
        "--target",
        required=True,
        help="One configured receptor target for this complete discovery run.",
    )
    discovery_plan = discovery_commands.add_parser(
        "plan", help="Freeze one production discovery run manifest."
    )
    _add_config_argument(discovery_plan)
    discovery_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the discovery run."
    )
    discovery_plan.add_argument(
        "--target",
        required=True,
        help="One configured receptor target for this complete discovery run.",
    )
    discovery_run = discovery_commands.add_parser(
        "run", help="Execute or resume the frozen CABS-dock sampling tasks."
    )
    _add_config_argument(discovery_run)
    discovery_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned run directory below outputs/runs.",
    )
    discovery_analyze = discovery_commands.add_parser(
        "analyze", help="Analyze normalized pose evidence from a completed run."
    )
    _add_config_argument(discovery_analyze)
    discovery_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed run directory below outputs/runs.",
    )

    validation_group = groups.add_parser(
        "validation", help="Prepare and qualify Stage 3 structural validation."
    )
    validation_commands = validation_group.add_subparsers(dest="command", required=True)
    validation_status = validation_commands.add_parser(
        "status", help="Inspect Stage 3 setup and production readiness."
    )
    _add_config_argument(validation_status)
    validation_prepare = validation_commands.add_parser(
        "prepare", help="Prepare paired and receptor-only bound-state assets."
    )
    _add_config_argument(validation_prepare)
    validation_tool_check = validation_commands.add_parser(
        "tool-check", help="Verify configured Stage 3 reconstruction and Rosetta tools."
    )
    _add_config_argument(validation_tool_check)
    validation_control_plan = validation_commands.add_parser(
        "control-plan", help="Freeze a configured local-recovery control run."
    )
    _add_config_argument(validation_control_plan)
    validation_control_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the control run."
    )
    validation_control_run = validation_commands.add_parser(
        "control-run", help="Execute or resume a frozen local-recovery control run."
    )
    _add_config_argument(validation_control_run)
    validation_control_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned control directory below outputs/validation/controls.",
    )
    validation_control_analyze = validation_commands.add_parser(
        "control-analyze",
        help="Analyze completed local-recovery tasks without overwriting prior reports.",
    )
    _add_config_argument(validation_control_analyze)
    validation_control_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed control directory below outputs/validation/controls.",
    )
    validation_replication_plan = validation_commands.add_parser(
        "replication-plan",
        help="Freeze stripped bound-state full-surface replication tasks.",
    )
    _add_config_argument(validation_replication_plan)
    validation_replication_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the replication run."
    )
    validation_replication_plan.add_argument(
        "--target",
        required=True,
        help="One configured receptor target for this replication run.",
    )
    validation_replication_run = validation_commands.add_parser(
        "replication-run",
        help="Execute or resume frozen bound-state replication tasks.",
    )
    _add_config_argument(validation_replication_run)
    validation_replication_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned run directory below outputs/validation/replications.",
    )
    validation_replication_analyze = validation_commands.add_parser(
        "replication-analyze",
        help="Analyze completed bound-state replication evidence separately.",
    )
    _add_config_argument(validation_replication_analyze)
    validation_replication_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed run directory below outputs/validation/replications.",
    )
    validation_replication_compare = validation_commands.add_parser(
        "replication-compare",
        help="Compare analyzed main-discovery and bound-state sites.",
    )
    _add_config_argument(validation_replication_compare)
    validation_replication_compare.add_argument(
        "--discovery-run",
        type=Path,
        required=True,
        help="Analyzed main-discovery directory below outputs/runs.",
    )
    validation_replication_compare.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Analyzed replication directory below outputs/validation/replications.",
    )
    validation_handoff_plan = validation_commands.add_parser(
        "handoff-plan", help="Freeze supported Stage 2 poses for all-atom handoff."
    )
    _add_config_argument(validation_handoff_plan)
    validation_handoff_plan.add_argument(
        "--discovery-run",
        type=Path,
        required=True,
        help="Analyzed discovery directory below outputs/runs.",
    )
    validation_handoff_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the handoff run."
    )
    validation_handoff_plan.add_argument(
        "--candidate-id",
        action="append",
        required=True,
        help="Supported candidate ID to include; repeat for an explicit reviewed set.",
    )
    validation_handoff_run = validation_commands.add_parser(
        "handoff-run", help="Execute or resume a frozen all-atom handoff."
    )
    _add_config_argument(validation_handoff_run)
    validation_handoff_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned handoff directory below outputs/validation/handoffs.",
    )
    qualification_handoff_plan = validation_commands.add_parser(
        "qualification-handoff-plan",
        help="Freeze a development-only Top-B handoff from qualification evidence.",
    )
    _add_config_argument(qualification_handoff_plan)
    qualification_handoff_plan.add_argument(
        "--qualification-run",
        type=Path,
        required=True,
        help="Analyzed directory below outputs/discovery/qualifications.",
    )
    qualification_handoff_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the development run."
    )
    qualification_handoff_plan.add_argument(
        "--site-budget",
        type=int,
        required=True,
        help="Native-free Top-B supported-site budget to freeze.",
    )
    qualification_handoff_run = validation_commands.add_parser(
        "qualification-handoff-run",
        help="Execute or resume a frozen development-only all-atom handoff.",
    )
    _add_config_argument(qualification_handoff_run)
    qualification_handoff_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help=("Planned directory below outputs/validation/qualification_handoffs."),
    )
    qualification_refinement_plan = validation_commands.add_parser(
        "qualification-refinement-plan",
        help="Freeze a native-aware Stage 2-to-3 recovery diagnostic.",
    )
    _add_config_argument(qualification_refinement_plan)
    qualification_refinement_plan.add_argument(
        "--source-run",
        type=Path,
        required=True,
        help="Completed Top-B qualification handoff directory.",
    )
    qualification_refinement_plan.add_argument(
        "--control-run",
        type=Path,
        required=True,
        help="Completed chemistry-aware FlexPepDock control directory.",
    )
    qualification_refinement_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the diagnostic run."
    )
    qualification_refinement_run = validation_commands.add_parser(
        "qualification-refinement-run",
        help="Execute or resume the frozen native-aware recovery diagnostic.",
    )
    _add_config_argument(qualification_refinement_run)
    qualification_refinement_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned directory below outputs/validation/qualification_refinements.",
    )
    qualification_refinement_analyze = validation_commands.add_parser(
        "qualification-refinement-analyze",
        help="Analyze a completed native-aware refinement diagnostic.",
    )
    _add_config_argument(qualification_refinement_analyze)
    qualification_refinement_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed directory below outputs/validation/qualification_refinements.",
    )
    validation_guided_plan = validation_commands.add_parser(
        "guided-plan",
        help="Freeze configured experimental-site ligand threading tasks.",
    )
    _add_config_argument(validation_guided_plan)
    validation_guided_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the guided handoff."
    )
    validation_guided_run = validation_commands.add_parser(
        "guided-run", help="Build configured guided all-atom ligand starts."
    )
    _add_config_argument(validation_guided_run)
    validation_guided_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned run directory below outputs/validation/guided.",
    )
    validation_refinement_plan = validation_commands.add_parser(
        "refinement-plan", help="Freeze qualified blind or guided refinement tasks."
    )
    _add_config_argument(validation_refinement_plan)
    validation_refinement_plan.add_argument(
        "--source-run",
        type=Path,
        required=True,
        help="Completed handoff or guided directory below outputs/validation.",
    )
    validation_refinement_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the refinement run."
    )
    validation_refinement_run = validation_commands.add_parser(
        "refinement-run", help="Execute or resume frozen local-refinement tasks."
    )
    _add_config_argument(validation_refinement_run)
    validation_refinement_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned refinement directory below outputs/validation/refinements.",
    )
    validation_refinement_analyze = validation_commands.add_parser(
        "refinement-analyze", help="QC and cluster a completed refinement run."
    )
    _add_config_argument(validation_refinement_analyze)
    validation_refinement_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed refinement directory below outputs/validation/refinements.",
    )
    validation_environment_map = validation_commands.add_parser(
        "environment-map",
        help="Map supported refined representatives into full-enzyme layouts.",
    )
    _add_config_argument(validation_environment_map)
    validation_environment_map.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Analyzed refinement directory below outputs/validation/refinements.",
    )
    validation_candidate_review = validation_commands.add_parser(
        "candidate-review",
        help="Aggregate blind candidate evidence without assigning final grades.",
    )
    _add_config_argument(validation_candidate_review)
    validation_candidate_review.add_argument(
        "--discovery-run",
        type=Path,
        required=True,
        help="Analyzed main-discovery directory below outputs/runs.",
    )
    validation_candidate_review.add_argument(
        "--replication-run",
        type=Path,
        required=True,
        help="Compared replication directory below outputs/validation/replications.",
    )
    validation_candidate_review.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Analyzed blind refinement directory below outputs/validation/refinements.",
    )

    design_group = groups.add_parser(
        "design", help="Plan and run Stage 4 multi-state sequence optimization."
    )
    design_commands = design_group.add_subparsers(dest="command", required=True)
    design_status = design_commands.add_parser(
        "status", help="Inspect Stage 4 setup and production readiness."
    )
    _add_config_argument(design_status)
    design_tool_check = design_commands.add_parser(
        "tool-check", help="Verify the configured RosettaScripts design backend."
    )
    _add_config_argument(design_tool_check)
    design_single_plan = design_commands.add_parser(
        "single-plan",
        help="Freeze an exhaustive single-mutation paired screen.",
    )
    _add_config_argument(design_single_plan)
    design_single_plan.add_argument(
        "--source-run",
        dest="design_source_run",
        type=Path,
        required=True,
        help="Reviewed blind refinement directory below outputs/validation/refinements.",
    )
    design_single_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the design screen."
    )
    design_single_plan.add_argument(
        "--target-cluster",
        action="append",
        required=True,
        help="Supported Stage 3 cluster for one target subtype; repeat as needed.",
    )
    design_combination_plan = design_commands.add_parser(
        "combination-plan",
        help="Build and freeze combinations from an analyzed single screen.",
    )
    _add_config_argument(design_combination_plan)
    design_combination_plan.add_argument(
        "--source-run",
        dest="design_source_run",
        type=Path,
        required=True,
        help="Analyzed single-screen directory below outputs/design/screens.",
    )
    design_combination_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the combination screen."
    )
    design_iteration_plan = design_commands.add_parser(
        "iteration-plan",
        help="Freeze a one-edit neighborhood around reviewed parent candidates.",
    )
    _add_config_argument(design_iteration_plan)
    design_iteration_plan.add_argument(
        "--source-run",
        dest="design_source_run",
        type=Path,
        required=True,
        help="Analyzed flexible-verification run below outputs/design/finalists.",
    )
    design_iteration_plan.add_argument(
        "--parent-candidate-id",
        dest="candidate_id",
        action="append",
        required=True,
        help="Eligible MD-queued parent candidate; repeat for a same-generation set.",
    )
    design_iteration_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for the iteration screen."
    )
    design_screen_run = design_commands.add_parser(
        "screen-run", help="Execute or resume a frozen paired design screen."
    )
    _add_config_argument(design_screen_run)
    design_screen_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned run directory below outputs/design/screens.",
    )
    design_screen_analyze = design_commands.add_parser(
        "screen-analyze", help="Aggregate WT-paired multi-state screen evidence."
    )
    _add_config_argument(design_screen_analyze)
    design_screen_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed run directory below outputs/design/screens.",
    )
    design_finalist_plan = design_commands.add_parser(
        "finalist-plan",
        help="Freeze a limited flexible-verification set from an analyzed screen.",
    )
    _add_config_argument(design_finalist_plan)
    design_finalist_plan.add_argument(
        "--source-run",
        dest="design_source_run",
        type=Path,
        required=True,
        help="Analyzed single or combination run below outputs/design/screens.",
    )
    design_finalist_plan.add_argument(
        "--run-id", required=True, help="Unique identifier for flexible verification."
    )
    design_finalist_run = design_commands.add_parser(
        "finalist-run", help="Execute or resume frozen flexible-verification tasks."
    )
    _add_config_argument(design_finalist_run)
    design_finalist_run.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Planned run directory below outputs/design/finalists.",
    )
    design_finalist_analyze = design_commands.add_parser(
        "finalist-analyze",
        help="Analyze flexible ensembles and write the finite Stage 5 queue.",
    )
    _add_config_argument(design_finalist_analyze)
    design_finalist_analyze.add_argument(
        "--run-dir",
        type=Path,
        required=True,
        help="Completed run directory below outputs/design/finalists.",
    )
    return parser
