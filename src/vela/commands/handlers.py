"""将 CLI 命令适配到阶段工作流。"""

from pathlib import Path

from vela.commands.stages import (
    design_combination_plan,
    design_finalist_analyze,
    design_finalist_plan,
    design_finalist_run,
    design_iteration_plan,
    design_screen_analyze,
    design_screen_run,
    design_single_plan,
    design_status,
    design_tool_check,
    discovery_analyze,
    discovery_plan,
    discovery_qualification_analyze,
    discovery_qualification_plan,
    discovery_qualification_run,
    discovery_run,
    discovery_status,
    discovery_topology_calibration_analyze,
    discovery_topology_calibration_plan,
    discovery_topology_calibration_run,
    preparation_run,
    preparation_status,
    validation_candidate_review,
    validation_control_plan,
    validation_control_run,
    validation_environment_map,
    validation_guided_plan,
    validation_guided_run,
    validation_handoff_plan,
    validation_handoff_run,
    validation_prepare,
    validation_refinement_analyze,
    validation_refinement_plan,
    validation_refinement_run,
    validation_replication_analyze,
    validation_replication_compare,
    validation_replication_plan,
    validation_replication_run,
    validation_status,
    validation_tool_check,
)
from vela.config import AppConfig
from vela.preparation.chemistry import (
    chemistry_record_relative_path,
    write_chemistry_record,
)
from vela.preparation.receptors.audit import audit_receptors
from vela.preparation.receptors.cleaning import prepare_receptors
from vela.preparation.receptors.download import download_receptors
from vela.preparation.receptors.models import ReceptorError


def _config_check(config: AppConfig) -> int:
    apo_blind = {
        item.receptor_id for item in config.receptors if "blind_discovery" in item.roles
    }
    bound_blind = {
        item.receptor_id
        for item in config.receptors
        if "bound_state_blind_replication" in item.roles
    }
    print(f"Config valid: {config.source_dir}")
    print(f"Data directory: {config.paths.data_dir}")
    print(
        f"Ligand: {config.chemistry.ligand_id}; "
        f"chemistry_id={config.chemistry.chemistry_id}"
    )
    print(
        "Download policy: "
        f"retries={config.download.retries}, "
        f"timeout_seconds={config.download.timeout_seconds:g}, "
        f"backoff_multiplier={config.download.backoff_multiplier:g}, "
        f"chunk_size_bytes={config.download.chunk_size_bytes}"
    )
    contacts = config.audit.crystal_contacts
    print(
        "Crystal contact audit: "
        f"distance_A={contacts.distance_A:g}, "
        f"min_occupancy={contacts.min_occupancy:g}, "
        f"include_hydrogens={str(contacts.include_hydrogens).lower()}"
    )
    print(f"Altloc preferred label: {config.preparation.altloc.preferred_label}")
    print(
        "Discovery ensemble: "
        f"targets={','.join(target.target_id for target in config.discovery.targets)}, "
        "min_receptors_per_target="
        f"{config.discovery.ensemble.min_receptors_per_target}, "
        "states="
        f"{','.join(config.discovery.ensemble.allowed_structure_states)}"
    )
    cabsdock = config.discovery.cabsdock
    print(f"CABS-dock executable: {cabsdock.executable}")
    print(
        "CABS-dock source: "
        f"{cabsdock.source_dir} @ {cabsdock.source_revision}; "
        f"patch={cabsdock.patch_file}"
    )
    print(
        "CABS-dock sampling: "
        f"formal_seeds={len(config.discovery.seeds)}, "
        f"seed_workers={cabsdock.seed_workers}, "
        f"annealing={cabsdock.mc_annealing}, "
        f"cycles={cabsdock.mc_cycles}, "
        f"steps={cabsdock.mc_steps}, "
        f"replicas={cabsdock.replicas}, "
        f"replicas_dtemp={cabsdock.replicas_dtemp:g}, "
        "temperature="
        f"{cabsdock.temperature_initial:g}->{cabsdock.temperature_final:g}, "
        f"filtered={cabsdock.filtering_count}, "
        "selection_pool=complete_traf, "
        "trajectory_contact_CA_A="
        f"{cabsdock.trajectory_contact_ca_threshold_A:g}, "
        f"topology_CA_threshold_A={cabsdock.max_disulfide_ca_distance_A:g}, "
        f"baseline_medoids={cabsdock.clustering_medoids}, "
        "selection=site_first_then_pose, "
        f"min_selection_models={cabsdock.min_models_for_selection}, "
        "pose_cluster_coverage=complete"
    )
    print("CABS-dock worker unit: seed batch; receptor conformations run sequentially")
    print(
        "Discovery coordinate frames: "
        + ", ".join(
            f"{target.target_id}={target.reference_receptor}"
            for target in config.discovery.targets
        )
    )
    qualification = config.discovery.qualification
    print(
        "Discovery qualification: "
        f"seeds={len(qualification.seeds)}, "
        f"control={qualification.control_bound_state_id}, "
        f"native_LRMSD_A={qualification.max_native_ligand_rmsd_A:g}, "
        "required_control_seeds="
        f"{qualification.min_successful_control_seeds}"
    )
    topology = qualification.topology_calibration
    print(
        "Discovery topology calibration: "
        f"status={qualification.topology_calibration_status}, "
        "candidate_CA_thresholds_A="
        f"{','.join(f'{value:g}' for value in topology.candidate_ca_thresholds_A)}, "
        "comparator_upper_A="
        f"{topology.above_threshold_comparator_upper_bound_A:g}, "
        "active_CA_threshold_A="
        f"{config.discovery.cabsdock.max_disulfide_ca_distance_A:g}, "
        f"models_per_stratum={topology.models_per_stratum}, "
        "min_success_fraction="
        f"{topology.min_success_fraction_per_stratum:g}, "
        "min_successful_seeds="
        f"{topology.min_successful_seeds_per_stratum}"
    )
    for target in config.discovery.targets:
        analysis = target.analysis
        print(
            f"Discovery target {target.target_id}: "
            f"pilot={target.pilot_receptor}, "
            f"qualification={target.qualification_status}, "
            "contact_jaccard_distance="
            f"{analysis.contact_jaccard_distance}, "
            f"position_distance_A={analysis.position_distance_A}, "
            f"min_seed_support={analysis.min_seed_support}, "
            f"min_receptor_support={analysis.min_receptor_support}"
        )
    print(f"Registered structures: {len(config.receptors)}")
    print(f"Stage 2 apo/apo-like blind discovery: {len(apo_blind)}")
    print(f"Stage 3 stripped bound-state blind replication: {len(bound_blind)}")
    print(
        "Total planned ligand full-surface receptor conformations: "
        f"{len(apo_blind | bound_blind)}"
    )
    print(
        "Currently prepared receptor-only bases: "
        f"{sum(item.prepare for item in config.receptors)}"
    )
    print(
        "Stage 3 bound-state registry: "
        f"{len(config.validation.bound_states)} states; "
        f"{len(config.validation.local_controls)} configured local recovery controls"
    )
    print(
        "Stage 3 guided templates and full-enzyme references: "
        f"{len(config.validation.guided_templates)} guided; "
        f"{len(config.validation.environment_references)} environments"
    )
    print(f"FlexPepDock executable: {config.validation.rosetta.executable}")
    print(
        "FlexPepDock sampling: "
        f"parallel_tasks={config.validation.rosetta.parallel_tasks}, "
        f"seeds_per_start={len(config.validation.seeds)}, "
        f"decoys_per_seed={config.validation.rosetta.decoys_per_seed}, "
        "total_decoys_per_start="
        f"{len(config.validation.seeds) * config.validation.rosetta.decoys_per_seed}, "
        f"score_function={config.validation.rosetta.score_function}, "
        "lowres_preoptimize="
        f"{str(config.validation.rosetta.lowres_preoptimize).lower()}"
    )
    print(
        "Stage 3 all-atom handoff: "
        f"cg2all={config.validation.cg2all.executable}, "
        f"representation={config.validation.cg2all.representation}, "
        "poses_per_receptor_site="
        f"{config.validation.handoff.poses_per_receptor_site}"
    )
    print(
        "Stage 3 local refinement: "
        f"prepack_seed={config.validation.refinement.prepack_seed}, "
        f"ranking_score={config.validation.refinement.ranking_score}, "
        "random_translation_A="
        f"{config.validation.refinement.random_translation_A:g}, "
        "random_rotation_degrees="
        f"{config.validation.refinement.random_rotation_degrees:g}"
    )
    print(
        "Stage 4 sequence space: "
        f"mutable_positions={len(config.design.sequence.mutable_positions)}, "
        f"allowed_amino_acids={len(config.design.sequence.allowed_amino_acids)}, "
        f"max_combinations={config.design.combination.max_candidates}"
    )
    print(
        "Stage 4 iterative neighborhood: "
        f"max_parents={config.design.iteration.max_parents}, "
        f"max_total_mutations={config.design.iteration.max_total_mutations}, "
        f"max_candidates={config.design.iteration.max_candidates}"
    )
    print(
        "Stage 4 flexible verification: "
        f"max_candidates={config.design.finalists.max_candidates}, "
        f"max_md_candidates={config.design.finalists.max_md_candidates}, "
        f"seeds={len(config.design.finalists.seeds)}, "
        f"ranking_score={config.design.finalists.ranking_score}, "
        f"interface_score={config.design.finalists.interface_score}"
    )
    return 0


def _chemistry_record(config: AppConfig) -> int:
    destination = config.paths.data_dir / chemistry_record_relative_path(
        config.chemistry
    )
    assessment = write_chemistry_record(
        definition=config.chemistry, destination=destination
    )
    print(f"Chemistry record written: {destination}")
    print(f"Production ready: {str(assessment.production_ready).lower()}")
    if assessment.unresolved_fields:
        print("Unresolved fields: " + ", ".join(assessment.unresolved_fields))
    return 0


def _receptors_download(config: AppConfig) -> int:
    files = download_receptors(
        definitions=config.receptors,
        settings=config.download,
        data_dir=config.paths.data_dir,
    )
    print(f"Receptor files verified: {len(files)}")
    print(
        f"Downloaded: {sum(item.status == 'downloaded' for item in files)}; "
        f"cached: {sum(item.status == 'cached' for item in files)}"
    )
    print(
        "Manifest: "
        + str(config.paths.data_dir / "receptors" / "raw" / "download_manifest.json")
    )
    return 0


def _receptors_audit(config: AppConfig) -> int:
    results = audit_receptors(
        definitions=config.receptors,
        settings=config.audit,
        data_dir=config.paths.data_dir,
    )
    failures = [
        item.receptor_id for item in results if item.identity_status != "passed"
    ]
    print(f"Receptors audited: {len(results)}")
    print(
        f"Manual review required: {sum(item.manual_review_required for item in results)}"
    )
    print(
        "Audit manifest: "
        + str(config.paths.data_dir / "receptors" / "audit" / "audit_manifest.json")
    )
    if failures:
        raise ReceptorError("receptor identity audit failed: " + ", ".join(failures))
    return 0


def _receptors_prepare(config: AppConfig) -> int:
    audit_results = audit_receptors(
        definitions=config.receptors,
        settings=config.audit,
        data_dir=config.paths.data_dir,
    )
    failures = [
        item.receptor_id for item in audit_results if item.identity_status != "passed"
    ]
    if failures:
        raise ReceptorError("receptor identity audit failed: " + ", ".join(failures))
    results = prepare_receptors(
        definitions=config.receptors,
        settings=config.preparation,
        data_dir=config.paths.data_dir,
    )
    print(f"Receptors prepared: {len(results)}")
    for result in results:
        print(f"Prepared {result.receptor_id}: {result.output_path}")
    return 0


def execute(
    *,
    group: str,
    command: str,
    config: AppConfig,
    run_id: str | None = None,
    run_dir: Path | None = None,
    source_run: Path | None = None,
    control_run: Path | None = None,
    replication_run: Path | None = None,
    refinement_source: Path | None = None,
    topology_source: Path | None = None,
    candidate_ids: tuple[str, ...] = (),
    design_source: Path | None = None,
    target_cluster_ids: tuple[str, ...] = (),
    target_id: str | None = None,
) -> int:
    """按解析后的稳定命令身份调用处理器。"""
    if (group, command) == ("preparation", "status"):
        return preparation_status(config)
    if (group, command) == ("preparation", "run"):
        return preparation_run(config)
    if (group, command) == ("discovery", "status"):
        if target_id is None:
            raise RuntimeError("discovery status requires target")
        return discovery_status(config=config, target_id=target_id)
    if (group, command) == ("discovery", "qualification-plan"):
        if run_id is None or target_id is None:
            raise RuntimeError(
                "discovery qualification-plan requires run_id and target"
            )
        return discovery_qualification_plan(
            config=config,
            run_id=run_id,
            target_id=target_id,
            control_run=control_run,
        )
    if (group, command) == ("discovery", "qualification-run"):
        if run_dir is None:
            raise RuntimeError("discovery qualification-run requires run_dir")
        return discovery_qualification_run(config=config, run_dir=run_dir)
    if (group, command) == ("discovery", "qualification-analyze"):
        if run_dir is None:
            raise RuntimeError("discovery qualification-analyze requires run_dir")
        return discovery_qualification_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("discovery", "topology-calibration-plan"):
        if run_id is None or topology_source is None:
            raise RuntimeError(
                "discovery topology-calibration-plan requires run_id and source_run"
            )
        return discovery_topology_calibration_plan(
            config=config, source_run=topology_source, run_id=run_id
        )
    if (group, command) == ("discovery", "topology-calibration-run"):
        if run_dir is None:
            raise RuntimeError("discovery topology-calibration-run requires run_dir")
        return discovery_topology_calibration_run(config=config, run_dir=run_dir)
    if (group, command) == ("discovery", "topology-calibration-analyze"):
        if run_dir is None:
            raise RuntimeError(
                "discovery topology-calibration-analyze requires run_dir"
            )
        return discovery_topology_calibration_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("discovery", "plan"):
        if run_id is None or target_id is None:
            raise RuntimeError("discovery plan requires run_id and target")
        return discovery_plan(config=config, run_id=run_id, target_id=target_id)
    if (group, command) == ("discovery", "analyze"):
        if run_dir is None:
            raise RuntimeError("discovery analyze requires run_dir")
        return discovery_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("discovery", "run"):
        if run_dir is None:
            raise RuntimeError("discovery run requires run_dir")
        return discovery_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "status"):
        return validation_status(config)
    if (group, command) == ("validation", "prepare"):
        return validation_prepare(config)
    if (group, command) == ("validation", "tool-check"):
        return validation_tool_check(config)
    if (group, command) == ("validation", "control-plan"):
        if run_id is None:
            raise RuntimeError("validation control-plan requires run_id")
        return validation_control_plan(config=config, run_id=run_id)
    if (group, command) == ("validation", "control-run"):
        if run_dir is None:
            raise RuntimeError("validation control-run requires run_dir")
        return validation_control_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "replication-plan"):
        if run_id is None or target_id is None:
            raise RuntimeError("validation replication-plan requires run_id and target")
        return validation_replication_plan(
            config=config, run_id=run_id, target_id=target_id
        )
    if (group, command) == ("validation", "replication-run"):
        if run_dir is None:
            raise RuntimeError("validation replication-run requires run_dir")
        return validation_replication_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "replication-analyze"):
        if run_dir is None:
            raise RuntimeError("validation replication-analyze requires run_dir")
        return validation_replication_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "replication-compare"):
        if run_dir is None or source_run is None:
            raise RuntimeError(
                "validation replication-compare requires discovery_run and run_dir"
            )
        return validation_replication_compare(
            config=config,
            discovery_run=source_run,
            replication_run=run_dir,
        )
    if (group, command) == ("validation", "handoff-plan"):
        if run_id is None or source_run is None:
            raise RuntimeError(
                "validation handoff-plan requires run_id and discovery_run"
            )
        return validation_handoff_plan(
            config=config,
            discovery_run=source_run,
            run_id=run_id,
            candidate_ids=candidate_ids,
        )
    if (group, command) == ("validation", "handoff-run"):
        if run_dir is None:
            raise RuntimeError("validation handoff-run requires run_dir")
        return validation_handoff_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "guided-plan"):
        if run_id is None:
            raise RuntimeError("validation guided-plan requires run_id")
        return validation_guided_plan(config=config, run_id=run_id)
    if (group, command) == ("validation", "guided-run"):
        if run_dir is None:
            raise RuntimeError("validation guided-run requires run_dir")
        return validation_guided_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "refinement-plan"):
        if run_id is None or refinement_source is None:
            raise RuntimeError(
                "validation refinement-plan requires run_id and source_run"
            )
        return validation_refinement_plan(
            config=config, source_run=refinement_source, run_id=run_id
        )
    if (group, command) == ("validation", "refinement-run"):
        if run_dir is None:
            raise RuntimeError("validation refinement-run requires run_dir")
        return validation_refinement_run(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "refinement-analyze"):
        if run_dir is None:
            raise RuntimeError("validation refinement-analyze requires run_dir")
        return validation_refinement_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "environment-map"):
        if run_dir is None:
            raise RuntimeError("validation environment-map requires run_dir")
        return validation_environment_map(config=config, run_dir=run_dir)
    if (group, command) == ("validation", "candidate-review"):
        if run_dir is None or source_run is None or replication_run is None:
            raise RuntimeError(
                "validation candidate-review requires discovery_run, replication_run, and run_dir"
            )
        return validation_candidate_review(
            config=config,
            discovery_run=source_run,
            replication_run=replication_run,
            refinement_run=run_dir,
        )
    if (group, command) == ("design", "status"):
        return design_status(config)
    if (group, command) == ("design", "tool-check"):
        return design_tool_check(config)
    if (group, command) == ("design", "single-plan"):
        if run_id is None or design_source is None:
            raise RuntimeError("design single-plan requires run_id and source_run")
        return design_single_plan(
            config=config,
            source_run=design_source,
            run_id=run_id,
            target_cluster_ids=target_cluster_ids,
        )
    if (group, command) == ("design", "combination-plan"):
        if run_id is None or design_source is None:
            raise RuntimeError("design combination-plan requires run_id and source_run")
        return design_combination_plan(
            config=config, source_run=design_source, run_id=run_id
        )
    if (group, command) == ("design", "iteration-plan"):
        if run_id is None or design_source is None:
            raise RuntimeError("design iteration-plan requires run_id and source_run")
        return design_iteration_plan(
            config=config,
            source_run=design_source,
            run_id=run_id,
            parent_candidate_ids=candidate_ids,
        )
    if (group, command) == ("design", "screen-run"):
        if run_dir is None:
            raise RuntimeError("design screen-run requires run_dir")
        return design_screen_run(config=config, run_dir=run_dir)
    if (group, command) == ("design", "screen-analyze"):
        if run_dir is None:
            raise RuntimeError("design screen-analyze requires run_dir")
        return design_screen_analyze(config=config, run_dir=run_dir)
    if (group, command) == ("design", "finalist-plan"):
        if run_id is None or design_source is None:
            raise RuntimeError("design finalist-plan requires run_id and source_run")
        return design_finalist_plan(
            config=config, source_run=design_source, run_id=run_id
        )
    if (group, command) == ("design", "finalist-run"):
        if run_dir is None:
            raise RuntimeError("design finalist-run requires run_dir")
        return design_finalist_run(config=config, run_dir=run_dir)
    if (group, command) == ("design", "finalist-analyze"):
        if run_dir is None:
            raise RuntimeError("design finalist-analyze requires run_dir")
        return design_finalist_analyze(config=config, run_dir=run_dir)
    handlers = {
        ("config", "check"): _config_check,
        ("chemistry", "record"): _chemistry_record,
        ("receptors", "download"): _receptors_download,
        ("receptors", "audit"): _receptors_audit,
        ("receptors", "prepare"): _receptors_prepare,
    }
    handler = handlers.get((group, command))
    if handler is None:
        raise RuntimeError(f"unhandled command: {group} {command}")
    return handler(config)
