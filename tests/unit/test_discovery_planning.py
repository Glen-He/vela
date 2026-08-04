import json
from dataclasses import replace
from inspect import signature
from pathlib import Path

import pytest

from vela.config import load_config
from vela.config.models import AppConfig, PathsConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
)
from vela.core.typed_data import object_mapping
from vela.discovery.models import (
    CabsDockSettings,
    DiscoveryError,
    DiscoveryQualificationSettings,
    DiscoverySettings,
    DiscoveryTargetSettings,
    ReceptorEnsembleSettings,
    SiteAnalysisSettings,
    TopologyCalibrationSettings,
)
from vela.discovery.qualification.topology import (
    topology_calibration_contract,
    validate_topology_source_evidence,
)
from vela.discovery.readiness import assess_discovery_readiness
from vela.discovery.sampling.cabsdock import build_cabsdock_command
from vela.discovery.sampling.evidence import (
    candidate_selection_contract,
    collect_cabsdock_evidence,
)
from vela.discovery.sampling.planning import build_tasks, write_discovery_plan
from vela.preparation.chemistry import (
    ChemistryDefinition,
    DisulfideBond,
    HistidineState,
    write_chemistry_record,
)
from vela.preparation.readiness import assess_preparation_readiness
from vela.preparation.receptors.models import (
    AltlocConfig,
    CrystalContactConfig,
    DownloadConfig,
    ReceptorAuditConfig,
    ReceptorDefinition,
    ReceptorPreparationConfig,
)
from vela.validation.bound_states.replication import (
    REPLICATION_EVIDENCE,
    build_replication_tasks,
    write_replication_plan,
)
from vela.validation.models import BoundStateDefinition


def test_topology_calibration_contract_is_native_free_and_all_atom() -> None:
    config = load_config(Path(__file__).resolve().parents[2] / "configs")

    contract = topology_calibration_contract(config)

    sampling = object_mapping(
        contract["stratified_sampling"], name="stratified sampling"
    )
    reconstruction = object_mapping(
        contract["all_atom_reconstruction"], name="all-atom reconstruction"
    )
    feasibility = object_mapping(
        contract["topology_feasibility"], name="topology feasibility"
    )
    assert sampling["native_information_used"] is False
    assert sampling["ca_strata_upper_bounds_A"] == [6.0, 7.0, 8.0, 10.0, 12.0]
    assert feasibility["candidate_thresholds_A"] == [6.0, 7.0, 8.0, 10.0]
    assert feasibility["threshold_selection"] == (
        "highest_contiguous_passing_candidate"
    )
    assert feasibility["above_threshold_comparator_upper_bound_A"] == 12.0
    assert feasibility["chemical_bond_claimed"] is False
    assert reconstruction["pipeline"] == (
        "cg2all_peptide_then_aligned_experimental_receptor_graft_then_"
        "RosettaScripts_ForceDisulfides_repack_then_FlexPepDock_prepack_"
        "then_site_coordinate_constrained_single_local_refine"
    )
    constraints = object_mapping(
        reconstruction["site_coordinate_constraints"], name="site constraints"
    )
    assert constraints == {
        "atoms": "all_peptide_CA",
        "reference_frame": "fixed_receptor_first_CA",
        "function": "FLAT_HARMONIC",
        "flat_width_A": 2.0,
        "standard_deviation_A": 1.0,
        "score_weight": 1.0,
    }
    assert reconstruction["terminal_chemistry_assessed"] is False


def test_topology_source_validation_rejects_changed_cabs_archive(
    tmp_path: Path,
) -> None:
    source = tmp_path / "qualification"
    task_id = "control__seed_1"
    plan_path = source / "qualification_plan.json"
    plan: dict[str, JsonValue] = {
        "schema": "vela.discovery-qualification-plan/6",
        "stage": "discovery_qualification",
        "target_id": "test_target",
        "tasks": [{"task_id": task_id}],
    }
    atomic_write_json(plan_path, plan)
    atomic_write_json(
        source / "qualification_sampling.json",
        {
            "schema": "vela.discovery-qualification-sampling/6",
            "stage": "discovery_qualification",
            "status": "sampling_completed",
            "target_id": "test_target",
            "qualification_plan_sha256": sha256_file(plan_path),
            "tasks": [{"task_id": task_id, "execution_status": "completed"}],
        },
    )
    task_dir = source / "tasks" / task_id
    archive = task_dir / "result.cbs"
    atomic_write_text(archive, "original archive\n")
    atomic_write_json(
        task_dir / "task_result.json",
        {
            "schema": "vela.cabsdock-task-result/5",
            "execution_status": "completed",
            "task_id": task_id,
            "run_manifest_sha256": sha256_file(plan_path),
            "outputs": {
                "cabs_archive": {
                    "path": archive.name,
                    "sha256": sha256_file(archive),
                }
            },
        },
    )
    validate_topology_source_evidence(source_run=source, source_plan=plan)

    atomic_write_text(archive, "changed archive\n")

    with pytest.raises(DiscoveryError, match="hash mismatch"):
        validate_topology_source_evidence(source_run=source, source_plan=plan)


def test_candidate_selection_contract_is_complete_and_has_no_native_input() -> None:
    config = load_config(Path(__file__).resolve().parents[2] / "configs")

    contract = candidate_selection_contract(config.discovery.cabsdock)
    pose_clustering = object_mapping(
        contract["pose_clustering"], name="pose clustering"
    )

    assert contract["native_information_used"] is False
    assert contract["minimum_input_models"] == 10
    assert pose_clustering["cluster_coverage"] == "complete"
    assert "max_clusters_per_site" not in pose_clustering
    assert set(contract) == {
        "mode",
        "native_information_used",
        "input_pool",
        "minimum_input_models",
        "topology_feasibility",
        "receptor_contact",
        "coordinate_frame",
        "site_clustering",
        "pose_clustering",
        "representatives",
    }
    assert not any(
        "native" in name for name in signature(collect_cabsdock_evidence).parameters
    )


PROJECT_CONFIG = Path(__file__).resolve().parents[2] / "configs"


def _receptor(
    receptor_id: str, *, pdb_id: str, target: str, chain: str, state: str
) -> ReceptorDefinition:
    return ReceptorDefinition(
        receptor_id=receptor_id,
        pdb_id=pdb_id,
        target=target,
        uniprot_accession="P68400" if target == "ck2_alpha" else "P19784",
        author_chain_id=chain,
        structure_state=state,
        roles=("blind_discovery",),
        prepare=True,
        water_policy="remove_all",
        remove_components=(),
        retain_components=(),
        selection_reason="Test ensemble member.",
    )


def _ready_config(tmp_path: Path) -> AppConfig:
    data_dir = tmp_path / "data"
    outputs_dir = tmp_path / "outputs"
    source_path = tmp_path / "configs" / "fixture.toml"
    atomic_write_text(source_path, "[test]\nfixture = true\n")
    chemistry = ChemistryDefinition(
        ligand_id="test-peptide",
        chemistry_id="p15-test",
        sequence="CWMSPRHLGTC",
        chirality="L",
        disulfide_bonds=(DisulfideBond(1, 11),),
        n_terminus="protonated",
        c_terminus="deprotonated",
        target_ph=7.4,
        net_charge=1,
        histidines=(HistidineState(7, "HIE"),),
        other_modifications_status="none",
        other_modifications=(),
        decision_sources=("test protocol",),
    )
    write_chemistry_record(
        definition=chemistry,
        destination=(data_dir / "chemistry" / "test-peptide" / "chemistry_record.json"),
    )
    receptors = (
        _receptor("3Q04_A", pdb_id="3Q04", target="ck2_alpha", chain="A", state="apo"),
        _receptor("3QA0_A", pdb_id="3QA0", target="ck2_alpha", chain="A", state="apo"),
        _receptor(
            "5YF9_X",
            pdb_id="5YF9",
            target="ck2_alpha_prime",
            chain="X",
            state="apo_like",
        ),
        _receptor(
            "5Y9M_A",
            pdb_id="5Y9M",
            target="ck2_alpha_prime",
            chain="A",
            state="apo_like",
        ),
    )

    def file_record(path: Path) -> dict[str, JsonValue]:
        return {
            "path": path.relative_to(data_dir).as_posix(),
            "sha256": sha256_file(path),
        }

    download_entries: list[dict[str, JsonValue]] = []
    preparation_entries: list[dict[str, JsonValue]] = []
    for receptor in receptors:
        raw_path = data_dir / "receptors" / "raw" / f"{receptor.pdb_id}.cif"
        metadata_path = data_dir / "receptors" / "raw" / f"{receptor.pdb_id}.entry.json"
        prepared_path = (
            data_dir / "receptors" / "prepared" / f"{receptor.receptor_id}.cif"
        )
        atomic_write_text(raw_path, f"data_{receptor.pdb_id}\n")
        atomic_write_text(metadata_path, f'{{"pdb_id":"{receptor.pdb_id}"}}\n')
        atomic_write_text(prepared_path, f"data_{receptor.receptor_id}\n")
        download_entries.append(
            {
                "receptor_id": receptor.receptor_id,
                "pdb_id": receptor.pdb_id,
                "files": [file_record(raw_path), file_record(metadata_path)],
            }
        )
        preparation_entries.append(
            {
                "receptor_id": receptor.receptor_id,
                "pdb_id": receptor.pdb_id,
                "author_chain_id": receptor.author_chain_id,
                "source": file_record(raw_path),
                "output": file_record(prepared_path),
            }
        )
    atomic_write_json(
        data_dir / "receptors" / "raw" / "download_manifest.json",
        {
            "schema": "vela.receptor-download-manifest/1",
            "entries": download_entries,
        },
    )

    audit_dir = data_dir / "receptors" / "audit"
    audit_outputs = (
        audit_dir / "structure_summary.tsv",
        audit_dir / "chain_summary.tsv",
        audit_dir / "component_summary.tsv",
        audit_dir / "missing_summary.tsv",
        audit_dir / "sequence_difference_summary.tsv",
    )
    atomic_write_text(
        audit_outputs[0],
        "receptor_id\n" + "".join(f"{item.receptor_id}\n" for item in receptors),
    )
    for path in audit_outputs[1:]:
        atomic_write_text(path, "fixture\n")
    atomic_write_json(
        audit_dir / "audit_manifest.json",
        {
            "schema": "vela.receptor-audit-manifest/1",
            "parameters": {
                "crystal_contacts": {
                    "distance_A": 4.5,
                    "min_occupancy": 0.01,
                    "include_hydrogens": False,
                }
            },
            "identity_failures": [],
            "outputs": [file_record(path) for path in audit_outputs],
        },
    )

    preparation_reports = (
        data_dir / "receptors" / "prepared" / "preparation_summary.tsv",
        data_dir / "receptors" / "prepared" / "altloc_decisions.tsv",
    )
    for path in preparation_reports:
        atomic_write_text(path, "fixture\n")
    atomic_write_json(
        data_dir / "receptors" / "prepared" / "preparation_manifest.json",
        {
            "schema": "vela.receptor-preparation-manifest/1",
            "parameters": {
                "altloc": {
                    "preferred_label": "A",
                }
            },
            "entries": preparation_entries,
            "reports": [file_record(path) for path in preparation_reports],
        },
    )

    def qualification_report(target_id: str) -> Path:
        report_path = tmp_path / "qualification" / f"{target_id}.json"
        atomic_write_json(
            report_path,
            {
                "schema": "vela.discovery-qualification-report/7",
                "status": "qualified",
                "target_id": target_id,
                "recommended_target_config": {
                    "qualification_status": "qualified",
                    "contact_jaccard_distance": 0.5,
                    "position_distance_A": 4.0,
                    "min_seed_support": 2,
                    "min_receptor_support": 2,
                },
            },
        )
        return report_path

    alpha_report = qualification_report("ck2_alpha")
    prime_report = qualification_report("ck2_alpha_prime")
    executable = tmp_path / "CABSdock"
    atomic_write_text(executable, "#!/bin/sh\n")
    source_dir = tmp_path / "cabs-source"
    source_dir.mkdir()
    for relative in (
        "CABS/analysis/restraints.py",
        "CABS/core/cabs.py",
        "CABS/core/job.py",
        "CABS/core/trajectory.py",
        "CABS/data/data0.dat",
        "CABS/io/config.json",
        "CABS/structures/atom.py",
        "CABS/utils/filter.py",
        "CABS/utils/utils.py",
    ):
        source = source_dir / relative
        atomic_write_text(source, f"# fixture {relative}\n")
    topology_report = tmp_path / "qualification" / "topology.json"
    atomic_write_json(
        topology_report,
        {"schema": "vela.topology-calibration-report/1", "status": "qualified"},
    )
    settings = DiscoverySettings(
        method_id="qualified-global-method",
        adapter_id="qualified-adapter",
        seeds=(11, 22),
        ensemble=ReceptorEnsembleSettings(
            min_receptors_per_target=2,
            allowed_structure_states=("apo", "apo_like"),
        ),
        cabsdock=CabsDockSettings(
            executable=executable,
            source_dir=source_dir,
            source_revision="1" * 40,
            seed_workers=2,
            peptide_secondary_structure="CCCCCCCCCCC",
            mc_annealing=20,
            mc_cycles=50,
            mc_steps=50,
            replicas=10,
            replicas_dtemp=0.5,
            temperature_initial=2.0,
            temperature_final=1.0,
            binding_interactions=1.0,
            protein_restraint_gap=5,
            protein_restraint_min_A=5.0,
            protein_restraint_max_A=15.0,
            filtering_count=1000,
            clustering_medoids=10,
            clustering_iterations=100,
            trajectory_contact_ca_threshold_A=10.0,
            disulfide_ca_restraint_distance_A=5.5,
            disulfide_ca_restraint_weight=1.0,
            max_reconstructable_disulfide_ca_distance_A=8.0,
            min_models_for_selection=10,
            selection_contact_jaccard_distance=0.8,
            selection_position_distance_A=12.0,
            pose_clustering_rmsd_A=4.0,
            max_sites_per_task=64,
            max_pose_clusters_per_site=4,
        ),
        qualification=DiscoveryQualificationSettings(
            seeds=(31, 32),
            control_bound_state_id="4IB5_A_D",
            control_receptor_id="3Q04_A",
            control_target_id="ck2_alpha",
            control_secondary_structure="CCCCCCCCCCCCC",
            max_native_ligand_rmsd_A=4.0,
            max_native_site_centroid_distance_A=4.0,
            min_native_receptor_contact_fraction=0.2,
            min_native_sampling_seed_support=1,
            min_native_site_seed_support=1,
            min_selection_native_seed_recall_fraction=1.0,
            topology_calibration_status="qualified",
            topology_calibration_report=topology_report,
            topology_calibration_report_sha256=sha256_file(topology_report),
            topology_calibration=TopologyCalibrationSettings(
                candidate_ca_thresholds_A=(6.0, 7.0, 8.0),
                above_threshold_comparator_upper_bound_A=10.0,
                models_per_stratum=8,
                min_success_fraction_per_stratum=0.75,
                min_successful_seeds_per_stratum=6,
                min_interchain_heavy_atom_distance_A=1.2,
                min_nonlocal_peptide_heavy_atom_distance_A=1.2,
                max_peptide_internal_ca_rmsd_A=4.0,
                max_ligand_centroid_displacement_A=4.0,
                contact_ca_threshold_A=10.0,
                min_receptor_contact_retention_fraction=0.5,
                site_coordinate_constraint_flat_width_A=2.0,
                site_coordinate_constraint_sd_A=1.0,
                site_coordinate_constraint_weight=1.0,
            ),
        ),
        targets=(
            DiscoveryTargetSettings(
                target_id="ck2_alpha",
                reference_receptor="3Q04_A",
                pilot_receptor="3Q04_A",
                qualification_status="qualified",
                qualification_report=alpha_report,
                qualification_report_sha256=sha256_file(alpha_report),
                analysis=SiteAnalysisSettings(
                    contact_jaccard_distance=0.5,
                    position_distance_A=4.0,
                    min_seed_support=2,
                    min_receptor_support=2,
                ),
            ),
            DiscoveryTargetSettings(
                target_id="ck2_alpha_prime",
                reference_receptor="5YF9_X",
                pilot_receptor="5YF9_X",
                qualification_status="qualified",
                qualification_report=prime_report,
                qualification_report_sha256=sha256_file(prime_report),
                analysis=SiteAnalysisSettings(
                    contact_jaccard_distance=0.5,
                    position_distance_A=4.0,
                    min_seed_support=2,
                    min_receptor_support=2,
                ),
            ),
        ),
    )
    return AppConfig(
        source_dir=source_path.parent,
        source_files=(source_path,),
        source_snapshot_text=source_path.read_text(encoding="utf-8"),
        source_snapshot_sha256=sha256_file(source_path),
        paths=PathsConfig(data_dir=data_dir, outputs_dir=outputs_dir),
        download=DownloadConfig(
            coordinate_base_url="https://example.test/coordinates",
            metadata_base_url="https://example.test/metadata",
            retries=1,
            timeout_seconds=1.0,
            backoff_initial_seconds=0.0,
            backoff_multiplier=2.0,
            chunk_size_bytes=1024,
            user_agent="vela-test",
        ),
        audit=ReceptorAuditConfig(
            crystal_contacts=CrystalContactConfig(
                distance_A=4.5,
                min_occupancy=0.01,
                include_hydrogens=False,
            )
        ),
        preparation=ReceptorPreparationConfig(altloc=AltlocConfig(preferred_label="A")),
        chemistry=chemistry,
        receptors=receptors,
        discovery=settings,
        validation=load_config(PROJECT_CONFIG).validation,
        design=load_config(PROJECT_CONFIG).design,
    )


def test_ready_target_expands_its_receptors_by_independent_seed(
    tmp_path: Path,
) -> None:
    config = _ready_config(tmp_path)

    readiness = assess_discovery_readiness(
        target_id="ck2_alpha",
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )
    plan = write_discovery_plan(
        config=config, run_id="production-001", target_id="ck2_alpha"
    )

    assert readiness.ready
    assert plan.target_id == "ck2_alpha"
    assert len(plan.tasks) == 4
    assert {task.receptor_id for task in plan.tasks} == {
        "3Q04_A",
        "3QA0_A",
    }
    assert {task.seed for task in plan.tasks} == {11, 22}
    raw_document: object = json.loads(
        (plan.run_dir / "run_manifest.json").read_text(encoding="utf-8")
    )
    document = object_mapping(raw_document, name="run_manifest")
    assert document["known_site_information_used"] is False
    assert document["target_id"] == "ck2_alpha"
    assert document["task_count"] == 4
    method_parameters = object_mapping(
        document.get("method_parameters"), name="method_parameters"
    )
    cabsdock = object_mapping(method_parameters.get("cabsdock"), name="cabsdock")
    assert cabsdock["executable_sha256"] == sha256_file(
        config.discovery.cabsdock.executable
    )


def test_stage_one_readiness_is_independent_of_discovery_method(
    tmp_path: Path,
) -> None:
    config = _ready_config(tmp_path)
    unresolved_discovery = replace(
        config.discovery,
        seeds=(),
        targets=(
            replace(
                config.discovery.target("ck2_alpha"),
                qualification_status="unresolved",
                qualification_report=None,
                qualification_report_sha256=None,
                analysis=SiteAnalysisSettings(
                    contact_jaccard_distance=None,
                    position_distance_A=None,
                    min_seed_support=None,
                    min_receptor_support=None,
                ),
            ),
        ),
    )

    readiness = assess_preparation_readiness(
        chemistry=config.chemistry,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )

    assert readiness.ready
    assert len(readiness.audited_receptor_ids) == 4
    assert len(readiness.prepared_receptor_ids) == 4
    assert unresolved_discovery.target("ck2_alpha").qualification_status == "unresolved"


def test_stage_one_readiness_rejects_tampered_download(tmp_path: Path) -> None:
    config = _ready_config(tmp_path)
    raw_path = config.paths.data_dir / "receptors" / "raw" / "3Q04.cif"
    atomic_write_text(raw_path, "tampered\n")

    readiness = assess_preparation_readiness(
        chemistry=config.chemistry,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )

    assert not readiness.ready
    codes = {issue.code for issue in readiness.issues}
    assert "download_hash_mismatch" in codes
    assert "preparation_source_hash_mismatch" in codes


def test_discovery_plan_never_overwrites_existing_run(tmp_path: Path) -> None:
    config = _ready_config(tmp_path)
    write_discovery_plan(config=config, run_id="production-001", target_id="ck2_alpha")

    with pytest.raises(DiscoveryError, match="already exists"):
        write_discovery_plan(
            config=config, run_id="production-001", target_id="ck2_alpha"
        )


def test_cabsdock_command_uses_ca_ring_restraint_and_no_known_site_restraint(
    tmp_path: Path,
) -> None:
    config = _ready_config(tmp_path)
    task = build_tasks(config, target_id="ck2_alpha")[0]

    command = build_cabsdock_command(
        task=task,
        settings=config.discovery.cabsdock,
        chemistry=config.chemistry,
        secondary_structure=config.discovery.cabsdock.peptide_secondary_structure,
        task_dir=tmp_path / "task",
    )

    assert command[0] == str(config.discovery.cabsdock.executable)
    assert command[command.index("-p") + 1] == "CWMSPRHLGTC:CCCCCCCCCCC"
    index = command.index("--ca-rest-add")
    assert command[index + 1 : index + 5] == (
        "1:PEP1",
        "11:PEP1",
        "5.5",
        "1.0",
    )
    assert command[command.index("-A") + 1] == "N"
    assert "--exclude" not in command
    assert "-F" not in command
    assert "--sc-rest-add" not in command


def test_readiness_rejects_support_threshold_above_frozen_seeds(
    tmp_path: Path,
) -> None:
    config = _ready_config(tmp_path)
    impossible = replace(config.discovery, seeds=(11,))

    readiness = assess_discovery_readiness(
        target_id="ck2_alpha",
        chemistry=config.chemistry,
        settings=impossible,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )

    assert not readiness.ready
    assert "seed_support_impossible" in {issue.code for issue in readiness.issues}


def test_readiness_rejects_stage_one_manifest_from_old_parameters(
    tmp_path: Path,
) -> None:
    config = _ready_config(tmp_path)
    changed_audit = replace(
        config.audit,
        crystal_contacts=replace(config.audit.crystal_contacts, distance_A=5.0),
    )

    readiness = assess_discovery_readiness(
        target_id="ck2_alpha",
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=changed_audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )

    assert not readiness.ready
    assert "audit_parameters_stale" in {issue.code for issue in readiness.issues}


def _replication_config(config: AppConfig) -> AppConfig:
    receptor = ReceptorDefinition(
        receptor_id="GENERIC_RECEPTOR",
        pdb_id="1ABC",
        target="ck2_alpha",
        uniprot_accession="P68400",
        author_chain_id="R",
        structure_state="ligand_bound",
        roles=("bound_state_review", "bound_state_blind_replication"),
        prepare=False,
        water_policy=None,
        remove_components=(),
        retain_components=(),
        selection_reason=None,
    )
    state = BoundStateDefinition(
        state_id="GENERIC_STATE",
        receptor_id=receptor.receptor_id,
        ligand_id="GENERIC_LIGAND",
        ligand_author_chain_id="L",
        local_control_kind="standard_cyclic_peptide",
        ligand_sequence="CAC",
        disulfide_bonds=(DisulfideBond(1, 3),),
        histidines=(),
        selection_reason="Synthetic replaceable validation state.",
    )
    control = replace(
        config.validation.local_controls[0],
        control_id="GENERIC_CONTROL",
        bound_state_id=state.state_id,
    )
    validation = replace(
        config.validation,
        bound_states=(state,),
        local_controls=(control,),
        guided_templates=(),
    )
    receptor_path = (
        config.paths.data_dir
        / "validation"
        / "bound_states"
        / state.state_id
        / "receptor_only.cif"
    )
    atomic_write_text(receptor_path, "data_generic_receptor\n")
    atomic_write_json(
        config.paths.data_dir
        / "validation"
        / "bound_states"
        / "preparation_manifest.json",
        {
            "schema": "vela.bound-state-preparation-manifest/1",
            "entries": [
                {
                    "state_id": state.state_id,
                    "receptor_id": state.receptor_id,
                    "ligand_id": state.ligand_id,
                    "source_ligand_chain": state.ligand_author_chain_id,
                    "local_control_kind": state.local_control_kind,
                    "outputs": {
                        "receptor_only": {
                            "path": receptor_path.relative_to(
                                config.paths.data_dir
                            ).as_posix(),
                            "sha256": sha256_file(receptor_path),
                        }
                    },
                }
            ],
        },
    )
    receptors = (*config.receptors, receptor)
    raw_dir = config.paths.data_dir / "receptors" / "raw"
    atomic_write_text(raw_dir / f"{receptor.pdb_id}.cif", "data_generic_raw\n")
    atomic_write_text(
        raw_dir / f"{receptor.pdb_id}.entry.json",
        f'{{"pdb_id":"{receptor.pdb_id}"}}\n',
    )
    download_entries: list[dict[str, JsonValue]] = []
    for item in receptors:
        coordinate = raw_dir / f"{item.pdb_id}.cif"
        metadata = raw_dir / f"{item.pdb_id}.entry.json"
        download_entries.append(
            {
                "receptor_id": item.receptor_id,
                "pdb_id": item.pdb_id,
                "files": [
                    {
                        "path": coordinate.relative_to(
                            config.paths.data_dir
                        ).as_posix(),
                        "sha256": sha256_file(coordinate),
                    },
                    {
                        "path": metadata.relative_to(config.paths.data_dir).as_posix(),
                        "sha256": sha256_file(metadata),
                    },
                ],
            }
        )
    atomic_write_json(
        raw_dir / "download_manifest.json",
        {
            "schema": "vela.receptor-download-manifest/1",
            "entries": download_entries,
        },
    )

    audit_dir = config.paths.data_dir / "receptors" / "audit"
    audit_outputs = (
        audit_dir / "structure_summary.tsv",
        audit_dir / "chain_summary.tsv",
        audit_dir / "component_summary.tsv",
        audit_dir / "missing_summary.tsv",
        audit_dir / "sequence_difference_summary.tsv",
    )
    atomic_write_text(
        audit_outputs[0],
        "receptor_id\n" + "".join(f"{item.receptor_id}\n" for item in receptors),
    )
    atomic_write_json(
        audit_dir / "audit_manifest.json",
        {
            "schema": "vela.receptor-audit-manifest/1",
            "parameters": {
                "crystal_contacts": {
                    "distance_A": config.audit.crystal_contacts.distance_A,
                    "min_occupancy": config.audit.crystal_contacts.min_occupancy,
                    "include_hydrogens": (
                        config.audit.crystal_contacts.include_hydrogens
                    ),
                }
            },
            "identity_failures": [],
            "outputs": [
                {
                    "path": path.relative_to(config.paths.data_dir).as_posix(),
                    "sha256": sha256_file(path),
                }
                for path in audit_outputs
            ],
        },
    )
    return replace(
        config,
        receptors=receptors,
        validation=validation,
    )


def test_bound_state_replication_is_derived_from_registry_and_frozen_seeds(
    tmp_path: Path,
) -> None:
    config = _replication_config(_ready_config(tmp_path))

    tasks = build_replication_tasks(config, target_id="ck2_alpha")
    plan = write_replication_plan(
        config=config, run_id="replication-001", target_id="ck2_alpha"
    )

    assert len(tasks) == 2
    assert {task.receptor_id for task in tasks} == {"GENERIC_STATE"}
    assert {task.seed for task in tasks} == {11, 22}
    assert {task.evidence_category for task in tasks} == {REPLICATION_EVIDENCE}
    assert plan.tasks == tasks
    raw_document: object = json.loads(
        (plan.run_dir / "run_manifest.json").read_text(encoding="utf-8")
    )
    document = object_mapping(raw_document, name="replication manifest")
    assert document["evidence_category"] == REPLICATION_EVIDENCE
    assert document["target_id"] == "ck2_alpha"
    assert document["known_site_information_used"] is False
    selection = object_mapping(
        document.get("receptor_selection"), name="receptor selection"
    )
    assert selection["role"] == "bound_state_blind_replication"
