import shutil
from pathlib import Path

import pytest

from vela.config import ConfigError, load_config
from vela.core.errors import VelaError
from vela.design.models import DesignError
from vela.discovery.models import DiscoveryError
from vela.validation.models import ValidationError

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROJECT_CONFIG_DIR = PROJECT_ROOT / "configs"


def _copy_project_config(destination: Path) -> None:
    """复制完整的四文件项目参数集, 供局部覆盖测试使用。"""
    destination.mkdir()
    for source in PROJECT_CONFIG_DIR.glob("*.toml"):
        shutil.copyfile(source, destination / source.name)


def test_project_config_resolves_paths_relative_to_config() -> None:
    config = load_config(PROJECT_CONFIG_DIR)
    assert config.paths.data_dir == PROJECT_ROOT / "data"
    assert config.paths.outputs_dir == PROJECT_ROOT / "outputs"
    assert len(config.receptors) == 10
    assert config.chemistry.ligand_id == "p15"
    assert config.chemistry.chemistry_id == "p15-free-n-amide-hie-ph7p4"
    assert config.chemistry.n_terminus == "NH3+"
    assert config.chemistry.c_terminus == "CONH2"
    assert config.chemistry.target_ph == 7.4
    assert config.chemistry.net_charge == 2
    assert config.chemistry.histidines[0].state == "HIE"
    assert [item.receptor_id for item in config.receptors if item.prepare] == [
        "3Q04_A",
        "3QA0_A",
        "3Q9X_A",
        "5YF9_X",
        "5Y9M_A",
    ]
    assert [
        item.receptor_id
        for item in config.receptors
        if "bound_state_blind_replication" in item.roles
    ] == ["4IB5_A", "9FBM_A", "9FBI_A"]
    assert [
        item.receptor_id
        for item in config.receptors
        if "full_enzyme_environment" in item.roles
    ] == ["4MD7_E", "1JWH_A"]
    assert (
        config.discovery.method_id
        == "cabsdock-global-cyclic-ca-restraint-site-first-cg"
    )
    assert config.discovery.adapter_id == "vela-cabsdock-cyclic-site-first-cg"
    assert config.discovery.seeds == tuple(range(120623, 120631))
    assert config.discovery.qualification.seeds == tuple(range(120679, 120687))
    assert not config.discovery.config_complete
    assert config.download.chunk_size_bytes == 1_048_576
    assert config.download.backoff_multiplier == 2.0
    assert config.audit.crystal_contacts.distance_A == 4.5
    assert config.audit.crystal_contacts.min_occupancy == 0.01
    assert not config.audit.crystal_contacts.include_hydrogens
    assert config.preparation.altloc.preferred_label == "A"
    assert tuple(target.target_id for target in config.discovery.targets) == (
        "ck2_alpha",
        "ck2_alpha_prime",
    )
    assert config.discovery.ensemble.min_receptors_per_target == 2
    assert config.discovery.target("ck2_alpha").reference_receptor == "3Q04_A"
    assert config.discovery.cabsdock.executable == Path(
        "/home/glen/apps/cabsflex/.venv-cabs/bin/CABSdock"
    )
    assert config.discovery.cabsdock.source_dir == Path("/home/glen/apps/cabsflex")
    assert (
        config.discovery.cabsdock.source_revision
        == "36ec10e681d0b6c5c991101bc0bdfc6e224c5e0b"
    )
    assert config.discovery.cabsdock.mc_annealing == 20
    assert config.discovery.cabsdock.seed_workers == 8
    assert config.discovery.cabsdock.replicas_dtemp == 0.5
    assert config.discovery.cabsdock.temperature_initial == 2.0
    assert config.discovery.cabsdock.temperature_final == 1.0
    assert config.discovery.target("ck2_alpha").analysis.min_seed_support == 2
    assert config.discovery.target("ck2_alpha").analysis.min_receptor_support == 2
    assert (
        config.discovery.target(
            "ck2_alpha"
        ).analysis.min_conformation_specific_seed_support
        == 4
    )
    assert config.discovery.target("ck2_alpha").analysis.ensemble_candidate_budget == 32
    assert (
        config.discovery.target(
            "ck2_alpha"
        ).analysis.conformation_specific_candidate_budget
        == 8
    )
    assert config.validation.rosetta.parallel_tasks == 8
    assert config.validation.rosetta.decoys_per_seed == 128
    assert config.validation.seeds == (120623, 120624, 120625, 120626)
    assert config.validation.refinement.seed_batch_sizes == (1, 1, 2)
    assert config.validation.funnel.ensemble_screening_budget == 10
    assert config.validation.funnel.confirmation_min_task_cells == 3
    assert config.design.objective == "single_supported_target"
    assert config.design.screen.parallel_tasks == 8
    assert config.design.finalists.parallel_tasks == 8
    assert config.discovery.cabsdock.filtering_count == 1000
    assert config.discovery.cabsdock.disulfide_ca_restraint_distance_A == 5.5
    assert config.discovery.cabsdock.disulfide_ca_restraint_weight == 1.0
    assert config.discovery.cabsdock.max_reconstructable_disulfide_ca_distance_A == 10.0
    assert config.discovery.cabsdock.min_models_for_selection == 10
    assert config.discovery.qualification.control_target_id == "ck2_alpha"
    assert config.discovery.qualification.control_receptor_ids == (
        "3Q04_A",
        "3QA0_A",
    )
    assert config.discovery.qualification.benchmark_receptor_id == "3Q9X_A"
    assert config.discovery.qualification.receptor_site_diagnostic_budget == 32
    assert config.discovery.qualification.topology_calibration_status == "qualified"
    assert config.discovery.qualification.topology_calibration_report == (
        PROJECT_ROOT / "outputs/discovery/topology_calibrations/"
        "ck2-alpha-ca-constrained-topology-20260804/"
        "topology_calibration_report.json"
    )
    assert config.validation.cg2all.expected_version == "1.2.0"
    assert config.validation.cg2all.representation == "CalphaSCModel"
    assert config.validation.cg2all.receptor_histidine_state == "HIE"
    assert config.validation.cg2all.processes == 2
    assert config.validation.handoff.poses_per_receptor_site == 4
    assert config.validation.refinement.prepack_seed == 3201
    assert config.validation.refinement.ranking_score == "reweighted_sc"
    assert config.validation.refinement.receptor_backbone_contact_A == 6.0
    assert config.validation.refinement.receptor_backbone_sequence_padding == 2
    assert config.design.sequence.mutable_positions == tuple(range(2, 11))
    assert config.design.combination.max_candidates == 96
    assert config.design.iteration.max_parents == 4
    assert config.design.iteration.max_total_mutations == 5
    assert config.design.iteration.max_candidates == 256
    assert config.design.finalists.max_candidates == 12
    assert config.design.finalists.max_md_candidates == 5
    assert config.design.finalists.seeds == ()
    assert tuple(path.name for path in config.source_files) == (
        "common.toml",
        "discovery.toml",
        "validation.toml",
        "design.toml",
    )
    assert "# --- common.toml ---" in config.source_snapshot_text


def test_environment_path_has_higher_priority(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    data_dir = tmp_path / "external-data"
    monkeypatch.setenv("VELA_DATA_DIR", str(data_dir))
    config = load_config(PROJECT_CONFIG_DIR)
    assert config.paths.data_dir == data_dir


def test_explicit_missing_config_fails(tmp_path: Path) -> None:
    with pytest.raises(ConfigError, match="config directory does not exist"):
        load_config(tmp_path / "missing")


@pytest.mark.parametrize(
    ("filename", "error_type"),
    [("validation.toml", ValidationError), ("design.toml", DesignError)],
)
def test_resolved_method_status_requires_report_and_hash(
    tmp_path: Path, filename: str, error_type: type[VelaError]
) -> None:
    config_dir = tmp_path / "configs"
    _copy_project_config(config_dir)
    path = config_dir / filename
    text = path.read_text(encoding="utf-8")
    if filename == "validation.toml":
        text = text.replace(
            'qualification_status = "qualified"',
            'qualification_status = "failed"',
            1,
        )
        text = text.replace(
            'qualification_report = "../outputs/validation/controls/4ib5-pc-chemistry-aware-local-control-20260804/qualification_analysis/qualification_report.json"',
            'qualification_report = "unresolved"',
            1,
        )
        text = text.replace(
            'qualification_report_sha256 = "69514733bc4b8b97bff0366ec29b3acb1491e5dc1e9d565b0f7072590d66e078"',
            'qualification_report_sha256 = "unresolved"',
            1,
        )
    else:
        text = text.replace(
            'qualification_status = "unresolved"',
            'qualification_status = "failed"',
            1,
        )
    path.write_text(text, encoding="utf-8")

    with pytest.raises(error_type, match="requires a report and SHA-256"):
        load_config(config_dir)


def test_project_file_recursively_overrides_selected_package_defaults(
    tmp_path: Path,
) -> None:
    config_dir = tmp_path / "configs"
    _copy_project_config(config_dir)
    discovery_path = config_dir / "discovery.toml"
    discovery_path.write_text(
        discovery_path.read_text(encoding="utf-8").replace(
            "[discovery.targets.ck2_alpha]",
            "[discovery.ensemble]\nmin_receptors_per_target = 1\n\n"
            "[discovery.targets.ck2_alpha]",
            1,
        ),
        encoding="utf-8",
    )

    config = load_config(config_dir)

    assert config.download.chunk_size_bytes == 1_048_576
    assert config.download.retries == 4
    assert config.audit.crystal_contacts.distance_A == 4.5
    assert config.audit.crystal_contacts.min_occupancy == 0.01
    assert config.preparation.altloc.preferred_label == "A"
    assert config.discovery.ensemble.min_receptors_per_target == 1
    assert config.discovery.ensemble.allowed_structure_states == ("apo", "apo_like")


def test_invalid_stage_default_override_is_rejected(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    _copy_project_config(config_dir)
    discovery_path = config_dir / "discovery.toml"
    discovery_path.write_text(
        discovery_path.read_text(encoding="utf-8").replace(
            "[discovery.targets.ck2_alpha]",
            "[discovery.ensemble]\nmin_receptors_per_target = 0\n\n"
            "[discovery.targets.ck2_alpha]",
            1,
        ),
        encoding="utf-8",
    )

    with pytest.raises(DiscoveryError, match="min_receptors_per_target must be"):
        load_config(config_dir)


def test_missing_required_project_file_is_rejected(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    _copy_project_config(config_dir)
    (config_dir / "design.toml").unlink()

    with pytest.raises(
        ConfigError, match=r"config directory is missing files: design.toml"
    ):
        load_config(config_dir)


def test_duplicate_project_field_is_rejected(tmp_path: Path) -> None:
    config_dir = tmp_path / "configs"
    _copy_project_config(config_dir)
    common_path = config_dir / "common.toml"
    common_path.write_text(
        common_path.read_text(encoding="utf-8")
        + """

[design]
method_id = "duplicate"
""",
        encoding="utf-8",
    )

    with pytest.raises(
        ConfigError,
        match=r"external config field is declared more than once: root.design.method_id",
    ):
        load_config(config_dir)
