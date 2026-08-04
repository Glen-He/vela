from pathlib import Path

import pytest

from vela.cli import run

PROJECT_CONFIG = Path(__file__).resolve().parents[2] / "configs"


def test_config_check_reports_all_planned_full_surface_receptors(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = run(["config", "check", "--config", str(PROJECT_CONFIG)])

    assert status == 0
    output = capsys.readouterr().out
    assert "Crystal contact audit: distance_A=4.5" in output
    assert "Altloc preferred label: A" in output
    assert "min_receptors_per_target=2" in output
    assert "replicas_dtemp=0.5, temperature=2->1" in output
    assert "CABS-dock sampling: formal_seeds=8, seed_workers=8" in output
    assert (
        "max_sites_per_task=64, max_pose_clusters_per_site=4, "
        "max_candidates_per_task=512" in output
    )
    assert "CABS-dock worker unit: seed batch" in output
    assert (
        "Discovery topology calibration: status=qualified, "
        "candidate_CA_thresholds_A=6,7,8,10, comparator_upper_A=12, "
        "active_CA_threshold_A=10, models_per_stratum=8" in output
    )
    assert "Stage 2 apo/apo-like blind discovery: 4" in output
    assert "Stage 3 stripped bound-state blind replication: 3" in output
    assert "Ligand: p15; chemistry_id=p15-free-n-amide-hie-ph7p4" in output
    assert "Total planned ligand full-surface receptor conformations: 7" in output
    assert "Currently prepared receptor-only bases: 5" in output
    assert "3 states; 1 configured local recovery controls" in output
    assert "1 guided; 2 environments" in output
    assert (
        "seeds_per_start=4, decoys_per_seed=128, total_decoys_per_start=512" in output
    )
    assert "Stage 4 sequence space: mutable_positions=9" in output
    assert "Stage 4 iterative neighborhood: max_parents=4" in output
    assert "Stage 4 flexible verification: max_candidates=12" in output


def test_design_status_separates_screen_and_flexible_gates(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = run(["design", "status", "--config", str(PROJECT_CONFIG)])

    assert status == 0
    output = capsys.readouterr().out
    assert "Stage 4 setup ready: true" in output
    assert "Stage 4 interface screen ready: false" in output
    assert "Stage 4 flexible finalist verification ready: false" in output
    assert "finalist_seeds_unresolved" in output


def test_preparation_status_reports_only_stage_one_data_boundary(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = run(["preparation", "status", "--config", str(PROJECT_CONFIG)])

    assert status == 0
    output = capsys.readouterr().out
    assert "Stage 1 preparation ready: true" in output
    assert "Receptors audited: 10" in output
    assert "Receptors prepared for Stage 2: 5" in output
    assert "method_not_qualified" not in output
    assert "seeds_unresolved" not in output
    assert "analysis_rules_unresolved" not in output


def test_discovery_status_exposes_unresolved_scientific_gates(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = run(
        [
            "discovery",
            "status",
            "--target",
            "ck2_alpha",
            "--config",
            str(PROJECT_CONFIG),
        ]
    )

    assert status == 0
    output = capsys.readouterr().out
    assert "Discovery target: ck2_alpha" in output
    assert "Discovery production ready: false" in output
    assert "chemistry_unresolved" not in output
    assert "method_not_qualified" in output
    assert "seeds_unresolved" not in output
    assert "analysis_rules_unresolved" not in output


def test_validation_status_separates_replication_and_refinement_gates(
    capsys: pytest.CaptureFixture[str],
) -> None:
    status = run(["validation", "status", "--config", str(PROJECT_CONFIG)])

    assert status == 0
    output = capsys.readouterr().out
    assert "Stage 3 independent setup ready: true" in output
    assert "Stage 3 bound-state blind replication ready: false" in output
    assert "Stage 3 candidate refinement ready: false" in output
    assert "global_method_not_qualified" in output
    assert "method_not_qualified" in output
