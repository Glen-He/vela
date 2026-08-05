from dataclasses import replace
from pathlib import Path

import pytest

from vela.config import load_config
from vela.design.finalists.evidence import PairedSummary, candidate_evidence
from vela.design.models import DesignError, DesignTemplate, ScreenTask
from vela.design.readiness import assess_design_readiness
from vela.design.scores import SCREEN_SCORE_COLUMNS, screen_metrics
from vela.design.screening.execution import task_histidine_pose_indices
from vela.design.screening.inputs import validate_objective
from vela.design.screening.planning import candidate_record
from vela.design.screening.records import candidate_from_record
from vela.design.sequence.library import (
    combination_library,
    first_generation_candidate,
    flexibility_reasons,
    sequence_facts,
    systematic_single_library,
)
from vela.design.sequence.neighborhood import iteration_library
from vela.validation.models import ValidationError
from vela.validation.scores import RosettaScoreRow

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROJECT_CONFIG = PROJECT_ROOT / "configs"


def test_systematic_single_library_covers_the_declared_sequence_space() -> None:
    config = load_config(PROJECT_CONFIG)

    candidates = systematic_single_library(
        chemistry=config.chemistry, settings=config.design
    )

    assert len(candidates) == 162
    assert {item.mutation_positions[0] for item in candidates} == set(range(2, 11))
    assert any(item.mutation_string == "P5A" for item in candidates)
    assert any(item.mutation_string == "G9A" for item in candidates)
    assert any(item.mutation_string == "W2H" for item in candidates)
    assert all(
        item.sequence[0] == "C" and item.sequence[-1] == "C" for item in candidates
    )
    assert all(item.sequence.count("C") == 2 for item in candidates)


def test_combination_library_rebuilds_and_caps_complete_sequences() -> None:
    config = load_config(PROJECT_CONFIG)

    candidates = combination_library(
        chemistry=config.chemistry,
        settings=config.design,
        substitutions={2: ("F", "H"), 5: ("A",), 9: ("S",)},
    )

    assert candidates
    assert len(candidates) <= config.design.combination.max_candidates
    assert {item.mutation_count for item in candidates} == {2}
    assert len({item.candidate_id for item in candidates}) == len(candidates)


def test_combination_resource_cap_balances_mutation_counts_and_positions() -> None:
    config = load_config(PROJECT_CONFIG)
    substitutions = {
        position: tuple(
            residue
            for residue in "AFK"
            if residue != config.chemistry.sequence[position - 1]
        )
        for position in config.design.sequence.mutable_positions
    }

    candidates = combination_library(
        chemistry=config.chemistry,
        settings=config.design,
        substitutions=substitutions,
    )

    assert len(candidates) == config.design.combination.max_candidates
    assert {item.mutation_count for item in candidates} == {2}
    covered = {position for item in candidates for position in item.mutation_positions}
    assert covered == set(config.design.sequence.mutable_positions)


def test_iteration_library_uses_one_edit_lineage_and_preserves_deferred_candidates() -> (
    None
):
    config = load_config(PROJECT_CONFIG)
    parents = (
        first_generation_candidate(
            chemistry=config.chemistry,
            settings=config.design,
            sequence="CFMSARHLGTC",
            design_round="combination",
            proposal_source="test",
        ),
        first_generation_candidate(
            chemistry=config.chemistry,
            settings=config.design,
            sequence="CHMSPRHLSTC",
            design_round="combination",
            proposal_source="test",
        ),
    )

    library = iteration_library(
        chemistry=config.chemistry,
        settings=config.design,
        parents=parents,
    )

    assert len(library.selected) == config.design.iteration.max_candidates
    assert library.deferred
    assert {item.generation for item in library.all_candidates} == {2}
    assert max(item.mutation_count for item in library.all_candidates) <= 5
    assert {
        parent.edit_type for item in library.all_candidates for parent in item.parents
    } == {
        "add",
        "revert",
        "swap",
    }
    assert {
        parent.candidate_id for item in library.selected for parent in item.parents
    } == {item.candidate_id for item in parents}
    assert all(
        parent.candidate_id in item.ancestor_candidate_ids
        for item in library.all_candidates
        for parent in item.parents
    )

    child = library.selected[0]
    assert candidate_from_record(raw=candidate_record(child), config=config) == child
    with pytest.raises(DesignError, match="same-generation"):
        iteration_library(
            chemistry=config.chemistry,
            settings=config.design,
            parents=(parents[0], child),
        )


def test_candidate_charge_and_transparent_sequence_facts_are_derived() -> None:
    config = load_config(PROJECT_CONFIG)
    candidate = first_generation_candidate(
        chemistry=config.chemistry,
        settings=config.design,
        sequence="CKMSPRHLGTC",
        design_round="single",
        proposal_source="test",
    )

    facts = sequence_facts(candidate.sequence)

    assert candidate.mutation_string == "W2K"
    assert candidate.net_charge == 3
    assert facts.methionine_count == 1
    assert facts.deamidation_motif_count == 0
    assert facts.aspartimide_motif_count == 0


def test_fixed_backbone_flexibility_reasons_are_deterministic() -> None:
    config = load_config(PROJECT_CONFIG)
    candidate = first_generation_candidate(
        chemistry=config.chemistry,
        settings=config.design,
        sequence="CKMSGRHLGTC",
        design_round="combination",
        proposal_source="test",
    )

    reasons = flexibility_reasons(
        reference=config.chemistry.sequence,
        candidate=candidate,
        disulfide_positions=frozenset({1, 11}),
    )

    assert reasons == tuple(sorted(reasons))
    assert "adjacent_to_disulfide" in reasons
    assert "glycine_proline_change" in reasons
    assert "charge_class_change" in reasons


def test_screen_metric_contract_requires_every_frozen_rosetta_column() -> None:
    scores = {
        column: float(index)
        for index, column in enumerate(SCREEN_SCORE_COLUMNS.values())
    }
    metrics = screen_metrics(RosettaScoreRow("model_1", scores))

    assert metrics.dG_separated == 0.0
    scores.pop("sc_value")
    with pytest.raises(ValidationError, match="sc_value"):
        screen_metrics(RosettaScoreRow("model_2", scores))


def test_histidine_indices_follow_each_paired_sequence_state() -> None:
    config = load_config(PROJECT_CONFIG)
    candidate = first_generation_candidate(
        chemistry=config.chemistry,
        settings=config.design,
        sequence="CHMSPRHLGTC",
        design_round="single",
        proposal_source="test",
    )
    template = DesignTemplate(
        "template_001",
        "positive",
        "cluster_001",
        "candidate_001",
        "3Q04_A",
        "ck2_alpha",
        PROJECT_ROOT / "unused.pdb",
        "0" * 64,
        331,
        (10, 338),
    )
    wt_task = ScreenTask(
        "screen_0000001",
        "pair_0000001",
        "wt",
        candidate,
        template,
        4001,
        PROJECT_ROOT / "unused.resfile",
        "0" * 64,
        "wt_context_001",
    )
    mutant_task = ScreenTask(
        "screen_0000002",
        "pair_0000001",
        "mutant",
        candidate,
        template,
        4001,
        PROJECT_ROOT / "unused.resfile",
        "0" * 64,
        None,
    )

    assert task_histidine_pose_indices(config=config, task=wt_task) == (10, 338)
    assert task_histidine_pose_indices(config=config, task=mutant_task) == (
        10,
        333,
        338,
    )


def test_stage4_status_separates_screen_and_flexible_gates() -> None:
    config = load_config(PROJECT_CONFIG)

    readiness = assess_design_readiness(config)
    issue_codes = {item.code for item in readiness.issues}

    assert readiness.setup_ready
    assert not readiness.screen_ready
    assert not readiness.finalist_ready
    assert "method_not_qualified" in issue_codes
    assert "seeds_unresolved" in issue_codes
    assert "finalist_seeds_unresolved" in issue_codes
    assert "finalist_rules_unresolved" in issue_codes


def test_single_supported_target_accepts_exactly_one_positive_subtype() -> None:
    config = load_config(PROJECT_CONFIG)
    alpha = DesignTemplate(
        "template_001",
        "positive",
        "cluster_001",
        "candidate_001",
        "3Q04_A",
        "ck2_alpha",
        PROJECT_ROOT / "unused.pdb",
        "0" * 64,
        331,
        (),
    )
    prime = replace(
        alpha,
        template_id="template_002",
        cluster_id="cluster_002",
        receptor_id="5YF9_X",
        target="ck2_alpha_prime",
    )

    validate_objective(config=config, templates=(alpha,))
    validate_objective(config=config, templates=(prime,))
    with pytest.raises(DesignError, match="exactly one CK2 subtype"):
        validate_objective(config=config, templates=(alpha, prime))
    with pytest.raises(DesignError, match="evidence_role must be positive"):
        replace(alpha, evidence_role="negative")
    with pytest.raises(DesignError, match="unsupported design objective"):
        replace(config.design, objective="pan_ck2")


def test_flexible_evidence_rejects_a_candidate_without_each_template_seed_support() -> (
    None
):
    config = load_config(PROJECT_CONFIG)
    finalists = replace(
        config.design.finalists,
        seeds=(5101, 5103),
        calibrated=True,
        min_passed_decoy_fraction=0.5,
        min_successful_seeds=2,
        max_positive_median_ranking_delta=0.0,
        max_positive_median_interface_delta=0.0,
    )
    configured = replace(config, design=replace(config.design, finalists=finalists))
    candidate = first_generation_candidate(
        chemistry=config.chemistry,
        settings=config.design,
        sequence="CFMSPRHLGTC",
        design_round="single",
        proposal_source="test",
    )
    pairs = (
        PairedSummary(
            "pair_1",
            candidate,
            "template_alpha",
            "ck2_alpha",
            "positive",
            5101,
            1.0,
            1.0,
            True,
            "matched_pose",
            "wt_pose_1",
            "mutant_pose_1",
            -1.0,
            -0.5,
            -0.2,
        ),
        PairedSummary(
            "pair_2",
            candidate,
            "template_alpha",
            "ck2_alpha",
            "positive",
            5103,
            1.0,
            1.0,
            True,
            "matched_pose",
            "wt_pose_2",
            "mutant_pose_2",
            -0.8,
            -0.3,
            -0.1,
        ),
        PairedSummary(
            "pair_3",
            candidate,
            "template_prime",
            "ck2_alpha_prime",
            "positive",
            5101,
            1.0,
            0.2,
            False,
            "unmatched_wt_pose",
            "wt_pose_3",
            None,
            None,
            None,
            None,
        ),
    )

    evidence = candidate_evidence(config=configured, pairs=pairs)

    assert len(evidence) == 1
    assert evidence[0].status == "rejected"
    assert "template_seed_support:template_prime" in evidence[0].failed_gates
