from dataclasses import replace
from pathlib import Path

import pytest

from vela.discovery.analysis.cluster_engine import complete_linkage
from vela.discovery.analysis.clustering import analyze_sites
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.models import DiscoveryError, SiteAnalysisSettings

SETTINGS = SiteAnalysisSettings(
    contact_jaccard_distance=0.5,
    position_distance_A=2.0,
    min_seed_support=2,
    min_receptor_support=2,
    min_conformation_specific_seed_support=2,
    ensemble_candidate_budget=1,
    conformation_specific_candidate_budget=1,
)


def test_complete_linkage_preserves_diameter_threshold_and_order() -> None:
    positions = {"a": 0.0, "b": 0.2, "c": 0.8, "d": 2.0, "e": 2.4}

    clusters = complete_linkage(
        tuple(reversed(tuple(positions))),
        distance=lambda first, second: abs(positions[first] - positions[second]),
        identity=lambda item: item,
    )

    assert clusters == [("a", "b", "c"), ("d", "e")]


def _pose(
    pose_id: str,
    *,
    receptor: str,
    target: str,
    seed: int,
    contacts: frozenset[str],
    x: float,
    frame: str,
    score: float = -10.0,
    score_name: str = "test_score",
) -> PoseEvidence:
    return PoseEvidence(
        task_id=f"{receptor}__seed_{seed}",
        pose_id=pose_id,
        receptor_id=receptor,
        target=target,
        seed=seed,
        model_path=Path(f"{pose_id}.pdb"),
        model_sha256="a" * 64,
        model_index=1,
        contact_residues=contacts,
        local_position=(x, 0.0, 0.0),
        coordinate_frame_id=frame,
        ranking_score=score,
        score_name=score_name,
        qc_status="passed",
    )


def test_site_support_counts_independent_seeds_and_receptor_conformations() -> None:
    contacts = frozenset({"A:10", "A:11"})
    poses = (
        _pose(
            "a_seed1_first",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=1,
            contacts=contacts,
            x=0.0,
            frame="ck2_alpha_reference_v1",
        ),
        _pose(
            "a_seed1_duplicate",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=1,
            contacts=contacts,
            x=0.1,
            frame="ck2_alpha_reference_v1",
        ),
        _pose(
            "a_seed2",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=2,
            contacts=frozenset({"A:10", "A:11", "A:12"}),
            x=0.3,
            frame="ck2_alpha_reference_v1",
        ),
        _pose(
            "a_single_seed_far",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=1,
            contacts=frozenset({"A:100"}),
            x=20.0,
            frame="ck2_alpha_reference_v1",
        ),
        _pose(
            "b_seed1",
            receptor="3QA0_A",
            target="ck2_alpha",
            seed=1,
            contacts=contacts,
            x=0.2,
            frame="ck2_alpha_reference_v1",
        ),
        _pose(
            "b_seed2",
            receptor="3QA0_A",
            target="ck2_alpha",
            seed=2,
            contacts=contacts,
            x=0.4,
            frame="ck2_alpha_reference_v1",
        ),
    )

    result = analyze_sites(poses=poses, settings=SETTINGS)

    assert len(result.receptor_sites) == 3
    supported_receptor_sites = [
        site for site in result.receptor_sites if site.supported
    ]
    assert len(supported_receptor_sites) == 2
    assert all(site.supporting_seeds == (1, 2) for site in supported_receptor_sites)
    assert len(result.candidate_sites) == 1
    candidate = result.candidate_sites[0]
    assert candidate.receptor_ids == ("3Q04_A", "3QA0_A")
    assert candidate.receptor_support == 2
    assert candidate.evidence_tier == "ensemble_consensus"
    assert candidate.rank_within_tier == 1
    assert candidate.minimum_seed_support == 2
    assert candidate.total_seed_support == 4
    assert candidate.handoff_eligible


def test_analysis_never_combines_subtypes() -> None:
    contacts = frozenset({"A:10", "A:11"})
    poses = tuple(
        _pose(
            f"{target}_{seed}",
            receptor=receptor,
            target=target,
            seed=seed,
            contacts=contacts,
            x=0.0,
            frame=frame,
        )
        for receptor, target, frame in (
            ("3Q04_A", "ck2_alpha", "ck2_alpha_reference_v1"),
            ("5YF9_X", "ck2_alpha_prime", "ck2_alpha_prime_reference_v1"),
        )
        for seed in (1, 2)
    )

    result = analyze_sites(poses=poses, settings=SETTINGS)

    assert {site.target for site in result.candidate_sites} == {
        "ck2_alpha",
        "ck2_alpha_prime",
    }
    assert all(site.receptor_support == 1 for site in result.candidate_sites)
    assert all(
        site.evidence_tier == "conformation_specific" for site in result.candidate_sites
    )
    assert all(site.handoff_eligible for site in result.candidate_sites)


def test_analysis_rejects_inconsistent_coordinate_frames_within_target() -> None:
    contacts = frozenset({"A:10"})
    poses = (
        _pose(
            "one",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=1,
            contacts=contacts,
            x=0.0,
            frame="frame_one",
        ),
        _pose(
            "two",
            receptor="3QA0_A",
            target="ck2_alpha",
            seed=1,
            contacts=contacts,
            x=0.0,
            frame="frame_two",
        ),
    )

    with pytest.raises(DiscoveryError, match="one aligned coordinate frame"):
        analyze_sites(poses=poses, settings=SETTINGS)


def test_candidate_ranking_enforces_tier_budget() -> None:
    frame = "ck2_alpha_reference_v1"
    poses = tuple(
        _pose(
            f"{receptor}_{site}_{seed}",
            receptor=receptor,
            target="ck2_alpha",
            seed=seed,
            contacts=contacts,
            x=x,
            frame=frame,
        )
        for receptor in ("3Q04_A", "3QA0_A")
        for site, seeds, contacts, x in (
            ("strong", (1, 2, 3), frozenset({"A:10", "A:11"}), 0.0),
            ("weak", (1, 2), frozenset({"A:100", "A:101"}), 10.0),
        )
        for seed in seeds
    )

    result = analyze_sites(poses=poses, settings=SETTINGS)
    ensemble = sorted(
        (
            candidate
            for candidate in result.candidate_sites
            if candidate.evidence_tier == "ensemble_consensus"
        ),
        key=lambda candidate: candidate.rank_within_tier,
    )

    assert [candidate.minimum_seed_support for candidate in ensemble] == [3, 2]
    assert [candidate.handoff_eligible for candidate in ensemble] == [True, False]


def test_candidate_score_tiebreak_is_invariant_to_receptor_score_scale() -> None:
    frame = "ck2_alpha_reference_v1"

    def poses(*, transform_second: bool) -> tuple[PoseEvidence, ...]:
        rows: list[PoseEvidence] = []
        for receptor in ("3Q04_A", "3QA0_A"):
            for site, contacts, x, score in (
                ("first", frozenset({"A:10"}), 0.0, -30.0),
                ("second", frozenset({"A:100"}), 20.0, -20.0),
                ("third", frozenset({"A:200"}), 40.0, -10.0),
            ):
                receptor_score = score
                if receptor == "3QA0_A" and transform_second:
                    receptor_score = score * 100.0 + 1_000_000.0
                for seed in (1, 2):
                    rows.append(
                        _pose(
                            f"{receptor}_{site}_{seed}",
                            receptor=receptor,
                            target="ck2_alpha",
                            seed=seed,
                            contacts=contacts,
                            x=x,
                            frame=frame,
                            score=receptor_score,
                        )
                    )
        return tuple(rows)

    original = analyze_sites(poses=poses(transform_second=False), settings=SETTINGS)
    transformed = analyze_sites(poses=poses(transform_second=True), settings=SETTINGS)

    original_metrics = {
        candidate.candidate_id: (
            candidate.rank_within_tier,
            candidate.median_receptor_score_quantile,
        )
        for candidate in original.candidate_sites
    }
    transformed_metrics = {
        candidate.candidate_id: (
            candidate.rank_within_tier,
            candidate.median_receptor_score_quantile,
        )
        for candidate in transformed.candidate_sites
    }
    assert transformed_metrics == original_metrics


def test_candidate_analysis_rejects_mixed_score_identities_per_receptor() -> None:
    poses = (
        _pose(
            "first",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=1,
            contacts=frozenset({"A:10"}),
            x=0.0,
            frame="ck2_alpha_reference_v1",
            score_name="first_score",
        ),
        _pose(
            "second",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=2,
            contacts=frozenset({"A:10"}),
            x=0.1,
            frame="ck2_alpha_reference_v1",
            score_name="second_score",
        ),
    )

    with pytest.raises(DiscoveryError, match="share one ranking score"):
        analyze_sites(poses=poses, settings=SETTINGS)


def test_conformation_specific_tier_requires_stronger_seed_support() -> None:
    settings = replace(SETTINGS, min_conformation_specific_seed_support=4)
    poses = tuple(
        _pose(
            f"site_{site}_{seed}",
            receptor="3Q04_A",
            target="ck2_alpha",
            seed=seed,
            contacts=contacts,
            x=x,
            frame="ck2_alpha_reference_v1",
        )
        for site, seeds, contacts, x in (
            ("strong", (1, 2, 3, 4), frozenset({"A:10", "A:11"}), 0.0),
            ("weak", (1, 2, 3), frozenset({"A:100", "A:101"}), 10.0),
        )
        for seed in seeds
    )

    result = analyze_sites(poses=poses, settings=settings)
    tiers = {candidate.evidence_tier for candidate in result.candidate_sites}

    assert tiers == {"conformation_specific", "insufficient_evidence"}
    assert sum(candidate.handoff_eligible for candidate in result.candidate_sites) == 1
