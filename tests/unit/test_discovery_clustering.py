from pathlib import Path

import pytest

from vela.discovery.analysis.clustering import analyze_sites
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.models import DiscoveryError, SiteAnalysisSettings
from vela.discovery.qualification.analysis import (
    best_native_coherent_site_evidence,
    selection_seed_recall,
)

SETTINGS = SiteAnalysisSettings(
    contact_jaccard_distance=0.5,
    position_distance_A=2.0,
    min_seed_support=2,
    min_receptor_support=2,
)


def _pose(
    pose_id: str,
    *,
    receptor: str,
    target: str,
    seed: int,
    contacts: frozenset[str],
    x: float,
    frame: str,
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
        ranking_score=-10.0,
        score_name="test_score",
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
    assert candidate.supported


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
    assert not any(site.supported for site in result.candidate_sites)


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


def test_native_site_support_requires_pairwise_compatible_poses() -> None:
    contacts = frozenset({"A:10", "A:11"})
    poses = {
        pose.pose_id: pose
        for pose in (
            _pose(
                "a_anchor",
                receptor="3Q04_A",
                target="ck2_alpha",
                seed=1,
                contacts=contacts,
                x=0.0,
                frame="3Q04_A",
            ),
            _pose(
                "b_left",
                receptor="3Q04_A",
                target="ck2_alpha",
                seed=2,
                contacts=contacts,
                x=-1.8,
                frame="3Q04_A",
            ),
            _pose(
                "c_right",
                receptor="3Q04_A",
                target="ck2_alpha",
                seed=3,
                contacts=contacts,
                x=1.8,
                frame="3Q04_A",
            ),
        )
    }

    seed_support, precision = best_native_coherent_site_evidence(
        recovered=set(poses),
        poses=poses,
        contact_limit=0.5,
        position_limit=2.0,
    )

    assert seed_support == 2
    assert precision == 1.0


def test_selection_seed_recall_allows_additional_site_only_seeds() -> None:
    recall = selection_seed_recall(
        sampling_seeds=(120666,),
        selection_site_seeds=(120664, 120666, 120667),
    )

    assert recall.fraction == 1.0
    assert recall.retained == (120666,)
    assert recall.missing == ()
    assert recall.site_only == (120664, 120667)


def test_selection_seed_recall_reports_missing_sampling_seeds() -> None:
    recall = selection_seed_recall(
        sampling_seeds=(120666, 120670),
        selection_site_seeds=(120664, 120666),
    )

    assert recall.fraction == 0.5
    assert recall.retained == (120666,)
    assert recall.missing == (120670,)
    assert recall.site_only == (120664,)


def test_selection_seed_recall_does_not_invent_success_without_sampling_hits() -> None:
    recall = selection_seed_recall(
        sampling_seeds=(),
        selection_site_seeds=(120664,),
    )

    assert recall.fraction == 0.0
    assert recall.retained == ()
    assert recall.missing == ()
    assert recall.site_only == (120664,)
