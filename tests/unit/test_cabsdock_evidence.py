import inspect
import io
import tarfile
from dataclasses import replace
from pathlib import Path

import pytest

from vela.core.provenance import atomic_write_json, atomic_write_text, sha256_file
from vela.discovery.models import CabsDockSettings, DiscoveryError, DiscoveryTask
from vela.discovery.qualification.evaluation import compare_poses_to_native
from vela.discovery.sampling.evidence import (
    CandidateSelectionSettings,
    collect_cabsdock_evidence,
)
from vela.discovery.sampling.execution import group_tasks_by_seed
from vela.discovery.sampling.materialization import CabsFrameIdentity
from vela.preparation.chemistry import (
    ChemistryDefinition,
    DisulfideBond,
    HistidineState,
)

P15_NAMES = (
    "CYS",
    "TRP",
    "MET",
    "SER",
    "PRO",
    "ARG",
    "HIS",
    "LEU",
    "GLY",
    "THR",
    "CYS",
)


def _atom_line(
    serial: int,
    *,
    atom: str,
    residue: str,
    chain: str,
    sequence_number: int,
    x: float,
    y: float,
) -> str:
    return (
        f"ATOM  {serial:5d} {atom:^4s} {residue:>3s} {chain}{sequence_number:4d}    "
        f"{x:8.3f}{y:8.3f}{0.0:8.3f}{1.0:6.2f}{0.0:6.2f}\n"
    )


def _model(number: int, *, ring_distance_A: float) -> str:
    lines = [f"MODEL     {number:4d}\n"]
    serial = 1
    for sequence_number in range(1, 61):
        x = float((sequence_number - 1) % 10) * 3.8
        y = float((sequence_number - 1) // 10) * 3.8
        for atom in ("CA", "SC"):
            lines.append(
                _atom_line(
                    serial,
                    atom=atom,
                    residue="ALA",
                    chain="A",
                    sequence_number=sequence_number,
                    x=x,
                    y=y,
                )
            )
            serial += 1
    for sequence_number, residue in enumerate(P15_NAMES, 1):
        x = (
            ring_distance_A
            if sequence_number == len(P15_NAMES)
            else 0.5 * (sequence_number - 1)
        )
        for atom in ("CA", "SC"):
            lines.append(
                _atom_line(
                    serial,
                    atom=atom,
                    residue=residue,
                    chain="P",
                    sequence_number=sequence_number,
                    x=x,
                    y=4.0,
                )
            )
            serial += 1
    lines.append("ENDMDL\n")
    return "".join(lines)


def _settings() -> CabsDockSettings:
    return CabsDockSettings(
        executable=Path("/test/CABSdock"),
        source_dir=Path("/test/cabs-source"),
        source_revision="1" * 40,
        patch_file=Path("/test/cabs.patch"),
        seed_workers=1,
        peptide_secondary_structure="CCCCCCCCCCC",
        mc_annealing=1,
        mc_cycles=1,
        mc_steps=1,
        replicas=1,
        replicas_dtemp=0.5,
        temperature_initial=2.0,
        temperature_final=1.0,
        binding_interactions=1.0,
        protein_restraint_gap=5,
        protein_restraint_min_A=5.0,
        protein_restraint_max_A=15.0,
        filtering_count=2,
        clustering_medoids=2,
        clustering_iterations=1,
        trajectory_contact_ca_threshold_A=10.0,
        max_disulfide_ca_distance_A=8.0,
        min_models_for_selection=2,
        selection_contact_jaccard_distance=0.8,
        selection_position_distance_A=12.0,
        pose_clustering_rmsd_A=4.0,
    )


def _chemistry() -> ChemistryDefinition:
    return ChemistryDefinition(
        ligand_id="test-peptide",
        chemistry_id="p15-test",
        sequence="CWMSPRHLGTC",
        chirality="L",
        disulfide_bonds=(DisulfideBond(1, 11),),
        n_terminus="NH3+",
        c_terminus="CONH2",
        target_ph=7.4,
        net_charge=2,
        histidines=(HistidineState(7, "HIE"),),
        other_modifications_status="none",
        other_modifications=(),
        decision_sources=("test",),
    )


def _selection() -> CandidateSelectionSettings:
    return CandidateSelectionSettings(0.8, 12.0, 4.0)


def _fake_materializer(monkeypatch: pytest.MonkeyPatch) -> None:
    def materialize(
        *,
        archive_path: Path,
        identities: tuple[CabsFrameIdentity, ...],
        task_dir: Path,
        settings: CabsDockSettings,
    ) -> tuple[Path, Path]:
        del archive_path, settings
        selection_path = task_dir / "output_data" / "trajectory_candidate_frames.json"
        output_path = task_dir / "output_pdbs" / "trajectory_candidates.pdb"
        atomic_write_json(
            selection_path,
            {
                "frames": [
                    {"replica": identity.replica, "model": identity.model}
                    for identity in identities
                ]
            },
        )
        atomic_write_text(
            output_path,
            "".join(
                _model(index, ring_distance_A=5.0)
                for index, _ in enumerate(identities, 1)
            )
            + "END\n",
        )
        return selection_path, output_path

    monkeypatch.setattr(
        "vela.discovery.sampling.evidence.materialize_cabs_frames", materialize
    )


def _task(reference_path: Path) -> DiscoveryTask:
    return DiscoveryTask(
        task_id="3Q04_A__seed_17",
        receptor_id="3Q04_A",
        target="ck2_alpha",
        receptor_path=reference_path,
        receptor_sha256=sha256_file(reference_path),
        chemistry_id="p15-test",
        method_id="cabsdock-test",
        adapter_id="vela-test",
        seed=17,
        evidence_category="test_blind_sampling",
    )


def _inputs(
    tmp_path: Path, *, filtered_distances: tuple[float, ...]
) -> tuple[Path, Path]:
    task_dir = tmp_path / "task"
    output_dir = task_dir / "output_pdbs"
    reference_path = tmp_path / "reference.pdb"
    atomic_write_text(reference_path, _model(1, ring_distance_A=5.0) + "END\n")
    atomic_write_text(
        output_dir / "top1000.pdb",
        "".join(
            _model(index, ring_distance_A=distance)
            for index, distance in enumerate(filtered_distances, 1)
        )
        + "END\n",
    )
    for index, distance in enumerate((5.0, 6.0)):
        atomic_write_text(
            output_dir / f"model_{index}.pdb",
            _model(1, ring_distance_A=distance) + "END\n",
        )
    _write_cabs_archive(task_dir / "test.cbs", filtered_distances=filtered_distances)
    return task_dir, reference_path


def _seq_line(number: int, *, residue: str, chain: str) -> str:
    return f"{number:5d}   {residue:>3s} {chain}  1  0.00\n"


def _traf_coordinates(*, count: int, endpoint_distance_A: float | None = None) -> str:
    values = [0, 0, 0]
    for index in range(count):
        x = index
        if endpoint_distance_A is not None and index == count - 1:
            x = round(endpoint_distance_A / 0.61)
        values.extend((x, 0, 0))
    values.extend((0, 0, 0))
    return " ".join(map(str, values)) + "\n"


def _write_cabs_archive(path: Path, *, filtered_distances: tuple[float, ...]) -> None:
    sequence = "".join(
        _seq_line(index, residue="ALA", chain="A") for index in range(1, 61)
    ) + "".join(
        _seq_line(index, residue=residue, chain="P")
        for index, residue in enumerate(P15_NAMES, 1)
    )
    trajectory: list[str] = []
    for model, distance in enumerate(filtered_distances, 1):
        interaction_energy = -13.0 + model
        trajectory.extend(
            (
                f"{model} 62 -100.0 {interaction_energy:.4f} 1.0 1\n",
                _traf_coordinates(count=60),
                f"{model} 13 {interaction_energy:.4f} -5.0 1.0 1\n",
                _traf_coordinates(count=11, endpoint_distance_A=distance),
            )
        )
    with tarfile.open(path, mode="w:gz") as archive:
        for name, content in (("SEQ", sequence), ("TRAF", "".join(trajectory))):
            payload = content.encode("utf-8")
            member = tarfile.TarInfo(name)
            member.size = len(payload)
            archive.addfile(member, io.BytesIO(payload))


def test_seed_batches_keep_receptors_sequential_within_each_worker(
    tmp_path: Path,
) -> None:
    reference = tmp_path / "reference.pdb"
    atomic_write_text(reference, "test\n")
    first = replace(_task(reference), task_id="receptor_a__seed_11", seed=11)
    second = replace(first, task_id="receptor_b__seed_11", receptor_id="receptor_b")
    third = replace(first, task_id="receptor_a__seed_22", seed=22)
    fourth = replace(third, task_id="receptor_b__seed_22", receptor_id="receptor_b")

    batches = group_tasks_by_seed((first, second, third, fourth))

    assert [[task.seed for task in batch] for batch in batches] == [
        [11, 11],
        [22, 22],
    ]
    assert [[task.receptor_id for task in batch] for batch in batches] == [
        ["3Q04_A", "receptor_b"],
        ["3Q04_A", "receptor_b"],
    ]
    with pytest.raises(DiscoveryError, match="exactly one target"):
        group_tasks_by_seed((first, replace(second, target="ck2_alpha_prime")))


def test_collects_site_first_candidates_and_keeps_top10_baseline(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _fake_materializer(monkeypatch)
    task_dir, reference_path = _inputs(tmp_path, filtered_distances=(5.0, 6.0))

    evidence = collect_cabsdock_evidence(
        task=_task(reference_path),
        task_dir=task_dir,
        reference_path=reference_path,
        reference_receptor_id="3Q04_A",
        chemistry=_chemistry(),
        settings=_settings(),
        selection=_selection(),
    )

    assert evidence.filtered_model_count == 2
    assert evidence.topology_feasible_model_count == 2
    assert evidence.topology_feasible_fraction == 1.0
    assert evidence.selection_status == "completed"
    assert evidence.selection_failure_reasons == ()
    assert evidence.sampling_model_count == 2
    assert evidence.contacting_topology_feasible_model_count == 2
    assert evidence.pose_cluster_count == evidence.selected_pose_cluster_count == 1
    assert len(evidence.poses) == 1
    assert len(evidence.baseline_poses) == 2
    assert evidence.poses[0].qc_status == "passed"
    assert evidence.poses[0].score_name == "cabsdock_interaction_energy"
    assert evidence.poses[0].contact_residues
    assert evidence.poses[0].model_path.name == "trajectory_candidates.pdb"
    assert evidence.baseline_poses[0].model_path.name == "model_0.pdb"


def test_small_selection_pool_is_skipped_as_a_technical_selection_limit(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _fake_materializer(monkeypatch)
    task_dir, reference_path = _inputs(tmp_path, filtered_distances=(5.0, 12.0))

    evidence = collect_cabsdock_evidence(
        task=_task(reference_path),
        task_dir=task_dir,
        reference_path=reference_path,
        reference_receptor_id="3Q04_A",
        chemistry=_chemistry(),
        settings=_settings(),
        selection=_selection(),
    )

    assert evidence.topology_feasible_model_count == 1
    assert evidence.selection_status == ("skipped_insufficient_models_for_selection")
    assert evidence.selection_failure_reasons == (
        "skipped_insufficient_models_for_selection",
    )
    assert evidence.poses == ()
    assert evidence.contacting_topology_feasible_model_count == 1


def test_candidate_collection_has_no_native_reference_input() -> None:
    parameters = inspect.signature(collect_cabsdock_evidence).parameters

    assert "native_pair_path" not in parameters
    assert "native" not in parameters


def test_rejects_output_count_that_differs_from_frozen_method(tmp_path: Path) -> None:
    task_dir, reference_path = _inputs(tmp_path, filtered_distances=(5.0,))

    with pytest.raises(DiscoveryError, match="filtered model count differs"):
        collect_cabsdock_evidence(
            task=_task(reference_path),
            task_dir=task_dir,
            reference_path=reference_path,
            reference_receptor_id="3Q04_A",
            chemistry=_chemistry(),
            settings=replace(_settings(), filtering_count=2),
            selection=_selection(),
        )


def test_native_comparison_uses_receptor_aligned_ligand_coordinates(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    _fake_materializer(monkeypatch)
    pair_path = tmp_path / "native_pair.pdb"
    pose_path = tmp_path / "sampled_pose.pdb"
    structure = _model(1, ring_distance_A=5.0) + "END\n"
    atomic_write_text(pair_path, structure)
    atomic_write_text(pose_path, structure)

    pose = collect_cabsdock_evidence(
        task=_task(reference_path=pose_path),
        task_dir=_inputs(tmp_path / "evidence", filtered_distances=(5.0, 6.0))[0],
        reference_path=pair_path,
        reference_receptor_id="3Q04_A",
        chemistry=_chemistry(),
        settings=_settings(),
        selection=_selection(),
    ).poses[0]
    metrics = compare_poses_to_native(
        poses=(pose,),
        native_pair_path=pair_path,
        peptide_sequence="CWMSPRHLGTC",
        contact_ca_threshold_A=10.0,
    )[pose.pose_id]

    assert metrics.ligand_ca_rmsd_A == pytest.approx(0.0)
    assert metrics.ligand_centroid_distance_A == pytest.approx(0.0)
    assert metrics.native_receptor_contact_fraction == pytest.approx(1.0)
