from dataclasses import replace
from pathlib import Path

import gemmi
import pytest

from vela.config import load_config
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    vela_software_identity,
)
from vela.discovery.analysis.clustering import candidate_analysis_contract
from vela.discovery.analysis.pose_table import POSE_FIELDS, read_pose_evidence
from vela.discovery.analysis.workflow import analyze_discovery_run
from vela.discovery.models import DiscoveryError, SiteAnalysisSettings
from vela.validation.models import ValidationError
from vela.validation.refinement.handoff_plan import (
    build_handoff_tasks,
    candidate_evidence_records,
)

PROJECT_CONFIG = Path(__file__).resolve().parents[2] / "configs"


def _write_pose_model(
    *, receptor_path: Path, destination: Path, peptide_offset: float
) -> None:
    structure = gemmi.read_structure(str(receptor_path))
    peptide = gemmi.Chain("P")
    residue_names = (
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
    for index, name in enumerate(residue_names, 1):
        residue = gemmi.Residue()
        residue.name = name
        residue.seqid = gemmi.SeqId(index, " ")
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(float(index), peptide_offset, 0.0)
        residue.add_atom(atom)
        peptide.add_residue(residue)
    structure[0].add_chain(peptide)
    destination.parent.mkdir(parents=True, exist_ok=True)
    atomic_write_text(destination, structure.make_pdb_string())


def test_pose_table_hashes_each_shared_model_once(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    run_dir = tmp_path / "run"
    model = run_dir / "models" / "shared.pdb"
    atomic_write_text(model, "MODEL shared\nEND\n")
    digest = sha256_file(model)
    rows: list[str] = []
    for index in (1, 2):
        values = {
            "task_id": "task_001",
            "pose_id": f"pose_{index:03d}",
            "receptor_id": "3Q04_A",
            "target": "ck2_alpha",
            "seed": "120623",
            "model_path": model.relative_to(run_dir).as_posix(),
            "model_sha256": digest,
            "model_index": str(index),
            "contact_residues": "A:10",
            "local_x_A": "0.0",
            "local_y_A": "0.0",
            "local_z_A": "0.0",
            "coordinate_frame_id": "3Q04_A",
            "ranking_score": "-1.0",
            "score_name": "interaction_energy",
            "qc_status": "passed",
        }
        rows.append("\t".join(values[field] for field in POSE_FIELDS))
    table = run_dir / "pose_evidence.tsv"
    atomic_write_text(
        table,
        "\t".join(POSE_FIELDS) + "\n" + "\n".join(rows) + "\n",
    )
    calls = 0

    def counting_sha256(path: Path) -> str:
        nonlocal calls
        calls += 1
        return digest

    monkeypatch.setattr(
        "vela.discovery.analysis.pose_table.sha256_file", counting_sha256
    )

    poses = read_pose_evidence(path=table, run_dir=run_dir)

    assert len(poses) == 2
    assert calls == 1


def test_completed_normalized_run_produces_supported_candidate(tmp_path: Path) -> None:
    run_dir = tmp_path / "outputs" / "runs" / "synthetic"
    models_dir = run_dir / "models"
    task_rows: list[dict[str, JsonValue]] = []
    evidence_rows: list[str] = []
    config = load_config(PROJECT_CONFIG)
    for receptor, offset in (("3Q04_A", 0.0), ("3QA0_A", 0.2)):
        for seed in (11, 22, 33, 44):
            task_id = f"{receptor}__seed_{seed}"
            pose_id = f"{task_id}__pose_001"
            model = models_dir / f"{pose_id}.pdb"
            _write_pose_model(
                receptor_path=(
                    config.paths.data_dir / "receptors" / "prepared" / f"{receptor}.cif"
                ),
                destination=model,
                peptide_offset=offset + seed / 100.0,
            )
            task_rows.append(
                {
                    "task_id": task_id,
                    "receptor_id": receptor,
                    "target": "ck2_alpha",
                    "seed": seed,
                    "chemistry_id": "p15-test",
                    "method_id": "synthetic-method",
                    "adapter_id": "synthetic-adapter",
                    "status": "completed",
                }
            )
            values = {
                "task_id": task_id,
                "pose_id": pose_id,
                "receptor_id": receptor,
                "target": "ck2_alpha",
                "seed": str(seed),
                "model_path": model.relative_to(run_dir).as_posix(),
                "model_sha256": sha256_file(model),
                "model_index": "1",
                "contact_residues": "A:10;A:11",
                "local_x_A": str(offset + seed / 100.0),
                "local_y_A": "0.0",
                "local_z_A": "0.0",
                "coordinate_frame_id": "ck2_alpha_reference_v1",
                "ranking_score": "-10.0",
                "score_name": "synthetic_score",
                "qc_status": "passed",
            }
            evidence_rows.append("\t".join(values[field] for field in POSE_FIELDS))
    run_manifest = run_dir / "run_manifest.json"
    settings = SiteAnalysisSettings(
        contact_jaccard_distance=0.5,
        position_distance_A=2.0,
        min_seed_support=2,
        min_receptor_support=2,
        min_conformation_specific_seed_support=2,
        ensemble_candidate_budget=8,
        conformation_specific_candidate_budget=2,
    )
    atomic_write_json(
        run_manifest,
        {
            "schema": "vela.discovery-run-manifest/8",
            "stage": "discovery",
            "target_id": "ck2_alpha",
            "status": "planned",
            "evidence_category": "main_discovery",
            "known_site_information_used": False,
            "software": vela_software_identity(),
            "analysis_contract": candidate_analysis_contract(settings),
            "tasks": [dict(task, status="planned") for task in task_rows],
        },
    )
    atomic_write_json(
        run_dir / "sampling_manifest.json",
        {
            "schema": "vela.discovery-sampling-manifest/7",
            "stage": "discovery",
            "target_id": "ck2_alpha",
            "status": "sampling_completed",
            "run_manifest_sha256": sha256_file(run_manifest),
            "software": {
                "method_version": "1.0",
                "adapter_version": "1.0",
            },
            "tasks": [
                dict(
                    task,
                    execution_status="completed",
                    selection_status="completed",
                )
                for task in task_rows
            ],
        },
    )
    atomic_write_text(
        run_dir / "pose_evidence.tsv",
        "\t".join(POSE_FIELDS) + "\n" + "\n".join(evidence_rows) + "\n",
    )

    with pytest.raises(DiscoveryError, match="differs from the frozen run"):
        analyze_discovery_run(
            run_dir=run_dir,
            settings=replace(settings, position_distance_A=3.0),
        )

    analyze_discovery_run(run_dir=run_dir, settings=settings)

    receptor_report = run_dir / "site_analysis" / "receptor_sites.tsv"
    candidate_report = run_dir / "site_analysis" / "candidate_sites.tsv"
    assert receptor_report.is_file()
    assert candidate_report.is_file()
    candidate_text = candidate_report.read_text(encoding="utf-8")
    assert "3Q04_A;3QA0_A" in candidate_text
    assert "ensemble_consensus" in candidate_text
    assert "pose_ids" in receptor_report.read_text(encoding="utf-8").splitlines()[0]

    tasks = build_handoff_tasks(
        config=config,
        discovery_run_dir=run_dir,
        candidate_ids=("ALPHA_C001",),
    )

    assert len(tasks) == 8
    assert {task.pose.seed for task in tasks} == {11, 22, 33, 44}
    assert {task.pose.receptor_id for task in tasks} == {"3Q04_A", "3QA0_A"}
    assert all(task.candidate_id == "ALPHA_C001" for task in tasks)
    assert candidate_evidence_records(
        discovery_run_dir=run_dir,
        candidate_ids=("ALPHA_C001",),
    ) == [
        {
            "candidate_id": "ALPHA_C001",
            "evidence_tier": "ensemble_consensus",
            "rank_within_tier": 1,
        }
    ]

    one_pose_config = replace(
        config,
        validation=replace(
            config.validation,
            handoff=replace(config.validation.handoff, poses_per_receptor_site=1),
        ),
    )
    one_pose_tasks = build_handoff_tasks(
        config=one_pose_config,
        discovery_run_dir=run_dir,
        candidate_ids=("ALPHA_C001",),
    )
    assert len(one_pose_tasks) == 2

    with pytest.raises(ValidationError, match="unknown candidate IDs"):
        build_handoff_tasks(
            config=config,
            discovery_run_dir=run_dir,
            candidate_ids=("OTHER_C001",),
        )
    with pytest.raises(ValidationError, match="explicit candidate ID"):
        build_handoff_tasks(
            config=config,
            discovery_run_dir=run_dir,
            candidate_ids=(),
        )
