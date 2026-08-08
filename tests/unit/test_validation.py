import json
from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import gemmi
import pytest

from vela.config import load_config
from vela.core.provenance import atomic_write_json, sha256_file
from vela.core.typed_data import object_mapping
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.analysis.reports import (
    ReportedCandidateSite,
    ReportedReceptorSite,
    SiteAnalysisReport,
)
from vela.preparation.chemistry import DisulfideBond, HistidineState
from vela.validation.bound_states.comparison import compare_sites
from vela.validation.bound_states.controls import ControlInput, ControlTask
from vela.validation.bound_states.recovery import analyze_control_recovery
from vela.validation.models import (
    BoundStateDefinition,
    LocalRecoveryControl,
    RosettaSettings,
    ValidationError,
)
from vela.validation.readiness import assess_validation_readiness
from vela.validation.refinement.analysis import (
    is_funnel_confirmation_hit,
    is_funnel_deep_hit,
)
from vela.validation.refinement.geometry import (
    RefinedCluster,
    RefinedDecoy,
    ResolvedAnalysisSettings,
    assess_complex_geometry,
    cluster_refined_decoys,
    read_complex_geometry,
)
from vela.validation.refinement.handoff_plan import (
    exploration_promotion_contract,
    funnel_screening_contract,
    select_ranked_pose_cluster_representatives,
    select_source_seed_representatives,
)
from vela.validation.refinement.planning import (
    build_refinement_tasks,
    select_exploration_candidates,
    select_funnel_screening_candidates,
)
from vela.validation.refinement.qualification_analysis import (
    DiagnosticDecoy,
    cluster_diagnostic_decoys,
)
from vela.validation.refinement.qualification_diagnostic import (
    DiagnosticStart,
    build_diagnostic_tasks,
    validated_handoff_task_counts,
)
from vela.validation.refinement.receptor_flexibility import (
    select_local_receptor_backbone,
    write_local_receptor_movemap,
)
from vela.validation.refinement.reconstruction import (
    assess_topology_reconstruction,
    build_cg2all_command,
    write_cg2all_input,
    write_chemistry_protocol,
    write_chemistry_refine_protocol,
    write_reference_receptor_complex,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    build_chemistry_refine_command,
    build_prepack_command,
    build_refine_command,
    rosetta_crash_log_dir,
    run_rosetta_command,
)
from vela.validation.scores import read_rosetta_scorefile

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROJECT_CONFIG = PROJECT_ROOT / "configs"


def _handoff_pose(*, pose_id: str, x: float, seed: int = 101) -> PoseEvidence:
    """构造只用于姿态代表选择测试的最小有效证据。"""
    return PoseEvidence(
        task_id="task",
        pose_id=pose_id,
        receptor_id="receptor",
        target="target",
        seed=seed,
        model_path=Path("models.pdb"),
        model_sha256="0" * 64,
        model_index=1,
        contact_residues=frozenset({"A:1"}),
        local_position=(x, 0.0, 0.0),
        coordinate_frame_id="frame",
        ranking_score=x,
        score_name="score",
        qc_status="passed",
    )


def _handoff_pose_distance(first: PoseEvidence, second: PoseEvidence) -> float:
    return abs(first.local_position[0] - second.local_position[0])


def test_qualification_diagnostic_accepts_an_explicit_single_start() -> None:
    start = DiagnosticStart(
        start_id="selected_start",
        receptor_site_id="site",
        receptor_id="receptor",
        target="target",
        source_seed=101,
        input_path=Path("input.pdb"),
        input_sha256="0" * 64,
        receptor_residue_count=10,
        fixed_histidine_pose_indices=(),
        direct_contact_receptor_pose_indices=(),
        flexible_receptor_pose_indices=(),
    )

    tasks = build_diagnostic_tasks((start,), (201, 203, 205, 207))

    assert len(tasks) == 4
    assert {task.start.start_id for task in tasks} == {"selected_start"}
    assert {task.refinement_seed for task in tasks} == {201, 203, 205, 207}


def test_qualification_diagnostic_requires_starts_and_seeds() -> None:
    with pytest.raises(ValidationError, match="requires starts and seeds"):
        build_diagnostic_tasks((), (201,))


def _write_local_receptor_selection_structure(path: Path) -> None:
    structure = gemmi.Structure()
    model = gemmi.Model(1)
    receptor = gemmi.Chain("A")
    for index, x in enumerate((0.0, 4.0, 8.0, 12.0, 16.0), 1):
        residue = gemmi.Residue()
        residue.name = "ALA"
        residue.seqid = gemmi.SeqId(index, " ")
        atom = gemmi.Atom()
        atom.name = "CA"
        atom.element = gemmi.Element("C")
        atom.pos = gemmi.Position(x, 0.0, 0.0)
        residue.add_atom(atom)
        receptor.add_residue(residue)
    peptide = gemmi.Chain("P")
    residue = gemmi.Residue()
    residue.name = "GLY"
    residue.seqid = gemmi.SeqId(1, " ")
    atom = gemmi.Atom()
    atom.name = "CA"
    atom.element = gemmi.Element("C")
    atom.pos = gemmi.Position(8.5, 0.0, 0.0)
    residue.add_atom(atom)
    peptide.add_residue(residue)
    model.add_chain(receptor)
    model.add_chain(peptide)
    structure.add_model(model)
    path.write_text(structure.make_minimal_pdb(), encoding="utf-8")


def test_local_receptor_backbone_selection_and_movemap_are_native_free(
    tmp_path: Path,
) -> None:
    complex_path = tmp_path / "start.pdb"
    _write_local_receptor_selection_structure(complex_path)

    selection = select_local_receptor_backbone(
        path=complex_path, contact_A=2.0, sequence_padding=1
    )
    movemap = tmp_path / "local.movemap"
    write_local_receptor_movemap(
        destination=movemap,
        receptor_residue_count=5,
        peptide_residue_count=1,
        flexible_receptor_pose_indices=selection.flexible_pose_indices,
    )

    assert selection.direct_contact_pose_indices == (3,)
    assert selection.flexible_pose_indices == (2, 3, 4)
    assert movemap.read_text(encoding="utf-8") == (
        "RESIDUE 1 5 CHI\nRESIDUE 2 4 BBCHI\nRESIDUE 6 6 BBCHI\nJUMP * NO\nJUMP 1 YES\n"
    )


def test_handoff_pose_selection_balances_cluster_centers_and_edges() -> None:
    first_cluster = tuple(
        _handoff_pose(pose_id=f"first_{x}", x=x) for x in (0.0, 1.0, 3.0)
    )
    second_cluster = tuple(
        _handoff_pose(pose_id=f"second_{x}", x=x) for x in (10.0, 11.0, 13.0)
    )

    selected = select_ranked_pose_cluster_representatives(
        ranked_clusters=(first_cluster, second_cluster),
        count=4,
        distance=_handoff_pose_distance,
    )

    assert [pose.pose_id for pose in selected] == [
        "first_1.0",
        "second_11.0",
        "first_3.0",
        "second_13.0",
    ]
    assert {pose.seed for pose in selected} == {101}


def test_handoff_pose_selection_rejects_insufficient_passed_poses() -> None:
    with pytest.raises(ValidationError, match="only 1 passed poses"):
        select_ranked_pose_cluster_representatives(
            ranked_clusters=((_handoff_pose(pose_id="only", x=0.0),),),
            count=2,
            distance=_handoff_pose_distance,
        )


def test_confirmation_pose_selection_requires_distinct_source_seeds() -> None:
    cluster = tuple(
        _handoff_pose(pose_id=f"seed_{seed}", x=float(index), seed=seed)
        for index, seed in enumerate((101, 102, 103, 104, 105))
    )

    selected = select_source_seed_representatives(
        ranked_clusters=(cluster,),
        count=4,
        distance=_handoff_pose_distance,
    )

    assert len({pose.seed for pose in selected}) == 4
    assert [pose.seed for pose in selected] == [103, 102, 104, 101]


def test_confirmation_pose_selection_rejects_insufficient_source_seeds() -> None:
    cluster = tuple(
        _handoff_pose(pose_id=f"pose_{index}", x=float(index), seed=101)
        for index in range(4)
    )

    with pytest.raises(ValidationError, match="only 1 distinct source seeds"):
        select_source_seed_representatives(
            ranked_clusters=(cluster,),
            count=4,
            distance=_handoff_pose_distance,
        )


def _rosetta_settings(*, lowres_preoptimize: bool = False) -> RosettaSettings:
    return RosettaSettings(
        executable=Path("/apps/rosetta/FlexPepDocking"),
        scripts_executable=Path("/apps/rosetta/rosetta_scripts"),
        database=Path("/apps/rosetta/database"),
        version_file=Path("/apps/rosetta/source/src/utility/version.hh"),
        expected_version="2025.37+release.test",
        parallel_tasks=4,
        decoys_per_seed=100,
        score_function="ref2015",
        lowres_preoptimize=lowres_preoptimize,
    )


def _local_control(
    *,
    ranking_score: str = "custom_rank",
    recovery_rmsd_score: str = "custom_rmsd",
) -> LocalRecoveryControl:
    return LocalRecoveryControl(
        control_id="replaceable_control",
        bound_state_id="replaceable_bound_state",
        prepack_seed=1901,
        seed_batches=((1901, 1903), (1905, 1907)),
        random_translation_A=1.5,
        random_rotation_degrees=12.0,
        ranking_score=ranking_score,
        recovery_rmsd_score=recovery_rmsd_score,
        top_clusters=2,
        max_recovery_rmsd_A=2.0,
        max_cluster_backbone_rmsd_A=2.0,
        min_cluster_seed_support=1,
        max_batch_pose_rmsd_A=2.0,
    )


@pytest.mark.parametrize(
    ("seed_batches", "message"),
    [
        (((1901, 1903),), "at least two"),
        (((1901, 1903), (1903, 1905)), "unique non-negative"),
    ],
)
def test_local_control_rejects_invalid_seed_batches(
    seed_batches: tuple[tuple[int, ...], ...], message: str
) -> None:
    with pytest.raises(ValidationError, match=message):
        replace(_local_control(), seed_batches=seed_batches)


def test_project_config_registers_distinct_bound_state_evidence() -> None:
    settings = load_config(PROJECT_CONFIG).validation

    assert [state.state_id for state in settings.bound_states] == [
        "4IB5_A_D",
        "9FBM_A_F",
        "9FBI_A_E",
    ]
    assert settings.bound_states[0].local_control_kind == "standard_cyclic_peptide"
    assert settings.bound_states[0].ligand_sequence == "GCRLYGFKIHGCG"
    assert settings.bound_states[0].ligand_n_terminus == "Acetyl"
    assert settings.bound_states[0].ligand_c_terminus == "CONH2"
    assert settings.bound_states[0].disulfide_bonds == (DisulfideBond(2, 12),)
    assert settings.bound_states[0].histidines == (HistidineState(10, "HIE"),)
    assert all(
        state.local_control_kind == "site_reference_only"
        for state in settings.bound_states[1:]
    )
    assert len(settings.local_controls) == 1
    assert settings.local_controls[0].bound_state_id == "4IB5_A_D"
    assert settings.local_controls[0].seed_batches == (
        (2701, 2703, 2705, 2707),
        (2711, 2713, 2715, 2717),
    )
    assert settings.guided_templates[0].ligand_positions == tuple(range(2, 13))
    assert [item.reference_id for item in settings.environment_references] == [
        "4MD7_assembly_1",
        "1JWH_assembly_1",
    ]
    assert settings.config_complete
    assert settings.qualification_status == "qualified"
    assert settings.analysis.complete


@pytest.mark.parametrize(
    ("sequence", "bond", "message"),
    [
        ("CAC", DisulfideBond(3, 1), "outside the ligand sequence"),
        ("CAA", DisulfideBond(1, 3), "position 3 is not cysteine"),
        ("CAC", DisulfideBond(1, 4), "outside the ligand sequence"),
    ],
)
def test_standard_control_rejects_invalid_disulfide_definition(
    sequence: str, bond: DisulfideBond, message: str
) -> None:
    with pytest.raises(ValidationError, match=message):
        BoundStateDefinition(
            state_id="control",
            receptor_id="4IB5_A",
            ligand_id="test",
            ligand_author_chain_id="D",
            local_control_kind="standard_cyclic_peptide",
            ligand_sequence=sequence,
            ligand_n_terminus="NH3+",
            ligand_c_terminus="CONH2",
            disulfide_bonds=(bond,),
            histidines=(),
            selection_reason="test",
        )


def test_prepack_command_freezes_seed_chemistry_and_chain_roles() -> None:
    command = build_prepack_command(
        settings=_rosetta_settings(),
        input_path=Path("input.pdb"),
        disulfide_path=Path("fix_disulfide.txt"),
        output_dir=Path("prepack"),
        seed=1701,
        fixed_histidine_pose_indices=(338,),
    )

    assert command[0] == "/apps/rosetta/FlexPepDocking"
    assert command[command.index("-flexPepDocking:receptor_chain") + 1] == "A"
    assert command[command.index("-flexPepDocking:peptide_chain") + 1] == "P"
    assert command[command.index("-in:fix_disulf") + 1] == "fix_disulfide.txt"
    assert command[command.index("-jran") + 1] == "1701"
    assert command[command.index("-packing:fix_his_tautomer") + 1] == "338"
    assert "-constant_seed" in command
    assert "-flexPepDocking:flexpep_prepack" in command
    assert "-flexPepDocking:pep_refine" not in command


def test_refine_command_is_single_process_and_only_enables_declared_lowres() -> None:
    high_resolution = build_refine_command(
        settings=_rosetta_settings(),
        input_path=Path("prepacked.pdb"),
        disulfide_path=Path("fix_disulfide.txt"),
        output_dir=Path("refine"),
        seed=1703,
        native_path=Path("native.pdb"),
        random_translation_A=1.25,
        random_rotation_degrees=9.5,
        fixed_histidine_pose_indices=(338,),
    )
    with_low_resolution = build_refine_command(
        settings=_rosetta_settings(lowres_preoptimize=True),
        input_path=Path("prepacked.pdb"),
        disulfide_path=Path("fix_disulfide.txt"),
        output_dir=Path("refine"),
        seed=1703,
    )

    assert high_resolution[0] == "/apps/rosetta/FlexPepDocking"
    assert "-np" not in high_resolution
    assert high_resolution[high_resolution.index("-nstruct") + 1] == "100"
    assert high_resolution[high_resolution.index("-jran") + 1] == "1703"
    assert high_resolution[high_resolution.index("-native") + 1] == "native.pdb"
    assert (
        high_resolution[high_resolution.index("-flexPepDocking:random_trans_start") + 1]
        == "1.25"
    )
    assert (
        high_resolution[high_resolution.index("-flexPepDocking:random_rot_start") + 1]
        == "9.5"
    )
    assert high_resolution[high_resolution.index("-packing:fix_his_tautomer") + 1] == (
        "338"
    )
    assert "-flexPepDocking:pep_refine" in high_resolution
    assert "-flexPepDocking:lowres_preoptimize" not in high_resolution
    assert "-flexPepDocking:lowres_preoptimize" in with_low_resolution


@pytest.mark.parametrize("builder", [build_prepack_command, build_refine_command])
def test_rosetta_commands_reject_negative_seed(
    builder: Callable[..., tuple[str, ...]],
) -> None:
    keyword_arguments = {
        "settings": _rosetta_settings(),
        "input_path": Path("input.pdb"),
        "disulfide_path": Path("fix_disulfide.txt"),
        "output_dir": Path("output"),
        "seed": -1,
    }
    with pytest.raises(ValidationError, match="seed must not be negative"):
        builder(**keyword_arguments)


def test_scorefile_uses_configured_columns(tmp_path: Path) -> None:
    scorefile = tmp_path / "refine.sc"
    scorefile.write_text(
        "SEQUENCE: CAC\n"
        "SCORE: total_score custom_rank custom_rmsd description\n"
        "SCORE: -10.0 -3.0 4.5 decoy_1\n"
        "SCORE: -9.0 -2.0 1.8 decoy_2\n"
        "SCORE: -8.0 -1.0 0.7 decoy_3\n",
        encoding="utf-8",
    )

    rows = read_rosetta_scorefile(scorefile)
    assert [row.score(_local_control().ranking_score) for row in rows] == [
        -3.0,
        -2.0,
        -1.0,
    ]


def test_score_row_rejects_missing_configured_column(tmp_path: Path) -> None:
    scorefile = tmp_path / "refine.sc"
    scorefile.write_text(
        "SCORE: total_score other_score description\n"
        "SCORE: -10.0 1.0 decoy_1\n"
        "SCORE: -9.0 2.0 decoy_2\n",
        encoding="utf-8",
    )

    rows = read_rosetta_scorefile(scorefile)
    with pytest.raises(ValidationError, match="lacks required column: custom_rank"):
        rows[0].score(_local_control().ranking_score)


def test_scorefile_keeps_decoy_with_undefined_optional_metric(tmp_path: Path) -> None:
    scorefile = tmp_path / "refine.sc"
    scorefile.write_text(
        "SCORE: total_score custom_rank optional_metric description\n"
        "SCORE: -10.0 -3.0 -nan decoy_1\n",
        encoding="utf-8",
    )

    row = read_rosetta_scorefile(scorefile)[0]
    assert row.score("custom_rank") == -3.0
    with pytest.raises(ValidationError, match="lacks required column"):
        row.score("optional_metric")


def _coarse_chain(
    *, chain_id: str, residue_names: tuple[str, ...], y_offset: float
) -> gemmi.Chain:
    chain = gemmi.Chain(chain_id)
    for index, residue_name in enumerate(residue_names, 1):
        residue = gemmi.Residue()
        residue.name = residue_name
        residue.seqid = gemmi.SeqId(index, " ")
        residue.entity_type = gemmi.EntityType.Polymer
        for atom_name, atom_y in (("CA", y_offset), ("SC", y_offset + 1.0)):
            if residue_name == "GLY" and atom_name == "SC":
                continue
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element("C")
            atom.pos = gemmi.Position(index * 3.8, atom_y, 0.0)
            atom.occ = 1.0
            residue.add_atom(atom)
        chain.add_residue(residue)
    return chain


def test_topology_reconstruction_requires_composite_all_atom_geometry(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    residue_names = {
        "C": "CYS",
        "G": "GLY",
        "H": "HIS",
        "L": "LEU",
        "M": "MET",
        "P": "PRO",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
    }
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model(1))
    receptor = _all_atom_chain(
        chain_id="A", residue_names=("ALA", "ALA", "ALA"), y_offset=0.0
    )
    for index, residue in enumerate(receptor):
        for atom in residue:
            atom.pos.y += float(index % 2)
            atom.pos.z += float(index // 2)
    structure[0].add_chain(receptor)
    peptide_chain = _all_atom_chain(
        chain_id="P",
        residue_names=tuple(residue_names[code] for code in config.chemistry.sequence),
        y_offset=8.0,
    )
    for index, residue in enumerate(peptide_chain):
        for atom in residue:
            if atom.name != "SG":
                atom.pos.y += float(index % 3) * 0.4
                atom.pos.z += float(index % 2) * 0.5
    structure[0].add_chain(peptide_chain)
    source = tmp_path / "source.pdb"
    source.write_text(structure.make_pdb_string(), encoding="utf-8")

    passed = assess_topology_reconstruction(
        input_path=source,
        output_path=source,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
        min_disulfide_sg_A=1.8,
        max_disulfide_sg_A=2.3,
        min_interchain_heavy_atom_distance_A=1.2,
        min_nonlocal_peptide_heavy_atom_distance_A=1.2,
        contact_ca_threshold_A=10.0,
        max_peptide_internal_ca_rmsd_A=4.0,
        max_ligand_centroid_displacement_A=4.0,
        min_receptor_contact_retention_fraction=0.5,
    )
    assert passed.successful
    assert passed.disulfide_sg_distances_A == pytest.approx((2.05,))

    peptide = structure[0][1]
    second_sg = next(atom for atom in peptide[-1] if atom.name == "SG")
    second_sg.pos = gemmi.Position(10.0, 8.0, 1.0)
    failed_path = tmp_path / "failed.pdb"
    failed_path.write_text(structure.make_pdb_string(), encoding="utf-8")
    failed = assess_topology_reconstruction(
        input_path=failed_path,
        output_path=failed_path,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
        min_disulfide_sg_A=1.8,
        max_disulfide_sg_A=2.3,
        min_interchain_heavy_atom_distance_A=1.2,
        min_nonlocal_peptide_heavy_atom_distance_A=1.2,
        contact_ca_threshold_A=10.0,
        max_peptide_internal_ca_rmsd_A=4.0,
        max_ligand_centroid_displacement_A=4.0,
        min_receptor_contact_retention_fraction=0.5,
    )
    assert not failed.successful
    assert "disulfide_sg_geometry_outside_limits" in failed.failures

    shifted = gemmi.read_structure(str(source))
    for residue in shifted[0][1]:
        for atom in residue:
            atom.pos.x += 5.0
    shifted_path = tmp_path / "shifted.pdb"
    shifted_path.write_text(shifted.make_pdb_string(), encoding="utf-8")
    displaced = assess_topology_reconstruction(
        input_path=source,
        output_path=shifted_path,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
        min_disulfide_sg_A=1.8,
        max_disulfide_sg_A=2.3,
        min_interchain_heavy_atom_distance_A=1.2,
        min_nonlocal_peptide_heavy_atom_distance_A=1.2,
        contact_ca_threshold_A=10.0,
        max_peptide_internal_ca_rmsd_A=4.0,
        max_ligand_centroid_displacement_A=4.0,
        min_receptor_contact_retention_fraction=0.5,
    )
    assert not displaced.successful
    assert "ligand_site_centroid_not_preserved" in displaced.failures


def test_cg2all_handoff_uses_configured_tool_and_normalizes_pose(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    source = tmp_path / "source.pdb"
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model(1))
    structure[0].add_chain(
        _coarse_chain(chain_id="R", residue_names=("ALA", "HIS"), y_offset=0.0)
    )
    residue_names = {
        "C": "CYS",
        "G": "GLY",
        "H": "HIS",
        "L": "LEU",
        "M": "MET",
        "P": "PRO",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
    }
    peptide_names = tuple(residue_names[code] for code in config.chemistry.sequence)
    structure[0].add_chain(
        _coarse_chain(chain_id="L", residue_names=peptide_names, y_offset=8.0)
    )
    structure.add_model(gemmi.Model(2))
    structure[1].add_chain(
        _coarse_chain(chain_id="R", residue_names=("ALA", "HIS"), y_offset=0.0)
    )
    structure[1].add_chain(
        _coarse_chain(chain_id="L", residue_names=peptide_names, y_offset=18.0)
    )
    source.write_text(structure.make_pdb_string(), encoding="utf-8")
    destination = tmp_path / "cg2all_input.pdb"

    asset = write_cg2all_input(
        source_path=source,
        model_index=2,
        destination=destination,
        chemistry=config.chemistry,
        settings=config.validation.cg2all,
    )
    command = build_cg2all_command(
        settings=config.validation.cg2all,
        input_path=destination,
        output_path=tmp_path / "all_atom.pdb",
    )

    output = gemmi.read_structure(str(destination))
    assert [chain.name for chain in output[0]] == ["A", "P"]
    assert output[0][0][1].name == "HSE"
    assert output[0][1][6].name == "HSE"
    assert output[0][1][0][0].pos.y == pytest.approx(18.0)
    assert asset.receptor_residue_count == 2
    assert asset.fixed_histidine_pose_indices == (9,)
    assert destination.read_text(encoding="utf-8").startswith("SSBOND")
    assert command[command.index("--ckpt") + 1] == str(
        config.validation.cg2all.checkpoint
    )
    assert command[command.index("--proc") + 1] == "2"
    assert "--standard-name" in command


def test_chemistry_restoration_uses_declared_score_and_seed(tmp_path: Path) -> None:
    config = load_config(PROJECT_CONFIG)
    protocol = tmp_path / "restore.xml"
    write_chemistry_protocol(
        destination=protocol,
        receptor_residue_count=331,
        chemistry=config.chemistry,
        score_function="custom_score",
    )
    command = build_chemistry_command(
        settings=config.validation.rosetta,
        input_path=Path("input.pdb"),
        protocol_path=protocol,
        disulfide_path=Path("fix_disulfide.txt"),
        output_dir=tmp_path,
        seed=3101,
        fixed_histidine_pose_indices=(338,),
    )

    protocol_text = protocol.read_text(encoding="utf-8")
    assert 'weights="custom_score"' in protocol_text
    assert 'resnums="342"' in protocol_text
    assert 'add_type="CTERM_AMIDATION"' in protocol_text
    assert '<ForceDisulfides name="form_disulfides"' in protocol_text
    assert 'disulfides="332:342"' in protocol_text
    assert 'repack="true"' in protocol_text
    assert 'new_res="HIS"' in protocol_text
    assert 'target="338"' in protocol_text
    assert command[0] == str(config.validation.rosetta.scripts_executable)
    assert command[command.index("-jran") + 1] == "3101"
    assert command[command.index("-packing:fix_his_tautomer") + 1] == "338"


def test_chemistry_restoration_supports_declared_n_terminal_acetylation(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    protocol = tmp_path / "restore_acetyl.xml"

    write_chemistry_protocol(
        destination=protocol,
        receptor_residue_count=331,
        chemistry=replace(config.chemistry, n_terminus="Acetyl"),
        score_function="ref2015",
    )

    text = protocol.read_text(encoding="utf-8")
    assert 'name="peptide_n_terminus" resnums="332"' in text
    assert 'add_type="N_ACETYLATION"' in text
    assert 'remove_type="LOWER_TERMINUS_VARIANT"' in text
    assert 'name="remove_n_terminal_charge"' in text
    assert 'name="add_n_terminal_acetyl"' in text
    assert text.index('<Add mover="add_c_terminal_amide" />') < text.index(
        '<Add mover="form_disulfides" />'
    )


def test_terminal_chemistry_refine_stays_in_one_rosetta_pose(tmp_path: Path) -> None:
    config = load_config(PROJECT_CONFIG)
    protocol = tmp_path / "restore_and_refine.xml"
    constraint = tmp_path / "site.cst"
    constraint.write_text("constraint\n", encoding="utf-8")
    write_chemistry_refine_protocol(
        destination=protocol,
        receptor_residue_count=331,
        chemistry=replace(config.chemistry, n_terminus="Acetyl"),
        score_function="ref2015",
    )

    command = build_chemistry_refine_command(
        settings=config.validation.rosetta,
        input_path=Path("input.pdb"),
        protocol_path=protocol,
        disulfide_path=Path("fix_disulfide.txt"),
        output_dir=tmp_path,
        seed=3101,
        fixed_histidine_pose_indices=(338,),
        site_constraint_path=constraint,
        site_constraint_weight=1.0,
    )

    text = protocol.read_text(encoding="utf-8")
    assert '<FlexPepDock name="run_flexpepdock"' in text
    assert text.index('<Add mover="form_disulfides" />') < text.index(
        '<Add mover="run_flexpepdock" />'
    )
    assert command[command.index("-flexPepDocking:receptor_chain") + 1] == "A"
    assert command[command.index("-flexPepDocking:peptide_chain") + 1] == "P"
    assert command[command.index("-constraints:cst_fa_file") + 1] == str(constraint)


def _all_atom_chain(
    *, chain_id: str, residue_names: tuple[str, ...], y_offset: float
) -> gemmi.Chain:
    chain = gemmi.Chain(chain_id)
    for index, residue_name in enumerate(residue_names, 1):
        residue = gemmi.Residue()
        residue.name = residue_name
        residue.seqid = gemmi.SeqId(index, " ")
        residue.entity_type = gemmi.EntityType.Polymer
        for atom_name, x_offset in (("N", -1.2), ("CA", 0.0), ("C", 1.2), ("O", 1.8)):
            atom = gemmi.Atom()
            atom.name = atom_name
            atom.element = gemmi.Element(atom_name[0])
            atom.pos = gemmi.Position(index * 3.8 + x_offset, y_offset, 0.0)
            residue.add_atom(atom)
        if chain_id == "P" and index == 1:
            for name in ("1H", "2H", "3H"):
                atom = gemmi.Atom()
                atom.name = name
                atom.element = gemmi.Element("H")
                atom.pos = gemmi.Position(index * 3.8 - 1.5, y_offset, 0.0)
                residue.add_atom(atom)
        if chain_id == "P" and index == len(residue_names):
            atom = gemmi.Atom()
            atom.name = "NT"
            atom.element = gemmi.Element("N")
            atom.pos = gemmi.Position(index * 3.8 + 2.0, y_offset, 0.0)
            residue.add_atom(atom)
        if residue_name == "HIS":
            atom = gemmi.Atom()
            atom.name = "HE2"
            atom.element = gemmi.Element("H")
            atom.pos = gemmi.Position(index * 3.8, y_offset + 1.0, 0.0)
            residue.add_atom(atom)
        if chain_id == "P" and residue_name == "CYS":
            atom = gemmi.Atom()
            atom.name = "SG"
            atom.element = gemmi.Element("S")
            atom.pos = gemmi.Position(0.0 if index == 1 else 2.05, y_offset, 2.0)
            residue.add_atom(atom)
        chain.add_residue(residue)
    return chain


def test_reference_receptor_graft_discards_reconstructed_receptor(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    residue_names = {
        "C": "CYS",
        "G": "GLY",
        "H": "HIS",
        "L": "LEU",
        "M": "MET",
        "P": "PRO",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
    }
    reference = _all_atom_chain(
        chain_id="R", residue_names=("ALA", "ALA", "ALA"), y_offset=0.0
    )
    for index, residue in enumerate(reference):
        for atom in residue:
            atom.pos.y += float(index % 2)
            atom.pos.z += float(index // 2)
    peptide = _all_atom_chain(
        chain_id="L",
        residue_names=tuple(residue_names[code] for code in config.chemistry.sequence),
        y_offset=8.0,
    )

    def write(path: Path, receptor: gemmi.Chain) -> None:
        structure = gemmi.Structure()
        structure.add_model(gemmi.Model(1))
        structure[0].add_chain(receptor)
        structure[0].add_chain(peptide.clone())
        path.write_text(structure.make_pdb_string(), encoding="utf-8")

    coarse_path = tmp_path / "coarse.pdb"
    write(coarse_path, reference.clone())
    reconstructed = reference.clone()
    for residue in reconstructed:
        for atom in residue:
            atom.pos.y += 20.0
    reconstructed_path = tmp_path / "reconstructed.pdb"
    write(reconstructed_path, reconstructed)
    reference_path = tmp_path / "reference.pdb"
    reference_structure = gemmi.Structure()
    reference_structure.add_model(gemmi.Model(1))
    reference_structure[0].add_chain(reference.clone())
    reference_path.write_text(reference_structure.make_pdb_string(), encoding="utf-8")
    output_path = tmp_path / "grafted.pdb"

    write_reference_receptor_complex(
        coarse_pose_path=coarse_path,
        reconstructed_path=reconstructed_path,
        reference_receptor_path=reference_path,
        destination=output_path,
        chemistry=config.chemistry,
    )

    output = gemmi.read_structure(str(output_path))
    output_ca = next(atom for atom in output[0][0][0] if atom.name == "CA")
    reference_ca = next(atom for atom in reference[0] if atom.name == "CA")
    assert output[0][0].name == "A"
    assert output[0][1].name == "P"
    assert output_ca.pos.dist(reference_ca.pos) < 1e-6


def _flexpepdock_input(path: Path, *, peptide_y: float = 4.0) -> None:
    structure = gemmi.Structure()
    structure.add_model(gemmi.Model(1))
    receptor = _all_atom_chain(
        chain_id="A", residue_names=("ALA", "ALA", "ALA"), y_offset=0.0
    )
    for index, residue in enumerate(receptor):
        for atom in residue:
            atom.pos.y += float(index % 2)
            atom.pos.z += float(index // 2)
    structure[0].add_chain(receptor)
    names = {
        "C": "CYS",
        "G": "GLY",
        "H": "HIS",
        "L": "LEU",
        "M": "MET",
        "P": "PRO",
        "R": "ARG",
        "S": "SER",
        "T": "THR",
        "W": "TRP",
    }
    config = load_config(PROJECT_CONFIG)
    structure[0].add_chain(
        _all_atom_chain(
            chain_id="P",
            residue_names=tuple(names[code] for code in config.chemistry.sequence),
            y_offset=peptide_y,
        )
    )
    path.write_text(structure.make_pdb_string(), encoding="utf-8")


def test_refinement_geometry_records_invalid_disulfide_as_decoy_qc(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    start_path = tmp_path / "start.pdb"
    invalid_path = tmp_path / "invalid_disulfide.pdb"
    _flexpepdock_input(start_path)
    structure = gemmi.read_structure(str(start_path))
    peptide = structure[0][1]
    first_sg = next(atom for atom in peptide[0] if atom.name == "SG")
    second_sg = next(atom for atom in peptide[-1] if atom.name == "SG")
    second_sg.pos = gemmi.Position(
        first_sg.pos.x + 2.667, first_sg.pos.y, first_sg.pos.z
    )
    invalid_path.write_text(structure.make_pdb_string(), encoding="utf-8")
    start = read_complex_geometry(
        path=start_path,
        interface_contact_A=config.validation.interface_contact_A,
    )
    settings = ResolvedAnalysisSettings(
        min_interface_contact_pairs=1,
        min_interface_receptor_residues=1,
        max_receptor_ca_rmsd_A=1.0,
        min_start_contact_overlap=0.0,
        max_start_site_displacement_A=10.0,
        max_cluster_backbone_rmsd_A=2.0,
        min_heavy_atom_distance_A=0.0,
        min_refinement_seed_support=2,
        min_refinement_start_support=2,
        min_refinement_source_seed_support=2,
    )

    assessment = assess_complex_geometry(
        path=invalid_path,
        chemistry=config.chemistry,
        start=start,
        cluster_reference=start,
        config=config,
        settings=settings,
    )

    assert not assessment.passed
    assert assessment.qc_failures == ("chemistry_invalid",)
    assert assessment.chemistry_failure is not None
    assert "disulfide SG distance is outside config" in assessment.chemistry_failure


def test_refinement_neighbor_search_matches_exhaustive_interface_geometry(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    path = tmp_path / "complex.pdb"
    _flexpepdock_input(path)
    structure = gemmi.read_structure(str(path))
    receptor_atoms = tuple(
        (residue, atom)
        for residue in structure[0][0]
        for atom in residue
        if atom.element.name != "H"
    )
    peptide_atoms = tuple(
        atom
        for residue in structure[0][1]
        for atom in residue
        if atom.element.name != "H"
    )
    expected_pairs = 0
    expected_contacts: set[str] = set()
    expected_minimum = float("inf")
    for residue, receptor_atom in receptor_atoms:
        for peptide_atom in peptide_atoms:
            distance = receptor_atom.pos.dist(peptide_atom.pos)
            expected_minimum = min(expected_minimum, distance)
            if distance <= config.validation.interface_contact_A:
                expected_pairs += 1
                expected_contacts.add(
                    f"{residue.seqid.num}{residue.seqid.icode.strip()}"
                )

    geometry = read_complex_geometry(
        path=path,
        interface_contact_A=config.validation.interface_contact_A,
    )

    assert geometry.interface_contact_pairs == expected_pairs
    assert geometry.receptor_contacts == expected_contacts
    assert geometry.minimum_interface_distance_A == pytest.approx(expected_minimum)


def test_control_recovery_pools_seed_batches_and_requires_repeatable_clusters(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    control = _local_control()
    native_path = tmp_path / "native.pdb"
    _flexpepdock_input(native_path)
    control_input = ControlInput(
        definition=config.validation.bound_states[0],
        complex_path=native_path,
        disulfide_path=tmp_path / "fix_disulfide.txt",
        chemistry=config.chemistry,
        receptor_residue_count=3,
        fixed_histidine_pose_indices=(),
    )
    tasks: list[ControlTask] = []
    for batch_index, batch in enumerate(control.seed_batches, 1):
        batch_id = f"batch_{batch_index:02d}"
        for seed in batch:
            task_id = f"{control.control_id}__{batch_id}__seed_{seed}"
            task = ControlTask(task_id, control, control_input, batch_id, seed)
            tasks.append(task)
            task_dir = tmp_path / "tasks" / task_id
            task_dir.mkdir(parents=True)
            _flexpepdock_input(task_dir / "decoy_1.pdb")
            _flexpepdock_input(task_dir / "decoy_2.pdb", peptide_y=10.0)
            (task_dir / "refine.sc").write_text(
                "SCORE: total_score custom_rank custom_rmsd description\n"
                "SCORE: -10.0 -2.0 0.5 decoy_1\n"
                "SCORE: -9.0 -1.0 6.0 decoy_2\n",
                encoding="utf-8",
            )

    passed, report = analyze_control_recovery(
        control=control, tasks=tuple(tasks), run_dir=tmp_path
    )

    assert passed
    assert report["all_batches_sampled"] is True
    assert report["all_batches_ranked"] is True
    assert report["batch_poses_consistent"] is True
    batches = report["batches"]
    assert isinstance(batches, list)
    assert len(batches) == 2
    assert all(
        isinstance(batch, dict)
        and batch["decoy_count"] == 4
        and batch["selected_cluster_seed_support"] == 2
        for batch in batches
    )


def test_control_recovery_uses_cluster_medoid_for_batch_pose_consistency(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    control = _local_control()
    native_path = tmp_path / "native.pdb"
    _flexpepdock_input(native_path)
    control_input = ControlInput(
        definition=config.validation.bound_states[0],
        complex_path=native_path,
        disulfide_path=tmp_path / "fix_disulfide.txt",
        chemistry=config.chemistry,
        receptor_residue_count=3,
        fixed_histidine_pose_indices=(),
    )
    tasks: list[ControlTask] = []
    for batch_index, batch in enumerate(control.seed_batches, 1):
        batch_id = f"batch_{batch_index:02d}"
        extreme_y = 2.9 if batch_index == 1 else 5.1
        for seed_index, seed in enumerate(batch):
            task_id = f"{control.control_id}__{batch_id}__seed_{seed}"
            tasks.append(ControlTask(task_id, control, control_input, batch_id, seed))
            task_dir = tmp_path / "tasks" / task_id
            task_dir.mkdir(parents=True)
            _flexpepdock_input(task_dir / "center.pdb")
            rows = ["SCORE: -9.0 -1.0 0.5 center\n"]
            if seed_index == 0:
                _flexpepdock_input(task_dir / "energy_best.pdb", peptide_y=extreme_y)
                rows.insert(0, "SCORE: -10.0 -2.0 0.5 energy_best\n")
            (task_dir / "refine.sc").write_text(
                "SCORE: total_score custom_rank custom_rmsd description\n"
                + "".join(rows),
                encoding="utf-8",
            )

    passed, report = analyze_control_recovery(
        control=control, tasks=tuple(tasks), run_dir=tmp_path
    )

    assert passed
    assert report["batch_poses_consistent"] is True
    assert report["maximum_observed_batch_pose_rmsd_A"] == pytest.approx(0.0)
    batches = report["batches"]
    assert isinstance(batches, list)
    assert all(
        isinstance(batch, dict)
        and batch["selected_energy_representative_decoy_id"]
        != batch["selected_geometric_medoid_decoy_id"]
        for batch in batches
    )


def test_refinement_tasks_are_derived_from_handoff_and_configured_seeds(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    config = replace(
        config,
        validation=replace(
            config.validation,
            seeds=(4101, 4103, 4105),
            refinement=replace(
                config.validation.refinement, seed_batch_sizes=(1, 1, 1)
            ),
        ),
    )
    run_dir = tmp_path / "handoff"
    input_path = run_dir / "tasks" / "source_001" / "flexpepdock_input.pdb"
    input_path.parent.mkdir(parents=True)
    _flexpepdock_input(input_path)
    result_path = input_path.parent / "task_result.json"
    atomic_write_json(result_path, {"status": "completed"})
    plan_path = run_dir / "handoff_plan.json"
    atomic_write_json(plan_path, {"status": "planned"})
    atomic_write_json(
        run_dir / "handoff_manifest.json",
        {
            "schema": "vela.validation-handoff-manifest/8",
            "status": "completed",
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "main_discovery_handoff",
            "source_evidence_category": "main_discovery",
            "known_site_information_used": False,
            "production_qualified": True,
            "handoff_plan": {
                "path": plan_path.relative_to(run_dir).as_posix(),
                "sha256": sha256_file(plan_path),
            },
            "tasks": [
                {
                    "task_id": "source_001",
                    "candidate_id": "CANDIDATE_001",
                    "receptor_site_id": "SITE_001",
                    "pose_id": "POSE_001",
                    "receptor_id": "RECEPTOR_001",
                    "target": "TARGET_001",
                    "source_seed": 2101,
                    "execution_status": "completed",
                    "reconstruction_status": "passed",
                    "task_result": {
                        "path": result_path.relative_to(run_dir).as_posix(),
                        "sha256": sha256_file(result_path),
                    },
                    "flexpepdock_input": {
                        "path": input_path.relative_to(run_dir).as_posix(),
                        "sha256": sha256_file(input_path),
                    },
                }
            ],
        },
    )

    tasks = build_refinement_tasks(config=config, source_run_dir=run_dir)

    assert [task.task_id for task in tasks] == [
        "refine_00001",
        "refine_00002",
        "refine_00003",
    ]
    assert [task.seed for task in tasks] == [4101, 4103, 4105]
    assert all(task.start.candidate_id == "CANDIDATE_001" for task in tasks)


def test_exploration_refinement_uses_frozen_arm_order_and_same_arm_failover() -> None:
    config = load_config(PROJECT_CONFIG)
    contract = exploration_promotion_contract(config)
    plan: dict[str, object] = {
        "selection": {
            "candidate_arms": {
                "blind_discovery_arm": ["BLIND_A", "BLIND_B", "BLIND_C"],
                "functional_annotation_arm": ["FUNCTIONAL_A", "FUNCTIONAL_B"],
            },
            "requested_candidate_ids": [
                "BLIND_A",
                "BLIND_B",
                "BLIND_C",
                "FUNCTIONAL_A",
                "FUNCTIONAL_B",
            ],
            "promotion_contract": contract,
        }
    }
    minimum = config.validation.analysis.min_refinement_start_support
    assert minimum is not None
    rows: list[dict[str, object]] = []
    for candidate_id in (
        "BLIND_A",
        "BLIND_B",
        "BLIND_C",
        "FUNCTIONAL_A",
        "FUNCTIONAL_B",
    ):
        for site_id in ("RECEPTOR_1", "RECEPTOR_2"):
            passed_count = minimum
            if candidate_id == "BLIND_A" and site_id == "RECEPTOR_2":
                passed_count = minimum - 1
            for _ in range(passed_count):
                rows.append(
                    {
                        "candidate_id": candidate_id,
                        "receptor_site_id": site_id,
                        "execution_status": "completed",
                        "reconstruction_status": "passed",
                    }
                )

    selected, record = select_exploration_candidates(
        config=config,
        handoff_plan=plan,
        rows=tuple(rows),
    )

    assert selected == ("BLIND_B", "BLIND_C", "FUNCTIONAL_A")
    assert record["selected_by_arm"] == {
        "blind_discovery_arm": ["BLIND_B", "BLIND_C"],
        "functional_annotation_arm": ["FUNCTIONAL_A"],
    }


def test_funnel_screening_excludes_candidates_with_incomplete_site_handoff() -> None:
    config = load_config(PROJECT_CONFIG)
    requested = ("COMPLETE", "INCOMPLETE")
    rows: list[dict[str, object]] = []
    for candidate_id in requested:
        for site_id in ("RECEPTOR_1", "RECEPTOR_2"):
            count = 2
            if candidate_id == "INCOMPLETE" and site_id == "RECEPTOR_2":
                count = 1
            for _ in range(count):
                rows.append(
                    {
                        "candidate_id": candidate_id,
                        "receptor_site_id": site_id,
                        "execution_status": "completed",
                        "reconstruction_status": "passed",
                    }
                )

    selected, record = select_funnel_screening_candidates(
        config=config,
        requested=requested,
        rows=tuple(rows),
        promotion_contract={"contract": "frozen"},
        funnel_audit={"audit": "frozen"},
    )

    assert selected == ("COMPLETE",)
    assert record["excluded_candidate_ids"] == ["INCOMPLETE"]
    support = object_mapping(record["candidate_support"], name="candidate support")
    incomplete = object_mapping(support["INCOMPLETE"], name="incomplete support")
    assert incomplete["eligible"] is False


def test_historical_local_qualification_remains_valid_after_source_changes(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    report = config.validation.qualification_report
    assert report is not None
    document = json.loads(report.read_text(encoding="utf-8"))
    document["analysis_software"]["vela_source_sha256"] = "0" * 64
    historical_report = tmp_path / "qualification_report.json"
    atomic_write_json(historical_report, document)
    config = replace(
        config,
        validation=replace(
            config.validation,
            qualification_report=historical_report,
            qualification_report_sha256=sha256_file(historical_report),
        ),
    )

    readiness = assess_validation_readiness(config)

    assert "qualification_report_mismatch" not in {
        issue.code for issue in readiness.issues
    }


def test_refinement_rejects_technically_invalid_handoff(tmp_path: Path) -> None:
    config = load_config(PROJECT_CONFIG)
    run_dir = tmp_path / "handoff"
    result_path = run_dir / "tasks" / "source_001" / "task_result.json"
    atomic_write_json(result_path, {"status": "completed"})
    plan_path = run_dir / "handoff_plan.json"
    atomic_write_json(plan_path, {"status": "planned"})
    atomic_write_json(
        run_dir / "handoff_manifest.json",
        {
            "schema": "vela.validation-handoff-manifest/8",
            "status": "completed",
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "main_discovery_handoff",
            "source_evidence_category": "main_discovery",
            "known_site_information_used": False,
            "production_qualified": True,
            "handoff_plan": {
                "path": plan_path.relative_to(run_dir).as_posix(),
                "sha256": sha256_file(plan_path),
            },
            "tasks": [
                {
                    "task_id": "source_001",
                    "candidate_id": "CANDIDATE_001",
                    "receptor_site_id": "SITE_001",
                    "pose_id": "POSE_001",
                    "receptor_id": "RECEPTOR_001",
                    "target": "TARGET_001",
                    "source_seed": 2101,
                    "execution_status": "invalid",
                    "reconstruction_status": "not_assessed",
                    "task_result": {
                        "path": result_path.relative_to(run_dir).as_posix(),
                        "sha256": sha256_file(result_path),
                    },
                    "flexpepdock_input": None,
                }
            ],
        },
    )

    with pytest.raises(ValidationError, match="invalid handoff task"):
        build_refinement_tasks(config=config, source_run_dir=run_dir)


def test_qualification_diagnostic_allows_candidate_qc_attrition() -> None:
    manifest: dict[str, object] = {
        "passed_task_count": 1,
        "failed_task_count": 1,
        "invalid_task_count": 0,
        "tasks": [
            {
                "execution_status": "completed",
                "reconstruction_status": "passed",
            },
            {
                "execution_status": "completed",
                "reconstruction_status": "failed",
            },
        ],
    }

    assert validated_handoff_task_counts(manifest) == (1, 1, 0)


def test_qualification_diagnostic_rejects_inconsistent_handoff_counts() -> None:
    manifest: dict[str, object] = {
        "passed_task_count": 2,
        "failed_task_count": 0,
        "invalid_task_count": 0,
        "tasks": [
            {
                "execution_status": "completed",
                "reconstruction_status": "passed",
            }
        ],
    }

    with pytest.raises(ValidationError, match="counts are inconsistent"):
        validated_handoff_task_counts(manifest)


def test_refinement_tasks_preserve_guided_source_identity(tmp_path: Path) -> None:
    config = load_config(PROJECT_CONFIG)
    config = replace(
        config,
        validation=replace(
            config.validation,
            seeds=(4201, 4203, 4205),
            refinement=replace(
                config.validation.refinement, seed_batch_sizes=(1, 1, 1)
            ),
        ),
    )
    run_dir = tmp_path / "guided"
    input_path = run_dir / "tasks" / "guided_001" / "flexpepdock_input.pdb"
    input_path.parent.mkdir(parents=True)
    _flexpepdock_input(input_path)
    result_path = input_path.parent / "task_result.json"
    atomic_write_json(result_path, {"status": "completed"})
    plan_path = run_dir / "guided_plan.json"
    atomic_write_json(plan_path, {"status": "planned"})
    atomic_write_json(
        run_dir / "guided_manifest.json",
        {
            "schema": "vela.validation-guided-manifest/1",
            "status": "completed",
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "guided_site_compatibility_handoff",
            "known_site_information_used": True,
            "guided_plan": {
                "path": plan_path.relative_to(run_dir).as_posix(),
                "sha256": sha256_file(plan_path),
            },
            "tasks": [
                {
                    "task_id": "guided_001",
                    "candidate_id": "guided__replaceable_template",
                    "receptor_site_id": "experimental__replaceable_state",
                    "pose_id": "replaceable_template",
                    "receptor_id": "RECEPTOR_001",
                    "target": "TARGET_001",
                    "source_seed": None,
                    "task_result": {
                        "path": result_path.relative_to(run_dir).as_posix(),
                        "sha256": sha256_file(result_path),
                    },
                    "flexpepdock_input": {
                        "path": input_path.relative_to(run_dir).as_posix(),
                        "sha256": sha256_file(input_path),
                    },
                }
            ],
        },
    )

    tasks = build_refinement_tasks(config=config, source_run_dir=run_dir)

    assert len(tasks) == 3
    assert all(task.start.source_seed is None for task in tasks)
    assert all(
        task.start.candidate_id == "guided__replaceable_template" for task in tasks
    )


def _refined_decoy(
    *,
    decoy_id: str,
    seed: int,
    start_id: str,
    x_offset: float,
    source_seed: int = 1001,
) -> RefinedDecoy:
    backbone = tuple(gemmi.Position(index + x_offset, 0.0, 0.0) for index in range(6))
    return RefinedDecoy(
        decoy_id=decoy_id,
        task_id=f"task_{seed}",
        start_id=start_id,
        candidate_id="replaceable_candidate",
        receptor_id="replaceable_receptor",
        source_seed=source_seed,
        refinement_seed=seed,
        path=Path(f"{decoy_id}.pdb"),
        sha256="0" * 64,
        ranking_score=-10.0,
        interface_contact_pairs=10,
        interface_receptor_residues=3,
        minimum_interface_distance_A=2.0,
        receptor_ca_rmsd_A=0.1,
        start_contact_overlap=1.0,
        start_site_displacement_A=0.1,
        qc_status="passed",
        qc_failures=(),
        chemistry_failure=None,
        cluster_backbone=backbone,
    )


def test_refinement_cluster_support_is_config_driven() -> None:
    decoys = (
        _refined_decoy(
            decoy_id="decoy_a",
            seed=5101,
            start_id="start_a",
            x_offset=0.0,
            source_seed=1001,
        ),
        _refined_decoy(
            decoy_id="decoy_b",
            seed=5103,
            start_id="start_b",
            x_offset=0.1,
            source_seed=1002,
        ),
    )
    settings = ResolvedAnalysisSettings(
        min_interface_contact_pairs=1,
        min_interface_receptor_residues=1,
        max_receptor_ca_rmsd_A=1.0,
        min_start_contact_overlap=0.5,
        max_start_site_displacement_A=1.0,
        max_cluster_backbone_rmsd_A=0.5,
        min_heavy_atom_distance_A=1.0,
        min_refinement_seed_support=2,
        min_refinement_start_support=2,
        min_refinement_source_seed_support=2,
    )

    clusters = cluster_refined_decoys(
        decoys=decoys, settings=settings, require_source_seed_support=True
    )

    assert len(clusters) == 1
    assert clusters[0].supported
    assert clusters[0].refinement_seeds == (5101, 5103)
    assert clusters[0].start_ids == ("start_a", "start_b")
    assert clusters[0].source_seeds == (1001, 1002)
    assert clusters[0].task_cells == ("start_a:5101", "start_b:5103")


def test_refinement_funnel_contract_freezes_incremental_budget() -> None:
    config = load_config(PROJECT_CONFIG)

    contract = funnel_screening_contract(config)

    assert contract["stage3a_screening"] == {
        "source_starts_per_receptor_site": 2,
        "rosetta_seed_count": 1,
        "promotion_status": "cross_source_screening_hit",
        "promotion_budget": 6,
    }
    assert contract["stage3b_confirmation"] == {
        "additional_rosetta_seed_count": 1,
        "minimum_source_starts": 2,
        "minimum_rosetta_seeds": 2,
        "minimum_source_seed_task_cells": 3,
        "promotion_budget": 3,
    }
    assert contract["stage3c_deep_confirmation"] == {
        "additional_rosetta_seed_count": 2,
        "minimum_source_starts": 2,
        "minimum_rosetta_seeds_per_source": 2,
        "source_pool_per_receptor_site": 4,
        "final_hypothesis_budget": 2,
    }


def test_funnel_confirmation_requires_three_distinct_source_seed_cells() -> None:
    settings = ResolvedAnalysisSettings(
        min_interface_contact_pairs=1,
        min_interface_receptor_residues=1,
        max_receptor_ca_rmsd_A=1.0,
        min_start_contact_overlap=0.2,
        max_start_site_displacement_A=6.0,
        max_cluster_backbone_rmsd_A=2.0,
        min_heavy_atom_distance_A=1.2,
        min_refinement_seed_support=2,
        min_refinement_start_support=2,
        min_refinement_source_seed_support=2,
    )

    def cluster(task_cells: tuple[str, ...]) -> RefinedCluster:
        return RefinedCluster(
            cluster_id="cluster",
            candidate_id="candidate",
            receptor_id="receptor",
            decoy_ids=("decoy",),
            refinement_seeds=(120623, 120624),
            start_ids=("start_a", "start_b"),
            source_seeds=(120628, 120629),
            task_cells=task_cells,
            representative_decoy_id="decoy",
            supported=True,
        )

    assert is_funnel_confirmation_hit(
        cluster=cluster(("a:1", "a:2", "b:1")),
        settings=settings,
        minimum_task_cells=3,
    )
    assert not is_funnel_confirmation_hit(
        cluster=cluster(("a:1", "b:2")),
        settings=settings,
        minimum_task_cells=3,
    )


def test_funnel_deep_requires_two_seeds_for_each_of_two_sources() -> None:
    def cluster(task_cells: tuple[str, ...]) -> RefinedCluster:
        return RefinedCluster(
            cluster_id="cluster",
            candidate_id="candidate",
            receptor_id="receptor",
            decoy_ids=("decoy",),
            refinement_seeds=(120623, 120624, 120625),
            start_ids=("start_a", "start_b"),
            source_seeds=(120628, 120629),
            task_cells=task_cells,
            representative_decoy_id="decoy",
            supported=True,
        )

    assert is_funnel_deep_hit(
        cluster=cluster(
            ("start_a:120623", "start_a:120624", "start_b:120623", "start_b:120625")
        ),
        minimum_source_starts=2,
        minimum_rosetta_seeds_per_source=2,
    )
    assert not is_funnel_deep_hit(
        cluster=cluster(("start_a:120623", "start_a:120624", "start_b:120625")),
        minimum_source_starts=2,
        minimum_rosetta_seeds_per_source=2,
    )


def test_refinement_cluster_rejects_same_cabs_source_seed() -> None:
    decoys = (
        _refined_decoy(decoy_id="decoy_a", seed=5101, start_id="start_a", x_offset=0.0),
        _refined_decoy(decoy_id="decoy_b", seed=5103, start_id="start_b", x_offset=0.1),
    )
    settings = ResolvedAnalysisSettings(
        min_interface_contact_pairs=1,
        min_interface_receptor_residues=1,
        max_receptor_ca_rmsd_A=1.0,
        min_start_contact_overlap=0.5,
        max_start_site_displacement_A=1.0,
        max_cluster_backbone_rmsd_A=0.5,
        min_heavy_atom_distance_A=1.0,
        min_refinement_seed_support=2,
        min_refinement_start_support=2,
        min_refinement_source_seed_support=2,
    )

    clusters = cluster_refined_decoys(
        decoys=decoys, settings=settings, require_source_seed_support=True
    )

    assert len(clusters) == 1
    assert not clusters[0].supported
    assert clusters[0].source_seeds == (1001,)


def _diagnostic_decoy(
    *,
    decoy_id: str,
    seed: int,
    start_id: str,
    ranking_score: float,
    x_offset: float,
    native_rmsd_A: float,
    chemistry_valid: bool = True,
) -> DiagnosticDecoy:
    return DiagnosticDecoy(
        decoy_id=decoy_id,
        description=decoy_id,
        task_id=f"task_{decoy_id}",
        start_id=start_id,
        receptor_site_id="replaceable_site",
        source_seed=1001,
        refinement_seed=seed,
        path=Path(f"{decoy_id}.pdb"),
        sha256="0" * 64,
        ranking_score=ranking_score,
        chemistry_valid=chemistry_valid,
        chemistry_failure=None if chemistry_valid else "invalid disulfide geometry",
        receptor_alignment_rmsd_A=0.1,
        native_backbone_rmsd_A=native_rmsd_A,
        aligned_backbone=tuple(
            gemmi.Position(index + x_offset, 0.0, 0.0) for index in range(6)
        ),
    )


def test_qualification_diagnostic_clusters_only_chemistry_valid_decoys() -> None:
    decoys = (
        _diagnostic_decoy(
            decoy_id="far_score_best",
            seed=5101,
            start_id="start_a",
            ranking_score=-20.0,
            x_offset=4.0,
            native_rmsd_A=4.0,
        ),
        _diagnostic_decoy(
            decoy_id="near_a",
            seed=5101,
            start_id="start_a",
            ranking_score=-10.0,
            x_offset=0.0,
            native_rmsd_A=0.2,
        ),
        _diagnostic_decoy(
            decoy_id="near_b",
            seed=5103,
            start_id="start_b",
            ranking_score=-9.0,
            x_offset=0.1,
            native_rmsd_A=0.3,
        ),
        _diagnostic_decoy(
            decoy_id="invalid_near",
            seed=5105,
            start_id="start_c",
            ranking_score=-30.0,
            x_offset=0.0,
            native_rmsd_A=0.1,
            chemistry_valid=False,
        ),
    )

    clusters = cluster_diagnostic_decoys(
        decoys=decoys,
        cluster_threshold_A=0.5,
        recovery_threshold_A=2.0,
        min_seed_support=2,
        min_start_support=2,
    )

    assert len(clusters) == 2
    assert clusters[0].energy_representative_id == "far_score_best"
    assert clusters[1].supported
    assert clusters[1].recovered_seeds == (5101, 5103)
    assert clusters[1].recovered_start_ids == ("start_a", "start_b")
    assert all(
        member.decoy_id != "invalid_near"
        for cluster in clusters
        for member in cluster.members
    )


def _reported_site(
    *, site_id: str, receptor_id: str, position: tuple[float, float, float]
) -> ReportedReceptorSite:
    return ReportedReceptorSite(
        site_id=site_id,
        target="replaceable_target",
        receptor_id=receptor_id,
        coordinate_frame_id="replaceable_frame",
        pose_count=1,
        pose_ids=(f"{site_id}_pose",),
        supporting_seeds=(11, 13),
        representative_pose_id=f"{site_id}_pose",
        representative_contacts=frozenset({"10", "11"}),
        representative_position=position,
        supported=True,
    )


def test_cross_state_comparison_uses_frozen_site_distance() -> None:
    main_site = _reported_site(
        site_id="MAIN_SITE", receptor_id="MAIN_RECEPTOR", position=(0.0, 0.0, 0.0)
    )
    matched = _reported_site(
        site_id="STATE_MATCH", receptor_id="STATE_A", position=(0.2, 0.0, 0.0)
    )
    state_only = _reported_site(
        site_id="STATE_ONLY", receptor_id="STATE_B", position=(8.0, 0.0, 0.0)
    )
    main = SiteAnalysisReport(
        evidence_category="main_discovery",
        pose_path=Path("main_poses.tsv"),
        receptor_sites={main_site.site_id: main_site},
        candidate_sites={
            "CANDIDATE": ReportedCandidateSite(
                candidate_id="CANDIDATE",
                target=main_site.target,
                coordinate_frame_id=main_site.coordinate_frame_id,
                receptor_ids=(main_site.receptor_id,),
                receptor_site_ids=(main_site.site_id,),
                representative_site_id=main_site.site_id,
                receptor_support=1,
                evidence_tier="conformation_specific",
                rank_within_tier=1,
                minimum_seed_support=2,
                total_seed_support=2,
                maximum_normalized_site_distance=0.0,
                minimum_selected_pose_fraction=1.0,
                total_selected_pose_fraction=1.0,
                median_receptor_score_quantile=0.5,
                handoff_eligible=True,
            )
        },
        ensemble_candidate_budget=8,
        conformation_specific_candidate_budget=2,
        manifest_path=Path("main_manifest.json"),
    )
    replication = SiteAnalysisReport(
        evidence_category="bound_state_blind_replication",
        pose_path=Path("replication_poses.tsv"),
        receptor_sites={matched.site_id: matched, state_only.site_id: state_only},
        candidate_sites={},
        ensemble_candidate_budget=8,
        conformation_specific_candidate_budget=2,
        manifest_path=Path("replication_manifest.json"),
    )

    comparisons, replication_only = compare_sites(
        main=main,
        replication=replication,
        contact_limit=0.5,
        position_limit=2.0,
    )

    assert comparisons[0].replication_site_ids == ("STATE_MATCH",)
    assert comparisons[0].replication_state_ids == ("STATE_A",)
    assert comparisons[0].minimum_normalized_distance == pytest.approx(0.1)
    assert [site.site_id for site in replication_only] == ["STATE_ONLY"]


def test_rosetta_crash_log_is_confined_and_renamed(tmp_path: Path) -> None:
    crash_dir = rosetta_crash_log_dir(outputs_dir=tmp_path / "outputs")
    log_path = tmp_path / "task" / "refine.log"
    log_path.parent.mkdir()
    with pytest.raises(ValidationError, match="crash report"):
        run_rosetta_command(
            command=(
                "/bin/sh",
                "-c",
                "echo crash-dump > ROSETTA_CRASH.log; exit 1",
            ),
            log_path=log_path,
            crash_dir=crash_dir,
        )
    assert not (crash_dir / "ROSETTA_CRASH.log").exists()
    renamed = tuple(crash_dir.glob("*-refine-ROSETTA_CRASH.log"))
    assert len(renamed) == 1
    assert renamed[0].read_text(encoding="utf-8").strip() == "crash-dump"


def test_rosetta_success_leaves_no_crash_log(tmp_path: Path) -> None:
    crash_dir = tmp_path / "logs" / "rosetta_crashes"
    log_path = tmp_path / "task" / "prepack.log"
    log_path.parent.mkdir()
    run_rosetta_command(
        command=("/bin/sh", "-c", "exit 0"),
        log_path=log_path,
        crash_dir=crash_dir,
    )
    assert log_path.is_file()
    assert not (crash_dir / "ROSETTA_CRASH.log").exists()
    assert tuple(crash_dir.glob("*ROSETTA_CRASH.log")) == ()


def test_exploration_promotion_contract_uses_frozen_support_thresholds() -> None:
    config = load_config(PROJECT_CONFIG)

    contract = exploration_promotion_contract(config)

    eligibility = object_mapping(
        contract["candidate_eligibility"], name="candidate eligibility"
    )
    selection = object_mapping(
        contract["deep_refinement_selection"], name="deep refinement selection"
    )
    success = object_mapping(
        contract["deep_refinement_success"], name="deep refinement success"
    )
    assert eligibility["minimum_passed_starts_per_receptor_site"] == 2
    assert eligibility["qc_metrics_used_for_candidate_reranking"] is False
    assert selection["blind_discovery_arm_budget"] == 2
    assert selection["functional_annotation_arm_budget"] == 1
    assert success["minimum_independent_cabs_starts"] == 2
    assert success["minimum_independent_rosetta_seeds"] == 2
