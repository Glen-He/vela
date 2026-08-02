from collections.abc import Callable
from dataclasses import replace
from pathlib import Path

import gemmi
import pytest

from vela.config import load_config
from vela.core.provenance import atomic_write_json, sha256_file
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
from vela.validation.refinement.geometry import (
    RefinedDecoy,
    ResolvedAnalysisSettings,
    cluster_refined_decoys,
)
from vela.validation.refinement.planning import build_refinement_tasks
from vela.validation.refinement.reconstruction import (
    assess_topology_reconstruction,
    build_cg2all_command,
    write_cg2all_input,
    write_chemistry_protocol,
    write_reference_receptor_complex,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    build_prepack_command,
    build_refine_command,
)
from vela.validation.scores import read_rosetta_scorefile

PROJECT_ROOT = Path(__file__).resolve().parents[2]
PROJECT_CONFIG = PROJECT_ROOT / "configs"


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
    assert not settings.config_complete


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
            atom.pos = gemmi.Position(0.0 if index == 1 else 2.05, y_offset, 1.0)
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


def test_refinement_tasks_are_derived_from_handoff_and_configured_seeds(
    tmp_path: Path,
) -> None:
    config = load_config(PROJECT_CONFIG)
    config = replace(
        config,
        validation=replace(config.validation, seeds=(4101, 4103)),
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
            "schema": "vela.validation-handoff-manifest/2",
            "status": "completed",
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": "main_discovery_handoff",
            "known_site_information_used": False,
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

    assert [task.task_id for task in tasks] == ["refine_00001", "refine_00002"]
    assert [task.seed for task in tasks] == [4101, 4103]
    assert all(task.start.candidate_id == "CANDIDATE_001" for task in tasks)


def test_refinement_tasks_preserve_guided_source_identity(tmp_path: Path) -> None:
    config = load_config(PROJECT_CONFIG)
    config = replace(
        config,
        validation=replace(config.validation, seeds=(4201,)),
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

    assert len(tasks) == 1
    assert tasks[0].start.source_seed is None
    assert tasks[0].start.candidate_id == "guided__replaceable_template"


def _refined_decoy(
    *, decoy_id: str, seed: int, start_id: str, x_offset: float
) -> RefinedDecoy:
    backbone = tuple(gemmi.Position(index + x_offset, 0.0, 0.0) for index in range(6))
    return RefinedDecoy(
        decoy_id=decoy_id,
        task_id=f"task_{seed}",
        start_id=start_id,
        candidate_id="replaceable_candidate",
        receptor_id="replaceable_receptor",
        source_seed=1001,
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
        cluster_backbone=backbone,
    )


def test_refinement_cluster_support_is_config_driven() -> None:
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
    )

    clusters = cluster_refined_decoys(decoys=decoys, settings=settings)

    assert len(clusters) == 1
    assert clusters[0].supported
    assert clusters[0].refinement_seeds == (5101, 5103)
    assert clusters[0].start_ids == ("start_a", "start_b")


def _reported_site(
    *, site_id: str, receptor_id: str, position: tuple[float, float, float]
) -> ReportedReceptorSite:
    return ReportedReceptorSite(
        site_id=site_id,
        target="replaceable_target",
        receptor_id=receptor_id,
        coordinate_frame_id="replaceable_frame",
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
                supported=True,
            )
        },
        manifest_path=Path("main_manifest.json"),
    )
    replication = SiteAnalysisReport(
        evidence_category="bound_state_blind_replication",
        pose_path=Path("replication_poses.tsv"),
        receptor_sites={matched.site_id: matched, state_only.site_id: state_only},
        candidate_sites={},
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
