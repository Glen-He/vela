"""CABS CA/SC pose 的已资格化全原子闭环恢复流水线。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import JsonValue, atomic_write_text
from vela.preparation.chemistry import ChemistryDefinition
from vela.validation.models import ValidationError
from vela.validation.records import file_record
from vela.validation.refinement.reconstruction import (
    TopologyReconstructionAssessment,
    assess_topology_reconstruction,
    build_cg2all_command,
    write_cg2all_input,
    write_disulfide_indices,
    write_peptide_site_coordinate_constraints,
    write_reference_receptor_complex,
    write_topology_rebuild_protocol,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    build_prepack_command,
    build_topology_refine_command,
    rosetta_crash_log_dir,
    run_rosetta_command,
    single_rosetta_pdb_output,
)
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class TopologyRecoveryResult:
    """一次粗粒化 pose 全原子闭环恢复的完整结果。"""

    execution_status: str
    reconstruction_status: str
    failure_reasons: tuple[str, ...]
    failure_detail: str | None
    metrics: dict[str, JsonValue] | None
    commands: dict[str, JsonValue]
    artifacts: dict[str, JsonValue]
    cg_input_path: Path
    final_path: Path | None
    receptor_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]


def topology_assessment_record(
    assessment: TopologyReconstructionAssessment,
) -> dict[str, JsonValue]:
    """将全原子拓扑评估转换为稳定、可审计的记录。"""
    return {
        "receptor_ca_rmsd_A": round(assessment.receptor_ca_rmsd_A, 6),
        "peptide_pose_ca_rmsd_A": round(assessment.peptide_pose_ca_rmsd_A, 6),
        "peptide_internal_ca_rmsd_A": round(assessment.peptide_internal_ca_rmsd_A, 6),
        "ligand_centroid_displacement_A": round(
            assessment.ligand_centroid_displacement_A, 6
        ),
        "receptor_contact_retention_fraction": round(
            assessment.receptor_contact_retention_fraction, 6
        ),
        "disulfide_sg_distances_A": [
            round(value, 6) for value in assessment.disulfide_sg_distances_A
        ],
        "min_interchain_heavy_atom_distance_A": round(
            assessment.min_interchain_heavy_atom_distance_A, 6
        ),
        "min_nonlocal_peptide_heavy_atom_distance_A": round(
            assessment.min_nonlocal_peptide_heavy_atom_distance_A, 6
        ),
    }


def assess_recovered_topology(
    *,
    config: AppConfig,
    input_path: Path,
    output_path: Path,
    chemistry: ChemistryDefinition,
) -> TopologyReconstructionAssessment:
    """以唯一的资格阈值评估粗粒化 pose 的全原子恢复。"""
    calibration = config.discovery.qualification.topology_calibration
    return assess_topology_reconstruction(
        input_path=input_path,
        output_path=output_path,
        chemistry=chemistry,
        settings=config.validation.cg2all,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        min_interchain_heavy_atom_distance_A=(
            calibration.min_interchain_heavy_atom_distance_A
        ),
        min_nonlocal_peptide_heavy_atom_distance_A=(
            calibration.min_nonlocal_peptide_heavy_atom_distance_A
        ),
        contact_ca_threshold_A=calibration.contact_ca_threshold_A,
        max_peptide_internal_ca_rmsd_A=calibration.max_peptide_internal_ca_rmsd_A,
        max_ligand_centroid_displacement_A=(
            calibration.max_ligand_centroid_displacement_A
        ),
        min_receptor_contact_retention_fraction=(
            calibration.min_receptor_contact_retention_fraction
        ),
    )


def _invalid(
    *,
    reason: str,
    detail: str,
    commands: dict[str, JsonValue],
    artifacts: dict[str, JsonValue],
    cg_input_path: Path,
    receptor_residue_count: int,
    fixed_histidines: tuple[int, ...],
) -> TopologyRecoveryResult:
    return TopologyRecoveryResult(
        execution_status="invalid",
        reconstruction_status="not_assessed",
        failure_reasons=(reason,),
        failure_detail=detail,
        metrics=None,
        commands=commands,
        artifacts=artifacts,
        cg_input_path=cg_input_path,
        final_path=None,
        receptor_residue_count=receptor_residue_count,
        fixed_histidine_pose_indices=fixed_histidines,
    )


def recover_topology(
    *,
    config: AppConfig,
    source_path: Path,
    model_index: int,
    reference_receptor_path: Path,
    chemistry: ChemistryDefinition,
    task_dir: Path,
    rosetta_seed: int,
) -> TopologyRecoveryResult:
    """执行校准过的 cg2all、二硫键重建、prepack 和受约束局部恢复。"""
    calibration = config.discovery.qualification.topology_calibration
    artifacts: dict[str, JsonValue] = {}
    commands: dict[str, JsonValue] = {}
    cg_input_path = task_dir / "cg2all_input.pdb"
    cg_input = write_cg2all_input(
        source_path=source_path,
        model_index=model_index,
        destination=cg_input_path,
        chemistry=chemistry,
        settings=config.validation.cg2all,
    )
    receptor_count = cg_input.receptor_residue_count
    fixed_histidines = cg_input.fixed_histidine_pose_indices
    artifacts["cg2all_input"] = file_record(cg_input_path, root=task_dir)
    cg_output_path = task_dir / "cg2all_raw.pdb"
    cg_command = build_cg2all_command(
        settings=config.validation.cg2all,
        input_path=cg_input_path,
        output_path=cg_output_path,
    )
    commands["cg2all"] = list(cg_command)
    cg_log_path = task_dir / "cg2all.log"
    try:
        run_rosetta_command(
            command=cg_command,
            log_path=cg_log_path,
            crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
            thread_count=config.validation.cg2all.processes,
        )
    except ValidationError as exc:
        artifacts["cg2all_log"] = file_record(cg_log_path, root=task_dir)
        return _invalid(
            reason="cg2all_reconstruction_failed",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    artifacts["cg2all_log"] = file_record(cg_log_path, root=task_dir)
    artifacts["cg2all_raw"] = file_record(cg_output_path, root=task_dir)
    grafted_path = task_dir / "all_atom_reference_receptor.pdb"
    try:
        write_reference_receptor_complex(
            coarse_pose_path=cg_input_path,
            reconstructed_path=cg_output_path,
            reference_receptor_path=reference_receptor_path,
            destination=grafted_path,
            chemistry=chemistry,
        )
    except ValidationError as exc:
        return _invalid(
            reason="cg2all_output_invalid",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    artifacts["all_atom_reference_receptor"] = file_record(grafted_path, root=task_dir)
    disulfide_path = task_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=receptor_count,
        chemistry=chemistry,
    )
    artifacts["fix_disulfide"] = file_record(disulfide_path, root=task_dir)
    protocol_path = task_dir / "rebuild_topology.xml"
    write_topology_rebuild_protocol(
        destination=protocol_path,
        receptor_residue_count=receptor_count,
        chemistry=chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    artifacts["topology_rebuild_protocol"] = file_record(protocol_path, root=task_dir)
    rebuild_dir = task_dir / "topology_rebuild"
    rebuild_dir.mkdir()
    rebuild_command = build_chemistry_command(
        settings=config.validation.rosetta,
        input_path=grafted_path,
        protocol_path=protocol_path,
        disulfide_path=disulfide_path,
        output_dir=rebuild_dir,
        seed=rosetta_seed,
        fixed_histidine_pose_indices=fixed_histidines,
    )
    commands["disulfide_rebuild"] = list(rebuild_command)
    rebuild_log = task_dir / "disulfide_rebuild.log"
    try:
        run_rosetta_command(
            command=rebuild_command,
            log_path=rebuild_log,
            crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        )
        rebuild_scores = read_rosetta_scorefile(rebuild_dir / "chemistry.sc")
        if len(rebuild_scores) != 1:
            raise ValidationError("disulfide rebuild must produce one score row")
        rebuilt_output = single_rosetta_pdb_output(rebuild_dir)
    except ValidationError as exc:
        artifacts["disulfide_rebuild_log"] = file_record(rebuild_log, root=task_dir)
        return _invalid(
            reason="rosetta_disulfide_rebuild_failed",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    artifacts["disulfide_rebuild_log"] = file_record(rebuild_log, root=task_dir)
    rebuilt_path = task_dir / "all_atom_disulfide_rebuilt.pdb"
    atomic_write_text(rebuilt_path, rebuilt_output.read_text(encoding="utf-8"))
    artifacts["all_atom_disulfide_rebuilt"] = file_record(rebuilt_path, root=task_dir)
    prepack_dir = task_dir / "prepack"
    prepack_dir.mkdir()
    prepack_command = build_prepack_command(
        settings=config.validation.rosetta,
        input_path=rebuilt_path,
        disulfide_path=disulfide_path,
        output_dir=prepack_dir,
        seed=rosetta_seed,
        fixed_histidine_pose_indices=fixed_histidines,
    )
    commands["disulfide_prepack"] = list(prepack_command)
    prepack_log = task_dir / "disulfide_prepack.log"
    try:
        run_rosetta_command(
            command=prepack_command,
            log_path=prepack_log,
            crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        )
        prepack_scores = read_rosetta_scorefile(prepack_dir / "prepack.sc")
        if len(prepack_scores) != 1:
            raise ValidationError("disulfide prepack must produce one score row")
        prepacked_output = single_rosetta_pdb_output(prepack_dir)
    except ValidationError as exc:
        artifacts["disulfide_prepack_log"] = file_record(prepack_log, root=task_dir)
        return _invalid(
            reason="rosetta_disulfide_prepack_failed",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    artifacts["disulfide_prepack_log"] = file_record(prepack_log, root=task_dir)
    prepacked_path = task_dir / "all_atom_prepacked.pdb"
    atomic_write_text(prepacked_path, prepacked_output.read_text(encoding="utf-8"))
    artifacts["all_atom_prepacked"] = file_record(prepacked_path, root=task_dir)
    constraint_path = task_dir / "site_coordinate_constraints.cst"
    constraint_count = write_peptide_site_coordinate_constraints(
        source_path=prepacked_path,
        destination=constraint_path,
        chemistry=chemistry,
        flat_width_A=calibration.site_coordinate_constraint_flat_width_A,
        standard_deviation_A=calibration.site_coordinate_constraint_sd_A,
    )
    if constraint_count != len(chemistry.sequence):
        raise ValidationError("topology refine site constraint count is invalid")
    artifacts["site_coordinate_constraints"] = file_record(
        constraint_path, root=task_dir
    )
    refine_dir = task_dir / "topology_refine"
    refine_dir.mkdir()
    refine_command = build_topology_refine_command(
        settings=config.validation.rosetta,
        input_path=prepacked_path,
        disulfide_path=disulfide_path,
        output_dir=refine_dir,
        seed=rosetta_seed,
        site_constraint_path=constraint_path,
        site_constraint_weight=calibration.site_coordinate_constraint_weight,
        fixed_histidine_pose_indices=fixed_histidines,
    )
    commands["topology_refine"] = list(refine_command)
    refine_log = task_dir / "topology_refine.log"
    try:
        run_rosetta_command(
            command=refine_command,
            log_path=refine_log,
            crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        )
        refine_scores = read_rosetta_scorefile(refine_dir / "refine.sc")
        if len(refine_scores) != 1:
            raise ValidationError("topology refine must produce one score row")
        refined_output = single_rosetta_pdb_output(refine_dir)
    except ValidationError as exc:
        artifacts["topology_refine_log"] = file_record(refine_log, root=task_dir)
        return _invalid(
            reason="flexpepdock_topology_refine_failed",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    artifacts["topology_refine_log"] = file_record(refine_log, root=task_dir)
    final_path = task_dir / "all_atom_topology.pdb"
    atomic_write_text(final_path, refined_output.read_text(encoding="utf-8"))
    artifacts["all_atom_topology"] = file_record(final_path, root=task_dir)
    try:
        assessment = assess_recovered_topology(
            config=config,
            input_path=cg_input_path,
            output_path=final_path,
            chemistry=chemistry,
        )
    except ValidationError as exc:
        return _invalid(
            reason="all_atom_output_invalid",
            detail=str(exc),
            commands=commands,
            artifacts=artifacts,
            cg_input_path=cg_input_path,
            receptor_residue_count=receptor_count,
            fixed_histidines=fixed_histidines,
        )
    refine_fa_rep = refine_scores[0].score("fa_rep")
    refine_backbone_strain = max(refine_scores[0].score("omega"), 0.0) + max(
        refine_scores[0].score("rama_prepro"), 0.0
    )
    metrics: dict[str, JsonValue] = {
        **topology_assessment_record(assessment),
        "rosetta_scores": {
            "rebuild_total_score": rebuild_scores[0].score("total_score"),
            "rebuild_dslf_fa13": rebuild_scores[0].score("dslf_fa13"),
            "prepack_total_score": prepack_scores[0].score("total_score"),
            "prepack_dslf_fa13": prepack_scores[0].score("dslf_fa13"),
            "prepack_fa_rep": prepack_scores[0].score("fa_rep"),
            "refine_total_score": refine_scores[0].score("total_score"),
            "refine_dslf_fa13": refine_scores[0].score("dslf_fa13"),
            "refine_fa_rep": refine_fa_rep,
            "refine_omega": refine_scores[0].score("omega"),
            "refine_rama_prepro": refine_scores[0].score("rama_prepro"),
            "refine_positive_backbone_strain": refine_backbone_strain,
        },
    }
    return TopologyRecoveryResult(
        execution_status="completed",
        reconstruction_status="passed" if not assessment.failures else "failed",
        failure_reasons=assessment.failures,
        failure_detail=None,
        metrics=metrics,
        commands=commands,
        artifacts=artifacts,
        cg_input_path=cg_input_path,
        final_path=final_path,
        receptor_residue_count=receptor_count,
        fixed_histidine_pose_indices=fixed_histidines,
    )
