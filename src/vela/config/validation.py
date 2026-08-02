"""阶段三结构复核与局部精修配置的严格组装。"""

from collections.abc import Mapping
from pathlib import Path

from vela.config.models import ConfigError
from vela.config.values import (
    assert_keys,
    boolean,
    document,
    integer,
    number,
    resolved_path,
    string,
    table,
)
from vela.core.typed_data import object_list
from vela.discovery.models import UNRESOLVED
from vela.preparation.chemistry import DisulfideBond, HistidineState
from vela.validation.models import (
    BoundStateDefinition,
    CandidateHandoffSettings,
    CandidateRefinementSettings,
    Cg2AllSettings,
    EnvironmentReference,
    GuidedTemplate,
    LocalRecoveryControl,
    RosettaSettings,
    ValidationAnalysisSettings,
    ValidationSettings,
)


def _optional_text(source: Mapping[str, object], key: str, *, path: str) -> str | None:
    value = string(source, key, path=path)
    return None if value == UNRESOLVED else value


def _optional_number(
    source: Mapping[str, object], key: str, *, path: str
) -> float | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    return number(source, key, path=path)


def _optional_integer(
    source: Mapping[str, object], key: str, *, path: str
) -> int | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    return integer(source, key, path=path)


def _integers(value: object, *, name: str) -> tuple[int, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of integers") from exc
    result: list[int] = []
    for item in values:
        if not isinstance(item, int) or isinstance(item, bool):
            raise ConfigError(f"{name} must be an array of integers")
        result.append(item)
    return tuple(result)


def _seed_batches(value: object, *, name: str) -> tuple[tuple[int, ...], ...]:
    try:
        batches = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of integer arrays") from exc
    return tuple(
        _integers(batch, name=f"{name}[{index}]") for index, batch in enumerate(batches)
    )


def _texts(value: object, *, name: str) -> tuple[str, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of strings") from exc
    if any(not isinstance(item, str) or not item.strip() for item in values):
        raise ConfigError(f"{name} must be an array of non-empty strings")
    return tuple(item for item in values if isinstance(item, str))


def _disulfides(value: object, *, path: str) -> tuple[DisulfideBond, ...]:
    try:
        values = object_list(value, name=path)
    except TypeError as exc:
        raise ConfigError(f"{path} must be an array of integer pairs") from exc
    result: list[DisulfideBond] = []
    for index, item in enumerate(values):
        try:
            pair = object_list(item, name=f"{path}[{index}]")
        except TypeError as exc:
            raise ConfigError(f"{path}[{index}] must contain two integers") from exc
        if len(pair) != 2:
            raise ConfigError(f"{path}[{index}] must contain two integers")
        first, second = pair
        if (
            not isinstance(first, int)
            or isinstance(first, bool)
            or not isinstance(second, int)
            or isinstance(second, bool)
        ):
            raise ConfigError(f"{path}[{index}] must contain two integers")
        result.append(DisulfideBond(first, second))
    return tuple(result)


def _histidines(value: object, *, path: str) -> tuple[HistidineState, ...]:
    section = document(value, name=path)
    result: list[HistidineState] = []
    for raw_position, raw_state in section.items():
        try:
            position = int(raw_position)
        except ValueError as exc:
            raise ConfigError(f"{path} keys must be residue positions") from exc
        if isinstance(raw_state, bool) or not isinstance(raw_state, str):
            raise ConfigError(f"{path}.{raw_position} must be text")
        result.append(HistidineState(position, raw_state))
    return tuple(sorted(result, key=lambda item: item.position))


def _bound_states(value: object) -> tuple[BoundStateDefinition, ...]:
    try:
        entries = object_list(value, name="validation.bound_states")
    except TypeError as exc:
        raise ConfigError("validation.bound_states must be an array of tables") from exc
    if not entries:
        raise ConfigError("validation.bound_states must not be empty")
    result: list[BoundStateDefinition] = []
    required = {
        "state_id",
        "receptor_id",
        "ligand_id",
        "ligand_author_chain_id",
        "local_control_kind",
        "ligand_sequence",
        "disulfide_bonds",
        "histidines",
        "selection_reason",
    }
    for index, item in enumerate(entries):
        path = f"validation.bound_states[{index}]"
        section = document(item, name=path)
        assert_keys(section, allowed=required, required=required, path=path)
        raw_sequence = string(section, "ligand_sequence", path=path)
        result.append(
            BoundStateDefinition(
                state_id=string(section, "state_id", path=path),
                receptor_id=string(section, "receptor_id", path=path),
                ligand_id=string(section, "ligand_id", path=path),
                ligand_author_chain_id=string(
                    section, "ligand_author_chain_id", path=path
                ),
                local_control_kind=string(section, "local_control_kind", path=path),
                ligand_sequence=(None if raw_sequence == UNRESOLVED else raw_sequence),
                disulfide_bonds=_disulfides(
                    section["disulfide_bonds"], path=f"{path}.disulfide_bonds"
                ),
                histidines=_histidines(
                    section["histidines"], path=f"{path}.histidines"
                ),
                selection_reason=string(section, "selection_reason", path=path),
            )
        )
    return tuple(result)


def _local_controls(value: object) -> tuple[LocalRecoveryControl, ...]:
    try:
        entries = object_list(value, name="validation.local_controls")
    except TypeError as exc:
        raise ConfigError(
            "validation.local_controls must be an array of tables"
        ) from exc
    result: list[LocalRecoveryControl] = []
    required = {
        "control_id",
        "bound_state_id",
        "prepack_seed",
        "seed_batches",
        "random_translation_A",
        "random_rotation_degrees",
        "ranking_score",
        "recovery_rmsd_score",
        "top_clusters",
        "max_recovery_rmsd_A",
        "max_cluster_backbone_rmsd_A",
        "min_cluster_seed_support",
        "max_batch_pose_rmsd_A",
    }
    for index, item in enumerate(entries):
        path = f"validation.local_controls[{index}]"
        section = document(item, name=path)
        assert_keys(section, allowed=required, required=required, path=path)
        result.append(
            LocalRecoveryControl(
                control_id=string(section, "control_id", path=path),
                bound_state_id=string(section, "bound_state_id", path=path),
                prepack_seed=integer(section, "prepack_seed", path=path),
                seed_batches=_seed_batches(
                    section["seed_batches"], name=f"{path}.seed_batches"
                ),
                random_translation_A=number(section, "random_translation_A", path=path),
                random_rotation_degrees=number(
                    section, "random_rotation_degrees", path=path
                ),
                ranking_score=string(section, "ranking_score", path=path),
                recovery_rmsd_score=string(section, "recovery_rmsd_score", path=path),
                top_clusters=integer(section, "top_clusters", path=path),
                max_recovery_rmsd_A=number(section, "max_recovery_rmsd_A", path=path),
                max_cluster_backbone_rmsd_A=number(
                    section, "max_cluster_backbone_rmsd_A", path=path
                ),
                min_cluster_seed_support=integer(
                    section, "min_cluster_seed_support", path=path
                ),
                max_batch_pose_rmsd_A=number(
                    section, "max_batch_pose_rmsd_A", path=path
                ),
            )
        )
    return tuple(result)


def _guided_templates(value: object) -> tuple[GuidedTemplate, ...]:
    try:
        entries = object_list(value, name="validation.guided_templates")
    except TypeError as exc:
        raise ConfigError(
            "validation.guided_templates must be an array of tables"
        ) from exc
    result: list[GuidedTemplate] = []
    required = {
        "template_id",
        "bound_state_id",
        "ligand_positions",
        "selection_reason",
    }
    for index, item in enumerate(entries):
        path = f"validation.guided_templates[{index}]"
        section = document(item, name=path)
        assert_keys(section, allowed=required, required=required, path=path)
        result.append(
            GuidedTemplate(
                template_id=string(section, "template_id", path=path),
                bound_state_id=string(section, "bound_state_id", path=path),
                ligand_positions=_integers(
                    section["ligand_positions"], name=f"{path}.ligand_positions"
                ),
                selection_reason=string(section, "selection_reason", path=path),
            )
        )
    return tuple(result)


def _environment_references(value: object) -> tuple[EnvironmentReference, ...]:
    try:
        entries = object_list(value, name="validation.environment_references")
    except TypeError as exc:
        raise ConfigError(
            "validation.environment_references must be an array of tables"
        ) from exc
    result: list[EnvironmentReference] = []
    required = {
        "reference_id",
        "receptor_id",
        "assembly_id",
        "beta_author_chain_ids",
        "other_catalytic_author_chain_ids",
        "evaluation_targets",
        "construct_note",
    }
    for index, item in enumerate(entries):
        path = f"validation.environment_references[{index}]"
        section = document(item, name=path)
        assert_keys(section, allowed=required, required=required, path=path)
        result.append(
            EnvironmentReference(
                reference_id=string(section, "reference_id", path=path),
                receptor_id=string(section, "receptor_id", path=path),
                assembly_id=string(section, "assembly_id", path=path),
                beta_author_chain_ids=_texts(
                    section["beta_author_chain_ids"],
                    name=f"{path}.beta_author_chain_ids",
                ),
                other_catalytic_author_chain_ids=_texts(
                    section["other_catalytic_author_chain_ids"],
                    name=f"{path}.other_catalytic_author_chain_ids",
                ),
                evaluation_targets=_texts(
                    section["evaluation_targets"],
                    name=f"{path}.evaluation_targets",
                ),
                construct_note=string(section, "construct_note", path=path),
            )
        )
    return tuple(result)


def parse_validation(
    source: Mapping[str, object], *, config_dir: Path
) -> ValidationSettings:
    """组装 `[validation]`、工具和实验结合态登记。"""
    section = table(source, "validation", path="")
    required = {
        "method_id",
        "qualification_status",
        "qualification_report",
        "qualification_report_sha256",
        "seeds",
        "min_disulfide_sg_A",
        "max_disulfide_sg_A",
        "interface_contact_A",
        "rosetta",
        "cg2all",
        "handoff",
        "refinement",
        "analysis",
        "bound_states",
        "local_controls",
        "guided_templates",
        "environment_references",
    }
    assert_keys(section, allowed=required, required=required, path="validation")
    rosetta = document(section["rosetta"], name="validation.rosetta")
    rosetta_keys = {
        "executable",
        "scripts_executable",
        "database",
        "version_file",
        "expected_version",
        "parallel_tasks",
        "decoys_per_seed",
        "score_function",
        "lowres_preoptimize",
    }
    assert_keys(
        rosetta,
        allowed=rosetta_keys,
        required=rosetta_keys,
        path="validation.rosetta",
    )
    cg2all = document(section["cg2all"], name="validation.cg2all")
    cg2all_keys = {
        "executable",
        "package_metadata",
        "checkpoint",
        "checkpoint_sha256",
        "expected_version",
        "representation",
        "receptor_histidine_state",
        "device",
        "processes",
        "batch_size",
        "chain_break_cutoff_A",
        "max_ca_rmsd_A",
    }
    assert_keys(
        cg2all,
        allowed=cg2all_keys,
        required=cg2all_keys,
        path="validation.cg2all",
    )
    handoff = document(section["handoff"], name="validation.handoff")
    handoff_keys = {"poses_per_receptor_site", "chemistry_seed"}
    assert_keys(
        handoff,
        allowed=handoff_keys,
        required=handoff_keys,
        path="validation.handoff",
    )
    refinement = document(section["refinement"], name="validation.refinement")
    refinement_keys = {
        "prepack_seed",
        "random_translation_A",
        "random_rotation_degrees",
        "ranking_score",
    }
    assert_keys(
        refinement,
        allowed=refinement_keys,
        required=refinement_keys,
        path="validation.refinement",
    )
    analysis = document(section["analysis"], name="validation.analysis")
    analysis_keys = {
        "min_interface_contact_pairs",
        "min_interface_receptor_residues",
        "max_receptor_ca_rmsd_A",
        "min_start_contact_overlap",
        "max_start_site_displacement_A",
        "max_cluster_backbone_rmsd_A",
        "min_heavy_atom_distance_A",
        "min_refinement_seed_support",
        "min_refinement_start_support",
    }
    assert_keys(
        analysis,
        allowed=analysis_keys,
        required=analysis_keys,
        path="validation.analysis",
    )
    report = _optional_text(section, "qualification_report", path="validation")
    return ValidationSettings(
        method_id=string(section, "method_id", path="validation"),
        qualification_status=string(section, "qualification_status", path="validation"),
        qualification_report=(
            resolved_path(report, config_dir=config_dir) if report else None
        ),
        qualification_report_sha256=_optional_text(
            section, "qualification_report_sha256", path="validation"
        ),
        seeds=_integers(section["seeds"], name="validation.seeds"),
        min_disulfide_sg_A=number(section, "min_disulfide_sg_A", path="validation"),
        max_disulfide_sg_A=number(section, "max_disulfide_sg_A", path="validation"),
        interface_contact_A=number(section, "interface_contact_A", path="validation"),
        rosetta=RosettaSettings(
            executable=resolved_path(
                string(rosetta, "executable", path="validation.rosetta"),
                config_dir=config_dir,
            ),
            scripts_executable=resolved_path(
                string(rosetta, "scripts_executable", path="validation.rosetta"),
                config_dir=config_dir,
            ),
            database=resolved_path(
                string(rosetta, "database", path="validation.rosetta"),
                config_dir=config_dir,
            ),
            version_file=resolved_path(
                string(rosetta, "version_file", path="validation.rosetta"),
                config_dir=config_dir,
            ),
            expected_version=string(
                rosetta, "expected_version", path="validation.rosetta"
            ),
            parallel_tasks=integer(
                rosetta, "parallel_tasks", path="validation.rosetta"
            ),
            decoys_per_seed=integer(
                rosetta, "decoys_per_seed", path="validation.rosetta"
            ),
            score_function=string(rosetta, "score_function", path="validation.rosetta"),
            lowres_preoptimize=boolean(
                rosetta, "lowres_preoptimize", path="validation.rosetta"
            ),
        ),
        cg2all=Cg2AllSettings(
            executable=resolved_path(
                string(cg2all, "executable", path="validation.cg2all"),
                config_dir=config_dir,
            ),
            package_metadata=resolved_path(
                string(cg2all, "package_metadata", path="validation.cg2all"),
                config_dir=config_dir,
            ),
            checkpoint=resolved_path(
                string(cg2all, "checkpoint", path="validation.cg2all"),
                config_dir=config_dir,
            ),
            checkpoint_sha256=string(
                cg2all, "checkpoint_sha256", path="validation.cg2all"
            ),
            expected_version=string(
                cg2all, "expected_version", path="validation.cg2all"
            ),
            representation=string(cg2all, "representation", path="validation.cg2all"),
            receptor_histidine_state=string(
                cg2all, "receptor_histidine_state", path="validation.cg2all"
            ),
            device=string(cg2all, "device", path="validation.cg2all"),
            processes=integer(cg2all, "processes", path="validation.cg2all"),
            batch_size=integer(cg2all, "batch_size", path="validation.cg2all"),
            chain_break_cutoff_A=number(
                cg2all, "chain_break_cutoff_A", path="validation.cg2all"
            ),
            max_ca_rmsd_A=number(cg2all, "max_ca_rmsd_A", path="validation.cg2all"),
        ),
        handoff=CandidateHandoffSettings(
            poses_per_receptor_site=integer(
                handoff, "poses_per_receptor_site", path="validation.handoff"
            ),
            chemistry_seed=integer(
                handoff, "chemistry_seed", path="validation.handoff"
            ),
        ),
        refinement=CandidateRefinementSettings(
            prepack_seed=integer(
                refinement, "prepack_seed", path="validation.refinement"
            ),
            random_translation_A=number(
                refinement, "random_translation_A", path="validation.refinement"
            ),
            random_rotation_degrees=number(
                refinement,
                "random_rotation_degrees",
                path="validation.refinement",
            ),
            ranking_score=string(
                refinement, "ranking_score", path="validation.refinement"
            ),
        ),
        analysis=ValidationAnalysisSettings(
            min_interface_contact_pairs=_optional_integer(
                analysis,
                "min_interface_contact_pairs",
                path="validation.analysis",
            ),
            min_interface_receptor_residues=_optional_integer(
                analysis,
                "min_interface_receptor_residues",
                path="validation.analysis",
            ),
            max_receptor_ca_rmsd_A=_optional_number(
                analysis, "max_receptor_ca_rmsd_A", path="validation.analysis"
            ),
            min_start_contact_overlap=_optional_number(
                analysis, "min_start_contact_overlap", path="validation.analysis"
            ),
            max_start_site_displacement_A=_optional_number(
                analysis,
                "max_start_site_displacement_A",
                path="validation.analysis",
            ),
            max_cluster_backbone_rmsd_A=_optional_number(
                analysis,
                "max_cluster_backbone_rmsd_A",
                path="validation.analysis",
            ),
            min_heavy_atom_distance_A=_optional_number(
                analysis, "min_heavy_atom_distance_A", path="validation.analysis"
            ),
            min_refinement_seed_support=_optional_integer(
                analysis,
                "min_refinement_seed_support",
                path="validation.analysis",
            ),
            min_refinement_start_support=_optional_integer(
                analysis,
                "min_refinement_start_support",
                path="validation.analysis",
            ),
        ),
        bound_states=_bound_states(section["bound_states"]),
        local_controls=_local_controls(section["local_controls"]),
        guided_templates=_guided_templates(section["guided_templates"]),
        environment_references=_environment_references(
            section["environment_references"]
        ),
    )
