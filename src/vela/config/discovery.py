"""阶段二全表面发现配置的严格组装。"""

from collections.abc import Mapping
from pathlib import Path

from vela.config.models import ConfigError
from vela.config.values import (
    assert_keys,
    document,
    integer,
    number,
    resolved_path,
    string,
    strings,
    table,
)
from vela.core.typed_data import object_list
from vela.discovery.models import (
    UNRESOLVED,
    CabsDockSettings,
    DiscoveryQualificationSettings,
    DiscoverySettings,
    DiscoveryTargetSettings,
    ReceptorEnsembleSettings,
    SiteAnalysisSettings,
    TopologyCalibrationSettings,
)


def _optional_text(source: Mapping[str, object], key: str, *, path: str) -> str | None:
    value = string(source, key, path=path)
    return None if value == UNRESOLVED else value


def _optional_float(
    source: Mapping[str, object], key: str, *, path: str
) -> float | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be a number or unresolved")
    return float(value)


def _optional_integer(
    source: Mapping[str, object], key: str, *, path: str
) -> int | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    if not isinstance(value, int) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be an integer or unresolved")
    return value


def _seeds(value: object, *, name: str = "discovery.seeds") -> tuple[int, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of integers") from exc
    seeds: list[int] = []
    for value in values:
        if not isinstance(value, int) or isinstance(value, bool):
            raise ConfigError(f"{name} must be an array of integers")
        seeds.append(value)
    return tuple(seeds)


def _number_array(value: object, *, name: str) -> tuple[float, ...]:
    try:
        values = object_list(value, name=name)
    except TypeError as exc:
        raise ConfigError(f"{name} must be an array of numbers") from exc
    result: list[float] = []
    for value in values:
        if not isinstance(value, (int, float)) or isinstance(value, bool):
            raise ConfigError(f"{name} must be an array of numbers")
        result.append(float(value))
    return tuple(result)


def parse_discovery(
    source: Mapping[str, object], *, config_dir: Path
) -> DiscoverySettings:
    """组装阶段二方法、资格规则和逐靶标放行配置。"""
    section = table(source, "discovery", path="")
    required = {
        "method_id",
        "adapter_id",
        "seeds",
        "ensemble",
        "cabsdock",
        "qualification",
        "targets",
    }
    assert_keys(section, allowed=required, required=required, path="discovery")
    ensemble = document(section["ensemble"], name="discovery.ensemble")
    cabsdock = document(section["cabsdock"], name="discovery.cabsdock")
    qualification = document(section["qualification"], name="discovery.qualification")
    topology_calibration = document(
        qualification["topology_calibration"],
        name="discovery.qualification.topology_calibration",
    )
    targets = document(section["targets"], name="discovery.targets")
    ensemble_keys = {
        "min_receptors_per_target",
        "allowed_structure_states",
    }
    assert_keys(
        ensemble,
        allowed=ensemble_keys,
        required=ensemble_keys,
        path="discovery.ensemble",
    )
    cabsdock_keys = {
        "executable",
        "source_dir",
        "source_revision",
        "patch_file",
        "seed_workers",
        "peptide_secondary_structure",
        "mc_annealing",
        "mc_cycles",
        "mc_steps",
        "replicas",
        "replicas_dtemp",
        "temperature_initial",
        "temperature_final",
        "binding_interactions",
        "protein_restraint_gap",
        "protein_restraint_min_A",
        "protein_restraint_max_A",
        "filtering_count",
        "clustering_medoids",
        "clustering_iterations",
        "trajectory_contact_ca_threshold_A",
        "max_disulfide_ca_distance_A",
        "min_models_for_selection",
        "selection_contact_jaccard_distance",
        "selection_position_distance_A",
        "pose_clustering_rmsd_A",
    }
    assert_keys(
        cabsdock,
        allowed=cabsdock_keys,
        required=cabsdock_keys,
        path="discovery.cabsdock",
    )
    qualification_keys = {
        "seeds",
        "control_bound_state_id",
        "control_target_id",
        "control_secondary_structure",
        "max_native_ligand_rmsd_A",
        "min_native_receptor_contact_fraction",
        "min_successful_control_seeds",
        "topology_calibration_status",
        "topology_calibration_report",
        "topology_calibration_report_sha256",
        "topology_calibration",
    }
    assert_keys(
        qualification,
        allowed=qualification_keys,
        required=qualification_keys,
        path="discovery.qualification",
    )
    topology_calibration_keys = {
        "candidate_ca_thresholds_A",
        "above_threshold_comparator_upper_bound_A",
        "models_per_stratum",
        "min_success_fraction_per_stratum",
        "min_successful_seeds_per_stratum",
        "min_interchain_heavy_atom_distance_A",
        "max_peptide_internal_ca_rmsd_A",
        "max_ligand_centroid_displacement_A",
        "contact_ca_threshold_A",
        "min_receptor_contact_retention_fraction",
        "max_refine_fa_rep_per_residue",
        "max_refine_backbone_strain_per_residue",
    }
    assert_keys(
        topology_calibration,
        allowed=topology_calibration_keys,
        required=topology_calibration_keys,
        path="discovery.qualification.topology_calibration",
    )
    target_keys = {
        "reference_receptor",
        "pilot_receptor",
        "qualification_status",
        "qualification_report",
        "qualification_report_sha256",
        "contact_jaccard_distance",
        "position_distance_A",
        "min_seed_support",
        "min_receptor_support",
    }
    parsed_targets: list[DiscoveryTargetSettings] = []
    for target_id, raw_target in sorted(targets.items()):
        target = document(raw_target, name=f"discovery.targets.{target_id}")
        target_path = f"discovery.targets.{target_id}"
        assert_keys(
            target,
            allowed=target_keys,
            required=target_keys,
            path=target_path,
        )
        report = _optional_text(target, "qualification_report", path=target_path)
        parsed_targets.append(
            DiscoveryTargetSettings(
                target_id=target_id,
                reference_receptor=string(
                    target, "reference_receptor", path=target_path
                ),
                pilot_receptor=string(target, "pilot_receptor", path=target_path),
                qualification_status=string(
                    target, "qualification_status", path=target_path
                ),
                qualification_report=(
                    resolved_path(report, config_dir=config_dir)
                    if report is not None
                    else None
                ),
                qualification_report_sha256=_optional_text(
                    target, "qualification_report_sha256", path=target_path
                ),
                analysis=SiteAnalysisSettings(
                    contact_jaccard_distance=_optional_float(
                        target, "contact_jaccard_distance", path=target_path
                    ),
                    position_distance_A=_optional_float(
                        target, "position_distance_A", path=target_path
                    ),
                    min_seed_support=_optional_integer(
                        target, "min_seed_support", path=target_path
                    ),
                    min_receptor_support=_optional_integer(
                        target, "min_receptor_support", path=target_path
                    ),
                ),
            )
        )
    topology_report = _optional_text(
        qualification,
        "topology_calibration_report",
        path="discovery.qualification",
    )
    return DiscoverySettings(
        method_id=_optional_text(section, "method_id", path="discovery"),
        adapter_id=_optional_text(section, "adapter_id", path="discovery"),
        seeds=_seeds(section["seeds"]),
        ensemble=ReceptorEnsembleSettings(
            min_receptors_per_target=integer(
                ensemble, "min_receptors_per_target", path="discovery.ensemble"
            ),
            allowed_structure_states=strings(
                ensemble, "allowed_structure_states", path="discovery.ensemble"
            ),
        ),
        cabsdock=CabsDockSettings(
            executable=resolved_path(
                string(cabsdock, "executable", path="discovery.cabsdock"),
                config_dir=config_dir,
            ),
            source_dir=resolved_path(
                string(cabsdock, "source_dir", path="discovery.cabsdock"),
                config_dir=config_dir,
            ),
            source_revision=string(
                cabsdock, "source_revision", path="discovery.cabsdock"
            ),
            patch_file=resolved_path(
                string(cabsdock, "patch_file", path="discovery.cabsdock"),
                config_dir=config_dir,
            ),
            seed_workers=integer(cabsdock, "seed_workers", path="discovery.cabsdock"),
            peptide_secondary_structure=string(
                cabsdock,
                "peptide_secondary_structure",
                path="discovery.cabsdock",
            ),
            mc_annealing=integer(cabsdock, "mc_annealing", path="discovery.cabsdock"),
            mc_cycles=integer(cabsdock, "mc_cycles", path="discovery.cabsdock"),
            mc_steps=integer(cabsdock, "mc_steps", path="discovery.cabsdock"),
            replicas=integer(cabsdock, "replicas", path="discovery.cabsdock"),
            replicas_dtemp=number(
                cabsdock, "replicas_dtemp", path="discovery.cabsdock"
            ),
            temperature_initial=number(
                cabsdock, "temperature_initial", path="discovery.cabsdock"
            ),
            temperature_final=number(
                cabsdock, "temperature_final", path="discovery.cabsdock"
            ),
            binding_interactions=number(
                cabsdock, "binding_interactions", path="discovery.cabsdock"
            ),
            protein_restraint_gap=integer(
                cabsdock, "protein_restraint_gap", path="discovery.cabsdock"
            ),
            protein_restraint_min_A=number(
                cabsdock, "protein_restraint_min_A", path="discovery.cabsdock"
            ),
            protein_restraint_max_A=number(
                cabsdock, "protein_restraint_max_A", path="discovery.cabsdock"
            ),
            filtering_count=integer(
                cabsdock, "filtering_count", path="discovery.cabsdock"
            ),
            clustering_medoids=integer(
                cabsdock, "clustering_medoids", path="discovery.cabsdock"
            ),
            clustering_iterations=integer(
                cabsdock, "clustering_iterations", path="discovery.cabsdock"
            ),
            trajectory_contact_ca_threshold_A=number(
                cabsdock,
                "trajectory_contact_ca_threshold_A",
                path="discovery.cabsdock",
            ),
            max_disulfide_ca_distance_A=number(
                cabsdock,
                "max_disulfide_ca_distance_A",
                path="discovery.cabsdock",
            ),
            min_models_for_selection=integer(
                cabsdock,
                "min_models_for_selection",
                path="discovery.cabsdock",
            ),
            selection_contact_jaccard_distance=number(
                cabsdock,
                "selection_contact_jaccard_distance",
                path="discovery.cabsdock",
            ),
            selection_position_distance_A=number(
                cabsdock,
                "selection_position_distance_A",
                path="discovery.cabsdock",
            ),
            pose_clustering_rmsd_A=number(
                cabsdock,
                "pose_clustering_rmsd_A",
                path="discovery.cabsdock",
            ),
        ),
        qualification=DiscoveryQualificationSettings(
            seeds=_seeds(qualification["seeds"], name="discovery.qualification.seeds"),
            control_bound_state_id=string(
                qualification,
                "control_bound_state_id",
                path="discovery.qualification",
            ),
            control_target_id=string(
                qualification,
                "control_target_id",
                path="discovery.qualification",
            ),
            control_secondary_structure=string(
                qualification,
                "control_secondary_structure",
                path="discovery.qualification",
            ),
            max_native_ligand_rmsd_A=number(
                qualification,
                "max_native_ligand_rmsd_A",
                path="discovery.qualification",
            ),
            min_native_receptor_contact_fraction=number(
                qualification,
                "min_native_receptor_contact_fraction",
                path="discovery.qualification",
            ),
            min_successful_control_seeds=integer(
                qualification,
                "min_successful_control_seeds",
                path="discovery.qualification",
            ),
            topology_calibration_status=string(
                qualification,
                "topology_calibration_status",
                path="discovery.qualification",
            ),
            topology_calibration_report=(
                resolved_path(topology_report, config_dir=config_dir)
                if topology_report is not None
                else None
            ),
            topology_calibration_report_sha256=_optional_text(
                qualification,
                "topology_calibration_report_sha256",
                path="discovery.qualification",
            ),
            topology_calibration=TopologyCalibrationSettings(
                candidate_ca_thresholds_A=_number_array(
                    topology_calibration["candidate_ca_thresholds_A"],
                    name=(
                        "discovery.qualification.topology_calibration."
                        "candidate_ca_thresholds_A"
                    ),
                ),
                above_threshold_comparator_upper_bound_A=number(
                    topology_calibration,
                    "above_threshold_comparator_upper_bound_A",
                    path="discovery.qualification.topology_calibration",
                ),
                models_per_stratum=integer(
                    topology_calibration,
                    "models_per_stratum",
                    path="discovery.qualification.topology_calibration",
                ),
                min_success_fraction_per_stratum=number(
                    topology_calibration,
                    "min_success_fraction_per_stratum",
                    path="discovery.qualification.topology_calibration",
                ),
                min_successful_seeds_per_stratum=integer(
                    topology_calibration,
                    "min_successful_seeds_per_stratum",
                    path="discovery.qualification.topology_calibration",
                ),
                min_interchain_heavy_atom_distance_A=number(
                    topology_calibration,
                    "min_interchain_heavy_atom_distance_A",
                    path="discovery.qualification.topology_calibration",
                ),
                max_peptide_internal_ca_rmsd_A=number(
                    topology_calibration,
                    "max_peptide_internal_ca_rmsd_A",
                    path="discovery.qualification.topology_calibration",
                ),
                max_ligand_centroid_displacement_A=number(
                    topology_calibration,
                    "max_ligand_centroid_displacement_A",
                    path="discovery.qualification.topology_calibration",
                ),
                contact_ca_threshold_A=number(
                    topology_calibration,
                    "contact_ca_threshold_A",
                    path="discovery.qualification.topology_calibration",
                ),
                min_receptor_contact_retention_fraction=number(
                    topology_calibration,
                    "min_receptor_contact_retention_fraction",
                    path="discovery.qualification.topology_calibration",
                ),
                max_refine_fa_rep_per_residue=number(
                    topology_calibration,
                    "max_refine_fa_rep_per_residue",
                    path="discovery.qualification.topology_calibration",
                ),
                max_refine_backbone_strain_per_residue=number(
                    topology_calibration,
                    "max_refine_backbone_strain_per_residue",
                    path="discovery.qualification.topology_calibration",
                ),
            ),
        ),
        targets=tuple(parsed_targets),
    )
