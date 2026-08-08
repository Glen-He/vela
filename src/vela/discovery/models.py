"""阶段二配置、任务和 site 证据领域模型。"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError

UNRESOLVED = "unresolved"
QUALIFICATION_STATUSES = frozenset(
    {UNRESOLVED, "unqualified", "transferability_unresolved", "qualified"}
)
SITE_EVIDENCE_TIERS = (
    "ensemble_consensus",
    "conformation_specific",
    "insufficient_evidence",
)
SHA256_PATTERN = re.compile(r"^[0-9a-f]{64}$")


class DiscoveryError(VelaError):
    """阶段二输入、规划或证据分析失败。"""


@dataclass(frozen=True, slots=True)
class SiteAnalysisSettings:
    """在独立控制中校准后才可用于 production 的 site 规则。"""

    contact_jaccard_distance: float | None
    position_distance_A: float | None
    min_seed_support: int | None
    min_receptor_support: int | None
    min_conformation_specific_seed_support: int | None
    ensemble_candidate_budget: int | None
    conformation_specific_candidate_budget: int | None

    def __post_init__(self) -> None:
        if (
            self.contact_jaccard_distance is not None
            and not 0.0 < self.contact_jaccard_distance <= 1.0
        ):
            raise DiscoveryError(
                "contact_jaccard_distance must be in (0, 1] or unresolved"
            )
        if self.position_distance_A is not None and self.position_distance_A <= 0:
            raise DiscoveryError("position_distance_A must be positive or unresolved")
        for name, value in (
            ("min_seed_support", self.min_seed_support),
            ("min_receptor_support", self.min_receptor_support),
            (
                "min_conformation_specific_seed_support",
                self.min_conformation_specific_seed_support,
            ),
            ("ensemble_candidate_budget", self.ensemble_candidate_budget),
            (
                "conformation_specific_candidate_budget",
                self.conformation_specific_candidate_budget,
            ),
        ):
            if value is not None and value < 1:
                raise DiscoveryError(f"{name} must be at least 1 or unresolved")
        if (
            self.min_seed_support is not None
            and self.min_conformation_specific_seed_support is not None
            and self.min_conformation_specific_seed_support < self.min_seed_support
        ):
            raise DiscoveryError(
                "min_conformation_specific_seed_support must be at least "
                "min_seed_support"
            )

    @property
    def complete(self) -> bool:
        return all(
            value is not None
            for value in (
                self.contact_jaccard_distance,
                self.position_distance_A,
                self.min_seed_support,
                self.min_receptor_support,
                self.min_conformation_specific_seed_support,
                self.ensemble_candidate_budget,
                self.conformation_specific_candidate_budget,
            )
        )


@dataclass(frozen=True, slots=True)
class DiscoveryTargetSettings:
    """一个阶段二靶标的坐标系、先导受体和独立放行记录。"""

    target_id: str
    reference_receptor: str
    pilot_receptor: str
    qualification_status: str
    qualification_report: Path | None
    qualification_report_sha256: str | None
    analysis: SiteAnalysisSettings

    def __post_init__(self) -> None:
        for name, value in (
            ("target_id", self.target_id),
            ("reference_receptor", self.reference_receptor),
            ("pilot_receptor", self.pilot_receptor),
        ):
            if not value.strip():
                raise DiscoveryError(f"{name} must not be empty")
        if self.qualification_status not in QUALIFICATION_STATUSES:
            raise DiscoveryError(
                "qualification_status must be unresolved, unqualified, "
                "transferability_unresolved, or qualified"
            )
        if (
            self.qualification_report_sha256 is not None
            and not SHA256_PATTERN.fullmatch(self.qualification_report_sha256)
        ):
            raise DiscoveryError(
                "qualification_report_sha256 must be a lowercase SHA-256 or unresolved"
            )
        report_complete = (
            self.qualification_report is not None
            and self.qualification_report_sha256 is not None
        )
        if self.qualification_status == UNRESOLVED and (
            self.qualification_report is not None
            or self.qualification_report_sha256 is not None
        ):
            raise DiscoveryError(
                "unresolved target qualification must not declare a report"
            )
        if self.qualification_status != UNRESOLVED and not report_complete:
            raise DiscoveryError(
                "resolved target qualification requires a report and SHA-256"
            )

    @property
    def qualified(self) -> bool:
        """返回该靶标是否已经具备完整且可核验的放行声明。"""
        return (
            self.qualification_status == "qualified"
            and self.qualification_report is not None
            and self.qualification_report_sha256 is not None
            and self.analysis.complete
        )


@dataclass(frozen=True, slots=True)
class ReceptorEnsembleSettings:
    """阶段二正式主受体集合的显式放行策略。"""

    min_receptors_per_target: int
    allowed_structure_states: tuple[str, ...]

    def __post_init__(self) -> None:
        if not self.allowed_structure_states or any(
            not value.strip() for value in self.allowed_structure_states
        ):
            raise DiscoveryError(
                "allowed_structure_states must contain non-empty values"
            )
        if len(self.allowed_structure_states) != len(
            set(self.allowed_structure_states)
        ):
            raise DiscoveryError("allowed_structure_states values must be unique")
        if self.min_receptors_per_target < 1:
            raise DiscoveryError("min_receptors_per_target must be at least 1")


@dataclass(frozen=True, slots=True)
class TopologyCalibrationSettings:
    """CABS 粗粒化闭环启发式的独立全原子校准设计。"""

    candidate_ca_thresholds_A: tuple[float, ...]
    above_threshold_comparator_upper_bound_A: float
    models_per_stratum: int
    min_success_fraction_per_stratum: float
    min_successful_seeds_per_stratum: int
    min_interchain_heavy_atom_distance_A: float
    min_nonlocal_peptide_heavy_atom_distance_A: float
    max_peptide_internal_ca_rmsd_A: float
    max_ligand_centroid_displacement_A: float
    contact_ca_threshold_A: float
    min_receptor_contact_retention_fraction: float
    site_coordinate_constraint_flat_width_A: float
    site_coordinate_constraint_sd_A: float
    site_coordinate_constraint_weight: float

    def __post_init__(self) -> None:
        if not self.candidate_ca_thresholds_A:
            raise DiscoveryError(
                "topology calibration requires at least one candidate CA threshold"
            )
        if any(value <= 0 for value in self.candidate_ca_thresholds_A) or any(
            left >= right
            for left, right in zip(
                self.candidate_ca_thresholds_A,
                self.candidate_ca_thresholds_A[1:],
                strict=False,
            )
        ):
            raise DiscoveryError(
                "topology calibration candidate CA thresholds must increase"
            )
        if (
            self.above_threshold_comparator_upper_bound_A
            <= (self.candidate_ca_thresholds_A[-1])
        ):
            raise DiscoveryError(
                "topology calibration comparator must exceed every candidate threshold"
            )
        if self.models_per_stratum < 1:
            raise DiscoveryError(
                "topology calibration models_per_stratum must be positive"
            )
        if not 0.0 < self.min_success_fraction_per_stratum <= 1.0:
            raise DiscoveryError(
                "topology calibration min_success_fraction_per_stratum must be in (0, 1]"
            )
        if not 1 <= self.min_successful_seeds_per_stratum <= (self.models_per_stratum):
            raise DiscoveryError(
                "topology calibration seed support must fit the stratum budget"
            )
        if self.min_interchain_heavy_atom_distance_A <= 0:
            raise DiscoveryError(
                "topology calibration interchain distance must be positive"
            )
        if self.min_nonlocal_peptide_heavy_atom_distance_A <= 0:
            raise DiscoveryError(
                "topology calibration nonlocal peptide distance must be positive"
            )
        if self.max_peptide_internal_ca_rmsd_A <= 0:
            raise DiscoveryError(
                "topology calibration peptide internal RMSD must be positive"
            )
        if self.max_ligand_centroid_displacement_A <= 0:
            raise DiscoveryError(
                "topology calibration ligand centroid displacement must be positive"
            )
        if self.contact_ca_threshold_A <= 0:
            raise DiscoveryError(
                "topology calibration CA contact threshold must be positive"
            )
        if not 0.0 < self.min_receptor_contact_retention_fraction <= 1.0:
            raise DiscoveryError(
                "topology calibration contact retention must be in (0, 1]"
            )
        if self.site_coordinate_constraint_flat_width_A < 0:
            raise DiscoveryError(
                "topology calibration site constraint flat width must not be negative"
            )
        if self.site_coordinate_constraint_sd_A <= 0:
            raise DiscoveryError(
                "topology calibration site constraint standard deviation must be positive"
            )
        if self.site_coordinate_constraint_weight <= 0:
            raise DiscoveryError(
                "topology calibration site constraint weight must be positive"
            )

    @property
    def strata_upper_bounds_A(self) -> tuple[float, ...]:
        """返回候选阈值层及一个明确不参与阈值选择的外部对照层。"""
        return (
            *self.candidate_ca_thresholds_A,
            self.above_threshold_comparator_upper_bound_A,
        )


@dataclass(frozen=True, slots=True)
class DiscoveryQualificationSettings:
    """阶段二靶标域 site 回收控制的冻结资格规则。"""

    seeds: tuple[int, ...]
    control_bound_state_id: str
    control_receptor_ids: tuple[str, ...]
    benchmark_receptor_id: str
    control_target_id: str
    control_secondary_structure: str
    receptor_site_diagnostic_budget: int
    max_native_ligand_rmsd_A: float
    max_native_site_centroid_distance_A: float
    min_native_receptor_contact_fraction: float
    min_native_site_seed_support: int
    min_native_receptor_support: int
    topology_calibration_status: str
    topology_calibration_report: Path | None
    topology_calibration_report_sha256: str | None
    topology_calibration: TopologyCalibrationSettings

    def __post_init__(self) -> None:
        if not self.seeds or len(self.seeds) != len(set(self.seeds)):
            raise DiscoveryError("qualification seeds must be non-empty and unique")
        if any(seed < 0 for seed in self.seeds):
            raise DiscoveryError("qualification seeds must not be negative")
        for name, value in (
            ("control_bound_state_id", self.control_bound_state_id),
            ("benchmark_receptor_id", self.benchmark_receptor_id),
            ("control_target_id", self.control_target_id),
        ):
            if not value.strip():
                raise DiscoveryError(f"{name} must not be empty")
        if not self.control_receptor_ids or any(
            not receptor_id.strip() for receptor_id in self.control_receptor_ids
        ):
            raise DiscoveryError("control_receptor_ids must contain receptor IDs")
        if len(self.control_receptor_ids) != len(set(self.control_receptor_ids)):
            raise DiscoveryError("control_receptor_ids must be unique")
        if self.benchmark_receptor_id in self.control_receptor_ids:
            raise DiscoveryError(
                "benchmark_receptor_id must be separate from production-domain controls"
            )
        if not self.control_secondary_structure or set(
            self.control_secondary_structure
        ) - {"C", "H", "E", "T"}:
            raise DiscoveryError("control_secondary_structure is invalid")
        if self.max_native_ligand_rmsd_A <= 0:
            raise DiscoveryError("max_native_ligand_rmsd_A must be positive")
        if self.receptor_site_diagnostic_budget < 1:
            raise DiscoveryError("receptor_site_diagnostic_budget must be positive")
        if self.max_native_site_centroid_distance_A <= 0:
            raise DiscoveryError("max_native_site_centroid_distance_A must be positive")
        if not 0.0 < self.min_native_receptor_contact_fraction <= 1.0:
            raise DiscoveryError(
                "min_native_receptor_contact_fraction must be in (0, 1]"
            )
        if not 1 <= self.min_native_site_seed_support <= len(self.seeds):
            raise DiscoveryError(
                "min_native_site_seed_support must fit the qualification seed count"
            )
        if not 1 <= self.min_native_receptor_support <= len(self.control_receptor_ids):
            raise DiscoveryError(
                "min_native_receptor_support must fit the control receptor count"
            )
        if self.topology_calibration_status not in {
            UNRESOLVED,
            "unqualified",
            "qualified",
        }:
            raise DiscoveryError(
                "topology_calibration_status must be unresolved, unqualified, or qualified"
            )
        if (
            self.topology_calibration_report_sha256 is not None
            and not SHA256_PATTERN.fullmatch(self.topology_calibration_report_sha256)
        ):
            raise DiscoveryError(
                "topology_calibration_report_sha256 must be a lowercase SHA-256 or unresolved"
            )
        calibration_report_complete = (
            self.topology_calibration_report is not None
            and self.topology_calibration_report_sha256 is not None
        )
        if self.topology_calibration_status == UNRESOLVED and (
            self.topology_calibration_report is not None
            or self.topology_calibration_report_sha256 is not None
        ):
            raise DiscoveryError(
                "unresolved topology calibration must not declare a report"
            )
        if (
            self.topology_calibration_status != UNRESOLVED
            and not calibration_report_complete
        ):
            raise DiscoveryError(
                "resolved topology calibration requires a report and SHA-256"
            )

    @property
    def topology_calibrated(self) -> bool:
        """返回粗粒化拓扑筛选是否有独立全原子校准证据。"""
        return (
            self.topology_calibration_status == "qualified"
            and self.topology_calibration_report is not None
            and self.topology_calibration_report_sha256 is not None
        )


@dataclass(frozen=True, slots=True)
class CabsDockSettings:
    """CABS-dock 全表面粗粒化采样和任务质控参数。"""

    executable: Path
    source_dir: Path
    source_revision: str
    seed_workers: int
    peptide_secondary_structure: str
    mc_annealing: int
    mc_cycles: int
    mc_steps: int
    replicas: int
    replicas_dtemp: float
    temperature_initial: float
    temperature_final: float
    binding_interactions: float
    protein_restraint_gap: int
    protein_restraint_min_A: float
    protein_restraint_max_A: float
    filtering_count: int
    clustering_medoids: int
    clustering_iterations: int
    trajectory_contact_ca_threshold_A: float
    disulfide_ca_restraint_distance_A: float
    disulfide_ca_restraint_weight: float
    max_reconstructable_disulfide_ca_distance_A: float
    min_models_for_selection: int
    selection_contact_jaccard_distance: float
    selection_position_distance_A: float
    pose_clustering_rmsd_A: float
    max_sites_per_task: int
    max_pose_clusters_per_site: int

    def __post_init__(self) -> None:
        if not re.fullmatch(r"[0-9a-f]{40,64}", self.source_revision):
            raise DiscoveryError(
                "source_revision must be a 40-64 character lowercase Git object ID"
            )
        positive_integers = (
            ("mc_annealing", self.mc_annealing),
            ("mc_cycles", self.mc_cycles),
            ("mc_steps", self.mc_steps),
            ("replicas", self.replicas),
            ("protein_restraint_gap", self.protein_restraint_gap),
            ("filtering_count", self.filtering_count),
            ("clustering_medoids", self.clustering_medoids),
            ("clustering_iterations", self.clustering_iterations),
            ("seed_workers", self.seed_workers),
            ("min_models_for_selection", self.min_models_for_selection),
            ("max_sites_per_task", self.max_sites_per_task),
            ("max_pose_clusters_per_site", self.max_pose_clusters_per_site),
        )
        for name, value in positive_integers:
            if value < 1:
                raise DiscoveryError(f"{name} must be positive")
        positive_numbers = (
            ("replicas_dtemp", self.replicas_dtemp),
            ("temperature_initial", self.temperature_initial),
            ("temperature_final", self.temperature_final),
            ("binding_interactions", self.binding_interactions),
            ("protein_restraint_min_A", self.protein_restraint_min_A),
            ("protein_restraint_max_A", self.protein_restraint_max_A),
            (
                "trajectory_contact_ca_threshold_A",
                self.trajectory_contact_ca_threshold_A,
            ),
            (
                "disulfide_ca_restraint_distance_A",
                self.disulfide_ca_restraint_distance_A,
            ),
            (
                "max_reconstructable_disulfide_ca_distance_A",
                self.max_reconstructable_disulfide_ca_distance_A,
            ),
            ("selection_position_distance_A", self.selection_position_distance_A),
            ("pose_clustering_rmsd_A", self.pose_clustering_rmsd_A),
        )
        for name, value in positive_numbers:
            if value <= 0:
                raise DiscoveryError(f"{name} must be positive")
        if not 0.0 < self.disulfide_ca_restraint_weight <= 1.0:
            raise DiscoveryError("disulfide_ca_restraint_weight must be in (0, 1]")
        if self.temperature_initial < self.temperature_final:
            raise DiscoveryError(
                "temperature_initial must not be below temperature_final"
            )
        if self.protein_restraint_min_A >= self.protein_restraint_max_A:
            raise DiscoveryError(
                "protein_restraint_min_A must be below protein_restraint_max_A"
            )
        if self.filtering_count < self.replicas:
            raise DiscoveryError("filtering_count must be at least the replica count")
        if self.filtering_count % self.replicas != 0:
            raise DiscoveryError("filtering_count must be divisible by replicas")
        if self.clustering_medoids > self.filtering_count:
            raise DiscoveryError("clustering_medoids must not exceed filtering_count")
        if self.min_models_for_selection > self.filtering_count:
            raise DiscoveryError(
                "min_models_for_selection must not exceed filtering_count"
            )
        if not 0.0 < self.selection_contact_jaccard_distance <= 1.0:
            raise DiscoveryError("selection_contact_jaccard_distance must be in (0, 1]")
        if not self.peptide_secondary_structure:
            raise DiscoveryError("peptide_secondary_structure must not be empty")
        if set(self.peptide_secondary_structure) - {"C", "H", "E", "T"}:
            raise DiscoveryError(
                "peptide_secondary_structure may contain only C, H, E, or T"
            )


@dataclass(frozen=True, slots=True)
class DiscoverySettings:
    """全表面方法资格、随机重复和分析规则声明。"""

    method_id: str | None
    adapter_id: str | None
    seeds: tuple[int, ...]
    ensemble: ReceptorEnsembleSettings
    cabsdock: CabsDockSettings
    qualification: DiscoveryQualificationSettings
    targets: tuple[DiscoveryTargetSettings, ...]

    def __post_init__(self) -> None:
        if len(self.seeds) != len(set(self.seeds)):
            raise DiscoveryError("discovery seeds must be unique")
        if any(seed < 0 for seed in self.seeds):
            raise DiscoveryError("discovery seeds must not be negative")
        target_ids = tuple(item.target_id for item in self.targets)
        if not target_ids or len(target_ids) != len(set(target_ids)):
            raise DiscoveryError("discovery target IDs must be non-empty and unique")

    def target(self, target_id: str) -> DiscoveryTargetSettings:
        """返回一个明确配置的阶段二靶标。"""
        for target in self.targets:
            if target.target_id == target_id:
                return target
        raise DiscoveryError(f"discovery target is not configured: {target_id}")

    @property
    def config_complete(self) -> bool:
        return (
            self.method_id is not None
            and self.adapter_id is not None
            and bool(self.seeds)
            and all(target.qualified for target in self.targets)
        )


@dataclass(frozen=True, slots=True)
class DiscoveryTask:
    """一个受体、chemistry、method 和独立 seed 的不可变任务。"""

    task_id: str
    receptor_id: str
    target: str
    receptor_path: Path
    receptor_sha256: str
    chemistry_id: str
    method_id: str
    adapter_id: str
    seed: int
    evidence_category: str

    def __post_init__(self) -> None:
        for name, value in (
            ("task_id", self.task_id),
            ("receptor_id", self.receptor_id),
            ("target", self.target),
            ("chemistry_id", self.chemistry_id),
            ("method_id", self.method_id),
            ("adapter_id", self.adapter_id),
            ("evidence_category", self.evidence_category),
        ):
            if not value.strip():
                raise DiscoveryError(f"{name} must not be empty")
        if self.seed < 0:
            raise DiscoveryError("task seed must not be negative")
        if not SHA256_PATTERN.fullmatch(self.receptor_sha256):
            raise DiscoveryError("receptor_sha256 must be a lowercase SHA-256")
