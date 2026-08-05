"""阶段三配置、结合态资产和科学放行模型。"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError
from vela.core.run_identity import RUN_ID_PATTERN
from vela.discovery.models import SHA256_PATTERN
from vela.preparation.chemistry import (
    HISTIDINE_STATES,
    STANDARD_AMINO_ACIDS,
    DisulfideBond,
    HistidineState,
)

LOCAL_CONTROL_KINDS = frozenset({"standard_cyclic_peptide", "site_reference_only"})
SCORE_TERM_PATTERN = re.compile(r"^[A-Za-z][A-Za-z0-9_.-]*$")
DEVICE_PATTERN = re.compile(r"^(?:cpu|cuda(?::[0-9]+)?)$")
METHOD_QUALIFICATION_STATUSES = frozenset({"unresolved", "failed", "qualified"})


class ValidationError(VelaError):
    """阶段三输入、结构资产或局部精修合同失败。"""


@dataclass(frozen=True, slots=True)
class BoundStateDefinition:
    """实验结合态中的明确受体链、配体链及允许用途。"""

    state_id: str
    receptor_id: str
    ligand_id: str
    ligand_author_chain_id: str
    local_control_kind: str
    ligand_sequence: str | None
    ligand_n_terminus: str | None
    ligand_c_terminus: str | None
    disulfide_bonds: tuple[DisulfideBond, ...]
    histidines: tuple[HistidineState, ...]
    selection_reason: str

    def __post_init__(self) -> None:
        for name, value in (
            ("state_id", self.state_id),
            ("receptor_id", self.receptor_id),
        ):
            if not RUN_ID_PATTERN.fullmatch(value):
                raise ValidationError(f"{name} must be a safe identifier")
        for name, value in (
            ("ligand_id", self.ligand_id),
            ("ligand_author_chain_id", self.ligand_author_chain_id),
            ("selection_reason", self.selection_reason),
        ):
            if not value.strip():
                raise ValidationError(f"{name} must not be empty")
        if self.local_control_kind not in LOCAL_CONTROL_KINDS:
            raise ValidationError(
                f"unsupported local_control_kind: {self.local_control_kind}"
            )
        if self.local_control_kind == "standard_cyclic_peptide":
            if not self.ligand_sequence:
                raise ValidationError(
                    f"{self.state_id} standard peptide control requires a sequence"
                )
            invalid = set(self.ligand_sequence) - STANDARD_AMINO_ACIDS
            if invalid:
                raise ValidationError(
                    f"{self.state_id} ligand sequence contains non-standard residues"
                )
            if not self.disulfide_bonds:
                raise ValidationError(
                    f"{self.state_id} cyclic peptide control requires a disulfide"
                )
            if self.ligand_n_terminus not in {"NH3+", "Acetyl"}:
                raise ValidationError(
                    f"{self.state_id} standard peptide N terminus is unsupported"
                )
            if self.ligand_c_terminus not in {"CONH2", "COO-"}:
                raise ValidationError(
                    f"{self.state_id} standard peptide C terminus is unsupported"
                )
            seen_positions: set[int] = set()
            for bond in self.disulfide_bonds:
                if not 1 <= bond.first < bond.second <= len(self.ligand_sequence):
                    raise ValidationError(
                        f"{self.state_id} disulfide positions are outside the ligand sequence"
                    )
                for position in (bond.first, bond.second):
                    if self.ligand_sequence[position - 1] != "C":
                        raise ValidationError(
                            f"{self.state_id} disulfide position {position} is not cysteine"
                        )
                    if position in seen_positions:
                        raise ValidationError(
                            f"{self.state_id} disulfide position {position} is reused"
                        )
                    seen_positions.add(position)
            expected_histidines = {
                index
                for index, residue in enumerate(self.ligand_sequence, 1)
                if residue == "H"
            }
            declared_histidines = {item.position for item in self.histidines}
            if len(declared_histidines) != len(self.histidines):
                raise ValidationError(
                    f"{self.state_id} ligand histidine positions must be unique"
                )
            if declared_histidines != expected_histidines:
                raise ValidationError(
                    f"{self.state_id} ligand histidine states do not match the sequence"
                )
            if any(item.state not in HISTIDINE_STATES for item in self.histidines):
                raise ValidationError(
                    f"{self.state_id} ligand histidine state must be HID, HIE, or HIP"
                )
        elif any(
            (
                self.ligand_sequence is not None,
                self.ligand_n_terminus is not None,
                self.ligand_c_terminus is not None,
                bool(self.disulfide_bonds),
                bool(self.histidines),
            )
        ):
            raise ValidationError(
                f"{self.state_id} site-only reference must not define standard peptide chemistry"
            )


@dataclass(frozen=True, slots=True)
class LocalRecoveryControl:
    """一个由配置选择的实验肽局部扰动恢复验证。"""

    control_id: str
    bound_state_id: str
    prepack_seed: int
    seed_batches: tuple[tuple[int, ...], ...]
    random_translation_A: float
    random_rotation_degrees: float
    ranking_score: str
    recovery_rmsd_score: str
    top_clusters: int
    max_recovery_rmsd_A: float
    max_cluster_backbone_rmsd_A: float
    min_cluster_seed_support: int
    max_batch_pose_rmsd_A: float

    def __post_init__(self) -> None:
        for name, value in (
            ("control_id", self.control_id),
            ("bound_state_id", self.bound_state_id),
        ):
            if not RUN_ID_PATTERN.fullmatch(value):
                raise ValidationError(f"local control {name} must be a safe identifier")
        if self.prepack_seed < 0:
            raise ValidationError("local control prepack_seed must not be negative")
        if len(self.seed_batches) < 2 or any(not batch for batch in self.seed_batches):
            raise ValidationError(
                "local control requires at least two non-empty seed batches"
            )
        seeds = tuple(seed for batch in self.seed_batches for seed in batch)
        if len(seeds) != len(set(seeds)) or any(seed < 0 for seed in seeds):
            raise ValidationError(
                "local control batch seeds must be unique non-negative integers"
            )
        if self.random_translation_A < 0 or self.random_rotation_degrees < 0:
            raise ValidationError("local control perturbations must not be negative")
        if self.random_translation_A == 0 and self.random_rotation_degrees == 0:
            raise ValidationError("local recovery control requires a perturbation")
        for name, value in (
            ("ranking_score", self.ranking_score),
            ("recovery_rmsd_score", self.recovery_rmsd_score),
        ):
            if not SCORE_TERM_PATTERN.fullmatch(value):
                raise ValidationError(f"local control {name} is invalid")
        if self.top_clusters < 1:
            raise ValidationError("local control top_clusters must be positive")
        if self.max_recovery_rmsd_A <= 0:
            raise ValidationError("local control max_recovery_rmsd_A must be positive")
        if self.max_cluster_backbone_rmsd_A <= 0:
            raise ValidationError(
                "local control max_cluster_backbone_rmsd_A must be positive"
            )
        if self.max_batch_pose_rmsd_A <= 0:
            raise ValidationError(
                "local control max_batch_pose_rmsd_A must be positive"
            )
        if (
            any(
                self.min_cluster_seed_support > len(batch)
                for batch in self.seed_batches
            )
            or self.min_cluster_seed_support < 1
        ):
            raise ValidationError(
                "local control cluster seed-support requirement is invalid"
            )


@dataclass(frozen=True, slots=True)
class GuidedTemplate:
    """把目标肽逐位映射到标准实验肽骨架的显式声明。"""

    template_id: str
    bound_state_id: str
    ligand_positions: tuple[int, ...]
    selection_reason: str

    def __post_init__(self) -> None:
        for name, value in (
            ("template_id", self.template_id),
            ("bound_state_id", self.bound_state_id),
        ):
            if not RUN_ID_PATTERN.fullmatch(value):
                raise ValidationError(
                    f"guided template {name} must be a safe identifier"
                )
        if not self.ligand_positions:
            raise ValidationError("guided template ligand_positions must not be empty")
        if len(self.ligand_positions) != len(set(self.ligand_positions)) or any(
            position < 1 for position in self.ligand_positions
        ):
            raise ValidationError(
                "guided template ligand_positions must be unique positive integers"
            )
        if tuple(sorted(self.ligand_positions)) != self.ligand_positions:
            raise ValidationError(
                "guided template ligand_positions must follow ligand sequence order"
            )
        if not self.selection_reason.strip():
            raise ValidationError("guided template selection_reason must not be empty")


@dataclass(frozen=True, slots=True)
class EnvironmentReference:
    """一个实验全酶装配及其催化链、CK2β 链和其余催化链。"""

    reference_id: str
    receptor_id: str
    assembly_id: str
    beta_author_chain_ids: tuple[str, ...]
    other_catalytic_author_chain_ids: tuple[str, ...]
    evaluation_targets: tuple[str, ...]
    construct_note: str

    def __post_init__(self) -> None:
        for name, value in (
            ("reference_id", self.reference_id),
            ("receptor_id", self.receptor_id),
        ):
            if not RUN_ID_PATTERN.fullmatch(value):
                raise ValidationError(
                    f"environment reference {name} must be a safe identifier"
                )
        if not self.assembly_id.strip():
            raise ValidationError("environment reference assembly_id must not be empty")
        chain_ids = self.beta_author_chain_ids + self.other_catalytic_author_chain_ids
        if (
            not self.beta_author_chain_ids
            or not self.other_catalytic_author_chain_ids
            or any(not chain_id.strip() for chain_id in chain_ids)
            or len(chain_ids) != len(set(chain_ids))
        ):
            raise ValidationError(
                "environment reference chains must be non-empty and unique"
            )
        if (
            not self.evaluation_targets
            or len(self.evaluation_targets) != len(set(self.evaluation_targets))
            or any(
                target not in {"ck2_alpha", "ck2_alpha_prime"}
                for target in self.evaluation_targets
            )
        ):
            raise ValidationError(
                "environment reference evaluation_targets are invalid"
            )
        if not self.construct_note.strip():
            raise ValidationError(
                "environment reference construct_note must not be empty"
            )


@dataclass(frozen=True, slots=True)
class RosettaSettings:
    """本机 FlexPepDock 工具身份和实际采样参数。"""

    executable: Path
    scripts_executable: Path
    database: Path
    version_file: Path
    expected_version: str
    parallel_tasks: int
    decoys_per_seed: int
    score_function: str
    lowres_preoptimize: bool

    def __post_init__(self) -> None:
        if not self.expected_version.strip():
            raise ValidationError("Rosetta expected_version must not be empty")
        if self.parallel_tasks < 1:
            raise ValidationError("Rosetta parallel_tasks must be positive")
        if self.decoys_per_seed < 1:
            raise ValidationError("Rosetta decoys_per_seed must be positive")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", self.score_function):
            raise ValidationError("Rosetta score_function contains unsafe characters")


@dataclass(frozen=True, slots=True)
class Cg2AllSettings:
    """CABS CA/SC pose 到全原子结构的重建工具与质量合同。"""

    executable: Path
    package_metadata: Path
    checkpoint: Path
    checkpoint_sha256: str
    expected_version: str
    representation: str
    receptor_histidine_state: str
    device: str
    processes: int
    batch_size: int
    chain_break_cutoff_A: float
    max_ca_rmsd_A: float

    def __post_init__(self) -> None:
        if not self.expected_version.strip():
            raise ValidationError("cg2all expected_version must not be empty")
        if self.representation != "CalphaSCModel":
            raise ValidationError("Stage 2 CA/SC handoff requires cg2all CalphaSCModel")
        if self.receptor_histidine_state not in {"HID", "HIE"}:
            raise ValidationError("cg2all receptor_histidine_state must be HID or HIE")
        if not DEVICE_PATTERN.fullmatch(self.device):
            raise ValidationError("cg2all device must be cpu or cuda[:index]")
        if self.processes < 1 or self.batch_size < 1:
            raise ValidationError("cg2all processes and batch_size must be positive")
        if self.chain_break_cutoff_A <= 0 or self.max_ca_rmsd_A <= 0:
            raise ValidationError("cg2all distance thresholds must be positive")
        if not SHA256_PATTERN.fullmatch(self.checkpoint_sha256):
            raise ValidationError(
                "cg2all checkpoint_sha256 must be a lowercase SHA-256"
            )


@dataclass(frozen=True, slots=True)
class CandidateHandoffSettings:
    """受支持 candidate site 到局部精修起点的选择预算。"""

    poses_per_receptor_site: int
    chemistry_seed: int

    def __post_init__(self) -> None:
        if self.poses_per_receptor_site < 1:
            raise ValidationError("poses_per_receptor_site must be positive")
        if self.chemistry_seed < 0:
            raise ValidationError("handoff chemistry_seed must not be negative")


@dataclass(frozen=True, slots=True)
class CandidateRefinementSettings:
    """正式 candidate 局部精修的起点处理和排名合同。"""

    prepack_seed: int
    random_translation_A: float
    random_rotation_degrees: float
    ranking_score: str

    def __post_init__(self) -> None:
        if self.prepack_seed < 0:
            raise ValidationError("refinement prepack_seed must not be negative")
        if self.random_translation_A < 0 or self.random_rotation_degrees < 0:
            raise ValidationError("refinement perturbations must not be negative")
        if not SCORE_TERM_PATTERN.fullmatch(self.ranking_score):
            raise ValidationError("refinement ranking_score is invalid")


@dataclass(frozen=True, slots=True)
class ValidationAnalysisSettings:
    """必须经独立控制校准后才能用于正式候选筛选的规则。"""

    min_interface_contact_pairs: int | None
    min_interface_receptor_residues: int | None
    max_receptor_ca_rmsd_A: float | None
    min_start_contact_overlap: float | None
    max_start_site_displacement_A: float | None
    max_cluster_backbone_rmsd_A: float | None
    min_heavy_atom_distance_A: float | None
    min_refinement_seed_support: int | None
    min_refinement_start_support: int | None

    def __post_init__(self) -> None:
        for name, value in (
            ("min_interface_contact_pairs", self.min_interface_contact_pairs),
            (
                "min_interface_receptor_residues",
                self.min_interface_receptor_residues,
            ),
        ):
            if value is not None and value < 1:
                raise ValidationError(f"{name} must be positive or unresolved")
        for name, value in (
            ("max_receptor_ca_rmsd_A", self.max_receptor_ca_rmsd_A),
            ("max_start_site_displacement_A", self.max_start_site_displacement_A),
            ("max_cluster_backbone_rmsd_A", self.max_cluster_backbone_rmsd_A),
            ("min_heavy_atom_distance_A", self.min_heavy_atom_distance_A),
        ):
            if value is not None and value <= 0:
                raise ValidationError(f"{name} must be positive or unresolved")
        for name, value in (
            ("min_refinement_seed_support", self.min_refinement_seed_support),
            ("min_refinement_start_support", self.min_refinement_start_support),
        ):
            if value is not None and value < 1:
                raise ValidationError(f"{name} must be positive or unresolved")
        if (
            self.min_start_contact_overlap is not None
            and not 0.0 <= self.min_start_contact_overlap <= 1.0
        ):
            raise ValidationError(
                "min_start_contact_overlap must be in [0, 1] or unresolved"
            )

    @property
    def complete(self) -> bool:
        return all(
            value is not None
            for value in (
                self.min_interface_contact_pairs,
                self.min_interface_receptor_residues,
                self.max_receptor_ca_rmsd_A,
                self.min_start_contact_overlap,
                self.max_start_site_displacement_A,
                self.max_cluster_backbone_rmsd_A,
                self.min_heavy_atom_distance_A,
                self.min_refinement_seed_support,
                self.min_refinement_start_support,
            )
        )


@dataclass(frozen=True, slots=True)
class ValidationSettings:
    """阶段三方法资格、随机重复、结构资产和分析规则。"""

    method_id: str
    qualification_status: str
    qualification_report: Path | None
    qualification_report_sha256: str | None
    seeds: tuple[int, ...]
    min_disulfide_sg_A: float
    max_disulfide_sg_A: float
    interface_contact_A: float
    rosetta: RosettaSettings
    cg2all: Cg2AllSettings
    handoff: CandidateHandoffSettings
    refinement: CandidateRefinementSettings
    analysis: ValidationAnalysisSettings
    bound_states: tuple[BoundStateDefinition, ...]
    local_controls: tuple[LocalRecoveryControl, ...]
    guided_templates: tuple[GuidedTemplate, ...]
    environment_references: tuple[EnvironmentReference, ...]

    def __post_init__(self) -> None:
        if not self.method_id.strip():
            raise ValidationError("validation method_id must not be empty")
        if self.qualification_status not in METHOD_QUALIFICATION_STATUSES:
            raise ValidationError(
                "validation qualification_status must be unresolved, failed, or qualified"
            )
        if (
            self.qualification_report_sha256 is not None
            and not SHA256_PATTERN.fullmatch(self.qualification_report_sha256)
        ):
            raise ValidationError(
                "validation qualification report hash must be a lowercase SHA-256"
            )
        report_complete = (
            self.qualification_report is not None
            and self.qualification_report_sha256 is not None
        )
        if self.qualification_status == "unresolved" and (
            self.qualification_report is not None
            or self.qualification_report_sha256 is not None
        ):
            raise ValidationError(
                "unresolved validation qualification must not declare a report"
            )
        if self.qualification_status != "unresolved" and not report_complete:
            raise ValidationError(
                "resolved validation qualification requires a report and SHA-256"
            )
        if len(self.seeds) != len(set(self.seeds)) or any(
            seed < 0 for seed in self.seeds
        ):
            raise ValidationError(
                "validation seeds must be unique non-negative integers"
            )
        if self.min_disulfide_sg_A <= 0:
            raise ValidationError("min_disulfide_sg_A must be positive")
        if self.min_disulfide_sg_A >= self.max_disulfide_sg_A:
            raise ValidationError("min_disulfide_sg_A must be below max_disulfide_sg_A")
        if self.interface_contact_A <= 0:
            raise ValidationError("interface_contact_A must be positive")
        identifiers = [item.state_id for item in self.bound_states]
        if len(identifiers) != len(set(identifiers)):
            raise ValidationError("bound-state IDs must be unique")
        controls = [item.control_id for item in self.local_controls]
        if not controls:
            raise ValidationError("at least one local recovery control is required")
        if len(controls) != len(set(controls)):
            raise ValidationError("local recovery control IDs must be unique")
        bound_state_by_id = {item.state_id: item for item in self.bound_states}
        for control in self.local_controls:
            state = bound_state_by_id.get(control.bound_state_id)
            if state is None:
                raise ValidationError(
                    f"unknown control bound state: {control.bound_state_id}"
                )
            if state.local_control_kind != "standard_cyclic_peptide":
                raise ValidationError(
                    f"{control.control_id} requires a standard cyclic-peptide bound state"
                )
        template_ids = [item.template_id for item in self.guided_templates]
        if len(template_ids) != len(set(template_ids)):
            raise ValidationError("guided template IDs must be unique")
        for template in self.guided_templates:
            state = bound_state_by_id.get(template.bound_state_id)
            if state is None:
                raise ValidationError(
                    f"unknown guided-template bound state: {template.bound_state_id}"
                )
            if state.local_control_kind != "standard_cyclic_peptide":
                raise ValidationError(
                    f"{template.template_id} requires a standard cyclic-peptide bound state"
                )
            ligand_length = len(state.ligand_sequence or "")
            if max(template.ligand_positions) > ligand_length:
                raise ValidationError(
                    f"{template.template_id} maps outside the experimental ligand"
                )
        reference_ids = [item.reference_id for item in self.environment_references]
        if not reference_ids:
            raise ValidationError("at least one full-enzyme reference is required")
        if len(reference_ids) != len(set(reference_ids)):
            raise ValidationError("environment reference IDs must be unique")

    @property
    def config_complete(self) -> bool:
        return (
            self.qualification_status == "qualified"
            and self.qualification_report is not None
            and self.qualification_report_sha256 is not None
            and bool(self.seeds)
            and self.analysis.complete
        )
