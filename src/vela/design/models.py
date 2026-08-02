"""阶段四配置、序列候选和多状态筛查领域模型。"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError
from vela.core.run_identity import RUN_ID_PATTERN
from vela.discovery.models import SHA256_PATTERN
from vela.preparation.chemistry import STANDARD_AMINO_ACIDS
from vela.validation.models import METHOD_QUALIFICATION_STATUSES, SCORE_TERM_PATTERN

DESIGN_OBJECTIVES = frozenset({"single_supported_target"})
DESIGN_ROUNDS = frozenset({"single", "combination", "iteration"})
CANDIDATE_EDIT_TYPES = frozenset({"origin", "add", "revert", "swap"})


class DesignError(VelaError):
    """阶段四输入、计划、执行或证据汇总失败。"""


@dataclass(frozen=True, slots=True)
class SequencePolicy:
    """配体变体允许改变的范围及统一化学解释。"""

    mutable_positions: tuple[int, ...]
    allowed_amino_acids: str
    candidate_histidine_state: str

    def __post_init__(self) -> None:
        if (
            not self.mutable_positions
            or len(self.mutable_positions) != len(set(self.mutable_positions))
            or tuple(sorted(self.mutable_positions)) != self.mutable_positions
            or any(position < 1 for position in self.mutable_positions)
        ):
            raise DesignError(
                "mutable_positions must be unique positive positions in sequence order"
            )
        if (
            not self.allowed_amino_acids
            or len(self.allowed_amino_acids) != len(set(self.allowed_amino_acids))
            or tuple(sorted(self.allowed_amino_acids))
            != tuple(self.allowed_amino_acids)
            or set(self.allowed_amino_acids) - STANDARD_AMINO_ACIDS
        ):
            raise DesignError(
                "allowed_amino_acids must be sorted unique standard one-letter codes"
            )
        if "C" in self.allowed_amino_acids:
            raise DesignError(
                "allowed_amino_acids must exclude cysteine; the disulfide pair is fixed"
            )
        if self.candidate_histidine_state != "HIE":
            raise DesignError(
                "candidate_histidine_state must be HIE for the current Rosetta chemistry contract"
            )


@dataclass(frozen=True, slots=True)
class InterfaceScreenSettings:
    """固定骨架成对界面筛查的运行参数。"""

    parallel_tasks: int
    neighbor_distance_A: float
    score_function: str
    ranking_score: str
    pack_separated: bool

    def __post_init__(self) -> None:
        if self.parallel_tasks < 1:
            raise DesignError("design screen parallel_tasks must be positive")
        if self.neighbor_distance_A <= 0:
            raise DesignError("design screen neighbor_distance_A must be positive")
        if not re.fullmatch(r"[A-Za-z0-9_.-]+", self.score_function):
            raise DesignError("design screen score_function contains unsafe characters")
        if not SCORE_TERM_PATTERN.fullmatch(self.ranking_score):
            raise DesignError("design screen ranking_score is invalid")


@dataclass(frozen=True, slots=True)
class CombinationSettings:
    """由单点证据提出有限多突变组合的资源边界。"""

    min_mutations: int
    max_mutations: int
    max_options_per_position: int
    max_candidates: int

    def __post_init__(self) -> None:
        if self.min_mutations < 2:
            raise DesignError("combination min_mutations must be at least 2")
        if self.max_mutations < self.min_mutations:
            raise DesignError(
                "combination max_mutations must not be below min_mutations"
            )
        if self.max_options_per_position < 1 or self.max_candidates < 1:
            raise DesignError("combination resource limits must be positive")


@dataclass(frozen=True, slots=True)
class IterationSettings:
    """父序列一步邻域搜索的代次和资源边界。"""

    max_parents: int
    max_total_mutations: int
    max_candidates: int

    def __post_init__(self) -> None:
        if self.max_parents < 1 or self.max_candidates < 1:
            raise DesignError("iteration resource limits must be positive")
        if self.max_candidates < self.max_parents:
            raise DesignError(
                "iteration max_candidates must cover every selected parent"
            )
        if self.max_total_mutations < 1:
            raise DesignError("iteration max_total_mutations must be positive")


@dataclass(frozen=True, slots=True)
class DesignAnalysisSettings:
    """必须通过独立校准冻结的候选判断门槛。"""

    calibrated: bool
    max_positive_median_delta: float | None
    max_positive_worst_delta: float | None

    @property
    def complete(self) -> bool:
        """返回门槛是否足以形成正式候选判定。"""
        return (
            self.calibrated
            and self.max_positive_median_delta is not None
            and self.max_positive_worst_delta is not None
        )


@dataclass(frozen=True, slots=True)
class FinalistSettings:
    """候选柔性复核、结构通过率和阶段五资源边界。"""

    parallel_tasks: int
    max_candidates: int
    max_md_candidates: int
    seeds: tuple[int, ...]
    ranking_score: str
    interface_score: str
    calibrated: bool
    min_passed_decoy_fraction: float | None
    min_successful_seeds: int | None
    max_positive_median_ranking_delta: float | None
    max_positive_worst_ranking_delta: float | None
    max_positive_median_interface_delta: float | None
    max_positive_worst_interface_delta: float | None

    def __post_init__(self) -> None:
        if self.parallel_tasks < 1:
            raise DesignError("design finalist parallel_tasks must be positive")
        if self.max_candidates < 1 or self.max_md_candidates < 1:
            raise DesignError("design finalist resource limits must be positive")
        if self.max_md_candidates > self.max_candidates:
            raise DesignError(
                "design finalist max_md_candidates must not exceed max_candidates"
            )
        if len(self.seeds) != len(set(self.seeds)) or any(
            seed < 0 for seed in self.seeds
        ):
            raise DesignError(
                "design finalist seeds must be unique non-negative integers"
            )
        for name, value in (
            ("ranking_score", self.ranking_score),
            ("interface_score", self.interface_score),
        ):
            if not SCORE_TERM_PATTERN.fullmatch(value):
                raise DesignError(f"design finalist {name} is invalid")
        if (
            self.min_passed_decoy_fraction is not None
            and not 0.0 <= self.min_passed_decoy_fraction <= 1.0
        ):
            raise DesignError(
                "design finalist min_passed_decoy_fraction must be in [0, 1]"
            )
        if self.min_successful_seeds is not None and self.min_successful_seeds < 1:
            raise DesignError("design finalist min_successful_seeds must be positive")
        if (
            self.min_successful_seeds is not None
            and self.seeds
            and self.min_successful_seeds > len(self.seeds)
        ):
            raise DesignError(
                "design finalist min_successful_seeds exceeds available seeds"
            )

    @property
    def complete(self) -> bool:
        """返回单亚型柔性复核合同是否完整。"""
        return (
            self.calibrated
            and bool(self.seeds)
            and self.min_passed_decoy_fraction is not None
            and self.min_successful_seeds is not None
            and self.max_positive_median_ranking_delta is not None
            and self.max_positive_worst_ranking_delta is not None
            and self.max_positive_median_interface_delta is not None
            and self.max_positive_worst_interface_delta is not None
        )


@dataclass(frozen=True, slots=True)
class DesignSettings:
    """阶段四方法、目标、序列空间和分析策略。"""

    method_id: str | None
    qualification_status: str
    qualification_report: Path | None
    qualification_report_sha256: str | None
    objective: str | None
    seeds: tuple[int, ...]
    sequence: SequencePolicy
    screen: InterfaceScreenSettings
    combination: CombinationSettings
    iteration: IterationSettings
    analysis: DesignAnalysisSettings
    finalists: FinalistSettings

    def __post_init__(self) -> None:
        if self.qualification_status not in METHOD_QUALIFICATION_STATUSES:
            raise DesignError(
                "design qualification_status must be unresolved, failed, or qualified"
            )
        if self.objective is not None and self.objective not in DESIGN_OBJECTIVES:
            raise DesignError(f"unsupported design objective: {self.objective}")
        if len(self.seeds) != len(set(self.seeds)) or any(
            seed < 0 for seed in self.seeds
        ):
            raise DesignError("design seeds must be unique non-negative integers")
        if (
            self.qualification_report_sha256 is not None
            and not SHA256_PATTERN.fullmatch(self.qualification_report_sha256)
        ):
            raise DesignError(
                "design qualification_report_sha256 must be a lowercase SHA-256 or unresolved"
            )
        report_complete = (
            self.qualification_report is not None
            and self.qualification_report_sha256 is not None
        )
        if self.qualification_status == "unresolved" and (
            self.qualification_report is not None
            or self.qualification_report_sha256 is not None
        ):
            raise DesignError(
                "unresolved design qualification must not declare a report"
            )
        if self.qualification_status != "unresolved" and not report_complete:
            raise DesignError(
                "resolved design qualification requires a report and SHA-256"
            )
        if self.iteration.max_total_mutations < self.combination.max_mutations:
            raise DesignError(
                "iteration max_total_mutations must cover first-generation combinations"
            )

    @property
    def screen_ready(self) -> bool:
        """返回单点和组合界面筛查所需的科学决定是否齐全。"""
        return (
            self.qualification_status == "qualified"
            and self.method_id is not None
            and self.qualification_report is not None
            and self.qualification_report_sha256 is not None
            and self.objective is not None
            and bool(self.seeds)
            and self.analysis.complete
        )

    @property
    def finalist_ready(self) -> bool:
        """返回柔性复核及 MD 队列所需的阶段四决定是否齐全。"""
        return self.screen_ready and self.finalists.complete

    @property
    def production_ready(self) -> bool:
        """返回阶段四全部正式步骤是否已经放行。"""
        return self.finalist_ready


@dataclass(frozen=True, slots=True)
class DesignTemplate:
    """从阶段三 blind supported cluster 冻结的一个设计模板。"""

    template_id: str
    evidence_role: str
    cluster_id: str
    candidate_id: str
    receptor_id: str
    target: str
    path: Path
    sha256: str
    receptor_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        for name, value in (
            ("template_id", self.template_id),
            ("cluster_id", self.cluster_id),
            ("candidate_id", self.candidate_id),
            ("receptor_id", self.receptor_id),
            ("target", self.target),
        ):
            if not RUN_ID_PATTERN.fullmatch(value):
                raise DesignError(f"{name} must be a safe identifier")
        if self.evidence_role != "positive":
            raise DesignError("single-target template evidence_role must be positive")
        if not SHA256_PATTERN.fullmatch(self.sha256):
            raise DesignError("template sha256 must be a lowercase SHA-256")
        if self.receptor_residue_count < 1:
            raise DesignError("template receptor_residue_count must be positive")


@dataclass(frozen=True, slots=True)
class CandidateParent:
    """候选与一个直接父序列之间的单代谱系。"""

    candidate_id: str
    edit: str
    edit_type: str

    def __post_init__(self) -> None:
        if not RUN_ID_PATTERN.fullmatch(self.candidate_id):
            raise DesignError("parent candidate ID must be a safe identifier")
        if not self.edit.strip():
            raise DesignError("parent edit must not be empty")
        if self.edit_type not in CANDIDATE_EDIT_TYPES:
            raise DesignError("parent edit type is invalid")


@dataclass(frozen=True, slots=True)
class SequenceCandidate:
    """满足固定环拓扑和本次设计轮次合同的一条序列。"""

    candidate_id: str
    sequence: str
    mutation_string: str
    mutation_positions: tuple[int, ...]
    mutation_count: int
    net_charge: int
    design_round: str
    generation: int
    parents: tuple[CandidateParent, ...]
    ancestor_candidate_ids: tuple[str, ...]
    proposal_source: str

    def __post_init__(self) -> None:
        if not RUN_ID_PATTERN.fullmatch(self.candidate_id):
            raise DesignError("candidate_id must be a safe identifier")
        if self.design_round not in DESIGN_ROUNDS:
            raise DesignError(f"unsupported design round: {self.design_round}")
        if self.generation < 1:
            raise DesignError("candidate generation must be positive")
        if not self.parents or len({item.candidate_id for item in self.parents}) != len(
            self.parents
        ):
            raise DesignError("candidate parents must be non-empty and unique")
        if (
            tuple(sorted(self.parents, key=lambda item: item.candidate_id))
            != self.parents
        ):
            raise DesignError("candidate parents must be sorted by candidate ID")
        if (
            not self.ancestor_candidate_ids
            or len(self.ancestor_candidate_ids) != len(set(self.ancestor_candidate_ids))
            or tuple(sorted(self.ancestor_candidate_ids)) != self.ancestor_candidate_ids
            or any(
                not RUN_ID_PATTERN.fullmatch(item)
                for item in self.ancestor_candidate_ids
            )
        ):
            raise DesignError("candidate ancestor IDs must be sorted and unique")
        if any(
            parent.candidate_id not in self.ancestor_candidate_ids
            for parent in self.parents
        ):
            raise DesignError("candidate ancestors must contain every direct parent")
        if self.design_round == "iteration" and self.generation < 2:
            raise DesignError("iteration candidates must be generation 2 or later")
        if self.design_round != "iteration" and self.generation != 1:
            raise DesignError("first-generation candidates must use generation 1")
        if self.generation == 1 and (
            self.parents
            != (CandidateParent("ligand_wt", self.mutation_string, "origin"),)
            or self.ancestor_candidate_ids != ("ligand_wt",)
        ):
            raise DesignError("first-generation candidate lineage is invalid")
        if self.generation > 1 and any(
            parent.edit_type == "origin" for parent in self.parents
        ):
            raise DesignError("iterative candidate cannot use an origin parent edit")
        if self.mutation_count != len(self.mutation_positions):
            raise DesignError("candidate mutation count is inconsistent")
        if not self.mutation_string or not self.proposal_source.strip():
            raise DesignError("mutant candidate provenance must not be empty")


@dataclass(frozen=True, slots=True)
class ScreenTask:
    """候选或配对 WT 在一个模板和 seed 上的界面筛查任务。"""

    task_id: str
    pair_id: str
    state: str
    candidate: SequenceCandidate
    template: DesignTemplate
    seed: int
    resfile_path: Path
    resfile_sha256: str

    def __post_init__(self) -> None:
        if not RUN_ID_PATTERN.fullmatch(self.task_id) or not RUN_ID_PATTERN.fullmatch(
            self.pair_id
        ):
            raise DesignError("screen task identifiers are invalid")
        if self.state not in {"wt", "mutant"}:
            raise DesignError("screen task state must be wt or mutant")
        if self.seed < 0:
            raise DesignError("screen task seed must not be negative")
        if not SHA256_PATTERN.fullmatch(self.resfile_sha256):
            raise DesignError("screen task resfile hash is invalid")


@dataclass(frozen=True, slots=True)
class FinalistStart:
    """由成对界面筛查冻结的一份 WT 或候选柔性精修起点。"""

    start_id: str
    pair_id: str
    state: str
    candidate: SequenceCandidate
    template: DesignTemplate
    path: Path
    sha256: str
    histidine_pose_indices: tuple[int, ...]

    def __post_init__(self) -> None:
        if not RUN_ID_PATTERN.fullmatch(self.start_id) or not RUN_ID_PATTERN.fullmatch(
            self.pair_id
        ):
            raise DesignError("finalist start identifiers are invalid")
        if self.state not in {"wt", "mutant"}:
            raise DesignError("finalist start state must be wt or mutant")
        if not SHA256_PATTERN.fullmatch(self.sha256):
            raise DesignError("finalist start hash is invalid")
        if any(index < 1 for index in self.histidine_pose_indices):
            raise DesignError("finalist histidine pose indices must be positive")


@dataclass(frozen=True, slots=True)
class FinalistTask:
    """一份冻结起点和一个独立 FlexPepDock seed 的柔性复核任务。"""

    task_id: str
    pair_id: str
    start: FinalistStart
    seed: int

    def __post_init__(self) -> None:
        if not RUN_ID_PATTERN.fullmatch(self.task_id) or not RUN_ID_PATTERN.fullmatch(
            self.pair_id
        ):
            raise DesignError("finalist task identifiers are invalid")
        if self.seed < 0:
            raise DesignError("finalist task seed must not be negative")
