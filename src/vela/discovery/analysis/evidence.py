"""阶段二分析使用的引擎无关 pose 与 site 数据模型。"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

from vela.discovery.models import SHA256_PATTERN, DiscoveryError

Point3D = tuple[float, float, float]


@dataclass(frozen=True, slots=True)
class PoseEvidence:
    """全局引擎输出的一个可追溯粗粒化候选或基线 pose。"""

    task_id: str
    pose_id: str
    receptor_id: str
    target: str
    seed: int
    model_path: Path
    model_sha256: str
    model_index: int
    contact_residues: frozenset[str]
    local_position: Point3D
    coordinate_frame_id: str
    ranking_score: float
    score_name: str
    qc_status: str

    def __post_init__(self) -> None:
        for name, value in (
            ("task_id", self.task_id),
            ("pose_id", self.pose_id),
            ("receptor_id", self.receptor_id),
            ("target", self.target),
            ("coordinate_frame_id", self.coordinate_frame_id),
            ("score_name", self.score_name),
        ):
            if not value.strip():
                raise DiscoveryError(f"{name} must not be empty")
        if self.seed < 0:
            raise DiscoveryError("pose seed must not be negative")
        if self.model_index < 1:
            raise DiscoveryError("model_index must be positive")
        if not SHA256_PATTERN.fullmatch(self.model_sha256):
            raise DiscoveryError("model_sha256 must be a lowercase SHA-256")
        if self.qc_status not in {"passed", "failed"}:
            raise DiscoveryError("qc_status must be passed or failed")
        if self.qc_status == "passed" and not self.contact_residues:
            raise DiscoveryError("passed pose must contain receptor contacts")
        if any(not value.strip() for value in self.contact_residues):
            raise DiscoveryError("contact residue identifiers must not be empty")
        if not all(math.isfinite(value) for value in self.local_position):
            raise DiscoveryError("local_position must contain finite coordinates")
        if not math.isfinite(self.ranking_score):
            raise DiscoveryError("ranking_score must be finite")


@dataclass(frozen=True, slots=True)
class ReceptorSite:
    """单个受体构象内由多个 pose 形成的 site。"""

    site_id: str
    receptor_id: str
    target: str
    coordinate_frame_id: str
    pose_ids: tuple[str, ...]
    supporting_seeds: tuple[int, ...]
    representative_pose_id: str
    representative_contacts: frozenset[str]
    representative_position: Point3D
    pose_count: int
    supported: bool


@dataclass(frozen=True, slots=True)
class CandidateSite:
    """同一亚型多个受体构象之间对应的 candidate site。"""

    candidate_id: str
    target: str
    coordinate_frame_id: str
    receptor_site_ids: tuple[str, ...]
    receptor_ids: tuple[str, ...]
    seed_support_by_receptor: tuple[str, ...]
    representative_site_id: str
    receptor_support: int
    evidence_tier: str
    rank_within_tier: int
    minimum_seed_support: int
    total_seed_support: int
    maximum_normalized_site_distance: float
    minimum_selected_pose_fraction: float
    total_selected_pose_fraction: float
    median_receptor_score_quantile: float
    handoff_eligible: bool


@dataclass(frozen=True, slots=True)
class SiteAnalysisResult:
    """阶段二按亚型隔离的完整 site 分析结果。"""

    receptor_sites: tuple[ReceptorSite, ...]
    candidate_sites: tuple[CandidateSite, ...]
