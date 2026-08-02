"""受体登记与下载配置的领域模型。"""

from __future__ import annotations

import re
from dataclasses import dataclass

from vela.core.errors import VelaError

PDB_ID_PATTERN = re.compile(r"^[0-9][A-Za-z0-9]{3}$")
TARGETS = frozenset({"ck2_alpha", "ck2_alpha_prime"})


class ReceptorError(VelaError):
    """受体登记、下载、审计或制备失败。"""


@dataclass(frozen=True, slots=True)
class ReceptorDefinition:
    """一个受体结构及其在项目中的预先声明角色。"""

    receptor_id: str
    pdb_id: str
    target: str
    uniprot_accession: str
    author_chain_id: str
    structure_state: str
    roles: tuple[str, ...]
    prepare: bool
    water_policy: str | None
    remove_components: tuple[str, ...]
    retain_components: tuple[str, ...]
    selection_reason: str | None

    def __post_init__(self) -> None:
        normalized = self.pdb_id.upper()
        object.__setattr__(self, "pdb_id", normalized)
        if not self.receptor_id.strip():
            raise ReceptorError("receptor_id must not be empty")
        if not PDB_ID_PATTERN.fullmatch(normalized):
            raise ReceptorError(f"invalid PDB ID: {self.pdb_id}")
        if self.target not in TARGETS:
            raise ReceptorError(f"unsupported receptor target: {self.target}")
        if not self.uniprot_accession.strip():
            raise ReceptorError("uniprot_accession must not be empty")
        if not self.author_chain_id.strip():
            raise ReceptorError("author_chain_id must not be empty")
        if not self.structure_state.strip():
            raise ReceptorError("structure_state must not be empty")
        if not self.roles or any(not role.strip() for role in self.roles):
            raise ReceptorError("roles must contain at least one non-empty value")
        if self.prepare:
            if self.water_policy not in {"remove_all", "retain_all"}:
                raise ReceptorError(
                    f"{self.receptor_id} water_policy must be remove_all or retain_all"
                )
            if not self.selection_reason:
                raise ReceptorError(
                    f"{self.receptor_id} selection_reason must not be empty"
                )
            overlap = set(self.remove_components) & set(self.retain_components)
            if overlap:
                raise ReceptorError(
                    f"{self.receptor_id} component policies overlap: {sorted(overlap)}"
                )
        elif any(
            (
                self.water_policy is not None,
                self.remove_components,
                self.retain_components,
                self.selection_reason is not None,
            )
        ):
            raise ReceptorError(
                f"{self.receptor_id} has preparation policy but prepare is false"
            )


@dataclass(frozen=True, slots=True)
class DownloadConfig:
    """RCSB 原始文件下载策略。"""

    coordinate_base_url: str
    metadata_base_url: str
    retries: int
    timeout_seconds: float
    backoff_initial_seconds: float
    backoff_multiplier: float
    chunk_size_bytes: int
    user_agent: str

    def __post_init__(self) -> None:
        for name, value in (
            ("coordinate_base_url", self.coordinate_base_url),
            ("metadata_base_url", self.metadata_base_url),
        ):
            if not value.startswith("https://"):
                raise ReceptorError(f"{name} must use HTTPS")
        if self.retries < 1:
            raise ReceptorError("download retries must be at least 1")
        if self.timeout_seconds <= 0:
            raise ReceptorError("download timeout_seconds must be positive")
        if self.backoff_initial_seconds < 0:
            raise ReceptorError("download backoff_initial_seconds must not be negative")
        if self.backoff_multiplier < 1:
            raise ReceptorError("download backoff_multiplier must be at least 1")
        if self.chunk_size_bytes < 1:
            raise ReceptorError("download chunk_size_bytes must be positive")
        if not self.user_agent.strip():
            raise ReceptorError("download user_agent must not be empty")


@dataclass(frozen=True, slots=True)
class CrystalContactConfig:
    """结构审计中的晶体接触搜索参数。"""

    distance_A: float
    min_occupancy: float
    include_hydrogens: bool

    def __post_init__(self) -> None:
        if self.distance_A <= 0:
            raise ReceptorError("crystal contact distance_A must be positive")
        if not 0.0 <= self.min_occupancy <= 1.0:
            raise ReceptorError("crystal contact min_occupancy must be in [0, 1]")


@dataclass(frozen=True, slots=True)
class ReceptorAuditConfig:
    """受体结构审计参数。"""

    crystal_contacts: CrystalContactConfig


@dataclass(frozen=True, slots=True)
class AltlocConfig:
    """替代构象选择参数。"""

    preferred_label: str

    def __post_init__(self) -> None:
        if (
            len(self.preferred_label) != 1
            or not self.preferred_label.isprintable()
            or self.preferred_label.isspace()
            or self.preferred_label == "."
        ):
            raise ReceptorError(
                "altloc preferred_label must be one printable non-blank label"
            )


@dataclass(frozen=True, slots=True)
class ReceptorPreparationConfig:
    """方法无关的基础受体制备参数。"""

    altloc: AltlocConfig
