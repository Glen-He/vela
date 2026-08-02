"""配置组合层的数据模型与错误类型。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError
from vela.design.models import DesignSettings
from vela.discovery.models import DiscoverySettings
from vela.preparation.chemistry import ChemistryDefinition
from vela.preparation.receptors.models import (
    DownloadConfig,
    ReceptorAuditConfig,
    ReceptorDefinition,
    ReceptorPreparationConfig,
)
from vela.validation.models import ValidationSettings


class ConfigError(VelaError):
    """配置文件缺失、字段错误或值冲突。"""


@dataclass(frozen=True, slots=True)
class PathsConfig:
    """项目数据和运行结果根目录。"""

    data_dir: Path
    outputs_dir: Path


@dataclass(frozen=True, slots=True)
class AppConfig:
    """已经合并并严格校验的项目配置。"""

    source_dir: Path
    source_files: tuple[Path, ...]
    source_snapshot_text: str
    source_snapshot_sha256: str
    paths: PathsConfig
    download: DownloadConfig
    audit: ReceptorAuditConfig
    preparation: ReceptorPreparationConfig
    chemistry: ChemistryDefinition
    receptors: tuple[ReceptorDefinition, ...]
    discovery: DiscoverySettings
    validation: ValidationSettings
    design: DesignSettings
