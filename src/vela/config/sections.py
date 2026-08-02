"""跨领域的路径和下载配置 section 组装。"""

from collections.abc import Mapping
from pathlib import Path

from vela.config.models import PathsConfig
from vela.config.values import (
    assert_keys,
    integer,
    number,
    resolved_path,
    string,
    table,
)
from vela.preparation.receptors.models import DownloadConfig


def parse_paths(source: Mapping[str, object], *, config_dir: Path) -> PathsConfig:
    section = table(source, "paths", path="")
    assert_keys(
        section,
        allowed={"data_dir", "outputs_dir"},
        required={"data_dir", "outputs_dir"},
        path="paths",
    )
    return PathsConfig(
        data_dir=resolved_path(
            string(section, "data_dir", path="paths"), config_dir=config_dir
        ),
        outputs_dir=resolved_path(
            string(section, "outputs_dir", path="paths"), config_dir=config_dir
        ),
    )


def parse_download(source: Mapping[str, object]) -> DownloadConfig:
    section = table(source, "download", path="")
    required = {
        "coordinate_base_url",
        "metadata_base_url",
        "retries",
        "timeout_seconds",
        "backoff_initial_seconds",
        "backoff_multiplier",
        "chunk_size_bytes",
        "user_agent",
    }
    assert_keys(section, allowed=required, required=required, path="download")
    return DownloadConfig(
        coordinate_base_url=string(
            section, "coordinate_base_url", path="download"
        ).rstrip("/"),
        metadata_base_url=string(section, "metadata_base_url", path="download").rstrip(
            "/"
        ),
        retries=integer(section, "retries", path="download"),
        timeout_seconds=number(section, "timeout_seconds", path="download"),
        backoff_initial_seconds=number(
            section, "backoff_initial_seconds", path="download"
        ),
        backoff_multiplier=number(section, "backoff_multiplier", path="download"),
        chunk_size_bytes=integer(section, "chunk_size_bytes", path="download"),
        user_agent=string(section, "user_agent", path="download"),
    )
