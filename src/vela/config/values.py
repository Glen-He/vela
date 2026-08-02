"""TOML、环境变量和基础值的严格解析。"""

from __future__ import annotations

import os
import tomllib
from collections.abc import Mapping
from importlib.resources import files
from pathlib import Path

from vela.config.models import ConfigError
from vela.core.typed_data import object_list, object_mapping


def document(value: object, *, name: str) -> dict[str, object]:
    """将外部映射收窄为字符串键字典。"""
    try:
        return object_mapping(value, name=name)
    except TypeError as exc:
        raise ConfigError(str(exc)) from exc


def table(source: Mapping[str, object], key: str, *, path: str) -> dict[str, object]:
    """读取必需 TOML table。"""
    if key not in source:
        raise ConfigError(
            f"missing [{path}.{key}] table" if path else f"missing [{key}] table"
        )
    name = f"{path}.{key}" if path else key
    return document(source[key], name=name)


def assert_keys(
    source: Mapping[str, object],
    *,
    allowed: set[str],
    required: set[str],
    path: str,
) -> None:
    """拒绝缺失字段和未声明扩展字段。"""
    missing = sorted(required - source.keys())
    extra = sorted(source.keys() - allowed)
    if missing:
        raise ConfigError(f"{path} is missing keys: {', '.join(missing)}")
    if extra:
        raise ConfigError(f"{path} contains unsupported keys: {', '.join(extra)}")


def string(source: Mapping[str, object], key: str, *, path: str) -> str:
    value = source.get(key)
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{path}.{key} must be a non-empty string")
    return value


def boolean(source: Mapping[str, object], key: str, *, path: str) -> bool:
    value = source.get(key)
    if not isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be a boolean")
    return value


def optional_string(source: Mapping[str, object], key: str, *, path: str) -> str | None:
    value = source.get(key)
    if value is None:
        return None
    if not isinstance(value, str) or not value.strip():
        raise ConfigError(f"{path}.{key} must be a non-empty string")
    return value


def optional_strings(
    source: Mapping[str, object], key: str, *, path: str
) -> tuple[str, ...]:
    if key not in source:
        return ()
    try:
        values = object_list(source[key], name=f"{path}.{key}")
    except TypeError as exc:
        raise ConfigError(f"{path}.{key} must be an array of strings") from exc
    result: list[str] = []
    for value in values:
        if not isinstance(value, str) or not value.strip():
            raise ConfigError(f"{path}.{key} must be an array of strings")
        result.append(value)
    return tuple(result)


def integer(source: Mapping[str, object], key: str, *, path: str) -> int:
    value = source.get(key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be an integer")
    return value


def number(source: Mapping[str, object], key: str, *, path: str) -> float:
    value = source.get(key)
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be a number")
    return float(value)


def strings(source: Mapping[str, object], key: str, *, path: str) -> tuple[str, ...]:
    value = source.get(key)
    try:
        values = object_list(value, name=f"{path}.{key}")
    except TypeError as exc:
        raise ConfigError(
            f"{path}.{key} must be an array of non-empty strings"
        ) from exc
    if not values and key not in {"other_modifications", "decision_sources"}:
        raise ConfigError(f"{path}.{key} must be an array of non-empty strings")
    result: list[str] = []
    for item in values:
        if not isinstance(item, str) or not item.strip():
            raise ConfigError(f"{path}.{key} must be an array of non-empty strings")
        result.append(item)
    return tuple(result)


def merge(
    base: Mapping[str, object], override: Mapping[str, object]
) -> dict[str, object]:
    """递归覆盖嵌套 table, 非 table 值整体替换。"""
    merged = dict(base)
    for key, value in override.items():
        previous = merged.get(key)
        try:
            previous_table = object_mapping(previous, name=key)
            value_table = object_mapping(value, name=key)
        except TypeError:
            merged[key] = value
        else:
            merged[key] = merge(previous_table, value_table)
    return merged


def merge_disjoint(
    base: Mapping[str, object], addition: Mapping[str, object], *, path: str = "root"
) -> dict[str, object]:
    """合并外部片段, 拒绝由多个文件声明同一个具体字段。"""
    merged = dict(base)
    for key, value in addition.items():
        field_path = f"{path}.{key}"
        if key not in merged:
            merged[key] = value
            continue
        previous = merged[key]
        try:
            previous_table = object_mapping(previous, name=field_path)
            addition_table = object_mapping(value, name=field_path)
        except TypeError as exc:
            raise ConfigError(
                f"external config field is declared more than once: {field_path}"
            ) from exc
        merged[key] = merge_disjoint(previous_table, addition_table, path=field_path)
    return merged


def read_toml(path: Path) -> dict[str, object]:
    try:
        with path.open("rb") as handle:
            content: object = tomllib.load(handle)
    except tomllib.TOMLDecodeError as exc:
        raise ConfigError(f"invalid TOML in {path}: {exc}") from exc
    return document(content, name=str(path))


def default_document() -> dict[str, object]:
    """按阶段合并随包安装的默认参数。"""
    defaults: dict[str, object] = {}
    resource_dir = files("vela.resources")
    for filename in (
        "preparation.toml",
        "discovery.toml",
        "validation.toml",
        "design.toml",
    ):
        resource = resource_dir.joinpath(filename)
        with resource.open("rb") as handle:
            content: object = tomllib.load(handle)
        defaults = merge_disjoint(
            defaults, document(content, name=f"package defaults {filename}")
        )
    return defaults


def apply_environment(source: dict[str, object]) -> dict[str, object]:
    """应用已声明的运行环境覆盖。"""
    result = dict(source)
    paths = table(result, "paths", path="")
    download = table(result, "download", path="")
    if value := os.environ.get("VELA_DATA_DIR"):
        paths["data_dir"] = str(Path(value).expanduser().resolve())
    if value := os.environ.get("VELA_OUTPUTS_DIR"):
        paths["outputs_dir"] = str(Path(value).expanduser().resolve())
    if value := os.environ.get("VELA_DOWNLOAD_RETRIES"):
        try:
            download["retries"] = int(value)
        except ValueError as exc:
            raise ConfigError("VELA_DOWNLOAD_RETRIES must be an integer") from exc
    if value := os.environ.get("VELA_DOWNLOAD_TIMEOUT_SECONDS"):
        try:
            download["timeout_seconds"] = float(value)
        except ValueError as exc:
            raise ConfigError("VELA_DOWNLOAD_TIMEOUT_SECONDS must be a number") from exc
    result["paths"] = paths
    result["download"] = download
    return result


def resolved_path(value: str, *, config_dir: Path) -> Path:
    candidate = Path(value).expanduser()
    if not candidate.is_absolute():
        candidate = config_dir / candidate
    return candidate.resolve()
