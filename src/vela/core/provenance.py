"""跨工作流使用的文件摘要和原子产物写入能力。"""

from __future__ import annotations

import hashlib
import json
import os
import tempfile
from collections.abc import Mapping, Sequence
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

JsonScalar = str | int | float | bool | None
JsonValue = JsonScalar | Mapping[str, "JsonValue"] | Sequence["JsonValue"]

_SOURCE_SUFFIXES = frozenset({".py", ".toml", ".xml"})


def vela_source_sha256() -> str:
    """计算当前 Vela 包源码和包内资源的可复现摘要。"""
    package_root = Path(__file__).resolve().parents[1]
    digest = hashlib.sha256()
    source_paths = sorted(
        path
        for path in package_root.rglob("*")
        if path.is_file()
        and path.suffix in _SOURCE_SUFFIXES
        and "__pycache__" not in path.parts
    )
    for path in source_paths:
        relative = path.relative_to(package_root).as_posix().encode("utf-8")
        digest.update(relative)
        digest.update(b"\0")
        digest.update(path.read_bytes())
        digest.update(b"\0")
    return digest.hexdigest()


def vela_software_identity() -> dict[str, JsonValue]:
    """返回计划文件必须冻结的 Vela 软件身份。"""
    return {
        "vela_version": version("vela"),
        "vela_source_sha256": vela_source_sha256(),
    }


def is_current_vela_software(raw: object) -> bool:
    """判断记录的软件身份是否与当前源码完全一致。"""
    current = vela_software_identity()
    match raw:
        case {
            "vela_version": str(recorded_version),
            "vela_source_sha256": str(recorded_source_sha256),
        }:
            return (
                recorded_version == current["vela_version"]
                and recorded_source_sha256 == current["vela_source_sha256"]
            )
        case _:
            return False


def sha256_file(path: Path, *, chunk_size: int = 1024 * 1024) -> str:
    """计算文件的 SHA-256; 只接受现有普通文件。"""
    if not path.is_file():
        raise FileNotFoundError(path)
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def sha256_text(value: str) -> str:
    """计算 UTF-8 文本的 SHA-256。"""
    return hashlib.sha256(value.encode("utf-8")).hexdigest()


def utc_now() -> str:
    """返回可排序且带时区的 UTC 时间。"""
    return datetime.now(tz=UTC).isoformat(timespec="seconds")


def atomic_write_json(path: Path, value: JsonValue) -> None:
    """先写同目录临时文件, 再原子替换 JSON 目标。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".part",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            json.dump(
                value,
                handle,
                ensure_ascii=False,
                indent=2,
                sort_keys=True,
                allow_nan=False,
            )
            handle.write("\n")
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise


def atomic_write_text(path: Path, value: str) -> None:
    """以 UTF-8 原子写入文本。"""
    path.parent.mkdir(parents=True, exist_ok=True)
    descriptor, temporary_name = tempfile.mkstemp(
        dir=path.parent,
        prefix=f".{path.name}.",
        suffix=".part",
    )
    temporary_path = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w", encoding="utf-8") as handle:
            handle.write(value)
            handle.flush()
            os.fsync(handle.fileno())
        temporary_path.replace(path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
