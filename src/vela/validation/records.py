"""跨阶段不可变 JSON 记录、文件哈希和安全身份边界。"""

import json
from collections.abc import Mapping
from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError
from vela.core.provenance import JsonScalar, JsonValue, sha256_file
from vela.core.run_identity import RUN_ID_PATTERN
from vela.core.typed_data import object_mapping
from vela.validation.models import ValidationError


@dataclass(frozen=True, slots=True)
class ResumedResult:
    """一个身份和全部文件记录均已核对的完成结果。"""

    path: Path
    document: dict[str, object]
    files: dict[str, Path]


def read_document(path: Path, *, name: str) -> dict[str, object]:
    """读取并收窄一个 UTF-8 JSON 对象。"""
    if not path.is_file():
        raise ValidationError(f"{name} does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name=name)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise ValidationError(f"invalid {name}: {path}") from exc


def file_record(path: Path, *, root: Path) -> dict[str, JsonValue]:
    """建立相对于声明根目录的文件哈希记录。"""
    return {
        "path": path.relative_to(root).as_posix(),
        "sha256": sha256_file(path),
    }


def validate_record(*, root: Path, raw: object, name: str) -> tuple[Path, str]:
    """核对记录路径没有越界且当前文件哈希未变化。"""
    try:
        record = object_mapping(raw, name=name)
    except TypeError as exc:
        raise ValidationError(f"invalid {name} record") from exc
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError(f"invalid {name} record")
    resolved_root = root.resolve()
    path = (resolved_root / relative).resolve()
    if not path.is_relative_to(resolved_root):
        raise ValidationError(f"{name} path escapes its source directory")
    if not path.is_file() or sha256_file(path) != expected_hash:
        raise ValidationError(f"{name} hash mismatch: {path}")
    return path, expected_hash


def resume_completed_result(
    *,
    directory: Path,
    filename: str,
    document_name: str,
    schema: str,
    identity: Mapping[str, JsonScalar],
    plan_hash_key: str,
    plan_hash: str,
    records: Mapping[str, str],
    stale_label: str,
    stale_error: type[VelaError] = ValidationError,
) -> ResumedResult | None:
    """恢复一个完成结果并统一核对身份、计划哈希和文件记录。"""
    result_path = directory / filename
    if not result_path.is_file():
        return None
    document = read_document(result_path, name=document_name)
    if (
        document.get("schema") != schema
        or document.get("status") != "completed"
        or document.get(plan_hash_key) != plan_hash
        or any(document.get(key) != value for key, value in identity.items())
    ):
        identity_value = next(iter(identity.values()), directory.name)
        raise stale_error(f"stale {stale_label}: {identity_value}")
    files = {
        key: validate_record(root=directory, raw=document.get(key), name=name)[0]
        for key, name in records.items()
    }
    return ResumedResult(result_path, document, files)


def safe_identifier(value: object, *, name: str) -> str:
    """把外部记录字段收窄为可安全用作单层路径的身份。"""
    if not isinstance(value, str) or not RUN_ID_PATTERN.fullmatch(value):
        raise ValidationError(f"{name} is not a safe identifier")
    return value


def nonnegative_integer(value: object, *, name: str) -> int:
    """把外部记录字段收窄为非负整数。"""
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValidationError(f"{name} must be a non-negative integer")
    return value
