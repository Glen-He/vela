"""带重试、校验和原子安装的单文件获取。"""

from __future__ import annotations

import http.client
import json
import os
import ssl
import tempfile
import time
import urllib.error
import urllib.request
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.core.provenance import sha256_file
from vela.core.typed_data import object_mapping
from vela.preparation.receptors.models import DownloadConfig, ReceptorError

RETRYABLE_HTTP_CODES = frozenset({408, 425, 429, 500, 502, 503, 504})


@dataclass(frozen=True, slots=True)
class DownloadedFile:
    """一个已经校验并安装到 raw/ 的文件。"""

    source_url: str
    path: Path
    status: str
    size_bytes: int
    sha256: str
    etag: str | None
    last_modified: str | None


def validate_mmcif(path: Path, *, pdb_id: str) -> None:
    """验证 mmCIF 语法和 data block 身份。"""
    try:
        document = gemmi.cif.read_file(str(path))
    except RuntimeError as exc:
        raise ReceptorError(f"downloaded file is not valid mmCIF: {pdb_id}") from exc
    if len(document) != 1:
        raise ReceptorError(f"mmCIF must contain one data block: {pdb_id}")
    block_name = document.sole_block().name.upper()
    accepted_names = {pdb_id.upper(), f"PDB_{pdb_id.upper()}"}
    if block_name not in accepted_names:
        raise ReceptorError(f"mmCIF data block does not match {pdb_id}: {block_name}")


def validate_metadata(path: Path, *, pdb_id: str) -> None:
    """验证 RCSB entry JSON 语法和条目标识。"""
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
    except (UnicodeDecodeError, json.JSONDecodeError) as exc:
        raise ReceptorError(f"downloaded metadata is not valid JSON: {pdb_id}") from exc
    try:
        metadata = object_mapping(value, name=f"metadata {pdb_id}")
    except TypeError as exc:
        raise ReceptorError(f"metadata must be a JSON object: {pdb_id}") from exc
    identifier = metadata.get("rcsb_id")
    if not isinstance(identifier, str) or identifier.upper() != pdb_id.upper():
        raise ReceptorError(f"metadata identifier does not match {pdb_id}")


def _retry_delay(*, error: urllib.error.HTTPError, fallback_seconds: float) -> float:
    value = error.headers.get("Retry-After")
    if value is not None:
        try:
            parsed = float(value)
        except ValueError:
            pass
        else:
            if parsed >= 0:
                return parsed
    return fallback_seconds


def download_file(
    *,
    url: str,
    destination: Path,
    settings: DownloadConfig,
    validator: Callable[[Path], None],
) -> DownloadedFile:
    """下载一个文件, 校验后原子安装; 现有缓存也必须重新校验。"""
    if destination.exists():
        if not destination.is_file():
            raise ReceptorError(f"download destination is not a file: {destination}")
        validator(destination)
        return DownloadedFile(
            source_url=url,
            path=destination,
            status="cached",
            size_bytes=destination.stat().st_size,
            sha256=sha256_file(destination),
            etag=None,
            last_modified=None,
        )

    destination.parent.mkdir(parents=True, exist_ok=True)
    request = urllib.request.Request(
        url,
        headers={"Accept": "*/*", "User-Agent": settings.user_agent},
        method="GET",
    )
    last_error: BaseException | None = None
    for attempt in range(1, settings.retries + 1):
        descriptor, temporary_name = tempfile.mkstemp(
            dir=destination.parent,
            prefix=f".{destination.name}.",
            suffix=".part",
        )
        temporary_path = Path(temporary_name)
        try:
            with (
                os.fdopen(descriptor, "wb") as output,
                urllib.request.urlopen(
                    request, timeout=settings.timeout_seconds
                ) as response,
            ):
                while chunk := response.read(settings.chunk_size_bytes):
                    output.write(chunk)
                output.flush()
                os.fsync(output.fileno())
                etag = response.headers.get("ETag")
                last_modified = response.headers.get("Last-Modified")
            validator(temporary_path)
            temporary_path.replace(destination)
            return DownloadedFile(
                source_url=url,
                path=destination,
                status="downloaded",
                size_bytes=destination.stat().st_size,
                sha256=sha256_file(destination),
                etag=etag,
                last_modified=last_modified,
            )
        except urllib.error.HTTPError as exc:
            last_error = exc
            if exc.code == 404:
                raise ReceptorError(f"remote structure does not exist: {url}") from exc
            if exc.code not in RETRYABLE_HTTP_CODES:
                raise ReceptorError(
                    f"download failed with HTTP {exc.code}: {url}"
                ) from exc
            if attempt < settings.retries:
                time.sleep(
                    _retry_delay(
                        error=exc,
                        fallback_seconds=settings.backoff_initial_seconds
                        * (settings.backoff_multiplier ** (attempt - 1)),
                    )
                )
        except (
            urllib.error.URLError,
            TimeoutError,
            http.client.HTTPException,
            ssl.SSLError,
            ReceptorError,
        ) as exc:
            last_error = exc
            if attempt < settings.retries:
                time.sleep(
                    settings.backoff_initial_seconds
                    * (settings.backoff_multiplier ** (attempt - 1))
                )
        finally:
            temporary_path.unlink(missing_ok=True)
    raise ReceptorError(
        f"download failed after {settings.retries} attempts: {url}"
    ) from last_error
