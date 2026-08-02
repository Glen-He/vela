"""受体登记表的 RCSB 获取与下载 manifest 编排。"""

from functools import partial
from pathlib import Path

from vela.core.provenance import JsonValue, atomic_write_json, utc_now
from vela.preparation.receptors.acquisition import (
    DownloadedFile,
    download_file,
    validate_metadata,
    validate_mmcif,
)
from vela.preparation.receptors.models import (
    DownloadConfig,
    ReceptorDefinition,
    ReceptorError,
)


def _file_record(file: DownloadedFile, *, data_dir: Path) -> dict[str, JsonValue]:
    try:
        relative_path = file.path.relative_to(data_dir)
    except ValueError as exc:
        raise ReceptorError(
            f"downloaded file is outside data_dir: {file.path}"
        ) from exc
    return {
        "source_url": file.source_url,
        "path": relative_path.as_posix(),
        "status": file.status,
        "size_bytes": file.size_bytes,
        "sha256": file.sha256,
        "etag": file.etag,
        "last_modified": file.last_modified,
    }


def download_receptors(
    *,
    definitions: tuple[ReceptorDefinition, ...],
    settings: DownloadConfig,
    data_dir: Path,
) -> tuple[DownloadedFile, ...]:
    """下载登记表中的全部原始对象, 并逐条更新 manifest。"""
    raw_dir = data_dir / "receptors" / "raw"
    manifest_path = raw_dir / "download_manifest.json"
    files: list[DownloadedFile] = []
    records: list[dict[str, JsonValue]] = []
    for definition in definitions:
        pdb_id = definition.pdb_id
        coordinate = download_file(
            url=f"{settings.coordinate_base_url}/{pdb_id}.cif",
            destination=raw_dir / f"{pdb_id}.cif",
            settings=settings,
            validator=partial(validate_mmcif, pdb_id=pdb_id),
        )
        metadata = download_file(
            url=f"{settings.metadata_base_url}/{pdb_id}",
            destination=raw_dir / f"{pdb_id}.entry.json",
            settings=settings,
            validator=partial(validate_metadata, pdb_id=pdb_id),
        )
        files.extend((coordinate, metadata))
        records.append(
            {
                "pdb_id": pdb_id,
                "receptor_id": definition.receptor_id,
                "files": [
                    _file_record(coordinate, data_dir=data_dir),
                    _file_record(metadata, data_dir=data_dir),
                ],
            }
        )
        atomic_write_json(
            manifest_path,
            {
                "schema": "vela.receptor-download-manifest/1",
                "verified_at": utc_now(),
                "entries": records,
            },
        )
    return tuple(files)
