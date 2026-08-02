"""调用 CABS 自身运行时恢复完整 TRAF 中选定帧的 CA/SC 表示。"""

from __future__ import annotations

import os
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import JsonValue, atomic_write_json, sha256_file
from vela.discovery.models import CabsDockSettings, DiscoveryError


@dataclass(frozen=True, slots=True, order=True)
class CabsFrameIdentity:
    """一个由 replica 与模型号共同确定的原始 TRAF 帧。"""

    replica: int
    model: int

    def __post_init__(self) -> None:
        if self.replica < 1 or self.model < 1:
            raise DiscoveryError("CABS frame identity values must be positive")


def materializer_script_path() -> Path:
    """返回受 Vela 计划哈希保护的 CABS 外部运行时适配脚本。"""
    path = Path(__file__).resolve().parents[4] / "scripts/materialize_cabs_frames.py"
    if not path.is_file():
        raise DiscoveryError(f"CABS materializer script is missing: {path}")
    return path


def materializer_record() -> dict[str, JsonValue]:
    """返回外部运行时适配脚本的可核验记录。"""
    path = materializer_script_path()
    return {"path": path.as_posix(), "sha256": sha256_file(path)}


def cabs_python_executable(settings: CabsDockSettings) -> Path:
    """返回与已验证 CABSdock 入口位于同一环境的 Python。"""
    executable = settings.executable.parent / "python"
    if not executable.is_file() or not os.access(executable, os.X_OK):
        raise DiscoveryError(f"CABS Python executable is not runnable: {executable}")
    return executable


def materialize_cabs_frames(
    *,
    archive_path: Path,
    identities: tuple[CabsFrameIdentity, ...],
    task_dir: Path,
    settings: CabsDockSettings,
) -> tuple[Path, Path]:
    """以冻结顺序写出候选帧清单和对应 CA/SC 多模型 PDB。"""
    if not identities or len(identities) != len(set(identities)):
        raise DiscoveryError(
            "CABS materialization identities must be non-empty and unique"
        )
    selection_path = task_dir / "output_data" / "trajectory_candidate_frames.json"
    output_path = task_dir / "output_pdbs" / "trajectory_candidates.pdb"
    atomic_write_json(
        selection_path,
        {
            "frames": [
                {"replica": identity.replica, "model": identity.model}
                for identity in identities
            ]
        },
    )
    descriptor, temporary_name = tempfile.mkstemp(
        dir=output_path.parent,
        prefix=f".{output_path.name}.",
        suffix=".part",
    )
    os.close(descriptor)
    temporary_path = Path(temporary_name)
    try:
        result = subprocess.run(
            [
                str(cabs_python_executable(settings)),
                str(materializer_script_path()),
                "--archive",
                str(archive_path),
                "--selection",
                str(selection_path),
                "--output",
                str(temporary_path),
            ],
            check=False,
            capture_output=True,
            text=True,
        )
        if result.returncode != 0 or not temporary_path.is_file():
            detail = (result.stderr or result.stdout).strip()
            raise DiscoveryError(
                "CABS trajectory candidate materialization failed: "
                + (detail or f"exit code {result.returncode}")
            )
        temporary_path.replace(output_path)
    except BaseException:
        temporary_path.unlink(missing_ok=True)
        raise
    return selection_path, output_path
