"""阶段二采样层与 site 分析层之间的规范 pose TSV 边界。"""

from __future__ import annotations

import csv
import io
from pathlib import Path

from vela.core.provenance import atomic_write_text, sha256_file
from vela.discovery.analysis.evidence import PoseEvidence
from vela.discovery.models import DiscoveryError

POSE_FIELDS = (
    "task_id",
    "pose_id",
    "receptor_id",
    "target",
    "seed",
    "model_path",
    "model_sha256",
    "model_index",
    "contact_residues",
    "local_x_A",
    "local_y_A",
    "local_z_A",
    "coordinate_frame_id",
    "ranking_score",
    "score_name",
    "qc_status",
)


def write_pose_evidence(
    *, poses: tuple[PoseEvidence, ...], path: Path, run_dir: Path
) -> None:
    """将粗粒化候选写入唯一的阶段二 TSV 合同。"""
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=POSE_FIELDS, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for pose in poses:
        try:
            model_path = pose.model_path.resolve().relative_to(run_dir.resolve())
        except ValueError as exc:
            raise DiscoveryError(
                f"pose model is outside run directory: {pose.model_path}"
            ) from exc
        writer.writerow(
            {
                "task_id": pose.task_id,
                "pose_id": pose.pose_id,
                "receptor_id": pose.receptor_id,
                "target": pose.target,
                "seed": str(pose.seed),
                "model_path": model_path.as_posix(),
                "model_sha256": pose.model_sha256,
                "model_index": str(pose.model_index),
                "contact_residues": ";".join(sorted(pose.contact_residues)),
                "local_x_A": f"{pose.local_position[0]:.6f}",
                "local_y_A": f"{pose.local_position[1]:.6f}",
                "local_z_A": f"{pose.local_position[2]:.6f}",
                "coordinate_frame_id": pose.coordinate_frame_id,
                "ranking_score": f"{pose.ranking_score:.6f}",
                "score_name": pose.score_name,
                "qc_status": pose.qc_status,
            }
        )
    atomic_write_text(path, output.getvalue())


def _rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file():
        raise DiscoveryError(f"pose evidence table does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            header = next(reader)
        except StopIteration as exc:
            raise DiscoveryError(f"pose evidence table is empty: {path}") from exc
        if len(header) != len(set(header)) or not set(POSE_FIELDS).issubset(header):
            raise DiscoveryError(
                "pose evidence table has missing or duplicate required columns"
            )
        rows: list[dict[str, str]] = []
        for line_number, values in enumerate(reader, 2):
            if len(values) != len(header):
                raise DiscoveryError(
                    f"pose evidence row {line_number} has {len(values)} columns; "
                    f"expected {len(header)}"
                )
            rows.append(dict(zip(header, values, strict=True)))
    return rows


def _model_path(*, raw: str, run_dir: Path) -> Path:
    path = Path(raw)
    resolved = path.resolve() if path.is_absolute() else (run_dir / path).resolve()
    try:
        resolved.relative_to(run_dir.resolve())
    except ValueError as exc:
        raise DiscoveryError(
            f"pose model is outside run directory: {resolved}"
        ) from exc
    if not resolved.is_file():
        raise DiscoveryError(f"pose model does not exist: {resolved}")
    return resolved


def _integer(value: str, *, field: str, pose_id: str) -> int:
    try:
        return int(value)
    except ValueError as exc:
        raise DiscoveryError(f"{pose_id} has invalid integer {field}: {value}") from exc


def _float(value: str, *, field: str, pose_id: str) -> float:
    try:
        return float(value)
    except ValueError as exc:
        raise DiscoveryError(f"{pose_id} has invalid number {field}: {value}") from exc


def read_pose_evidence(*, path: Path, run_dir: Path) -> tuple[PoseEvidence, ...]:
    """读取、收窄并复核所有规范粗粒化候选记录。"""
    poses: list[PoseEvidence] = []
    for row in _rows(path):
        pose_id = row["pose_id"]
        model_path = _model_path(raw=row["model_path"], run_dir=run_dir)
        expected_hash = row["model_sha256"]
        if sha256_file(model_path) != expected_hash:
            raise DiscoveryError(f"pose model hash mismatch: {model_path}")
        poses.append(
            PoseEvidence(
                task_id=row["task_id"],
                pose_id=pose_id,
                receptor_id=row["receptor_id"],
                target=row["target"],
                seed=_integer(row["seed"], field="seed", pose_id=pose_id),
                model_path=model_path,
                model_sha256=expected_hash,
                model_index=_integer(
                    row["model_index"], field="model_index", pose_id=pose_id
                ),
                contact_residues=frozenset(
                    value for value in row["contact_residues"].split(";") if value
                ),
                local_position=(
                    _float(row["local_x_A"], field="local_x_A", pose_id=pose_id),
                    _float(row["local_y_A"], field="local_y_A", pose_id=pose_id),
                    _float(row["local_z_A"], field="local_z_A", pose_id=pose_id),
                ),
                coordinate_frame_id=row["coordinate_frame_id"],
                ranking_score=_float(
                    row["ranking_score"], field="ranking_score", pose_id=pose_id
                ),
                score_name=row["score_name"],
                qc_status=row["qc_status"],
            )
        )
    return tuple(poses)
