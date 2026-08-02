"""Rosetta scorefile 的严格解析和局部恢复判定。"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import atomic_write_text, sha256_file
from vela.validation.models import LocalRecoveryControl, ValidationError

DESCRIPTION_PATTERN = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_.-]*$")


@dataclass(frozen=True, slots=True)
class RosettaScoreRow:
    """一个具有唯一 decoy 身份的 Rosetta 分数行。"""

    description: str
    scores: dict[str, float]

    def score(self, name: str) -> float:
        try:
            return self.scores[name]
        except KeyError as exc:
            raise ValidationError(
                f"Rosetta scorefile lacks required column: {name}"
            ) from exc


def index_rosetta_pdb_outputs(directory: Path) -> dict[str, Path]:
    """按 Rosetta description 索引一个任务目录中的 PDB 输出。"""
    return {path.stem: path for path in directory.glob("*.pdb")}


def read_rosetta_scorefile(path: Path) -> tuple[RosettaScoreRow, ...]:
    """读取一个完整、有限值且 decoy 身份唯一的 scorefile。"""
    if not path.is_file():
        raise ValidationError(f"Rosetta scorefile does not exist: {path}")
    header: tuple[str, ...] | None = None
    rows: list[RosettaScoreRow] = []
    descriptions: set[str] = set()
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except UnicodeDecodeError as exc:
        raise ValidationError(f"Rosetta scorefile is not UTF-8: {path}") from exc
    for line_number, line in enumerate(lines, 1):
        if not line.startswith("SCORE:"):
            continue
        fields = tuple(line.split()[1:])
        if not fields:
            continue
        if fields[0] == "total_score":
            if not fields or fields[-1] != "description":
                raise ValidationError(
                    f"invalid Rosetta score header at {path}:{line_number}"
                )
            if header is not None and fields != header:
                raise ValidationError(f"inconsistent Rosetta score headers: {path}")
            header = fields
            continue
        if header is None:
            raise ValidationError(
                f"Rosetta score row appears before its header: {path}:{line_number}"
            )
        if len(fields) != len(header):
            raise ValidationError(
                f"Rosetta score row width mismatch at {path}:{line_number}"
            )
        description = fields[-1]
        if not DESCRIPTION_PATTERN.fullmatch(description):
            raise ValidationError(
                f"unsafe Rosetta decoy description at {path}:{line_number}"
            )
        if description in descriptions:
            raise ValidationError(f"duplicate Rosetta decoy description: {description}")
        scores: dict[str, float] = {}
        for name, raw_value in zip(header[:-1], fields[:-1], strict=True):
            try:
                value = float(raw_value)
            except ValueError as exc:
                raise ValidationError(
                    f"invalid Rosetta score value at {path}:{line_number}"
                ) from exc
            if not math.isfinite(value):
                raise ValidationError(
                    f"non-finite Rosetta score value at {path}:{line_number}"
                )
            scores[name] = value
        descriptions.add(description)
        rows.append(RosettaScoreRow(description, scores))
    if header is None or not rows:
        raise ValidationError(f"Rosetta scorefile contains no decoys: {path}")
    return tuple(rows)


def write_decoy_manifest(
    *,
    rows: tuple[RosettaScoreRow, ...],
    control: LocalRecoveryControl,
    decoy_paths: dict[str, Path],
    task_dir: Path,
) -> Path:
    """记录全部 decoy 的关键分数、恢复判定、相对路径和哈希。"""
    expected = {row.description for row in rows}
    if set(decoy_paths) != expected:
        missing = sorted(expected - set(decoy_paths))
        extra = sorted(set(decoy_paths) - expected)
        raise ValidationError(
            "Rosetta score/PDB identities differ: "
            f"missing={missing[:3]}, extra={extra[:3]}"
        )
    lines = ["description\tranking_score\trecovery_rmsd_A\trecovered\tpath\tsha256\n"]
    for row in rows:
        path = decoy_paths[row.description]
        relative = path.relative_to(task_dir).as_posix()
        rmsd = row.score(control.recovery_rmsd_score)
        lines.append(
            "\t".join(
                (
                    row.description,
                    f"{row.score(control.ranking_score):.6f}",
                    f"{rmsd:.6f}",
                    str(rmsd <= control.max_recovery_rmsd_A).lower(),
                    relative,
                    sha256_file(path),
                )
            )
            + "\n"
        )
    path = task_dir / "decoy_manifest.tsv"
    atomic_write_text(path, "".join(lines))
    return path


def write_refinement_decoy_manifest(
    *,
    rows: tuple[RosettaScoreRow, ...],
    ranking_score: str,
    decoy_paths: dict[str, Path],
    task_dir: Path,
) -> Path:
    """记录正式候选精修的全部 decoy、配置排名分数、路径和哈希。"""
    expected = {row.description for row in rows}
    if set(decoy_paths) != expected:
        missing = sorted(expected - set(decoy_paths))
        extra = sorted(set(decoy_paths) - expected)
        raise ValidationError(
            "Rosetta score/PDB identities differ: "
            f"missing={missing[:3]}, extra={extra[:3]}"
        )
    lines = ["description\tranking_score\tpath\tsha256\n"]
    for row in rows:
        path = decoy_paths[row.description]
        lines.append(
            "\t".join(
                (
                    row.description,
                    f"{row.score(ranking_score):.6f}",
                    path.relative_to(task_dir).as_posix(),
                    sha256_file(path),
                )
            )
            + "\n"
        )
    path = task_dir / "decoy_manifest.tsv"
    atomic_write_text(path, "".join(lines))
    return path
