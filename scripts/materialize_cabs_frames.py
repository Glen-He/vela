"""在 CABS 自身运行时中把指定 TRAF 帧恢复为 CA/SC 多模型 PDB。"""

from __future__ import annotations

import argparse
import json
import tarfile
import tempfile
from collections.abc import Iterable, Iterator
from copy import deepcopy
from pathlib import Path
from typing import Protocol, runtime_checkable

from CABS.core.trajectory import Trajectory
from CABS.structures.atom import Atoms


@runtime_checkable
class _ObjectMapping(Protocol):
    def items(self) -> Iterable[tuple[object, object]]: ...


@runtime_checkable
class _ObjectIterable(Protocol):
    def __iter__(self) -> Iterator[object]: ...


def _mapping(value: object, *, name: str) -> dict[str, object]:
    if not isinstance(value, _ObjectMapping):
        raise ValueError(f"{name} must be an object")
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise ValueError(f"{name} keys must be strings")
        result[key] = item
    return result


def _list(value: object, *, name: str) -> list[object]:
    if not isinstance(value, _ObjectIterable) or value.__class__ is not list:
        raise ValueError(f"{name} must be a list")
    return list(value)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--archive", required=True, type=Path)
    parser.add_argument("--selection", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def _selection(path: Path) -> tuple[tuple[int, int], ...]:
    value: object = json.loads(path.read_text(encoding="utf-8"))
    document = _mapping(value, name="selection")
    if set(document) != {"frames"}:
        raise ValueError("selection must contain only frames")
    rows = _list(document["frames"], name="selection frames")
    if not rows:
        raise ValueError("selection frames must be a non-empty list")
    identities: list[tuple[int, int]] = []
    for raw_row in rows:
        row = _mapping(raw_row, name="selection frame identity")
        replica = row.get("replica")
        model = row.get("model")
        if (
            set(row) != {"replica", "model"}
            or not isinstance(replica, int)
            or isinstance(replica, bool)
            or not isinstance(model, int)
            or isinstance(model, bool)
            or replica < 1
            or model < 1
        ):
            raise ValueError("selection frame identity is invalid")
        identities.append((replica, model))
    if len(identities) != len(set(identities)):
        raise ValueError("selection contains duplicate frame identities")
    return tuple(identities)


def _archive_member(archive: tarfile.TarFile, name: str) -> bytes:
    matches = [member for member in archive.getmembers() if member.name == name]
    if len(matches) != 1 or not matches[0].isfile():
        raise ValueError(f"archive must contain one regular {name} member")
    handle = archive.extractfile(matches[0])
    if handle is None:
        raise ValueError(f"archive member is unreadable: {name}")
    return handle.read()


def _trajectory(archive_path: Path) -> Trajectory:
    with tarfile.open(archive_path, mode="r:*") as archive:
        traf = _archive_member(archive, "TRAF")
        seq = _archive_member(archive, "SEQ")
    with tempfile.TemporaryDirectory(prefix="vela-cabs-materialize-") as directory:
        root = Path(directory)
        traf_path = root / "TRAF"
        seq_path = root / "SEQ"
        traf_path.write_bytes(traf)
        seq_path.write_bytes(seq)
        return Trajectory.read_trajectory(str(traf_path), str(seq_path))


def main() -> None:
    arguments = _arguments()
    identities = _selection(arguments.selection)
    trajectory = _trajectory(arguments.archive)
    headers = {
        (header.replica, header.model): index
        for index, header in enumerate(trajectory.headers)
    }
    if len(headers) != len(trajectory.headers):
        raise ValueError("CABS trajectory contains duplicate frame identities")
    coordinates = trajectory.coordinates.reshape(-1, len(trajectory.template), 3)
    output = Atoms()
    for output_index, identity in enumerate(identities, 1):
        source_index = headers.get(identity)
        if source_index is None:
            raise ValueError(f"selected frame is absent: {identity}")
        atoms = deepcopy(trajectory.template)
        atoms.set_model_number(output_index)
        atoms.from_numpy(coordinates[source_index])
        atoms.add_side_chain_centers()
        output.extend(atoms)
    arguments.output.parent.mkdir(parents=True, exist_ok=True)
    output.save_to_pdb(str(arguments.output))


if __name__ == "__main__":
    main()
