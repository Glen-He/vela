"""FlexPepDock 局部受体主链自由度的 native-free 选择与 MoveMap 生成。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Literal

import gemmi

from vela.core.provenance import atomic_write_text
from vela.validation.bound_states.assets import PEPTIDE_CHAIN, RECEPTOR_CHAIN
from vela.validation.models import ValidationError

type ReceptorBackboneMode = Literal["fixed", "local_constrained"]


@dataclass(frozen=True, slots=True)
class ReceptorBackboneSelection:
    """由起始复合物界面几何确定的受体主链可动范围。"""

    direct_contact_pose_indices: tuple[int, ...]
    flexible_pose_indices: tuple[int, ...]


def resolve_receptor_backbone_mode(value: object) -> ReceptorBackboneMode:
    """把外部模式值收窄为诊断协议允许的两种受体主链合同。"""
    if value == "fixed":
        return "fixed"
    if value == "local_constrained":
        return "local_constrained"
    raise ValidationError("receptor backbone mode must be fixed or local_constrained")


def _amino_acids(chain: gemmi.Chain) -> tuple[gemmi.Residue, ...]:
    return tuple(
        residue for residue in chain if any(atom.name == "CA" for atom in residue)
    )


def select_local_receptor_backbone(
    *, path: Path, contact_A: float, sequence_padding: int
) -> ReceptorBackboneSelection:
    """按起始肽重原子接触壳层选择受体残基, 并沿序列扩展连续自由度。"""
    if contact_A <= 0:
        raise ValidationError("local receptor contact distance must be positive")
    if sequence_padding < 0:
        raise ValidationError("local receptor sequence padding must not be negative")
    try:
        structure = gemmi.read_structure(str(path))
    except RuntimeError as exc:
        raise ValidationError(
            f"invalid local receptor selection structure: {path}"
        ) from exc
    if len(structure) != 1:
        raise ValidationError(
            f"local receptor selection requires one structure model: {path}"
        )
    chains = {chain.name: chain for chain in structure[0]}
    if set(chains) != {RECEPTOR_CHAIN, PEPTIDE_CHAIN}:
        raise ValidationError(
            f"local receptor selection requires only A/P chains: {path}"
        )
    receptor = _amino_acids(chains[RECEPTOR_CHAIN])
    peptide = _amino_acids(chains[PEPTIDE_CHAIN])
    peptide_atoms = tuple(
        atom for residue in peptide for atom in residue if atom.element.name != "H"
    )
    if not receptor or not peptide_atoms:
        raise ValidationError(f"local receptor selection lacks interface atoms: {path}")
    direct = tuple(
        index
        for index, residue in enumerate(receptor, 1)
        if any(
            receptor_atom.element.name != "H"
            and any(
                receptor_atom.pos.dist(peptide_atom.pos) <= contact_A
                for peptide_atom in peptide_atoms
            )
            for receptor_atom in residue
        )
    )
    if not direct:
        raise ValidationError(
            "local receptor selection found no receptor residue in the contact shell"
        )
    flexible = tuple(
        sorted(
            {
                neighbor
                for index in direct
                for neighbor in range(
                    max(1, index - sequence_padding),
                    min(len(receptor), index + sequence_padding) + 1,
                )
            }
        )
    )
    return ReceptorBackboneSelection(direct, flexible)


def _contiguous_ranges(indices: tuple[int, ...]) -> tuple[tuple[int, int], ...]:
    if not indices or tuple(sorted(set(indices))) != indices:
        raise ValidationError(
            "local receptor flexible indices must be non-empty, unique, and sorted"
        )
    ranges: list[tuple[int, int]] = []
    start = previous = indices[0]
    for index in indices[1:]:
        if index == previous + 1:
            previous = index
            continue
        ranges.append((start, previous))
        start = previous = index
    ranges.append((start, previous))
    return tuple(ranges)


def write_local_receptor_movemap(
    *,
    destination: Path,
    receptor_residue_count: int,
    peptide_residue_count: int,
    flexible_receptor_pose_indices: tuple[int, ...],
) -> None:
    """只开放界面受体主链、全部受体侧链、肽自由度和唯一对接 jump。"""
    if receptor_residue_count < 1 or peptide_residue_count < 1:
        raise ValidationError("FlexPepDock MoveMap residue counts must be positive")
    ranges = _contiguous_ranges(flexible_receptor_pose_indices)
    if ranges[0][0] < 1 or ranges[-1][1] > receptor_residue_count:
        raise ValidationError("local receptor flexible index is outside the receptor")
    peptide_first = receptor_residue_count + 1
    peptide_last = receptor_residue_count + peptide_residue_count
    lines = [f"RESIDUE 1 {receptor_residue_count} CHI"]
    lines.extend(f"RESIDUE {first} {last} BBCHI" for first, last in ranges)
    lines.extend(
        (
            f"RESIDUE {peptide_first} {peptide_last} BBCHI",
            "JUMP * NO",
            "JUMP 1 YES",
        )
    )
    atomic_write_text(destination, "\n".join(lines) + "\n")
