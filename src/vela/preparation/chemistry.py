"""配体化学身份及 production 就绪性检查。"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from vela.core.errors import VelaError
from vela.core.provenance import JsonValue, atomic_write_json, utc_now

STANDARD_AMINO_ACIDS = frozenset("ACDEFGHIKLMNPQRSTVWY")
HISTIDINE_STATES = frozenset({"HID", "HIE", "HIP"})
UNRESOLVED = "unresolved"
LIGAND_ID_PATTERN = re.compile(r"^[a-z0-9][a-z0-9_-]*$")


class ChemistryError(VelaError):
    """配体化学定义不满足结构合同。"""


@dataclass(frozen=True, slots=True)
class DisulfideBond:
    """以一基编号表示的二硫键。"""

    first: int
    second: int


@dataclass(frozen=True, slots=True)
class HistidineState:
    """一个组氨酸位点的微状态。"""

    position: int
    state: str


@dataclass(frozen=True, slots=True)
class ChemistryDefinition:
    """贯穿对接、设计和 MD 的单一配体化学声明。"""

    ligand_id: str
    chemistry_id: str
    sequence: str
    chirality: str
    disulfide_bonds: tuple[DisulfideBond, ...]
    n_terminus: str
    c_terminus: str
    target_ph: float | None
    net_charge: int | None
    histidines: tuple[HistidineState, ...]
    other_modifications_status: str
    other_modifications: tuple[str, ...]
    decision_sources: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class ChemistryAssessment:
    """化学定义的结构错误和仍待关闭的科学决定。"""

    errors: tuple[str, ...]
    unresolved_fields: tuple[str, ...]

    @property
    def schema_valid(self) -> bool:
        return not self.errors

    @property
    def production_ready(self) -> bool:
        return self.schema_valid and not self.unresolved_fields


def chemistry_record_relative_path(definition: ChemistryDefinition) -> Path:
    """返回配体化学记录在 data 目录下的稳定路径。"""
    if not LIGAND_ID_PATTERN.fullmatch(definition.ligand_id):
        raise ChemistryError(
            f"invalid ligand_id for record path: {definition.ligand_id}"
        )
    return Path("chemistry") / definition.ligand_id / "chemistry_record.json"


def assess_chemistry(definition: ChemistryDefinition) -> ChemistryAssessment:
    """区分无效输入与仍待实验依据关闭的字段。"""
    errors: list[str] = []
    unresolved: list[str] = []
    sequence = definition.sequence

    if not LIGAND_ID_PATTERN.fullmatch(definition.ligand_id):
        errors.append(
            "ligand_id must start with a lowercase letter or digit and contain only lowercase letters, digits, underscores, or hyphens"
        )
    if not definition.chemistry_id.strip():
        errors.append("chemistry_id must not be empty")
    if not sequence:
        errors.append("sequence must not be empty")
    invalid_residues = sorted(set(sequence) - STANDARD_AMINO_ACIDS)
    if invalid_residues:
        errors.append(
            "sequence contains non-standard one-letter residues: "
            + ", ".join(invalid_residues)
        )
    if definition.chirality != "L":
        errors.append("chirality must be L for the current ligand contract")
    if not definition.disulfide_bonds:
        errors.append("at least one disulfide bond is required")

    seen_positions: set[int] = set()
    for bond in definition.disulfide_bonds:
        if bond.first >= bond.second:
            errors.append(
                f"disulfide bond positions must be increasing: {bond.first}, {bond.second}"
            )
            continue
        for position in (bond.first, bond.second):
            if position < 1 or position > len(sequence):
                errors.append(f"disulfide position is outside sequence: {position}")
            elif sequence[position - 1] != "C":
                errors.append(f"disulfide position {position} is not cysteine")
            if position in seen_positions:
                errors.append(
                    f"cysteine position occurs in multiple disulfides: {position}"
                )
            seen_positions.add(position)

    histidine_positions = {
        index for index, residue in enumerate(sequence, 1) if residue == "H"
    }
    declared_histidines = {item.position for item in definition.histidines}
    if histidine_positions != declared_histidines:
        errors.append(
            "histidine state positions do not match sequence: "
            f"expected {sorted(histidine_positions)}, got {sorted(declared_histidines)}"
        )
    for histidine in definition.histidines:
        if histidine.state == UNRESOLVED:
            unresolved.append(f"histidines.{histidine.position}")
        elif histidine.state not in HISTIDINE_STATES:
            errors.append(
                f"histidine {histidine.position} state must be HID, HIE, HIP, or unresolved"
            )

    if definition.n_terminus == UNRESOLVED:
        unresolved.append("n_terminus")
    elif not definition.n_terminus.strip():
        errors.append("n_terminus must not be empty")
    if definition.c_terminus == UNRESOLVED:
        unresolved.append("c_terminus")
    elif not definition.c_terminus.strip():
        errors.append("c_terminus must not be empty")
    if definition.target_ph is None:
        unresolved.append("target_ph")
    elif not 0.0 <= definition.target_ph <= 14.0:
        errors.append("target_ph must be between 0 and 14")
    if definition.net_charge is None:
        unresolved.append("net_charge")
    if definition.other_modifications_status == UNRESOLVED:
        unresolved.append("other_modifications_status")
    elif definition.other_modifications_status not in {"none", "defined"}:
        errors.append("other_modifications_status must be none, defined, or unresolved")
    elif (
        definition.other_modifications_status == "none"
        and definition.other_modifications
    ):
        errors.append("other_modifications must be empty when status is none")
    elif (
        definition.other_modifications_status == "defined"
        and not definition.other_modifications
    ):
        errors.append("other_modifications must not be empty when status is defined")
    if not definition.decision_sources:
        unresolved.append("decision_sources")

    return ChemistryAssessment(tuple(errors), tuple(sorted(set(unresolved))))


def validate_chemistry(definition: ChemistryDefinition) -> ChemistryAssessment:
    """校验结构合同; 科学上未关闭的字段保留在返回值中。"""
    assessment = assess_chemistry(definition)
    if assessment.errors:
        raise ChemistryError("; ".join(assessment.errors))
    return assessment


def chemistry_record(
    definition: ChemistryDefinition, assessment: ChemistryAssessment
) -> dict[str, JsonValue]:
    """建立可供后续阶段引用的显式化学记录。"""
    return {
        "schema": "vela.chemistry-record/1",
        "generated_at": utc_now(),
        "chemistry": {
            "ligand_id": definition.ligand_id,
            "chemistry_id": definition.chemistry_id,
            "sequence": definition.sequence,
            "chirality": definition.chirality,
            "disulfide_bonds": [
                [bond.first, bond.second] for bond in definition.disulfide_bonds
            ],
            "n_terminus": definition.n_terminus,
            "c_terminus": definition.c_terminus,
            "target_ph": definition.target_ph,
            "net_charge": definition.net_charge,
            "histidines": {
                str(item.position): item.state for item in definition.histidines
            },
            "other_modifications_status": definition.other_modifications_status,
            "other_modifications": list(definition.other_modifications),
            "decision_sources": list(definition.decision_sources),
        },
        "assessment": {
            "schema_valid": assessment.schema_valid,
            "production_ready": assessment.production_ready,
            "errors": list(assessment.errors),
            "unresolved_fields": list(assessment.unresolved_fields),
        },
    }


def chemistry_identity_record(definition: ChemistryDefinition) -> dict[str, JsonValue]:
    """返回跨阶段计划中使用的稳定化学身份; 不包含运行时状态。"""
    return {
        "chemistry_id": definition.chemistry_id,
        "ligand_id": definition.ligand_id,
        "sequence": definition.sequence,
        "chirality": definition.chirality,
        "n_terminus": definition.n_terminus,
        "c_terminus": definition.c_terminus,
        "target_ph": definition.target_ph,
        "histidines": {
            str(item.position): item.state for item in definition.histidines
        },
        "disulfide_bonds": [
            [bond.first, bond.second] for bond in definition.disulfide_bonds
        ],
        "other_modifications_status": definition.other_modifications_status,
    }


def write_chemistry_record(
    *, definition: ChemistryDefinition, destination: Path
) -> ChemistryAssessment:
    """校验并原子写入当前配体化学记录。"""
    assessment = validate_chemistry(definition)
    atomic_write_json(destination, chemistry_record(definition, assessment))
    return assessment
