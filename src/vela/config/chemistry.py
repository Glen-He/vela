"""配体化学配置 section 的严格组装。"""

from collections.abc import Mapping

from vela.config.models import ConfigError
from vela.config.values import assert_keys, document, string, strings, table
from vela.core.typed_data import object_list
from vela.preparation.chemistry import (
    UNRESOLVED,
    ChemistryDefinition,
    DisulfideBond,
    HistidineState,
)


def _unresolved_number(
    source: Mapping[str, object], key: str, *, path: str
) -> float | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    if not isinstance(value, (int, float)) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be a number or unresolved")
    return float(value)


def _unresolved_integer(
    source: Mapping[str, object], key: str, *, path: str
) -> int | None:
    value = source.get(key)
    if value == UNRESOLVED:
        return None
    if not isinstance(value, int) or isinstance(value, bool):
        raise ConfigError(f"{path}.{key} must be an integer or unresolved")
    return value


def _disulfides(value: object) -> tuple[DisulfideBond, ...]:
    try:
        values = object_list(value, name="chemistry.ligand.disulfide_bonds")
    except TypeError as exc:
        raise ConfigError("chemistry.ligand.disulfide_bonds must be an array") from exc
    bonds: list[DisulfideBond] = []
    for index, item in enumerate(values):
        message = f"chemistry.ligand.disulfide_bonds[{index}] must contain two integers"
        try:
            positions = object_list(
                item, name=f"chemistry.ligand.disulfide_bonds[{index}]"
            )
        except TypeError as exc:
            raise ConfigError(message) from exc
        if len(positions) != 2:
            raise ConfigError(message)
        first, second = positions
        if (
            not isinstance(first, int)
            or isinstance(first, bool)
            or not isinstance(second, int)
            or isinstance(second, bool)
        ):
            raise ConfigError(message)
        bonds.append(DisulfideBond(first, second))
    return tuple(bonds)


def parse_chemistry(source: Mapping[str, object]) -> ChemistryDefinition:
    """将 `[chemistry.ligand]` 组装为配体化学身份。"""
    chemistry = table(source, "chemistry", path="")
    section = table(chemistry, "ligand", path="chemistry")
    required = {
        "ligand_id",
        "chemistry_id",
        "sequence",
        "chirality",
        "disulfide_bonds",
        "n_terminus",
        "c_terminus",
        "target_ph",
        "net_charge",
        "histidines",
        "other_modifications_status",
        "other_modifications",
        "decision_sources",
    }
    assert_keys(section, allowed=required, required=required, path="chemistry.ligand")
    histidines = document(section["histidines"], name="chemistry.ligand.histidines")
    histidine_states: list[HistidineState] = []
    for raw_position, value in histidines.items():
        try:
            position = int(raw_position)
        except ValueError as exc:
            raise ConfigError(
                f"chemistry.ligand.histidines key is not a residue number: {raw_position}"
            ) from exc
        if not isinstance(value, str):
            raise ConfigError(
                f"chemistry.ligand.histidines.{raw_position} must be a string"
            )
        histidine_states.append(HistidineState(position, value))
    return ChemistryDefinition(
        ligand_id=string(section, "ligand_id", path="chemistry.ligand"),
        chemistry_id=string(section, "chemistry_id", path="chemistry.ligand"),
        sequence=string(section, "sequence", path="chemistry.ligand"),
        chirality=string(section, "chirality", path="chemistry.ligand"),
        disulfide_bonds=_disulfides(section["disulfide_bonds"]),
        n_terminus=string(section, "n_terminus", path="chemistry.ligand"),
        c_terminus=string(section, "c_terminus", path="chemistry.ligand"),
        target_ph=_unresolved_number(section, "target_ph", path="chemistry.ligand"),
        net_charge=_unresolved_integer(section, "net_charge", path="chemistry.ligand"),
        histidines=tuple(sorted(histidine_states, key=lambda item: item.position)),
        other_modifications_status=string(
            section, "other_modifications_status", path="chemistry.ligand"
        ),
        other_modifications=strings(
            section, "other_modifications", path="chemistry.ligand"
        ),
        decision_sources=strings(section, "decision_sources", path="chemistry.ligand"),
    )
