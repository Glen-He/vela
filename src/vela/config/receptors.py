"""受体登记表配置的严格组装。"""

from collections.abc import Mapping

from vela.config.models import ConfigError
from vela.config.values import (
    assert_keys,
    boolean,
    document,
    optional_string,
    optional_strings,
    string,
    strings,
)
from vela.core.typed_data import object_list
from vela.preparation.receptors.models import ReceptorDefinition


def parse_receptors(
    source: Mapping[str, object],
) -> tuple[ReceptorDefinition, ...]:
    """组装并校验全部 `[[receptors]]` 条目。"""
    try:
        entries = object_list(source.get("receptors"), name="receptors")
    except TypeError as exc:
        raise ConfigError("[[receptors]] must contain at least one entry") from exc
    if not entries:
        raise ConfigError("[[receptors]] must contain at least one entry")
    required = {
        "receptor_id",
        "pdb_id",
        "target",
        "uniprot_accession",
        "author_chain_id",
        "structure_state",
        "roles",
        "prepare",
    }
    optional = {
        "water_policy",
        "remove_components",
        "retain_components",
        "selection_reason",
    }
    receptors: list[ReceptorDefinition] = []
    for index, value in enumerate(entries):
        path = f"receptors[{index}]"
        section = document(value, name=path)
        assert_keys(
            section,
            allowed=required | optional,
            required=required,
            path=path,
        )
        receptors.append(
            ReceptorDefinition(
                receptor_id=string(section, "receptor_id", path=path),
                pdb_id=string(section, "pdb_id", path=path),
                target=string(section, "target", path=path),
                uniprot_accession=string(section, "uniprot_accession", path=path),
                author_chain_id=string(section, "author_chain_id", path=path),
                structure_state=string(section, "structure_state", path=path),
                roles=strings(section, "roles", path=path),
                prepare=boolean(section, "prepare", path=path),
                water_policy=optional_string(section, "water_policy", path=path),
                remove_components=optional_strings(
                    section, "remove_components", path=path
                ),
                retain_components=optional_strings(
                    section, "retain_components", path=path
                ),
                selection_reason=optional_string(
                    section, "selection_reason", path=path
                ),
            )
        )
    identifiers = [item.receptor_id for item in receptors]
    if len(identifiers) != len(set(identifiers)):
        raise ConfigError("receptor_id values must be unique")
    pdb_ids = [item.pdb_id for item in receptors]
    if len(pdb_ids) != len(set(pdb_ids)):
        raise ConfigError(
            "pdb_id values must be unique in the current receptor registry"
        )
    return tuple(receptors)
