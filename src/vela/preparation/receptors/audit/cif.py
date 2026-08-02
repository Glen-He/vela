"""mmCIF category、条目元数据和审计行提取。"""

from __future__ import annotations

import json
from collections import defaultdict
from collections.abc import Mapping
from pathlib import Path
from typing import Protocol

import gemmi

from vela.core.typed_data import object_list, object_mapping
from vela.preparation.receptors.models import ReceptorDefinition, ReceptorError


class CategoryBlock(Protocol):
    def get_mmcif_category(self, name: str, raw: bool = False) -> object:
        """返回一个 mmCIF category。"""


def clean(value: object) -> str:
    if value is None or isinstance(value, bool):
        return ""
    text = str(value).strip()
    if text in {"?", "."}:
        return ""
    if len(text) >= 2 and text[0] == text[-1] and text[0] in {"'", '"'}:
        text = text[1:-1]
    return text.replace("\n", " ")


def category_rows(block: CategoryBlock, category: str) -> list[dict[str, str]]:
    try:
        raw = object_mapping(block.get_mmcif_category(f"_{category}."), name=category)
    except TypeError as exc:
        raise ReceptorError(f"invalid mmCIF category: {category}") from exc
    if not raw:
        return []
    columns: dict[str, list[object]] = {}
    for key, values in raw.items():
        try:
            columns[key] = object_list(values, name=f"{category}.{key}")
        except TypeError as exc:
            raise ReceptorError(f"invalid mmCIF column: {category}.{key}") from exc
    size = max((len(values) for values in columns.values()), default=0)
    return [
        {
            key: clean(values[index]) if index < len(values) else ""
            for key, values in columns.items()
        }
        for index in range(size)
    ]


def first_value(block: gemmi.cif.Block, tags: tuple[str, ...]) -> str:
    for tag in tags:
        value = clean(block.find_value(tag))
        if value:
            return value
    return ""


def read_metadata(path: Path) -> dict[str, object]:
    value: object = json.loads(path.read_text(encoding="utf-8"))
    try:
        return object_mapping(value, name=str(path))
    except TypeError as exc:
        raise ReceptorError(f"metadata must contain an object: {path}") from exc


def joined(values: set[str]) -> str:
    return ";".join(sorted(value for value in values if value))


def optional_float(value: str) -> float | None:
    if not value:
        return None
    try:
        return float(value)
    except ValueError:
        return None


def mapping_by_author_chain(
    block: gemmi.cif.Block,
) -> dict[str, list[dict[str, str]]]:
    references = {row.get("id", ""): row for row in category_rows(block, "struct_ref")}
    mappings: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in category_rows(block, "struct_ref_seq"):
        reference = references.get(row.get("ref_id", ""), {})
        mapping = {
            "accession": row.get("pdbx_db_accession", "")
            or reference.get("pdbx_db_accession", ""),
            "database": reference.get("db_name", ""),
            "database_code": reference.get("db_code", ""),
            "database_start": row.get("db_align_beg", ""),
            "database_end": row.get("db_align_end", ""),
            "author_start": row.get("pdbx_auth_seq_align_beg", ""),
            "author_end": row.get("pdbx_auth_seq_align_end", ""),
        }
        for chain in row.get("pdbx_strand_id", "").split(","):
            if chain.strip():
                mappings[chain.strip()].append(mapping)
    return mappings


def revision_date(metadata: Mapping[str, object]) -> str:
    try:
        accession = object_mapping(
            metadata.get("rcsb_accession_info"), name="rcsb_accession_info"
        )
    except TypeError:
        return ""
    value = accession.get("revision_date")
    return value if isinstance(value, str) else ""


def missing_rows(
    *, block: gemmi.cif.Block, definition: ReceptorDefinition
) -> list[dict[str, str]]:
    return [
        {
            "receptor_id": definition.receptor_id,
            "pdb_id": definition.pdb_id,
            "auth_chain_id": row.get("auth_asym_id", ""),
            "label_chain_id": row.get("label_asym_id", ""),
            "residue_name": row.get("auth_comp_id", ""),
            "auth_seq_id": row.get("auth_seq_id", ""),
            "label_seq_id": row.get("label_seq_id", ""),
            "polymer_flag": row.get("polymer_flag", ""),
            "occupancy_flag": row.get("occupancy_flag", ""),
            "configured_chain": str(
                row.get("auth_asym_id", "") == definition.author_chain_id
            ).lower(),
        }
        for row in category_rows(block, "pdbx_unobs_or_zero_occ_residues")
    ]


def _difference_class(row: Mapping[str, str]) -> str:
    details = row.get("details", "").lower()
    if "expression tag" in details:
        return "expression_tag"
    pdb_residue = row.get("mon_id", "")
    reference_residue = row.get("db_mon_id", "")
    if pdb_residue and reference_residue and pdb_residue != reference_residue:
        return "sequence_substitution"
    return "other_sequence_difference"


def difference_rows(
    *, block: gemmi.cif.Block, definition: ReceptorDefinition
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for source in category_rows(block, "struct_ref_seq_dif"):
        auth_chain = source.get("pdbx_pdb_strand_id", "") or source.get(
            "pdbx_auth_asym_id", ""
        )
        rows.append(
            {
                "receptor_id": definition.receptor_id,
                "pdb_id": definition.pdb_id,
                "auth_chain_id": auth_chain,
                "auth_seq_id": source.get("pdbx_auth_seq_num", ""),
                "pdb_residue": source.get("mon_id", ""),
                "reference_residue": source.get("db_mon_id", ""),
                "reference_seq_id": source.get("pdbx_seq_db_seq_num", ""),
                "details": source.get("details", ""),
                "difference_class": _difference_class(source),
                "configured_chain": str(
                    auth_chain == definition.author_chain_id
                ).lower(),
            }
        )
    return rows
