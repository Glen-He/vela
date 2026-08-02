"""结构审计阶段内部传递的稳定数据模型。"""

from dataclasses import dataclass


@dataclass(frozen=True, slots=True)
class AuditResult:
    """一个受体登记项的审计状态。"""

    receptor_id: str
    pdb_id: str
    identity_status: str
    manual_review_required: bool


@dataclass(frozen=True, slots=True)
class EntryAudit:
    """一个结构条目的全部表格行。"""

    summary: dict[str, str]
    chains: list[dict[str, str]]
    components: list[dict[str, str]]
    missing: list[dict[str, str]]
    differences: list[dict[str, str]]
    result: AuditResult


@dataclass(frozen=True, slots=True)
class CrystalContacts:
    """配置链参与的晶体接触计数。"""

    asu_other_polymer_atom_pairs: int
    asu_target_residues: int
    symmetry_polymer_atom_pairs: int
    symmetry_target_residues: int
    minimum_symmetry_distance: float | None
