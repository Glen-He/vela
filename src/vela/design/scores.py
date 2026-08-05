"""阶段四 Rosetta 分数字段的唯一语义映射和严格提取。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.core.provenance import atomic_write_text, sha256_file
from vela.design.models import DesignError
from vela.validation.scores import RosettaScoreRow

SCREEN_SCORE_COLUMNS = {
    "dG_separated": "dG_separated",
    "dSASA_int_A2": "dSASA_int",
    "dG_separated_per_dSASAx100": "dG_separated/dSASAx100",
    "delta_unsat_hbonds": "delta_unsatHbonds",
    "interface_hbonds": "hbonds_int",
    "shape_complementarity": "sc_value",
    "interface_residue_count": "nres_int",
    "receptor_separated_score": "side1_score",
    "peptide_separated_score": "side2_score",
}

FINALIST_SCORE_COLUMNS = {
    "reweighted_sc": "reweighted_sc",
    "I_sc": "I_sc",
    "pep_sc": "pep_sc",
    "pep_sc_noref": "pep_sc_noref",
    "I_bsa_A2": "I_bsa",
    "I_hb": "I_hb",
    "I_pack": "I_pack",
    "I_unsat": "I_unsat",
    "disulfide_score": "dslf_fa13",
}


@dataclass(frozen=True, slots=True)
class ScreenMetrics:
    """固定骨架 InterfaceAnalyzer 的必需原始指标。"""

    dG_separated: float
    dSASA_int_A2: float
    dG_separated_per_dSASAx100: float
    delta_unsat_hbonds: float
    interface_hbonds: float
    shape_complementarity: float
    interface_residue_count: float
    receptor_separated_score: float
    peptide_separated_score: float

    def as_dict(self) -> dict[str, float]:
        """返回可直接序列化且字段名安全的指标。"""
        return {name: getattr(self, name) for name in SCREEN_SCORE_COLUMNS}


@dataclass(frozen=True, slots=True)
class FinalistMetrics:
    """FlexPepDock 已产生且阶段四合同要求保留的指标。"""

    reweighted_sc: float
    I_sc: float
    pep_sc: float
    pep_sc_noref: float
    I_bsa_A2: float
    I_hb: float
    I_pack: float
    I_unsat: float
    disulfide_score: float

    def as_dict(self) -> dict[str, float]:
        """返回可直接序列化的全部 finalist 指标。"""
        return {name: getattr(self, name) for name in FINALIST_SCORE_COLUMNS}


def screen_metrics(row: RosettaScoreRow) -> ScreenMetrics:
    """按冻结列名严格读取 InterfaceAnalyzer 指标。"""
    return ScreenMetrics(
        dG_separated=row.score("dG_separated"),
        dSASA_int_A2=row.score("dSASA_int"),
        dG_separated_per_dSASAx100=row.score("dG_separated/dSASAx100"),
        delta_unsat_hbonds=row.score("delta_unsatHbonds"),
        interface_hbonds=row.score("hbonds_int"),
        shape_complementarity=row.score("sc_value"),
        interface_residue_count=row.score("nres_int"),
        receptor_separated_score=row.score("side1_score"),
        peptide_separated_score=row.score("side2_score"),
    )


def finalist_metrics(row: RosettaScoreRow) -> FinalistMetrics:
    """按冻结列名严格读取 FlexPepDock 指标。"""
    return FinalistMetrics(
        reweighted_sc=row.score("reweighted_sc"),
        I_sc=row.score("I_sc"),
        pep_sc=row.score("pep_sc"),
        pep_sc_noref=row.score("pep_sc_noref"),
        I_bsa_A2=row.score("I_bsa"),
        I_hb=row.score("I_hb"),
        I_pack=row.score("I_pack"),
        I_unsat=row.score("I_unsat"),
        disulfide_score=row.score("dslf_fa13"),
    )


def write_finalist_decoy_manifest(
    *,
    rows: tuple[RosettaScoreRow, ...],
    decoy_paths: dict[str, Path],
    task_dir: Path,
) -> Path:
    """无损记录 finalist 合同要求的分数、结构身份和来源。"""
    if set(decoy_paths) != {row.description for row in rows}:
        raise DesignError("finalist score/PDB identities differ")
    metric_names = tuple(FINALIST_SCORE_COLUMNS)
    lines = ["\t".join(("description", *metric_names, "path", "sha256")) + "\n"]
    for row in rows:
        path = decoy_paths[row.description]
        metrics = finalist_metrics(row).as_dict()
        lines.append(
            "\t".join(
                (
                    row.description,
                    *(f"{metrics[name]:.6f}" for name in metric_names),
                    path.relative_to(task_dir).as_posix(),
                    sha256_file(path),
                )
            )
            + "\n"
        )
    path = task_dir / "decoy_manifest.tsv"
    atomic_write_text(path, "".join(lines))
    return path
