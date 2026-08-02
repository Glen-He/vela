"""阶段四配体变体合同、单点库和有限组合提案。"""

from __future__ import annotations

import csv
import io
from collections import defaultdict
from collections.abc import Iterable
from dataclasses import dataclass
from hashlib import sha256
from itertools import combinations, product

from vela.design.models import (
    CandidateParent,
    DesignError,
    DesignSettings,
    SequenceCandidate,
)
from vela.preparation.chemistry import STANDARD_AMINO_ACIDS, ChemistryDefinition

POSITIVE_RESIDUES = frozenset({"K", "R"})
NEGATIVE_RESIDUES = frozenset({"D", "E"})
HYDROPHOBIC_RESIDUES = frozenset({"A", "F", "I", "L", "M", "V", "W", "Y"})
AROMATIC_RESIDUES = frozenset({"F", "W", "Y"})


@dataclass(frozen=True, slots=True)
class SequenceFacts:
    """不作为自动药物性结论的可核查序列描述量。"""

    hydrophobic_count: int
    max_hydrophobic_run: int
    aromatic_count: int
    basic_count: int
    acidic_count: int
    glycine_proline_count: int
    methionine_count: int
    deamidation_motif_count: int
    aspartimide_motif_count: int


def sequence_facts(sequence: str) -> SequenceFacts:
    """计算透明序列计数; 用于后续实验开发性审查而非硬性淘汰。"""
    runs: list[int] = []
    current = 0
    for residue in sequence:
        if residue in HYDROPHOBIC_RESIDUES:
            current += 1
            runs.append(current)
        else:
            current = 0
    return SequenceFacts(
        hydrophobic_count=sum(item in HYDROPHOBIC_RESIDUES for item in sequence),
        max_hydrophobic_run=max(runs, default=0),
        aromatic_count=sum(item in AROMATIC_RESIDUES for item in sequence),
        basic_count=sum(item in POSITIVE_RESIDUES for item in sequence),
        acidic_count=sum(item in NEGATIVE_RESIDUES for item in sequence),
        glycine_proline_count=sum(item in {"G", "P"} for item in sequence),
        methionine_count=sequence.count("M"),
        deamidation_motif_count=sum(
            sequence[index] == "N" and sequence[index + 1] in {"G", "S"}
            for index in range(len(sequence) - 1)
        ),
        aspartimide_motif_count=sum(
            sequence[index] == "D" and sequence[index + 1] in {"G", "S", "N"}
            for index in range(len(sequence) - 1)
        ),
    )


def mutation_string(*, reference: str, sequence: str) -> str:
    """返回按配体位置排序的稳定突变描述。"""
    if len(reference) != len(sequence):
        raise DesignError("candidate sequence length differs from the ligand")
    return ";".join(
        f"{wild_type}{position}{candidate}"
        for position, (wild_type, candidate) in enumerate(
            zip(reference, sequence, strict=True), 1
        )
        if wild_type != candidate
    )


def sequence_id(sequence: str) -> str:
    """根据完整序列建立跨模板稳定的候选身份。"""
    return f"seq_{sha256(sequence.encode('ascii')).hexdigest()[:12]}"


def _histidine_charge(state: str) -> int:
    return 1 if state == "HIP" else 0


def _side_chain_charge(
    *, sequence: str, histidine_states: dict[int, str], fallback_histidine_state: str
) -> int:
    charge = 0
    for position, residue in enumerate(sequence, 1):
        if residue in POSITIVE_RESIDUES:
            charge += 1
        elif residue in NEGATIVE_RESIDUES:
            charge -= 1
        elif residue == "H":
            charge += _histidine_charge(
                histidine_states.get(position, fallback_histidine_state)
            )
    return charge


def candidate_net_charge(
    *, chemistry: ChemistryDefinition, settings: DesignSettings, sequence: str
) -> int:
    """在端基不变的前提下由 WT 形式净电荷推导候选形式净电荷。"""
    if chemistry.net_charge is None:
        raise DesignError("ligand net charge must be resolved before sequence design")
    wt_histidines = {item.position: item.state for item in chemistry.histidines}
    wt_side = _side_chain_charge(
        sequence=chemistry.sequence,
        histidine_states=wt_histidines,
        fallback_histidine_state=settings.sequence.candidate_histidine_state,
    )
    candidate_side = _side_chain_charge(
        sequence=sequence,
        histidine_states={},
        fallback_histidine_state=settings.sequence.candidate_histidine_state,
    )
    return chemistry.net_charge - wt_side + candidate_side


def build_candidate(
    *,
    chemistry: ChemistryDefinition,
    settings: DesignSettings,
    sequence: str,
    design_round: str,
    generation: int,
    parents: tuple[CandidateParent, ...],
    ancestor_candidate_ids: tuple[str, ...],
    proposal_source: str,
) -> SequenceCandidate:
    """校验固定拓扑、变化范围和轮次后建立候选。"""
    reference = chemistry.sequence
    if len(sequence) != len(reference):
        raise DesignError("candidate sequence length differs from the ligand")
    if set(sequence) - STANDARD_AMINO_ACIDS:
        raise DesignError("candidate sequence contains a non-standard amino acid")
    disulfide_positions = {
        position
        for bond in chemistry.disulfide_bonds
        for position in (bond.first, bond.second)
    }
    if any(sequence[position - 1] != "C" for position in disulfide_positions):
        raise DesignError("candidate changes a disulfide cysteine")
    if sequence.count("C") != len(disulfide_positions):
        raise DesignError("candidate introduces an additional cysteine")
    mutations = tuple(
        position
        for position, (wild_type, candidate) in enumerate(
            zip(reference, sequence, strict=True), 1
        )
        if wild_type != candidate
    )
    if not mutations:
        raise DesignError("candidate sequence is identical to the ligand")
    if any(
        position not in settings.sequence.mutable_positions for position in mutations
    ):
        raise DesignError("candidate changes a fixed ligand position")
    if any(
        sequence[position - 1] not in settings.sequence.allowed_amino_acids
        for position in mutations
    ):
        raise DesignError("candidate uses an amino acid outside the allowed design set")
    if design_round == "single" and len(mutations) != 1:
        raise DesignError("single design candidate must contain one mutation")
    if design_round == "combination" and not (
        settings.combination.min_mutations
        <= len(mutations)
        <= settings.combination.max_mutations
    ):
        raise DesignError("combination candidate mutation count is outside its policy")
    if (
        design_round == "iteration"
        and len(mutations) > settings.iteration.max_total_mutations
    ):
        raise DesignError("iteration candidate exceeds its total mutation limit")
    notation = mutation_string(reference=reference, sequence=sequence)
    return SequenceCandidate(
        candidate_id=sequence_id(sequence),
        sequence=sequence,
        mutation_string=notation,
        mutation_positions=mutations,
        mutation_count=len(mutations),
        net_charge=candidate_net_charge(
            chemistry=chemistry, settings=settings, sequence=sequence
        ),
        design_round=design_round,
        generation=generation,
        parents=parents,
        ancestor_candidate_ids=ancestor_candidate_ids,
        proposal_source=proposal_source,
    )


def first_generation_candidate(
    *,
    chemistry: ChemistryDefinition,
    settings: DesignSettings,
    sequence: str,
    design_round: str,
    proposal_source: str,
) -> SequenceCandidate:
    """建立以原始配体为父序列的第一代单点或组合候选。"""
    edit = mutation_string(reference=chemistry.sequence, sequence=sequence)
    return build_candidate(
        chemistry=chemistry,
        settings=settings,
        sequence=sequence,
        design_round=design_round,
        generation=1,
        parents=(CandidateParent("ligand_wt", edit, "origin"),),
        ancestor_candidate_ids=("ligand_wt",),
        proposal_source=proposal_source,
    )


def systematic_single_library(
    *, chemistry: ChemistryDefinition, settings: DesignSettings
) -> tuple[SequenceCandidate, ...]:
    """完整枚举全部配置允许的单点替换, 不依赖生成模型。"""
    reference = chemistry.sequence
    candidates = [
        first_generation_candidate(
            chemistry=chemistry,
            settings=settings,
            sequence=(reference[: position - 1] + residue + reference[position:]),
            design_round="single",
            proposal_source="systematic_single_scan",
        )
        for position in settings.sequence.mutable_positions
        for residue in settings.sequence.allowed_amino_acids
        if residue != reference[position - 1]
    ]
    if not candidates or len({item.candidate_id for item in candidates}) != len(
        candidates
    ):
        raise DesignError("systematic single-mutation library is empty or duplicated")
    return tuple(
        sorted(candidates, key=lambda item: (item.mutation_positions, item.sequence))
    )


def combination_library(
    *,
    chemistry: ChemistryDefinition,
    settings: DesignSettings,
    substitutions: dict[int, tuple[str, ...]],
) -> tuple[SequenceCandidate, ...]:
    """从已筛选单点替换提出有限组合, 不把单点分数当组合结论。"""
    positions = tuple(sorted(substitutions))
    if any(
        not residues
        or len(residues) != len(set(residues))
        or position not in settings.sequence.mutable_positions
        for position, residues in substitutions.items()
    ):
        raise DesignError("combination substitution map is invalid")
    candidates: list[SequenceCandidate] = []
    for mutation_count in range(
        settings.combination.min_mutations,
        settings.combination.max_mutations + 1,
    ):
        for selected_positions in combinations(positions, mutation_count):
            residue_sets = tuple(
                substitutions[position] for position in selected_positions
            )
            for residues in product(*residue_sets):
                sequence = list(chemistry.sequence)
                for position, residue in zip(selected_positions, residues, strict=True):
                    sequence[position - 1] = residue
                candidates.append(
                    first_generation_candidate(
                        chemistry=chemistry,
                        settings=settings,
                        sequence="".join(sequence),
                        design_round="combination",
                        proposal_source="cross_template_single_evidence",
                    )
                )
    unique = {item.candidate_id: item for item in candidates}
    ordered = tuple(
        sorted(
            unique.values(),
            key=lambda item: (
                item.mutation_count,
                item.mutation_positions,
                item.sequence,
            ),
        )
    )
    if len(ordered) <= settings.combination.max_candidates:
        return ordered
    by_count: dict[int, list[SequenceCandidate]] = {}
    for candidate in ordered:
        by_count.setdefault(candidate.mutation_count, []).append(candidate)
    counts = tuple(sorted(by_count))
    base, remainder = divmod(settings.combination.max_candidates, len(counts))
    selected: list[SequenceCandidate] = []
    for index, mutation_count in enumerate(counts):
        budget = base + (1 if index < remainder else 0)
        selected.extend(
            _balanced_position_selection(by_count[mutation_count], budget=budget)
        )
    return tuple(
        sorted(
            selected,
            key=lambda item: (
                item.mutation_count,
                item.mutation_positions,
                item.sequence,
            ),
        )
    )


def _balanced_position_selection(
    candidates: list[SequenceCandidate], *, budget: int
) -> tuple[SequenceCandidate, ...]:
    """在资源截断时平衡位置组合覆盖, 避免按序号偏向肽 N 端。"""
    buckets: dict[tuple[int, ...], list[SequenceCandidate]] = {}
    for candidate in candidates:
        buckets.setdefault(candidate.mutation_positions, []).append(candidate)
    for values in buckets.values():
        values.sort(key=lambda item: item.sequence)
    coverage: dict[int, int] = defaultdict(int)
    selected: list[SequenceCandidate] = []
    while len(selected) < budget:
        available = [positions for positions, values in buckets.items() if values]
        if not available:
            break
        positions = min(
            available,
            key=lambda item: (
                max(coverage[position] for position in item),
                sum(coverage[position] for position in item),
                tuple(coverage[position] for position in item),
                item,
            ),
        )
        selected.append(buckets[positions].pop(0))
        for position in positions:
            coverage[position] += 1
    return tuple(selected)


def candidate_table(candidates: Iterable[SequenceCandidate]) -> str:
    """建立稳定、可审计的候选 TSV。"""
    fields = (
        "candidate_id",
        "sequence",
        "mutation_string",
        "mutation_positions",
        "mutation_count",
        "net_charge",
        "design_round",
        "generation",
        "parent_candidate_ids",
        "parent_edits",
        "ancestor_candidate_ids",
        "proposal_source",
        "hydrophobic_count",
        "max_hydrophobic_run",
        "aromatic_count",
        "methionine_count",
        "deamidation_motif_count",
        "aspartimide_motif_count",
    )
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for candidate in candidates:
        facts = sequence_facts(candidate.sequence)
        writer.writerow(
            {
                "candidate_id": candidate.candidate_id,
                "sequence": candidate.sequence,
                "mutation_string": candidate.mutation_string,
                "mutation_positions": ";".join(map(str, candidate.mutation_positions)),
                "mutation_count": str(candidate.mutation_count),
                "net_charge": str(candidate.net_charge),
                "design_round": candidate.design_round,
                "generation": str(candidate.generation),
                "parent_candidate_ids": ";".join(
                    item.candidate_id for item in candidate.parents
                ),
                "parent_edits": ";".join(
                    f"{item.candidate_id}:{item.edit}:{item.edit_type}"
                    for item in candidate.parents
                ),
                "ancestor_candidate_ids": ";".join(candidate.ancestor_candidate_ids),
                "proposal_source": candidate.proposal_source,
                "hydrophobic_count": str(facts.hydrophobic_count),
                "max_hydrophobic_run": str(facts.max_hydrophobic_run),
                "aromatic_count": str(facts.aromatic_count),
                "methionine_count": str(facts.methionine_count),
                "deamidation_motif_count": str(facts.deamidation_motif_count),
                "aspartimide_motif_count": str(facts.aspartimide_motif_count),
            }
        )
    return output.getvalue()
