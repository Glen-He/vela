"""阶段四父序列一步邻域、谱系和资源平衡。"""

from __future__ import annotations

import csv
import io
from collections import defaultdict
from dataclasses import dataclass

from vela.design.models import (
    CandidateParent,
    DesignError,
    DesignSettings,
    SequenceCandidate,
)
from vela.design.sequence.library import build_candidate, candidate_table, sequence_id
from vela.preparation.chemistry import ChemistryDefinition


@dataclass(frozen=True, slots=True)
class IterationLibrary:
    """完整一步邻域及资源预算内实际进入初筛的候选。"""

    selected: tuple[SequenceCandidate, ...]
    deferred: tuple[SequenceCandidate, ...]

    @property
    def all_candidates(self) -> tuple[SequenceCandidate, ...]:
        """返回按稳定身份排序的全部已选和延期提案。"""
        return tuple(
            sorted(
                (*self.selected, *self.deferred),
                key=lambda item: item.candidate_id,
            )
        )


def candidate_parent_edit(
    *, reference: str, parent: SequenceCandidate, child_sequence: str
) -> CandidateParent:
    """建立并分类父序列到一步后代的唯一编辑。"""
    differences = [
        position
        for position, (old, new) in enumerate(
            zip(parent.sequence, child_sequence, strict=True), 1
        )
        if old != new
    ]
    if len(differences) != 1:
        raise DesignError("iteration child must be one edit from its parent")
    position = differences[0]
    old = parent.sequence[position - 1]
    new = child_sequence[position - 1]
    if old == reference[position - 1]:
        edit_type = "add"
    elif new == reference[position - 1]:
        edit_type = "revert"
    else:
        edit_type = "swap"
    return CandidateParent(parent.candidate_id, f"{old}{position}{new}", edit_type)


def _iteration_dimensions(
    candidate: SequenceCandidate,
) -> tuple[tuple[str, ...], tuple[str, ...], tuple[str, ...]]:
    parents: set[str] = set()
    positions: set[str] = set()
    edit_types: set[str] = set()
    for parent in candidate.parents:
        position = "".join(
            character for character in parent.edit if character.isdigit()
        )
        parents.add(parent.candidate_id)
        positions.add(position)
        edit_types.add(parent.edit_type)
    return tuple(sorted(parents)), tuple(sorted(positions)), tuple(sorted(edit_types))


def _balanced_iteration_selection(
    candidates: tuple[SequenceCandidate, ...], *, budget: int
) -> tuple[SequenceCandidate, ...]:
    """在父序列、位置和编辑类型间平衡有限的一步邻域预算。"""
    remaining = list(candidates)
    parent_coverage: dict[str, int] = defaultdict(int)
    position_coverage: dict[str, int] = defaultdict(int)
    type_coverage: dict[str, int] = defaultdict(int)
    selected: list[SequenceCandidate] = []
    while remaining and len(selected) < budget:

        def priority(
            item: SequenceCandidate,
        ) -> tuple[int, int, int, int, int, int, int, tuple[int, ...], str]:
            parents, positions, edit_types = _iteration_dimensions(item)
            return (
                min(parent_coverage[value] for value in parents),
                sum(parent_coverage[value] for value in parents),
                min(position_coverage[value] for value in positions),
                sum(position_coverage[value] for value in positions),
                min(type_coverage[value] for value in edit_types),
                sum(type_coverage[value] for value in edit_types),
                item.mutation_count,
                item.mutation_positions,
                item.sequence,
            )

        candidate = min(remaining, key=priority)
        remaining.remove(candidate)
        selected.append(candidate)
        parents, positions, edit_types = _iteration_dimensions(candidate)
        for value in parents:
            parent_coverage[value] += 1
        for value in positions:
            position_coverage[value] += 1
        for value in edit_types:
            type_coverage[value] += 1
    return tuple(selected)


def iteration_library(
    *,
    chemistry: ChemistryDefinition,
    settings: DesignSettings,
    parents: tuple[SequenceCandidate, ...],
) -> IterationLibrary:
    """生成多个同代父序列的完整一步邻域并执行可审计资源截断。"""
    if (
        not parents
        or len(parents) > settings.iteration.max_parents
        or len({item.candidate_id for item in parents}) != len(parents)
        or len({item.generation for item in parents}) != 1
    ):
        raise DesignError(
            "iteration parents must be unique, same-generation candidates"
        )
    parent_sequences = {item.sequence for item in parents}
    proposals: dict[str, dict[str, CandidateParent]] = defaultdict(dict)
    ancestors: dict[str, set[str]] = defaultdict(set)
    for parent in sorted(parents, key=lambda item: item.candidate_id):
        for position in settings.sequence.mutable_positions:
            for residue in settings.sequence.allowed_amino_acids:
                if residue == parent.sequence[position - 1]:
                    continue
                sequence = (
                    parent.sequence[: position - 1]
                    + residue
                    + parent.sequence[position:]
                )
                if sequence == chemistry.sequence or sequence in parent_sequences:
                    continue
                total_mutations = sum(
                    old != new
                    for old, new in zip(chemistry.sequence, sequence, strict=True)
                )
                if total_mutations > settings.iteration.max_total_mutations:
                    continue
                candidate_id = sequence_id(sequence)
                if candidate_id in parent.ancestor_candidate_ids:
                    continue
                proposals[sequence][parent.candidate_id] = candidate_parent_edit(
                    reference=chemistry.sequence,
                    parent=parent,
                    child_sequence=sequence,
                )
                ancestors[sequence].update(parent.ancestor_candidate_ids)
                ancestors[sequence].add(parent.candidate_id)
    candidates = tuple(
        sorted(
            (
                build_candidate(
                    chemistry=chemistry,
                    settings=settings,
                    sequence=sequence,
                    design_round="iteration",
                    generation=parents[0].generation + 1,
                    parents=tuple(
                        sorted(values.values(), key=lambda item: item.candidate_id)
                    ),
                    ancestor_candidate_ids=tuple(sorted(ancestors[sequence])),
                    proposal_source="parent_one_edit_neighborhood",
                )
                for sequence, values in proposals.items()
            ),
            key=lambda item: (
                item.mutation_count,
                item.mutation_positions,
                item.sequence,
            ),
        )
    )
    if not candidates:
        raise DesignError("iteration parents produced no legal one-edit candidates")
    selected = _balanced_iteration_selection(
        candidates,
        budget=min(settings.iteration.max_candidates, len(candidates)),
    )
    selected_ids = {item.candidate_id for item in selected}
    deferred = tuple(
        item for item in candidates if item.candidate_id not in selected_ids
    )
    return IterationLibrary(selected=selected, deferred=deferred)


def iteration_candidate_table(library: IterationLibrary) -> str:
    """记录完整邻域中实际进入初筛和仅因资源延期的提案。"""
    base = candidate_table(library.all_candidates)
    rows = list(csv.DictReader(io.StringIO(base), delimiter="\t"))
    if not rows:
        raise DesignError("iteration candidate table cannot be empty")
    selected_ids = {item.candidate_id for item in library.selected}
    fields = (*tuple(rows[0]), "selection_status")
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    for row in rows:
        row["selection_status"] = (
            "selected" if row["candidate_id"] in selected_ids else "resource_deferred"
        )
        writer.writerow(row)
    return output.getvalue()
