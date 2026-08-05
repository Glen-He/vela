"""阶段四候选柔性复核报表与阶段五输入队列。"""

from __future__ import annotations

import csv
import io
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.design.finalists.evidence import (
    CandidateEvidence,
    FinalistDecoy,
    PairedSummary,
    StateSummary,
    candidate_evidence,
    manifest_results,
    paired_summaries,
    state_summary,
    task_decoys,
)
from vela.design.finalists.records import FinalistPlan, read_finalist_plan
from vela.design.models import DesignError
from vela.design.sequence.library import sequence_facts


@dataclass(frozen=True, slots=True)
class FinalistAnalysisOutcome:
    """柔性复核报告和阶段五队列规模。"""

    manifest_path: Path
    candidate_count: int
    eligible_count: int
    md_system_count: int


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _number(value: float | None) -> str:
    return "" if value is None else f"{value:.6f}"


def _representative(
    decoys: list[FinalistDecoy], *, ranking_median: float
) -> FinalistDecoy:
    passed = [item for item in decoys if item.geometry.passed]
    if not passed:
        raise DesignError("cannot select a representative from failed decoys")
    return min(
        passed,
        key=lambda item: (
            abs(item.metrics.reweighted_sc - ranking_median),
            item.metrics.I_sc,
            item.decoy_id,
        ),
    )


def _dominates(first: CandidateEvidence, second: CandidateEvidence) -> bool:
    first_values = (
        first.positive_ranking_median,
        first.positive_peptide_median,
        -first.successful_pair_fraction,
    )
    second_values = (
        second.positive_ranking_median,
        second.positive_peptide_median,
        -second.successful_pair_fraction,
    )
    if any(value is None for value in (*first_values, *second_values)):
        raise DesignError("eligible finalist lacks a Pareto objective")
    return all(
        left <= right
        for left, right in zip(first_values, second_values, strict=True)
        if left is not None and right is not None
    ) and any(
        left < right
        for left, right in zip(first_values, second_values, strict=True)
        if left is not None and right is not None
    )


def _pareto_resource_selection(
    candidates: list[CandidateEvidence], *, budget: int
) -> list[CandidateEvidence]:
    """逐层选择非支配候选; 同层按位置覆盖和稳定身份截断。"""
    remaining = list(candidates)
    selected: list[CandidateEvidence] = []
    position_coverage: dict[int, int] = defaultdict(int)
    while remaining and len(selected) < budget:
        front = [
            candidate
            for candidate in remaining
            if not any(
                _dominates(other, candidate)
                for other in remaining
                if other is not candidate
            )
        ]
        front.sort(
            key=lambda item: (
                max(
                    (
                        position_coverage[position]
                        for position in item.candidate.mutation_positions
                    ),
                    default=0,
                ),
                sum(
                    position_coverage[position]
                    for position in item.candidate.mutation_positions
                ),
                item.candidate.mutation_positions,
                item.candidate.candidate_id,
            )
        )
        for candidate in front[: budget - len(selected)]:
            selected.append(candidate)
            for position in candidate.candidate.mutation_positions:
                position_coverage[position] += 1
        front_ids = {item.candidate.candidate_id for item in front}
        remaining = [
            item for item in remaining if item.candidate.candidate_id not in front_ids
        ]
    return selected


def _md_queue(
    *,
    config: AppConfig,
    plan: FinalistPlan,
    evidence: tuple[CandidateEvidence, ...],
    states: tuple[StateSummary, ...],
    pairs: tuple[PairedSummary, ...],
) -> tuple[list[dict[str, str]], dict[str, str]]:
    eligible = [item for item in evidence if item.status == "eligible"]
    selected = _pareto_resource_selection(
        eligible, budget=config.design.finalists.max_md_candidates
    )
    selected_ids = {item.candidate.candidate_id for item in selected}
    resource_status = {
        item.candidate.candidate_id: (
            "md_selected"
            if item.candidate.candidate_id in selected_ids
            else "md_resource_deferred"
        )
        for item in eligible
    }
    cluster_by_id = {
        cluster.cluster_id: cluster
        for state in states
        for cluster in state.pose_clusters
    }
    matched_clusters: dict[tuple[str, str], list[FinalistDecoy]] = defaultdict(list)
    for pair in pairs:
        if pair.mutant_pose_cluster_id is None:
            continue
        cluster = cluster_by_id.get(pair.mutant_pose_cluster_id)
        if cluster is None:
            raise DesignError("matched finalist pose cluster is missing")
        matched_clusters[(pair.candidate.candidate_id, pair.template_id)].append(
            cluster.medoid
        )
    rows: list[dict[str, str]] = []
    for template in plan.templates:
        rows.append(
            {
                "system_id": f"md_wt__{template.template_id}",
                "candidate_id": "ligand_wt",
                "sequence": config.chemistry.sequence,
                "mutation_string": "WT",
                "generation": "0",
                "parent_candidate_ids": "",
                "system_role": "wild_type_control",
                "template_id": template.template_id,
                "target": template.target,
                "evidence_role": template.evidence_role,
                "structure_path": template.path.relative_to(plan.run_dir).as_posix(),
                "structure_sha256": template.sha256,
            }
        )
    for item in selected:
        for template in plan.templates:
            passed = matched_clusters[
                (item.candidate.candidate_id, template.template_id)
            ]
            if not passed:
                raise DesignError("MD-selected candidate lacks a matched template pose")
            median = float(
                statistics.median(member.metrics.reweighted_sc for member in passed)
            )
            representative = _representative(passed, ranking_median=median)
            rows.append(
                {
                    "system_id": (
                        f"md_{item.candidate.candidate_id}__{template.template_id}"
                    ),
                    "candidate_id": item.candidate.candidate_id,
                    "sequence": item.candidate.sequence,
                    "mutation_string": item.candidate.mutation_string,
                    "generation": str(item.candidate.generation),
                    "parent_candidate_ids": ";".join(
                        parent.candidate_id for parent in item.candidate.parents
                    ),
                    "system_role": "optimized_candidate",
                    "template_id": template.template_id,
                    "target": template.target,
                    "evidence_role": template.evidence_role,
                    "structure_path": representative.path.relative_to(
                        plan.run_dir
                    ).as_posix(),
                    "structure_sha256": representative.sha256,
                }
            )
    return rows, resource_status


def _decoy_rows(
    *, decoys: tuple[FinalistDecoy, ...], run_dir: Path
) -> list[dict[str, str]]:
    return [
        {
            "decoy_id": item.decoy_id,
            "task_id": item.task.task_id,
            "pair_id": item.task.pair_id,
            "state": item.task.start.state,
            "candidate_id": item.task.start.candidate.candidate_id,
            "template_id": item.task.start.template.template_id,
            "target": item.task.start.template.target,
            "evidence_role": item.task.start.template.evidence_role,
            "seed": str(item.task.seed),
            **{name: f"{value:.6f}" for name, value in item.metrics.as_dict().items()},
            "qc_status": "passed" if item.geometry.passed else "failed",
            "interface_contact_pairs": str(item.geometry.interface_contact_pairs),
            "interface_receptor_residues": str(
                item.geometry.interface_receptor_residues
            ),
            "minimum_interface_distance_A": (
                f"{item.geometry.minimum_interface_distance_A:.6f}"
            ),
            "receptor_ca_rmsd_A": f"{item.geometry.receptor_ca_rmsd_A:.6f}",
            "start_contact_overlap": f"{item.geometry.start_contact_overlap:.6f}",
            "start_site_displacement_A": (
                f"{item.geometry.start_site_displacement_A:.6f}"
            ),
            "path": item.path.relative_to(run_dir).as_posix(),
            "sha256": item.sha256,
        }
        for item in decoys
    ]


def _pair_rows(pairs: tuple[PairedSummary, ...]) -> list[dict[str, str]]:
    return [
        {
            "pair_id": item.pair_id,
            "candidate_id": item.candidate.candidate_id,
            "template_id": item.template_id,
            "target": item.target,
            "evidence_role": item.evidence_role,
            "seed": str(item.seed),
            "wt_passed_fraction": f"{item.wt_passed_fraction:.6f}",
            "mutant_passed_fraction": f"{item.mutant_passed_fraction:.6f}",
            "successful": str(item.successful).lower(),
            "pose_relation": item.pose_relation,
            "wt_pose_cluster_id": item.wt_pose_cluster_id or "",
            "mutant_pose_cluster_id": item.mutant_pose_cluster_id or "",
            "paired_reweighted_sc_delta_REU": _number(item.paired_reweighted_sc_delta),
            "paired_I_sc_delta_REU": _number(item.paired_I_sc_delta),
            "paired_pep_sc_delta_REU": _number(item.paired_pep_sc_delta),
        }
        for item in pairs
    ]


def _candidate_rows(
    *, evidence: tuple[CandidateEvidence, ...], resource_status: dict[str, str]
) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for item in evidence:
        facts = sequence_facts(item.candidate.sequence)
        rows.append(
            {
                "candidate_id": item.candidate.candidate_id,
                "sequence": item.candidate.sequence,
                "mutation_string": item.candidate.mutation_string,
                "mutation_count": str(item.candidate.mutation_count),
                "generation": str(item.candidate.generation),
                "parent_candidate_ids": ";".join(
                    parent.candidate_id for parent in item.candidate.parents
                ),
                "parent_edits": ";".join(
                    f"{parent.candidate_id}:{parent.edit}:{parent.edit_type}"
                    for parent in item.candidate.parents
                ),
                "net_charge": str(item.candidate.net_charge),
                "hydrophobic_count": str(facts.hydrophobic_count),
                "max_hydrophobic_run": str(facts.max_hydrophobic_run),
                "aromatic_count": str(facts.aromatic_count),
                "methionine_count": str(facts.methionine_count),
                "deamidation_motif_count": str(facts.deamidation_motif_count),
                "aspartimide_motif_count": str(facts.aspartimide_motif_count),
                "positive_ranking_median": _number(item.positive_ranking_median),
                "positive_ranking_worst": _number(item.positive_ranking_worst),
                "positive_interface_median": _number(item.positive_interface_median),
                "positive_interface_worst": _number(item.positive_interface_worst),
                "positive_peptide_median": _number(item.positive_peptide_median),
                "successful_pair_fraction": f"{item.successful_pair_fraction:.6f}",
                "candidate_status": item.status,
                "failed_gates": ";".join(item.failed_gates),
                "resource_status": resource_status.get(
                    item.candidate.candidate_id, "not_eligible_for_md"
                ),
            }
        )
    return rows


def _pose_cluster_rows(states: tuple[StateSummary, ...]) -> list[dict[str, str]]:
    return [
        {
            "pose_cluster_id": cluster.cluster_id,
            "task_id": state.task.task_id,
            "state": state.task.start.state,
            "candidate_id": state.task.start.candidate.candidate_id,
            "template_id": state.task.start.template.template_id,
            "seed": str(state.task.seed),
            "cluster_size": str(len(cluster.members)),
            "cluster_fraction_of_qc_passed": (
                f"{len(cluster.members) / config_count:.6f}"
            ),
            "medoid_decoy_id": cluster.medoid.decoy_id,
            "median_reweighted_sc_REU": f"{cluster.ranking_median:.6f}",
            "median_I_sc_REU": f"{cluster.interface_median:.6f}",
            "median_pep_sc_REU": f"{cluster.peptide_median:.6f}",
            "primary": str(cluster.cluster_id == state.primary_pose_cluster_id).lower(),
        }
        for state in states
        for config_count in (sum(len(item.members) for item in state.pose_clusters),)
        for cluster in state.pose_clusters
    ]


def _write_reports(
    *,
    run_dir: Path,
    output_dir: Path,
    decoys: tuple[FinalistDecoy, ...],
    states: tuple[StateSummary, ...],
    pairs: tuple[PairedSummary, ...],
    evidence: tuple[CandidateEvidence, ...],
    queue_rows: list[dict[str, str]],
    resource_status: dict[str, str],
) -> tuple[Path, Path, Path, Path, Path]:
    decoy_path = output_dir / "decoys.tsv"
    cluster_path = output_dir / "pose_clusters.tsv"
    pair_path = output_dir / "paired_seed_summary.tsv"
    candidate_path = output_dir / "candidate_summary.tsv"
    queue_path = output_dir / "md_queue.tsv"
    atomic_write_text(
        decoy_path,
        _tsv(
            (
                "decoy_id",
                "task_id",
                "pair_id",
                "state",
                "candidate_id",
                "template_id",
                "target",
                "evidence_role",
                "seed",
                "reweighted_sc",
                "I_sc",
                "pep_sc",
                "pep_sc_noref",
                "I_bsa_A2",
                "I_hb",
                "I_pack",
                "I_unsat",
                "disulfide_score",
                "qc_status",
                "interface_contact_pairs",
                "interface_receptor_residues",
                "minimum_interface_distance_A",
                "receptor_ca_rmsd_A",
                "start_contact_overlap",
                "start_site_displacement_A",
                "path",
                "sha256",
            ),
            _decoy_rows(decoys=decoys, run_dir=run_dir),
        ),
    )
    atomic_write_text(
        cluster_path,
        _tsv(
            (
                "pose_cluster_id",
                "task_id",
                "state",
                "candidate_id",
                "template_id",
                "seed",
                "cluster_size",
                "cluster_fraction_of_qc_passed",
                "medoid_decoy_id",
                "median_reweighted_sc_REU",
                "median_I_sc_REU",
                "median_pep_sc_REU",
                "primary",
            ),
            _pose_cluster_rows(states),
        ),
    )
    atomic_write_text(
        pair_path,
        _tsv(
            (
                "pair_id",
                "candidate_id",
                "template_id",
                "target",
                "evidence_role",
                "seed",
                "wt_passed_fraction",
                "mutant_passed_fraction",
                "successful",
                "pose_relation",
                "wt_pose_cluster_id",
                "mutant_pose_cluster_id",
                "paired_reweighted_sc_delta_REU",
                "paired_I_sc_delta_REU",
                "paired_pep_sc_delta_REU",
            ),
            _pair_rows(pairs),
        ),
    )
    atomic_write_text(
        candidate_path,
        _tsv(
            (
                "candidate_id",
                "sequence",
                "mutation_string",
                "mutation_count",
                "generation",
                "parent_candidate_ids",
                "parent_edits",
                "net_charge",
                "hydrophobic_count",
                "max_hydrophobic_run",
                "aromatic_count",
                "methionine_count",
                "deamidation_motif_count",
                "aspartimide_motif_count",
                "positive_ranking_median",
                "positive_ranking_worst",
                "positive_interface_median",
                "positive_interface_worst",
                "positive_peptide_median",
                "successful_pair_fraction",
                "candidate_status",
                "failed_gates",
                "resource_status",
            ),
            _candidate_rows(evidence=evidence, resource_status=resource_status),
        ),
    )
    atomic_write_text(
        queue_path,
        _tsv(
            (
                "system_id",
                "candidate_id",
                "sequence",
                "mutation_string",
                "generation",
                "parent_candidate_ids",
                "system_role",
                "template_id",
                "target",
                "evidence_role",
                "structure_path",
                "structure_sha256",
            ),
            queue_rows,
        ),
    )
    return decoy_path, cluster_path, pair_path, candidate_path, queue_path


def analyze_finalists(*, config: AppConfig, run_dir: Path) -> FinalistAnalysisOutcome:
    """形成柔性复核证据、开发性事实和有限阶段五输入队列。"""
    output_dir = run_dir / "finalist_analysis"
    if output_dir.exists():
        raise DesignError(f"finalist analysis already exists: {output_dir}")
    plan = read_finalist_plan(config=config, run_dir=run_dir)
    result_paths = manifest_results(config=config, plan=plan)
    decoys = tuple(
        decoy
        for task in plan.tasks
        for decoy in task_decoys(
            config=config, task=task, result_path=result_paths[task.task_id]
        )
    )
    by_task: dict[str, list[FinalistDecoy]] = defaultdict(list)
    for decoy in decoys:
        by_task[decoy.task.task_id].append(decoy)
    states = tuple(
        state_summary(config=config, task=task, decoys=tuple(by_task[task.task_id]))
        for task in plan.tasks
    )
    pairs = paired_summaries(config=config, states=states)
    evidence = candidate_evidence(config=config, pairs=pairs)
    queue_rows, resource_status = _md_queue(
        config=config,
        plan=plan,
        evidence=evidence,
        states=states,
        pairs=pairs,
    )
    decoy_path, cluster_path, pair_path, candidate_path, queue_path = _write_reports(
        run_dir=run_dir,
        output_dir=output_dir,
        decoys=decoys,
        states=states,
        pairs=pairs,
        evidence=evidence,
        queue_rows=queue_rows,
        resource_status=resource_status,
    )
    manifest_path = output_dir / "analysis_manifest.json"
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.design-finalist-analysis-manifest/1",
            "stage": "design_flexible_verification_analysis",
            "status": "completed",
            "generated_at": utc_now(),
            "method_id": config.design.method_id,
            "flexpepdock_method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "interpretation": (
                "Paired flexible Rosetta ensembles and geometry prioritize candidates; "
                "they do not establish experimental affinity, developability, or efficacy."
            ),
            "sequence_fact_policy": (
                "Transparent counts are reported for experimental review and are not "
                "uncalibrated automatic rejection gates."
            ),
            "finalist_manifest": {
                "path": "../finalist_manifest.json",
                "sha256": sha256_file(run_dir / "finalist_manifest.json"),
            },
            "decoys": {
                "path": decoy_path.name,
                "sha256": sha256_file(decoy_path),
                "count": len(decoys),
                "passed_count": sum(item.geometry.passed for item in decoys),
            },
            "paired_seed_summary": {
                "path": pair_path.name,
                "sha256": sha256_file(pair_path),
                "count": len(pairs),
                "successful_count": sum(item.successful for item in pairs),
            },
            "pose_clusters": {
                "path": cluster_path.name,
                "sha256": sha256_file(cluster_path),
                "count": sum(len(state.pose_clusters) for state in states),
            },
            "candidate_summary": {
                "path": candidate_path.name,
                "sha256": sha256_file(candidate_path),
                "count": len(evidence),
                "eligible_count": sum(item.status == "eligible" for item in evidence),
            },
            "md_queue": {
                "path": queue_path.name,
                "sha256": sha256_file(queue_path),
                "system_count": len(queue_rows),
                "candidate_count": len(
                    {
                        row["candidate_id"]
                        for row in queue_rows
                        if row["system_role"] == "optimized_candidate"
                    }
                ),
            },
        },
    )
    eligible_count = sum(item.status == "eligible" for item in evidence)
    return FinalistAnalysisOutcome(
        manifest_path, len(evidence), eligible_count, len(queue_rows)
    )
