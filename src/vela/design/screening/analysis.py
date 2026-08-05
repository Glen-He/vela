"""阶段四初筛的 WT 配对、多模板聚合和门槛判定。"""

from __future__ import annotations

import csv
import io
import math
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
)
from vela.core.typed_data import object_list, object_mapping
from vela.design.models import DesignError, ScreenTask, SequenceCandidate
from vela.design.scores import SCREEN_SCORE_COLUMNS, ScreenMetrics, screen_metrics
from vela.design.screening.execution import candidate_chemistry
from vela.design.screening.records import ScreenPlan, read_screen_plan
from vela.design.sequence.library import flexibility_reasons, sequence_id
from vela.validation.records import read_document, validate_record
from vela.validation.refinement.reconstruction import validate_flexpepdock_input
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class PairScore:
    """同一模板和 seed 下候选相对 WT 的界面分数差。"""

    pair_id: str
    candidate: SequenceCandidate
    template_id: str
    evidence_role: str
    target: str
    seed: int
    wt_metrics: ScreenMetrics
    mutant_metrics: ScreenMetrics

    @property
    def paired_dG_separated_delta(self) -> float:
        """候选减 WT; 负值只表示当前 Rosetta 协议下更有利。"""
        return self.mutant_metrics.dG_separated - self.wt_metrics.dG_separated


@dataclass(frozen=True, slots=True)
class ScreenAnalysisOutcome:
    """成对任务、模板和候选三级报告。"""

    manifest_path: Path
    candidate_count: int
    eligible_count: int


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _manifest_results(*, config: AppConfig, plan: ScreenPlan) -> dict[str, Path]:
    manifest = read_document(
        plan.run_dir / "screen_manifest.json", name="design screen manifest"
    )
    if (
        manifest.get("schema") != "vela.design-screen-manifest/2"
        or manifest.get("stage") != "design_interface_screen"
        or manifest.get("status") != "completed"
        or manifest.get("design_round") != plan.design_round
        or manifest.get("method_id") != config.design.method_id
        or manifest.get("chemistry_id") != config.chemistry.chemistry_id
        or manifest.get("objective") != config.design.objective
    ):
        raise DesignError("design screen manifest identity is invalid")
    validate_record(
        root=plan.run_dir,
        raw=manifest.get("screen_plan"),
        name="design screen plan",
    )
    try:
        rows = object_list(manifest.get("tasks"), name="design screen tasks")
    except TypeError as exc:
        raise DesignError("design screen manifest tasks are invalid") from exc
    results: dict[str, Path] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="design screen task")
        except TypeError as exc:
            raise DesignError("design screen manifest task is invalid") from exc
        task_id = row.get("task_id")
        if not isinstance(task_id, str) or task_id in results:
            raise DesignError("design screen task identity is invalid or duplicated")
        result, _ = validate_record(
            root=plan.run_dir,
            raw=row.get("task_result"),
            name=f"{task_id} result",
        )
        results[task_id] = result
    expected = {task.task_id for task in plan.tasks}
    if set(results) != expected:
        raise DesignError("design screen manifest task coverage is incomplete")
    expected_wt_contexts = {
        task.wt_context_id for task in plan.tasks if task.wt_context_id is not None
    }
    expected_executions = sum(task.state == "mutant" for task in plan.tasks) + len(
        expected_wt_contexts
    )
    if (
        manifest.get("logical_pair_count") != len(plan.tasks) // 2
        or manifest.get("executed_task_count") != expected_executions
        or manifest.get("shared_wt_context_count") != len(expected_wt_contexts)
    ):
        raise DesignError("design screen execution accounting is invalid")
    return results


def _task_metrics(
    *, config: AppConfig, task: ScreenTask, result_path: Path
) -> ScreenMetrics:
    task_dir = result_path.parent
    result = read_document(result_path, name="design screen task result")
    if (
        result.get("schema") != "vela.design-screen-task-result/2"
        or result.get("status") != "completed"
        or result.get("state") != task.state
        or result.get("template_id") != task.template.template_id
        or result.get("seed") != task.seed
        or result.get("score_columns") != SCREEN_SCORE_COLUMNS
    ):
        raise DesignError(f"design task result identity is invalid: {task.task_id}")
    if task.state == "mutant" and (
        result.get("task_id") != task.task_id
        or result.get("pair_id") != task.pair_id
        or result.get("candidate_id") != task.candidate.candidate_id
        or result.get("wt_context_id") is not None
    ):
        raise DesignError(f"mutant task result identity is invalid: {task.task_id}")
    if task.state == "wt" and result.get("wt_context_id") != task.wt_context_id:
        raise DesignError(f"shared WT result identity is invalid: {task.task_id}")
    score_path, _ = validate_record(
        root=task_dir, raw=result.get("scorefile"), name="design scorefile"
    )
    output_path, _ = validate_record(
        root=task_dir, raw=result.get("output"), name="design output"
    )
    for key in ("fix_disulfide", "log"):
        validate_record(root=task_dir, raw=result.get(key), name=f"design {key}")
    rows = read_rosetta_scorefile(score_path)
    if len(rows) != 1:
        raise DesignError("design task scorefile must contain one row")
    metrics = screen_metrics(rows[0])
    if result.get("interface_metrics") != metrics.as_dict():
        raise DesignError("design task recorded metrics differ from its scorefile")
    receptor_count, _ = validate_flexpepdock_input(
        path=output_path,
        chemistry=candidate_chemistry(config=config, task=task),
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    if receptor_count != task.template.receptor_residue_count:
        raise DesignError("design task output receptor length changed")
    return metrics


def _pair_scores(
    *, config: AppConfig, plan: ScreenPlan, result_paths: dict[str, Path]
) -> tuple[PairScore, ...]:
    grouped: dict[str, list[tuple[ScreenTask, ScreenMetrics]]] = defaultdict(list)
    for task in plan.tasks:
        grouped[task.pair_id].append(
            (
                task,
                _task_metrics(
                    config=config,
                    task=task,
                    result_path=result_paths[task.task_id],
                ),
            )
        )
    results: list[PairScore] = []
    for pair_id, items in sorted(grouped.items()):
        by_state = {task.state: (task, score) for task, score in items}
        if set(by_state) != {"wt", "mutant"}:
            raise DesignError(f"design pair is incomplete: {pair_id}")
        mutant, mutant_metrics = by_state["mutant"]
        wt, wt_metrics = by_state["wt"]
        if (
            mutant.candidate != wt.candidate
            or mutant.template != wt.template
            or mutant.seed != wt.seed
        ):
            raise DesignError(f"design pair inputs differ: {pair_id}")
        results.append(
            PairScore(
                pair_id=pair_id,
                candidate=mutant.candidate,
                template_id=mutant.template.template_id,
                evidence_role=mutant.template.evidence_role,
                target=mutant.template.target,
                seed=mutant.seed,
                wt_metrics=wt_metrics,
                mutant_metrics=mutant_metrics,
            )
        )
    return tuple(results)


def _write_pairs(path: Path, pairs: tuple[PairScore, ...]) -> None:
    fields = (
        "pair_id",
        "candidate_id",
        "sequence",
        "mutation_string",
        "mutation_positions",
        "mutation_count",
        "net_charge",
        "template_id",
        "evidence_role",
        "target",
        "seed",
        "wt_dG_separated_REU",
        "mutant_dG_separated_REU",
        "paired_dG_separated_delta_REU",
        "wt_dSASA_int_A2",
        "mutant_dSASA_int_A2",
        "wt_dG_separated_per_dSASAx100",
        "mutant_dG_separated_per_dSASAx100",
        "wt_delta_unsat_hbonds",
        "mutant_delta_unsat_hbonds",
        "wt_interface_hbonds",
        "mutant_interface_hbonds",
        "wt_shape_complementarity",
        "mutant_shape_complementarity",
        "wt_interface_residue_count",
        "mutant_interface_residue_count",
        "wt_receptor_separated_score_REU",
        "mutant_receptor_separated_score_REU",
        "wt_peptide_separated_score_REU",
        "mutant_peptide_separated_score_REU",
    )
    rows = [
        {
            "pair_id": pair.pair_id,
            "candidate_id": pair.candidate.candidate_id,
            "sequence": pair.candidate.sequence,
            "mutation_string": pair.candidate.mutation_string,
            "mutation_positions": ";".join(map(str, pair.candidate.mutation_positions)),
            "mutation_count": str(pair.candidate.mutation_count),
            "net_charge": str(pair.candidate.net_charge),
            "template_id": pair.template_id,
            "evidence_role": pair.evidence_role,
            "target": pair.target,
            "seed": str(pair.seed),
            "wt_dG_separated_REU": f"{pair.wt_metrics.dG_separated:.6f}",
            "mutant_dG_separated_REU": (f"{pair.mutant_metrics.dG_separated:.6f}"),
            "paired_dG_separated_delta_REU": (f"{pair.paired_dG_separated_delta:.6f}"),
            "wt_dSASA_int_A2": f"{pair.wt_metrics.dSASA_int_A2:.6f}",
            "mutant_dSASA_int_A2": f"{pair.mutant_metrics.dSASA_int_A2:.6f}",
            "wt_dG_separated_per_dSASAx100": (
                f"{pair.wt_metrics.dG_separated_per_dSASAx100:.6f}"
            ),
            "mutant_dG_separated_per_dSASAx100": (
                f"{pair.mutant_metrics.dG_separated_per_dSASAx100:.6f}"
            ),
            "wt_delta_unsat_hbonds": (f"{pair.wt_metrics.delta_unsat_hbonds:.6f}"),
            "mutant_delta_unsat_hbonds": (
                f"{pair.mutant_metrics.delta_unsat_hbonds:.6f}"
            ),
            "wt_interface_hbonds": f"{pair.wt_metrics.interface_hbonds:.6f}",
            "mutant_interface_hbonds": (f"{pair.mutant_metrics.interface_hbonds:.6f}"),
            "wt_shape_complementarity": (
                f"{pair.wt_metrics.shape_complementarity:.6f}"
            ),
            "mutant_shape_complementarity": (
                f"{pair.mutant_metrics.shape_complementarity:.6f}"
            ),
            "wt_interface_residue_count": (
                f"{pair.wt_metrics.interface_residue_count:.6f}"
            ),
            "mutant_interface_residue_count": (
                f"{pair.mutant_metrics.interface_residue_count:.6f}"
            ),
            "wt_receptor_separated_score_REU": (
                f"{pair.wt_metrics.receptor_separated_score:.6f}"
            ),
            "mutant_receptor_separated_score_REU": (
                f"{pair.mutant_metrics.receptor_separated_score:.6f}"
            ),
            "wt_peptide_separated_score_REU": (
                f"{pair.wt_metrics.peptide_separated_score:.6f}"
            ),
            "mutant_peptide_separated_score_REU": (
                f"{pair.mutant_metrics.peptide_separated_score:.6f}"
            ),
        }
        for pair in pairs
    ]
    atomic_write_text(path, _tsv(fields, rows))


def _template_rows(pairs: tuple[PairScore, ...]) -> list[dict[str, str]]:
    grouped: dict[tuple[str, str], list[PairScore]] = defaultdict(list)
    for pair in pairs:
        grouped[(pair.candidate.candidate_id, pair.template_id)].append(pair)
    rows: list[dict[str, str]] = []
    for (_, _), values in sorted(grouped.items()):
        first = values[0]
        deltas = [item.paired_dG_separated_delta for item in values]
        rows.append(
            {
                "candidate_id": first.candidate.candidate_id,
                "template_id": first.template_id,
                "evidence_role": first.evidence_role,
                "target": first.target,
                "seed_count": str(len(values)),
                "median_paired_dG_separated_delta_REU": (
                    f"{statistics.median(deltas):.6f}"
                ),
                "best_paired_dG_separated_delta_REU": f"{min(deltas):.6f}",
                "worst_paired_dG_separated_delta_REU": f"{max(deltas):.6f}",
            }
        )
    return rows


def _candidate_status(
    *,
    config: AppConfig,
    median_delta: float,
    favorable_seed_fraction: float,
) -> tuple[str, tuple[str, ...]]:
    settings = config.design.analysis
    if not settings.complete:
        return "calibration_required", ("analysis_not_calibrated",)
    failed: list[str] = []
    median_gate = settings.max_median_paired_dG_separated_delta_REU
    favorable_gate = settings.min_favorable_seed_fraction
    if median_gate is not None and median_delta > median_gate:
        failed.append("median_paired_dG_separated_delta")
    if favorable_gate is not None and favorable_seed_fraction < favorable_gate:
        failed.append("favorable_seed_fraction")
    return (
        ("screen_deprioritized", tuple(failed)) if failed else ("screen_supported", ())
    )


def _candidate_rows(
    *, config: AppConfig, pairs: tuple[PairScore, ...]
) -> list[dict[str, str]]:
    grouped: dict[str, list[PairScore]] = defaultdict(list)
    for pair in pairs:
        grouped[pair.candidate.candidate_id].append(pair)
    rows: list[dict[str, str]] = []
    for _, values in sorted(grouped.items()):
        first = values[0]
        if any(item.evidence_role != "positive" for item in values):
            raise DesignError("single-target design contains a non-target template")
        positive = [item.paired_dG_separated_delta for item in values]
        positive_median = float(statistics.median(positive))
        positive_worst = max(positive)
        favorable_threshold = (
            config.design.analysis.max_median_paired_dG_separated_delta_REU
        )
        if favorable_threshold is None:
            raise DesignError("design screen favorable threshold is unresolved")
        favorable_seed_fraction = sum(
            value <= favorable_threshold for value in positive
        ) / len(positive)
        status, failed = _candidate_status(
            config=config,
            median_delta=positive_median,
            favorable_seed_fraction=favorable_seed_fraction,
        )
        disulfide_positions = frozenset(
            position
            for bond in config.chemistry.disulfide_bonds
            for position in (bond.first, bond.second)
        )
        flexibility = flexibility_reasons(
            reference=config.chemistry.sequence,
            candidate=first.candidate,
            disulfide_positions=disulfide_positions,
        )
        if flexibility:
            status = "flexibility_required"
        rows.append(
            {
                "candidate_id": first.candidate.candidate_id,
                "sequence": first.candidate.sequence,
                "mutation_string": first.candidate.mutation_string,
                "mutation_positions": ";".join(
                    map(str, first.candidate.mutation_positions)
                ),
                "mutation_count": str(first.candidate.mutation_count),
                "net_charge": str(first.candidate.net_charge),
                "design_round": first.candidate.design_round,
                "generation": str(first.candidate.generation),
                "parent_candidate_ids": ";".join(
                    item.candidate_id for item in first.candidate.parents
                ),
                "parent_edits": ";".join(
                    f"{item.candidate_id}:{item.edit}:{item.edit_type}"
                    for item in first.candidate.parents
                ),
                "positive_pair_count": str(len(positive)),
                "positive_median_paired_dG_separated_delta_REU": (
                    f"{positive_median:.6f}"
                ),
                "positive_worst_paired_dG_separated_delta_REU": (
                    f"{positive_worst:.6f}"
                ),
                "favorable_seed_fraction": f"{favorable_seed_fraction:.6f}",
                "candidate_status": status,
                "failed_gates": ";".join(failed),
                "flexibility_reasons": ";".join(flexibility),
            }
        )
    return rows


def _combination_epistasis_rows(
    *, config: AppConfig, plan: ScreenPlan, pairs: tuple[PairScore, ...]
) -> list[dict[str, str]]:
    """在同模板、同 seed 的固定骨架上下文中计算双突变非加和性。"""
    if plan.design_round != "combination":
        return []
    plan_document = read_document(
        plan.run_dir / "screen_plan.json", name="combination screen plan"
    )
    inputs = object_mapping(plan_document.get("inputs"), name="combination inputs")
    parent = object_mapping(inputs.get("parent_single_screen"), name="parent screen")
    relative = parent.get("path")
    if not isinstance(relative, str):
        raise DesignError("parent single-screen path is invalid")
    parent_dir = (config.paths.outputs_dir / relative).resolve()
    screen_root = (config.paths.outputs_dir / "design" / "screens").resolve()
    if not parent_dir.is_relative_to(screen_root):
        raise DesignError("parent single-screen path escapes screen outputs")
    parent_manifest = read_document(
        parent_dir / "screen_analysis" / "analysis_manifest.json",
        name="parent single-screen analysis",
    )
    if (
        parent_manifest.get("schema") != "vela.design-screen-analysis-manifest/2"
        or parent_manifest.get("design_round") != "single"
        or parent_manifest.get("status") != "completed"
    ):
        raise DesignError("parent single-screen analysis identity is invalid")
    paired_path, _ = validate_record(
        root=parent_dir / "screen_analysis",
        raw=parent_manifest.get("paired_scores"),
        name="parent paired scores",
    )
    with paired_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {
            "candidate_id",
            "template_id",
            "seed",
            "paired_dG_separated_delta_REU",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise DesignError("parent paired-score columns are invalid")
        parent_rows = tuple(dict(row) for row in reader)
    singles: dict[tuple[str, str, int], float] = {}
    for row in parent_rows:
        try:
            key = (row["candidate_id"], row["template_id"], int(row["seed"]))
            value = float(row["paired_dG_separated_delta_REU"])
        except (ValueError, OverflowError) as exc:
            raise DesignError("parent paired score is invalid") from exc
        if key in singles or not math.isfinite(value):
            raise DesignError("parent paired-score identity is duplicated or invalid")
        singles[key] = value
    rows: list[dict[str, str]] = []
    for pair in pairs:
        if pair.candidate.mutation_count != 2:
            raise DesignError("combination epistasis requires a double mutant")
        single_ids: list[str] = []
        single_deltas: list[float] = []
        for position in pair.candidate.mutation_positions:
            sequence = list(config.chemistry.sequence)
            sequence[position - 1] = pair.candidate.sequence[position - 1]
            candidate_id = sequence_id("".join(sequence))
            single_ids.append(candidate_id)
            value = singles.get((candidate_id, pair.template_id, pair.seed))
            if value is not None:
                single_deltas.append(value)
        computable = len(single_deltas) == 2
        epsilon = (
            pair.paired_dG_separated_delta - sum(single_deltas) if computable else None
        )
        rows.append(
            {
                "candidate_id": pair.candidate.candidate_id,
                "template_id": pair.template_id,
                "seed": str(pair.seed),
                "fixed_pose_family": pair.template_id,
                "single_candidate_ids": ";".join(single_ids),
                "paired_double_delta_REU": (f"{pair.paired_dG_separated_delta:.6f}"),
                "single_delta_sum_REU": (
                    f"{sum(single_deltas):.6f}" if computable else ""
                ),
                "epistasis_REU": f"{epsilon:.6f}" if epsilon is not None else "",
                "status": "computed" if computable else "not_computable",
            }
        )
    return rows


def analyze_screen(*, config: AppConfig, run_dir: Path) -> ScreenAnalysisOutcome:
    """验证完整运行并写出配对、模板和候选三级不可变报告。"""
    output_dir = run_dir / "screen_analysis"
    if output_dir.exists():
        raise DesignError(f"design screen analysis already exists: {output_dir}")
    plan = read_screen_plan(config=config, run_dir=run_dir)
    results = _manifest_results(config=config, plan=plan)
    pairs = _pair_scores(config=config, plan=plan, result_paths=results)
    pair_path = output_dir / "paired_scores.tsv"
    template_path = output_dir / "template_summary.tsv"
    candidate_path = output_dir / "candidate_summary.tsv"
    _write_pairs(pair_path, pairs)
    template_rows = _template_rows(pairs)
    template_fields = (
        "candidate_id",
        "template_id",
        "evidence_role",
        "target",
        "seed_count",
        "median_paired_dG_separated_delta_REU",
        "best_paired_dG_separated_delta_REU",
        "worst_paired_dG_separated_delta_REU",
    )
    atomic_write_text(template_path, _tsv(template_fields, template_rows))
    candidate_rows = _candidate_rows(config=config, pairs=pairs)
    candidate_fields = (
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
        "positive_pair_count",
        "positive_median_paired_dG_separated_delta_REU",
        "positive_worst_paired_dG_separated_delta_REU",
        "favorable_seed_fraction",
        "candidate_status",
        "failed_gates",
        "flexibility_reasons",
    )
    atomic_write_text(candidate_path, _tsv(candidate_fields, candidate_rows))
    epistasis_rows = _combination_epistasis_rows(config=config, plan=plan, pairs=pairs)
    epistasis_path = output_dir / "epistasis.tsv"
    if plan.design_round == "combination":
        atomic_write_text(
            epistasis_path,
            _tsv(
                (
                    "candidate_id",
                    "template_id",
                    "seed",
                    "fixed_pose_family",
                    "single_candidate_ids",
                    "paired_double_delta_REU",
                    "single_delta_sum_REU",
                    "epistasis_REU",
                    "status",
                ),
                epistasis_rows,
            ),
        )
    manifest_path = output_dir / "analysis_manifest.json"
    parameters: dict[str, JsonValue] = {
        "primary_metric": "dG_separated",
        "score_units": "Rosetta score units (REU)",
        "delta_definition": "mutant_score - wt_score",
        "score_direction": "negative_delta_favors_mutant_in_this_protocol",
        "score_columns": SCREEN_SCORE_COLUMNS,
        "thresholds_calibrated": config.design.analysis.calibrated,
        "max_median_paired_dG_separated_delta_REU": (
            config.design.analysis.max_median_paired_dG_separated_delta_REU
        ),
        "min_favorable_seed_fraction": (
            config.design.analysis.min_favorable_seed_fraction
        ),
    }
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.design-screen-analysis-manifest/2",
            "stage": "design_interface_screen_analysis",
            "status": "completed",
            "generated_at": utc_now(),
            "design_round": plan.design_round,
            "method_id": config.design.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "interpretation": (
                "Paired Rosetta score deltas prioritize sequences under the frozen "
                "screen protocol; they are not experimental affinity or efficacy."
            ),
            "parameters": parameters,
            "screen_manifest": {
                "path": "../screen_manifest.json",
                "sha256": sha256_file(run_dir / "screen_manifest.json"),
            },
            "paired_scores": {
                "path": pair_path.name,
                "sha256": sha256_file(pair_path),
                "count": len(pairs),
            },
            "template_summary": {
                "path": template_path.name,
                "sha256": sha256_file(template_path),
                "count": len(template_rows),
            },
            "candidate_summary": {
                "path": candidate_path.name,
                "sha256": sha256_file(candidate_path),
                "count": len(candidate_rows),
                "eligible_count": sum(
                    row["candidate_status"]
                    in {"screen_supported", "flexibility_required"}
                    for row in candidate_rows
                ),
            },
            "epistasis": (
                {
                    "path": epistasis_path.name,
                    "sha256": sha256_file(epistasis_path),
                    "count": len(epistasis_rows),
                    "interpretation": (
                        "Protocol-level non-additivity in matched fixed-backbone "
                        "template/seed contexts; not experimental thermodynamic synergy."
                    ),
                }
                if plan.design_round == "combination"
                else None
            ),
        },
    )
    return ScreenAnalysisOutcome(
        manifest_path,
        len(candidate_rows),
        sum(row["candidate_status"] == "eligible" for row in candidate_rows),
    )
