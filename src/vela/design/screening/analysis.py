"""阶段四初筛的 WT 配对、多模板聚合和门槛判定。"""

from __future__ import annotations

import csv
import io
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
from vela.design.screening.execution import candidate_chemistry
from vela.design.screening.records import ScreenPlan, read_screen_plan
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
    wt_score: float
    mutant_score: float

    @property
    def delta(self) -> float:
        """候选减 WT; 负值只表示当前 Rosetta 协议下更有利。"""
        return self.mutant_score - self.wt_score


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
        manifest.get("schema") != "vela.design-screen-manifest/1"
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
    return results


def _task_score(*, config: AppConfig, task: ScreenTask, result_path: Path) -> float:
    task_dir = result_path.parent
    result = read_document(result_path, name="design screen task result")
    if (
        result.get("schema") != "vela.design-screen-task-result/1"
        or result.get("status") != "completed"
        or result.get("task_id") != task.task_id
        or result.get("pair_id") != task.pair_id
        or result.get("state") != task.state
        or result.get("candidate_id") != task.candidate.candidate_id
        or result.get("template_id") != task.template.template_id
        or result.get("seed") != task.seed
        or result.get("ranking_score_name") != config.design.screen.ranking_score
    ):
        raise DesignError(f"design task result identity is invalid: {task.task_id}")
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
    score = rows[0].score(config.design.screen.ranking_score)
    if result.get("ranking_score") != score:
        raise DesignError("design task recorded score differs from its scorefile")
    receptor_count, _ = validate_flexpepdock_input(
        path=output_path,
        chemistry=candidate_chemistry(config=config, task=task),
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    if receptor_count != task.template.receptor_residue_count:
        raise DesignError("design task output receptor length changed")
    return score


def _pair_scores(
    *, config: AppConfig, plan: ScreenPlan, result_paths: dict[str, Path]
) -> tuple[PairScore, ...]:
    grouped: dict[str, list[tuple[ScreenTask, float]]] = defaultdict(list)
    for task in plan.tasks:
        grouped[task.pair_id].append(
            (
                task,
                _task_score(
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
        mutant, mutant_score = by_state["mutant"]
        wt, wt_score = by_state["wt"]
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
                wt_score=wt_score,
                mutant_score=mutant_score,
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
        "wt_score",
        "mutant_score",
        "delta_score",
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
            "wt_score": f"{pair.wt_score:.6f}",
            "mutant_score": f"{pair.mutant_score:.6f}",
            "delta_score": f"{pair.delta:.6f}",
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
        deltas = [item.delta for item in values]
        rows.append(
            {
                "candidate_id": first.candidate.candidate_id,
                "template_id": first.template_id,
                "evidence_role": first.evidence_role,
                "target": first.target,
                "seed_count": str(len(values)),
                "median_delta_score": f"{statistics.median(deltas):.6f}",
                "best_delta_score": f"{min(deltas):.6f}",
                "worst_delta_score": f"{max(deltas):.6f}",
            }
        )
    return rows


def _candidate_status(
    *,
    config: AppConfig,
    positive_median: float,
    positive_worst: float,
) -> tuple[str, tuple[str, ...]]:
    settings = config.design.analysis
    if not settings.complete:
        return "calibration_required", ("analysis_not_calibrated",)
    failed: list[str] = []
    if (
        settings.max_positive_median_delta is not None
        and positive_median > settings.max_positive_median_delta
    ):
        failed.append("positive_median")
    if (
        settings.max_positive_worst_delta is not None
        and positive_worst > settings.max_positive_worst_delta
    ):
        failed.append("positive_worst")
    return ("rejected", tuple(failed)) if failed else ("eligible", ())


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
        positive = [item.delta for item in values]
        positive_median = float(statistics.median(positive))
        positive_worst = max(positive)
        status, failed = _candidate_status(
            config=config,
            positive_median=positive_median,
            positive_worst=positive_worst,
        )
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
                "positive_median_delta_score": f"{positive_median:.6f}",
                "positive_worst_delta_score": f"{positive_worst:.6f}",
                "candidate_status": status,
                "failed_gates": ";".join(failed),
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
        "median_delta_score",
        "best_delta_score",
        "worst_delta_score",
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
        "positive_median_delta_score",
        "positive_worst_delta_score",
        "candidate_status",
        "failed_gates",
    )
    atomic_write_text(candidate_path, _tsv(candidate_fields, candidate_rows))
    manifest_path = output_dir / "analysis_manifest.json"
    parameters: dict[str, JsonValue] = {
        "ranking_score": config.design.screen.ranking_score,
        "score_direction": "negative_delta_favors_mutant_in_this_protocol",
        "thresholds_calibrated": config.design.analysis.calibrated,
        "max_positive_median_delta": (config.design.analysis.max_positive_median_delta),
        "max_positive_worst_delta": (config.design.analysis.max_positive_worst_delta),
    }
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.design-screen-analysis-manifest/1",
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
                    row["candidate_status"] == "eligible" for row in candidate_rows
                ),
            },
        },
    )
    return ScreenAnalysisOutcome(
        manifest_path,
        len(candidate_rows),
        sum(row["candidate_status"] == "eligible" for row in candidate_rows),
    )
