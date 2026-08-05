"""阶段四初筛计划的数据合同、读取和完整性复核。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import JsonValue, is_current_vela_software, sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.design.models import (
    DESIGN_ROUNDS,
    CandidateParent,
    DesignError,
    DesignTemplate,
    ScreenTask,
    SequenceCandidate,
)
from vela.design.scores import SCREEN_SCORE_COLUMNS
from vela.design.sequence.library import build_candidate
from vela.design.sequence.neighborhood import candidate_parent_edit
from vela.validation.records import (
    nonnegative_integer,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.rosetta import verify_rosetta_scripts_tool

SCREEN_PLAN_NAME = "screen_plan.json"


@dataclass(frozen=True, slots=True)
class ScreenPlan:
    """一个已冻结的单点或组合成对界面筛查计划。"""

    run_dir: Path
    design_round: str
    candidates: tuple[SequenceCandidate, ...]
    templates: tuple[DesignTemplate, ...]
    tasks: tuple[ScreenTask, ...]
    protocol_path: Path


def design_parameters(config: AppConfig) -> dict[str, JsonValue]:
    """返回必须随计划冻结的阶段四参数。"""
    settings = config.design
    return {
        "objective": settings.objective,
        "seeds": list(settings.seeds),
        "sequence": {
            "mutable_positions": list(settings.sequence.mutable_positions),
            "allowed_amino_acids": settings.sequence.allowed_amino_acids,
            "candidate_histidine_state": settings.sequence.candidate_histidine_state,
        },
        "screen": {
            "parallel_tasks": settings.screen.parallel_tasks,
            "neighbor_distance_A": settings.screen.neighbor_distance_A,
            "score_function": settings.screen.score_function,
            "ranking_score": settings.screen.ranking_score,
            "pack_separated": settings.screen.pack_separated,
            "ligand_chain": "P",
            "pack_input": False,
            "packstat": False,
            "interface_sc": True,
            "score_columns": SCREEN_SCORE_COLUMNS,
        },
        "combination": {
            "max_options_per_position": settings.combination.max_options_per_position,
            "max_candidates": settings.combination.max_candidates,
        },
        "iteration": {
            "max_parents": settings.iteration.max_parents,
            "max_total_mutations": settings.iteration.max_total_mutations,
            "max_candidates": settings.iteration.max_candidates,
        },
        "analysis": {
            "calibrated": settings.analysis.calibrated,
            "max_median_paired_dG_separated_delta_REU": (
                settings.analysis.max_median_paired_dG_separated_delta_REU
            ),
            "min_favorable_seed_fraction": (
                settings.analysis.min_favorable_seed_fraction
            ),
        },
        "finalists": {
            "parallel_tasks": settings.finalists.parallel_tasks,
            "max_candidates": settings.finalists.max_candidates,
            "max_flexibility_required_candidates": (
                settings.finalists.max_flexibility_required_candidates
            ),
            "max_md_candidates": settings.finalists.max_md_candidates,
            "seeds": list(settings.finalists.seeds),
            "ranking_score": settings.finalists.ranking_score,
            "interface_score": settings.finalists.interface_score,
            "calibrated": settings.finalists.calibrated,
            "min_passed_decoy_fraction": (settings.finalists.min_passed_decoy_fraction),
            "min_successful_seeds": settings.finalists.min_successful_seeds,
            "max_positive_median_ranking_delta": (
                settings.finalists.max_positive_median_ranking_delta
            ),
            "max_positive_median_interface_delta": (
                settings.finalists.max_positive_median_interface_delta
            ),
        },
    }


def _integers(value: object, *, name: str) -> tuple[int, ...]:
    try:
        items = object_list(value, name=name)
    except TypeError as exc:
        raise DesignError(f"{name} is invalid") from exc
    result: list[int] = []
    for item in items:
        if not isinstance(item, int) or isinstance(item, bool):
            raise DesignError(f"{name} is invalid")
        result.append(item)
    return tuple(result)


def candidate_from_record(*, raw: object, config: AppConfig) -> SequenceCandidate:
    """读取候选记录并从完整序列重建其派生身份。"""
    try:
        row = object_mapping(raw, name="design candidate")
    except TypeError as exc:
        raise DesignError("design candidate record is invalid") from exc
    sequence = row.get("sequence")
    mutation_string = row.get("mutation_string")
    proposal_source = row.get("proposal_source")
    design_round = row.get("design_round")
    generation = row.get("generation")
    mutation_count = row.get("mutation_count")
    net_charge = row.get("net_charge")
    if (
        not isinstance(sequence, str)
        or not isinstance(mutation_string, str)
        or not isinstance(proposal_source, str)
        or not isinstance(design_round, str)
        or not isinstance(generation, int)
        or isinstance(generation, bool)
        or not isinstance(mutation_count, int)
        or isinstance(mutation_count, bool)
        or not isinstance(net_charge, int)
        or isinstance(net_charge, bool)
    ):
        raise DesignError("design candidate fields are invalid")
    try:
        parent_rows = object_list(row.get("parents"), name="candidate parents")
    except TypeError as exc:
        raise DesignError("design candidate parents are invalid") from exc
    parents: list[CandidateParent] = []
    for raw_parent in parent_rows:
        try:
            parent = object_mapping(raw_parent, name="candidate parent")
        except TypeError as exc:
            raise DesignError("design candidate parent is invalid") from exc
        parent_id = parent.get("candidate_id")
        edit = parent.get("edit")
        edit_type = parent.get("edit_type")
        if (
            not isinstance(parent_id, str)
            or not isinstance(edit, str)
            or not isinstance(edit_type, str)
        ):
            raise DesignError("design candidate parent fields are invalid")
        parents.append(CandidateParent(parent_id, edit, edit_type))
    try:
        ancestor_rows = object_list(
            row.get("ancestor_candidate_ids"), name="candidate ancestor IDs"
        )
    except TypeError as exc:
        raise DesignError("design candidate ancestor IDs are invalid") from exc
    ancestor_candidate_ids = tuple(
        safe_identifier(item, name="ancestor candidate ID") for item in ancestor_rows
    )
    candidate = SequenceCandidate(
        candidate_id=safe_identifier(row.get("candidate_id"), name="candidate ID"),
        sequence=sequence,
        mutation_string=mutation_string,
        mutation_positions=_integers(
            row.get("mutation_positions"), name="candidate mutation positions"
        ),
        mutation_count=mutation_count,
        net_charge=net_charge,
        design_round=design_round,
        generation=generation,
        parents=tuple(parents),
        ancestor_candidate_ids=ancestor_candidate_ids,
        proposal_source=proposal_source,
    )
    rebuilt = build_candidate(
        chemistry=config.chemistry,
        settings=config.design,
        sequence=sequence,
        design_round=design_round,
        generation=generation,
        parents=tuple(parents),
        ancestor_candidate_ids=ancestor_candidate_ids,
        proposal_source=proposal_source,
    )
    if candidate != rebuilt:
        raise DesignError("design candidate record differs from its derived identity")
    return candidate


def _validate_source_inputs(
    *, config: AppConfig, run_dir: Path, inputs: dict[str, object], design_round: str
) -> None:
    validate_record(
        root=config.paths.outputs_dir,
        raw=inputs.get("qualification_report"),
        name="design qualification report",
    )
    validate_record(
        root=run_dir,
        raw=inputs.get("candidate_library"),
        name="design candidate library",
    )
    if design_round == "single":
        source_key = "stage3_refinement_run"
    elif design_round == "combination":
        source_key = "parent_single_screen"
    else:
        source_key = "parent_finalist"
    try:
        source = object_mapping(inputs.get(source_key), name=source_key)
    except TypeError as exc:
        raise DesignError(f"{source_key} record is invalid") from exc
    relative = source.get("path")
    if not isinstance(relative, str):
        raise DesignError(f"{source_key} path is invalid")
    source_dir = (config.paths.outputs_dir / relative).resolve()
    if design_round == "single":
        expected_root = (
            config.paths.outputs_dir / "validation" / "refinements"
        ).resolve()
    elif design_round == "combination":
        expected_root = (config.paths.outputs_dir / "design" / "screens").resolve()
    else:
        expected_root = (config.paths.outputs_dir / "design" / "finalists").resolve()
    try:
        source_dir.relative_to(expected_root)
    except ValueError as exc:
        raise DesignError(f"{source_key} escapes its expected output area") from exc
    if design_round == "single":
        checks = (
            ("refinement_manifest_sha256", source_dir / "refinement_manifest.json"),
            (
                "analysis_manifest_sha256",
                source_dir / "refinement_analysis" / "analysis_manifest.json",
            ),
            (
                "candidate_review_manifest_sha256",
                source_dir / "candidate_review" / "review_manifest.json",
            ),
        )
    elif design_round == "combination":
        checks = (
            ("screen_plan_sha256", source_dir / SCREEN_PLAN_NAME),
            (
                "analysis_manifest_sha256",
                source_dir / "screen_analysis" / "analysis_manifest.json",
            ),
        )
    else:
        checks = (
            ("finalist_plan_sha256", source_dir / "finalist_plan.json"),
            ("finalist_manifest_sha256", source_dir / "finalist_manifest.json"),
            (
                "analysis_manifest_sha256",
                source_dir / "finalist_analysis" / "analysis_manifest.json",
            ),
        )
    for key, path in checks:
        expected_hash = source.get(key)
        if (
            not isinstance(expected_hash, str)
            or not path.is_file()
            or sha256_file(path) != expected_hash
        ):
            raise DesignError(f"{source_key} hash mismatch: {path}")


def _validate_iteration_lineage(
    *,
    config: AppConfig,
    inputs: dict[str, object],
    candidates: tuple[SequenceCandidate, ...],
) -> None:
    try:
        source = object_mapping(inputs.get("parent_finalist"), name="parent finalist")
        selected_rows = object_list(
            source.get("selected_parent_candidate_ids"),
            name="selected parent candidate IDs",
        )
    except TypeError as exc:
        raise DesignError("iteration parent lineage record is invalid") from exc
    selected_ids = tuple(
        safe_identifier(item, name="selected parent candidate ID")
        for item in selected_rows
    )
    if (
        not selected_ids
        or len(selected_ids) != len(set(selected_ids))
        or selected_ids != tuple(sorted(selected_ids))
    ):
        raise DesignError("selected parent candidate IDs are invalid")
    relative = source.get("path")
    if not isinstance(relative, str):
        raise DesignError("iteration parent finalist path is invalid")
    parent_dir = (config.paths.outputs_dir / relative).resolve()
    parent_document = read_document(
        parent_dir / "finalist_plan.json", name="iteration parent finalist plan"
    )
    try:
        parent_rows = object_list(
            parent_document.get("candidates"), name="iteration parent candidates"
        )
    except TypeError as exc:
        raise DesignError("iteration parent candidate records are invalid") from exc
    parents = tuple(
        candidate_from_record(raw=item, config=config) for item in parent_rows
    )
    parent_by_id = {item.candidate_id: item for item in parents}
    if set(selected_ids) != {
        parent.candidate_id for candidate in candidates for parent in candidate.parents
    }:
        raise DesignError("iteration plan does not cover its selected parent set")
    for candidate in candidates:
        direct_parents: list[SequenceCandidate] = []
        for link in candidate.parents:
            parent = parent_by_id.get(link.candidate_id)
            if parent is None or parent.candidate_id not in selected_ids:
                raise DesignError("iteration candidate references an unknown parent")
            if link != candidate_parent_edit(
                reference=config.chemistry.sequence,
                parent=parent,
                child_sequence=candidate.sequence,
            ):
                raise DesignError("iteration candidate parent edit is inconsistent")
            direct_parents.append(parent)
        if any(
            parent.generation + 1 != candidate.generation for parent in direct_parents
        ):
            raise DesignError("iteration candidate generation is inconsistent")
        expected_ancestors = {
            ancestor
            for parent in direct_parents
            for ancestor in (*parent.ancestor_candidate_ids, parent.candidate_id)
        }
        if set(candidate.ancestor_candidate_ids) != expected_ancestors:
            raise DesignError("iteration candidate ancestor lineage is inconsistent")


def read_screen_plan(*, config: AppConfig, run_dir: Path) -> ScreenPlan:
    """复核当前配置、工具、文件哈希和完整成对任务覆盖。"""
    plan = read_document(run_dir / SCREEN_PLAN_NAME, name="design screen plan")
    if (
        plan.get("schema") != "vela.design-screen-plan/2"
        or plan.get("stage") != "design_interface_screen"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("method_id") != config.design.method_id
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
        or plan.get("objective") != config.design.objective
        or plan.get("parameters") != design_parameters(config)
    ):
        raise DesignError("design screen plan identity or parameters changed")
    raw_round = plan.get("design_round")
    if not isinstance(raw_round, str) or raw_round not in DESIGN_ROUNDS:
        raise DesignError("design screen round is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="design screen inputs")
        candidate_rows = object_list(plan.get("candidates"), name="design candidates")
        template_rows = object_list(plan.get("templates"), name="design templates")
        task_rows = object_list(plan.get("tasks"), name="design tasks")
    except TypeError as exc:
        raise DesignError("design screen plan structure is invalid") from exc
    snapshot, _ = validate_record(
        root=run_dir, raw=inputs.get("config_snapshot"), name="config snapshot"
    )
    if sha256_file(snapshot) != config.source_snapshot_sha256:
        raise DesignError("current project config differs from the design plan")
    protocol_path, _ = validate_record(
        root=run_dir, raw=inputs.get("protocol"), name="design protocol"
    )
    _validate_source_inputs(
        config=config, run_dir=run_dir, inputs=inputs, design_round=raw_round
    )
    candidates = tuple(
        candidate_from_record(raw=raw, config=config) for raw in candidate_rows
    )
    candidate_by_id = {item.candidate_id: item for item in candidates}
    if not candidates or len(candidate_by_id) != len(candidates):
        raise DesignError("design candidate identities are empty or duplicated")
    if raw_round == "iteration":
        _validate_iteration_lineage(
            config=config,
            inputs=inputs,
            candidates=candidates,
        )
    templates: list[DesignTemplate] = []
    for raw in template_rows:
        try:
            row = object_mapping(raw, name="design template")
        except TypeError as exc:
            raise DesignError("design template record is invalid") from exc
        path, digest = validate_record(
            root=run_dir, raw=row.get("structure"), name="design template structure"
        )
        role = row.get("evidence_role")
        target = row.get("target")
        if not isinstance(role, str) or not isinstance(target, str):
            raise DesignError("design template role or target is invalid")
        templates.append(
            DesignTemplate(
                safe_identifier(row.get("template_id"), name="template ID"),
                role,
                safe_identifier(row.get("cluster_id"), name="cluster ID"),
                safe_identifier(row.get("candidate_id"), name="Stage 3 candidate ID"),
                safe_identifier(row.get("receptor_id"), name="receptor ID"),
                safe_identifier(target, name="target ID"),
                path,
                digest,
                nonnegative_integer(
                    row.get("receptor_residue_count"), name="receptor residue count"
                ),
                _integers(
                    row.get("fixed_histidine_pose_indices"),
                    name="fixed histidine indices",
                ),
            )
        )
    template_by_id = {item.template_id: item for item in templates}
    if not templates or len(template_by_id) != len(templates):
        raise DesignError("design template identities are empty or duplicated")
    tasks: list[ScreenTask] = []
    for raw in task_rows:
        try:
            row = object_mapping(raw, name="design task")
        except TypeError as exc:
            raise DesignError("design task record is invalid") from exc
        if row.get("status") != "planned":
            raise DesignError("design task status is invalid")
        candidate_id = safe_identifier(row.get("candidate_id"), name="candidate ID")
        template_id = safe_identifier(row.get("template_id"), name="template ID")
        state = row.get("state")
        if not isinstance(state, str):
            raise DesignError("design task state is invalid")
        raw_wt_context_id = row.get("wt_context_id")
        if state == "mutant" and raw_wt_context_id is not None:
            raise DesignError("mutant design task must not declare a WT context")
        resfile_path, resfile_hash = validate_record(
            root=run_dir, raw=row.get("resfile"), name="design task resfile"
        )
        try:
            candidate = candidate_by_id[candidate_id]
            template = template_by_id[template_id]
        except KeyError as exc:
            raise DesignError("design task references an unknown input") from exc
        tasks.append(
            ScreenTask(
                safe_identifier(row.get("task_id"), name="task ID"),
                safe_identifier(row.get("pair_id"), name="pair ID"),
                state,
                candidate,
                template,
                nonnegative_integer(row.get("seed"), name="design seed"),
                resfile_path,
                resfile_hash,
                (
                    safe_identifier(row.get("wt_context_id"), name="WT context ID")
                    if state == "wt"
                    else None
                ),
            )
        )
    expected_count = len(candidates) * len(templates) * len(config.design.seeds) * 2
    if len(tasks) != expected_count or len({item.task_id for item in tasks}) != len(
        tasks
    ):
        raise DesignError("design task coverage is incomplete or duplicated")
    pairs: dict[str, list[ScreenTask]] = {}
    for task in tasks:
        pairs.setdefault(task.pair_id, []).append(task)
    if any(
        len(items) != 2
        or {item.state for item in items} != {"wt", "mutant"}
        or len({item.candidate.candidate_id for item in items}) != 1
        or len({item.template.template_id for item in items}) != 1
        or len({item.seed for item in items}) != 1
        for items in pairs.values()
    ):
        raise DesignError("design WT/mutant pair coverage is invalid")
    wt_contexts: dict[str, list[ScreenTask]] = {}
    for task in tasks:
        if task.wt_context_id is not None:
            wt_contexts.setdefault(task.wt_context_id, []).append(task)
    if any(
        len({item.template.template_id for item in items}) != 1
        or len({item.seed for item in items}) != 1
        or len({item.resfile_sha256 for item in items}) != 1
        for items in wt_contexts.values()
    ):
        raise DesignError("shared WT context inputs differ")
    verify_rosetta_scripts_tool(config.validation.rosetta)
    return ScreenPlan(
        run_dir,
        raw_round,
        candidates,
        tuple(templates),
        tuple(tasks),
        protocol_path,
    )
