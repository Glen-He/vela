"""阶段三正式 candidate 局部精修来源与不可变计划。"""

from __future__ import annotations

import csv
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    is_vela_software_identity,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import json_value, object_list, object_mapping
from vela.validation.models import ValidationError
from vela.validation.readiness import assess_validation_readiness
from vela.validation.records import (
    file_record,
    nonnegative_integer,
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.guided import GUIDED_EVIDENCE
from vela.validation.refinement.handoff_plan import (
    EXPLORATORY_DISCOVERY_EVIDENCE,
    EXPLORATORY_HANDOFF_EVIDENCE,
    FUNNEL_SCREENING_HANDOFF_EVIDENCE,
    HANDOFF_PLAN_SCHEMA,
    MAIN_DISCOVERY_EVIDENCE,
    MAIN_HANDOFF_EVIDENCE,
    SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE,
    exploration_promotion_contract,
    funnel_screening_contract,
    source_seed_confirmation_contract,
)
from vela.validation.refinement.reconstruction import validate_flexpepdock_input
from vela.validation.rosetta import (
    verify_flexpepdock_tool,
    verify_rosetta_scripts_tool,
)

REFINEMENT_PLAN_NAME = "refinement_plan.json"
BLIND_REFINEMENT_EVIDENCE = "main_discovery_local_refinement"
EXPLORATORY_REFINEMENT_EVIDENCE = "exploratory_discovery_local_refinement"
SOURCE_SEED_CONFIRMATION_REFINEMENT_EVIDENCE = (
    "exploratory_source_seed_confirmation_refinement"
)
FUNNEL_SCREENING_REFINEMENT_EVIDENCE = "exploratory_funnel_screening_refinement"
FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE = "exploratory_funnel_confirmation_refinement"
FUNNEL_DEEP_REFINEMENT_EVIDENCE = "exploratory_funnel_deep_confirmation_refinement"
GUIDED_REFINEMENT_EVIDENCE = "guided_site_compatibility_refinement"
REFINEMENT_PLAN_SCHEMA = "vela.validation-refinement-plan/4"


@dataclass(frozen=True, slots=True)
class RefinementSource:
    """局部精修起点清单及其不可改变的证据身份。"""

    kind: str
    category: str
    manifest_path: Path
    source_evidence_category: str
    evidence_category: str
    known_site_information_used: bool
    production_qualified: bool
    selected_candidate_ids: tuple[str, ...]
    selection: dict[str, JsonValue]


@dataclass(frozen=True, slots=True)
class RefinementStart:
    """一个通过阶段二来源和全原子化学核对的局部精修起点。"""

    start_id: str
    candidate_id: str
    receptor_site_id: str
    pose_id: str
    receptor_id: str
    target: str
    source_seed: int | None
    input_path: Path
    input_sha256: str
    receptor_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class RefinementTask:
    """一个全原子起点与独立局部精修 seed 的组合。"""

    task_id: str
    start: RefinementStart
    seed: int


@dataclass(frozen=True, slots=True)
class RefinementPlan:
    """冻结在独立目录中的正式 candidate 局部精修计划。"""

    run_dir: Path
    tasks: tuple[RefinementTask, ...]
    selected_candidate_ids: tuple[str, ...]
    start_count: int
    total_decoy_count: int


def _optional_seed(value: object) -> int | None:
    if value is None:
        return None
    if not isinstance(value, int) or isinstance(value, bool) or value < 0:
        raise ValidationError("source seed must be a non-negative integer or null")
    return value


def _candidate_ids(raw: object, *, name: str) -> tuple[str, ...]:
    """读取保持顺序且不重复的候选 ID 列表。"""
    try:
        values = object_list(raw, name=name)
    except TypeError as exc:
        raise ValidationError(f"{name} is invalid") from exc
    identifiers = tuple(safe_identifier(value, name=name) for value in values)
    if not identifiers or len(set(identifiers)) != len(identifiers):
        raise ValidationError(f"{name} must contain unique candidate IDs")
    return identifiers


def select_exploration_candidates(
    *,
    config: AppConfig,
    handoff_plan: dict[str, object],
    rows: tuple[dict[str, object], ...],
) -> tuple[tuple[str, ...], dict[str, JsonValue]]:
    """按冻结的双臂规则选择进入完整局部精修的探索候选。"""
    try:
        selection = object_mapping(
            handoff_plan.get("selection"), name="exploration handoff selection"
        )
        arms = object_mapping(
            selection.get("candidate_arms"), name="exploration candidate arms"
        )
        recorded_contract = object_mapping(
            selection.get("promotion_contract"), name="exploration promotion contract"
        )
    except TypeError as exc:
        raise ValidationError("exploration handoff selection is invalid") from exc
    expected_contract = exploration_promotion_contract(config)
    if recorded_contract != expected_contract:
        raise ValidationError("exploration promotion contract has changed")
    blind = _candidate_ids(arms.get("blind_discovery_arm"), name="blind discovery arm")
    functional = _candidate_ids(
        arms.get("functional_annotation_arm"), name="functional annotation arm"
    )
    if set(blind).intersection(functional):
        raise ValidationError("exploration candidate arms overlap")
    requested = _candidate_ids(
        selection.get("requested_candidate_ids"), name="requested candidate IDs"
    )
    arm_candidates = blind + functional
    if requested != arm_candidates:
        raise ValidationError("requested candidates differ from frozen arm order")

    sites_by_candidate: dict[str, set[str]] = {item: set() for item in requested}
    passed_by_site: dict[tuple[str, str], int] = {}
    for row in rows:
        candidate_id = safe_identifier(row.get("candidate_id"), name="candidate ID")
        if candidate_id not in sites_by_candidate:
            raise ValidationError("handoff task candidate is outside frozen arms")
        receptor_site_id = safe_identifier(
            row.get("receptor_site_id"), name="receptor site ID"
        )
        sites_by_candidate[candidate_id].add(receptor_site_id)
        if (
            row.get("execution_status") == "completed"
            and row.get("reconstruction_status") == "passed"
        ):
            key = (candidate_id, receptor_site_id)
            passed_by_site[key] = passed_by_site.get(key, 0) + 1

    contract_eligibility = object_mapping(
        expected_contract["candidate_eligibility"], name="candidate eligibility"
    )
    minimum = nonnegative_integer(
        contract_eligibility.get("minimum_passed_starts_per_receptor_site"),
        name="minimum passed starts per receptor site",
    )
    eligible: list[str] = []
    support: dict[str, JsonValue] = {}
    for candidate_id in requested:
        sites = tuple(sorted(sites_by_candidate[candidate_id]))
        if not sites:
            raise ValidationError("frozen candidate has no handoff receptor sites")
        site_counts = {
            site_id: passed_by_site.get((candidate_id, site_id), 0) for site_id in sites
        }
        is_eligible = all(count >= minimum for count in site_counts.values())
        if is_eligible:
            eligible.append(candidate_id)
        support[candidate_id] = {
            "eligible": is_eligible,
            "passed_starts_by_receptor_site": site_counts,
        }

    deep_selection = object_mapping(
        expected_contract["deep_refinement_selection"],
        name="deep refinement selection",
    )
    blind_budget = nonnegative_integer(
        deep_selection.get("blind_discovery_arm_budget"),
        name="blind discovery arm budget",
    )
    functional_budget = nonnegative_integer(
        deep_selection.get("functional_annotation_arm_budget"),
        name="functional annotation arm budget",
    )
    eligible_set = set(eligible)
    selected_blind = tuple(item for item in blind if item in eligible_set)[
        :blind_budget
    ]
    selected_functional = tuple(item for item in functional if item in eligible_set)[
        :functional_budget
    ]
    selected = selected_blind + selected_functional
    if not selected:
        raise ValidationError(
            "no exploration candidate passed the frozen promotion rule"
        )
    return selected, {
        "mode": "frozen_exploration_promotion",
        "selected_candidate_ids": list(selected),
        "selected_by_arm": {
            "blind_discovery_arm": list(selected_blind),
            "functional_annotation_arm": list(selected_functional),
        },
        "arm_budget": {
            "blind_discovery_arm": blind_budget,
            "functional_annotation_arm": functional_budget,
        },
        "eligible_candidate_ids": eligible,
        "candidate_support": support,
        "promotion_contract": expected_contract,
    }


def select_funnel_screening_candidates(
    *,
    config: AppConfig,
    requested: tuple[str, ...],
    rows: tuple[dict[str, object], ...],
    promotion_contract: JsonValue,
    funnel_audit: JsonValue,
) -> tuple[tuple[str, ...], dict[str, JsonValue]]:
    """只保留每个受体site均有完整双起点证据的Stage 3A候选。"""
    requested_set = set(requested)
    sites_by_candidate: dict[str, set[str]] = {
        candidate_id: set() for candidate_id in requested
    }
    passed_by_site: dict[tuple[str, str], int] = {}
    for row in rows:
        candidate_id = safe_identifier(row.get("candidate_id"), name="candidate ID")
        if candidate_id not in requested_set:
            raise ValidationError("handoff task candidate is outside frozen funnel")
        receptor_site_id = safe_identifier(
            row.get("receptor_site_id"), name="receptor site ID"
        )
        sites_by_candidate[candidate_id].add(receptor_site_id)
        if (
            row.get("execution_status") == "completed"
            and row.get("reconstruction_status") == "passed"
        ):
            key = (candidate_id, receptor_site_id)
            passed_by_site[key] = passed_by_site.get(key, 0) + 1

    minimum = config.validation.funnel.screening_starts_per_receptor_site
    selected: list[str] = []
    excluded: list[str] = []
    support: dict[str, JsonValue] = {}
    for candidate_id in requested:
        sites = tuple(sorted(sites_by_candidate[candidate_id]))
        if not sites:
            raise ValidationError("frozen funnel candidate has no receptor sites")
        site_counts = {
            site_id: passed_by_site.get((candidate_id, site_id), 0) for site_id in sites
        }
        eligible = all(count >= minimum for count in site_counts.values())
        (selected if eligible else excluded).append(candidate_id)
        support[candidate_id] = {
            "eligible": eligible,
            "passed_starts_by_receptor_site": site_counts,
        }
    if not selected:
        raise ValidationError("no funnel candidate passed all-atom handoff QC")
    return tuple(selected), {
        "mode": "stage3a_cross_source_screening",
        "selected_candidate_ids": selected,
        "excluded_candidate_ids": excluded,
        "minimum_passed_starts_per_receptor_site": minimum,
        "candidate_support": support,
        "promotion_contract": promotion_contract,
        "funnel_audit": funnel_audit,
    }


def read_refinement_source(
    *, config: AppConfig, source_run_dir: Path
) -> tuple[RefinementSource, tuple[RefinementStart, ...]]:
    """读取 blind 或 guided 起点; 并由清单决定证据类别。"""
    handoff_path = source_run_dir / "handoff_manifest.json"
    guided_path = source_run_dir / "guided_manifest.json"
    existing = tuple(path for path in (handoff_path, guided_path) if path.is_file())
    if len(existing) != 1:
        raise ValidationError(
            "refinement source must contain exactly one supported completed manifest"
        )
    manifest_path = existing[0]
    manifest = read_document(manifest_path, name="refinement source manifest")
    if manifest_path == handoff_path:
        expected_schema = "vela.validation-handoff-manifest/8"
        source_plan_key = "handoff_plan"
        recorded_source_evidence = manifest.get("evidence_category")
        if recorded_source_evidence == MAIN_HANDOFF_EVIDENCE:
            source_evidence = MAIN_HANDOFF_EVIDENCE
            expected_discovery_evidence = MAIN_DISCOVERY_EVIDENCE
            refinement_evidence = BLIND_REFINEMENT_EVIDENCE
            production_qualified = True
        elif recorded_source_evidence == EXPLORATORY_HANDOFF_EVIDENCE:
            source_evidence = EXPLORATORY_HANDOFF_EVIDENCE
            expected_discovery_evidence = EXPLORATORY_DISCOVERY_EVIDENCE
            refinement_evidence = EXPLORATORY_REFINEMENT_EVIDENCE
            production_qualified = False
        elif recorded_source_evidence == SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE:
            source_evidence = SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE
            expected_discovery_evidence = EXPLORATORY_DISCOVERY_EVIDENCE
            refinement_evidence = SOURCE_SEED_CONFIRMATION_REFINEMENT_EVIDENCE
            production_qualified = False
        elif recorded_source_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE:
            source_evidence = FUNNEL_SCREENING_HANDOFF_EVIDENCE
            expected_discovery_evidence = EXPLORATORY_DISCOVERY_EVIDENCE
            refinement_evidence = FUNNEL_SCREENING_REFINEMENT_EVIDENCE
            production_qualified = False
        else:
            raise ValidationError("unsupported blind handoff evidence category")
    else:
        expected_schema = "vela.validation-guided-manifest/1"
        source_plan_key = "guided_plan"
        source_evidence = GUIDED_EVIDENCE
        expected_discovery_evidence = GUIDED_EVIDENCE
        refinement_evidence = GUIDED_REFINEMENT_EVIDENCE
        production_qualified = False
    if (
        manifest.get("schema") != expected_schema
        or manifest.get("status") != "completed"
        or manifest.get("chemistry_id") != config.chemistry.chemistry_id
        or manifest.get("known_site_information_used")
        is not (manifest_path == guided_path)
        or manifest.get("evidence_category") != source_evidence
    ):
        raise ValidationError("refinement source manifest identity is invalid")
    if manifest_path == handoff_path and (
        manifest.get("source_evidence_category") != expected_discovery_evidence
        or manifest.get("production_qualified") is not production_qualified
    ):
        raise ValidationError("blind handoff authorization is invalid")
    try:
        raw_rows = object_list(manifest.get("tasks"), name="handoff manifest tasks")
        rows = tuple(
            object_mapping(raw, name="handoff manifest task") for raw in raw_rows
        )
    except TypeError as exc:
        raise ValidationError("handoff manifest structure is invalid") from exc
    source_plan_path, _ = validate_record(
        root=source_run_dir,
        raw=manifest.get(source_plan_key),
        name="refinement source plan",
    )
    source_plan = read_document(source_plan_path, name="refinement source plan")
    selected_candidate_ids: tuple[str, ...] = ()
    selection: dict[str, JsonValue] = {
        "mode": "all_passed_starts",
        "selected_candidate_ids": [],
    }
    if source_evidence in {
        EXPLORATORY_HANDOFF_EVIDENCE,
        SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE,
        FUNNEL_SCREENING_HANDOFF_EVIDENCE,
    }:
        if (
            source_plan.get("schema") != HANDOFF_PLAN_SCHEMA
            or source_plan.get("evidence_category") != source_evidence
            or source_plan.get("source_evidence_category")
            != expected_discovery_evidence
            or source_plan.get("production_qualified") is not False
            or source_plan.get("known_site_information_used") is not False
            or source_plan.get("chemistry_id") != config.chemistry.chemistry_id
        ):
            raise ValidationError("exploration handoff plan identity is invalid")
        if source_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE:
            try:
                source_selection = object_mapping(
                    source_plan.get("selection"), name="funnel handoff selection"
                )
            except TypeError as exc:
                raise ValidationError("funnel handoff selection is invalid") from exc
            requested_candidate_ids = _candidate_ids(
                source_selection.get("requested_candidate_ids"),
                name="funnel candidate IDs",
            )
            if source_selection.get("promotion_contract") != funnel_screening_contract(
                config
            ):
                raise ValidationError("funnel screening contract has changed")
            try:
                promotion_contract = json_value(
                    source_selection.get("promotion_contract"),
                    name="funnel promotion contract",
                )
                funnel_audit = json_value(
                    source_selection.get("funnel_audit"), name="funnel audit"
                )
            except TypeError as exc:
                raise ValidationError("funnel screening evidence is invalid") from exc
            selected_candidate_ids, selection = select_funnel_screening_candidates(
                config=config,
                requested=requested_candidate_ids,
                rows=rows,
                promotion_contract=promotion_contract,
                funnel_audit=funnel_audit,
            )
        elif source_evidence == EXPLORATORY_HANDOFF_EVIDENCE:
            selected_candidate_ids, selection = select_exploration_candidates(
                config=config,
                handoff_plan=source_plan,
                rows=rows,
            )
        else:
            try:
                source_selection = object_mapping(
                    source_plan.get("selection"), name="confirmation handoff selection"
                )
            except TypeError as exc:
                raise ValidationError(
                    "source-seed confirmation selection is invalid"
                ) from exc
            if source_selection.get(
                "promotion_contract"
            ) != source_seed_confirmation_contract(config):
                raise ValidationError("source-seed confirmation contract has changed")
            selected_candidate_ids = _candidate_ids(
                source_selection.get("requested_candidate_ids"),
                name="confirmation candidate IDs",
            )
            selection = {
                "mode": "all_passed_source_seed_confirmation_starts",
                "selected_candidate_ids": list(selected_candidate_ids),
                "promotion_contract": source_seed_confirmation_contract(config),
            }

    starts: list[RefinementStart] = []
    for row in rows:
        start_id = safe_identifier(row.get("task_id"), name="handoff task ID")
        validate_record(
            root=source_run_dir,
            raw=row.get("task_result"),
            name=f"{start_id} task result",
        )
        if manifest_path == handoff_path:
            execution_status = row.get("execution_status")
            reconstruction_status = row.get("reconstruction_status")
            if execution_status == "invalid":
                raise ValidationError("invalid handoff task cannot be refined")
            if execution_status != "completed":
                raise ValidationError("handoff execution status is invalid")
            if reconstruction_status == "failed":
                if row.get("flexpepdock_input") is not None:
                    raise ValidationError(
                        "failed handoff task must not declare a FlexPepDock input"
                    )
                continue
            if reconstruction_status != "passed":
                raise ValidationError("handoff reconstruction status is invalid")
            candidate_id = safe_identifier(row.get("candidate_id"), name="candidate ID")
            if selected_candidate_ids and candidate_id not in selected_candidate_ids:
                continue
        input_path, input_hash = validate_record(
            root=source_run_dir,
            raw=row.get("flexpepdock_input"),
            name=f"{start_id} FlexPepDock input",
        )
        receptor_count, histidines = validate_flexpepdock_input(
            path=input_path,
            chemistry=config.chemistry,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        )
        starts.append(
            RefinementStart(
                start_id=start_id,
                candidate_id=safe_identifier(
                    row.get("candidate_id"), name="candidate ID"
                ),
                receptor_site_id=safe_identifier(
                    row.get("receptor_site_id"), name="receptor site ID"
                ),
                pose_id=safe_identifier(row.get("pose_id"), name="pose ID"),
                receptor_id=safe_identifier(row.get("receptor_id"), name="receptor ID"),
                target=safe_identifier(row.get("target"), name="target ID"),
                source_seed=_optional_seed(row.get("source_seed")),
                input_path=input_path,
                input_sha256=input_hash,
                receptor_residue_count=receptor_count,
                fixed_histidine_pose_indices=histidines,
            )
        )
    if not starts or len({item.start_id for item in starts}) != len(starts):
        raise ValidationError("handoff manifest contains no unique starts")
    if not selected_candidate_ids:
        selected_candidate_ids = tuple(
            dict.fromkeys(item.candidate_id for item in starts)
        )
        selection["selected_candidate_ids"] = list(selected_candidate_ids)
    source = RefinementSource(
        kind="guided_handoff" if manifest_path == guided_path else "blind_handoff",
        category="guided" if manifest_path == guided_path else "handoffs",
        manifest_path=manifest_path,
        source_evidence_category=source_evidence,
        evidence_category=refinement_evidence,
        known_site_information_used=manifest_path == guided_path,
        production_qualified=production_qualified,
        selected_candidate_ids=selected_candidate_ids,
        selection=selection,
    )
    return source, tuple(starts)


def build_refinement_tasks(
    *, config: AppConfig, source_run_dir: Path
) -> tuple[RefinementTask, ...]:
    """将每个合格起点与全部已冻结局部精修 seed 组合。"""
    source, starts = read_refinement_source(
        config=config, source_run_dir=source_run_dir
    )
    return _tasks_from_starts(
        starts=starts,
        seeds=_refinement_seeds(config=config, source=source),
    )


def _tasks_from_starts(
    *,
    starts: tuple[RefinementStart, ...],
    seeds: tuple[int, ...],
    task_prefix: str = "refine",
) -> tuple[RefinementTask, ...]:
    """从已经校验的起点构造确定性的 seed 任务笛卡尔积。"""
    tasks: list[RefinementTask] = []
    index = 1
    for start in starts:
        for seed in seeds:
            tasks.append(RefinementTask(f"{task_prefix}_{index:05d}", start, seed))
            index += 1
    if not tasks:
        raise ValidationError("Stage 3 production seeds are not frozen")
    return tuple(tasks)


def _refinement_seeds(
    *, config: AppConfig, source: RefinementSource
) -> tuple[int, ...]:
    """按证据来源选择完整协议或Stage 3A首批随机流。"""
    if source.source_evidence_category != FUNNEL_SCREENING_HANDOFF_EVIDENCE:
        return config.validation.seeds
    count = config.validation.refinement.seed_batch_sizes[0]
    return config.validation.seeds[:count]


def refinement_parameters(
    config: AppConfig, *, seeds: tuple[int, ...]
) -> dict[str, JsonValue]:
    """返回写入计划并在运行前复核的全部局部精修参数。"""
    settings = config.validation
    return {
        "prepack_seed": settings.refinement.prepack_seed,
        "random_translation_A": settings.refinement.random_translation_A,
        "random_rotation_degrees": settings.refinement.random_rotation_degrees,
        "ranking_score": settings.refinement.ranking_score,
        "seeds": list(seeds),
        "parallel_tasks": settings.rosetta.parallel_tasks,
        "decoys_per_seed": settings.rosetta.decoys_per_seed,
        "score_function": settings.rosetta.score_function,
        "lowres_preoptimize": settings.rosetta.lowres_preoptimize,
        "application": "rosetta_scripts_with_in_pose_chemistry_restoration",
    }


def refinement_identity(plan: dict[str, object]) -> tuple[str, bool]:
    """从已验证计划中收窄证据类别和已知位点标记。"""
    evidence = plan.get("evidence_category")
    known = plan.get("known_site_information_used")
    if not isinstance(evidence, str) or not isinstance(known, bool):
        raise ValidationError("refinement evidence identity is invalid")
    return evidence, known


def refinement_authorization(plan: dict[str, object]) -> tuple[str, bool]:
    """从已验证计划中读取上游证据和生产资格状态。"""
    source_evidence = plan.get("source_evidence_category")
    production_qualified = plan.get("production_qualified")
    if not isinstance(source_evidence, str) or not isinstance(
        production_qualified, bool
    ):
        raise ValidationError("refinement authorization is invalid")
    return source_evidence, production_qualified


def _relative_output(path: Path, *, config: AppConfig, name: str) -> str:
    try:
        return path.resolve().relative_to(config.paths.outputs_dir.resolve()).as_posix()
    except ValueError as exc:
        raise ValidationError(
            f"{name} is outside the configured outputs directory"
        ) from exc


def _require_local_refinement_ready(config: AppConfig) -> None:
    readiness = assess_validation_readiness(config)
    if readiness.production_ready:
        return
    local_issues = tuple(
        issue.code for issue in readiness.issues if not issue.code.startswith("global_")
    )
    raise ValidationError(
        "Stage 3 candidate refinement is not ready: " + "; ".join(local_issues)
    )


def _write_refinement_plan(
    *,
    config: AppConfig,
    source: RefinementSource,
    starts: tuple[RefinementStart, ...],
    seeds: tuple[int, ...],
    run_id: str,
    task_prefix: str,
    additional_inputs: dict[str, JsonValue] | None = None,
) -> RefinementPlan:
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    _require_local_refinement_ready(config)
    tasks = _tasks_from_starts(
        starts=starts,
        seeds=seeds,
        task_prefix=task_prefix,
    )
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    run_dir = config.paths.outputs_dir / "validation" / "refinements" / run_id
    if run_dir.exists():
        raise ValidationError(f"refinement run directory already exists: {run_dir}")
    report = config.validation.qualification_report
    if report is None:
        raise ValidationError("qualified method lacks a qualification report")
    report_relative = _relative_output(
        report, config=config, name="qualification report"
    )
    source_relative = _relative_output(
        source.manifest_path, config=config, name="refinement source manifest"
    )
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    inputs: dict[str, JsonValue] = {
        "config_snapshot": file_record(snapshot, root=run_dir),
        "qualification_report": {
            "path": report_relative,
            "sha256": sha256_file(report),
        },
        "source_manifest": {
            "kind": source.kind,
            "category": source.category,
            "path": source_relative,
            "sha256": sha256_file(source.manifest_path),
        },
    }
    if additional_inputs is not None:
        inputs.update(additional_inputs)
    atomic_write_json(
        run_dir / REFINEMENT_PLAN_NAME,
        {
            "schema": REFINEMENT_PLAN_SCHEMA,
            "stage": "validation_local_refinement",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "source_evidence_category": source.source_evidence_category,
            "evidence_category": source.evidence_category,
            "known_site_information_used": source.known_site_information_used,
            "production_qualified": source.production_qualified,
            "software": {
                **vela_software_identity(),
                "rosetta_version": flexpepdock.version,
                "flexpepdock_sha256": flexpepdock.executable_sha256,
                "rosetta_scripts_sha256": scripts.executable_sha256,
            },
            "inputs": inputs,
            "parameters": refinement_parameters(config, seeds=seeds),
            "selection": source.selection,
            "budget": {
                "candidate_count": len(source.selected_candidate_ids),
                "start_count": len(starts),
                "seed_count": len(seeds),
                "task_count": len(tasks),
                "decoys_per_task": config.validation.rosetta.decoys_per_seed,
                "total_decoy_count": (
                    len(tasks) * config.validation.rosetta.decoys_per_seed
                ),
            },
            "tasks": [
                {
                    "task_id": task.task_id,
                    "start_id": task.start.start_id,
                    "candidate_id": task.start.candidate_id,
                    "receptor_site_id": task.start.receptor_site_id,
                    "pose_id": task.start.pose_id,
                    "receptor_id": task.start.receptor_id,
                    "target": task.start.target,
                    "source_seed": task.start.source_seed,
                    "refinement_seed": task.seed,
                    "input_sha256": task.start.input_sha256,
                    "status": "planned",
                }
                for task in tasks
            ],
        },
    )
    return RefinementPlan(
        run_dir=run_dir,
        tasks=tasks,
        selected_candidate_ids=source.selected_candidate_ids,
        start_count=len(starts),
        total_decoy_count=len(tasks) * config.validation.rosetta.decoys_per_seed,
    )


def write_refinement_plan(
    *, config: AppConfig, source_run_dir: Path, run_id: str
) -> RefinementPlan:
    """冻结已放行方法、显式起点来源、seed 和局部精修参数。"""
    source, starts = read_refinement_source(
        config=config, source_run_dir=source_run_dir
    )
    seeds = _refinement_seeds(config=config, source=source)
    return _write_refinement_plan(
        config=config,
        source=source,
        starts=starts,
        seeds=seeds,
        run_id=run_id,
        task_prefix="refine",
    )


def _source_manifest(
    *, config: AppConfig, record: dict[str, object], category: str, name: str
) -> Path:
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError(f"invalid {name} record")
    path = (config.paths.outputs_dir / relative).resolve()
    root = (config.paths.outputs_dir / "validation" / category).resolve()
    try:
        path.relative_to(root)
    except ValueError as exc:
        raise ValidationError(f"{name} is outside {category} runs") from exc
    if not path.is_file() or sha256_file(path) != expected_hash:
        raise ValidationError(f"{name} hash mismatch")
    return path


def _funnel_confirmation_source(
    *, config: AppConfig, screening_run_dir: Path
) -> tuple[RefinementSource, tuple[RefinementStart, ...], Path]:
    """从完整Stage 3A分析中确定性读取Stage 3B候选和原起点。"""
    screening_plan, screening_tasks = verify_refinement_plan(
        config=config,
        run_dir=screening_run_dir,
        require_current_software=False,
    )
    if (
        screening_plan.get("source_evidence_category")
        != FUNNEL_SCREENING_HANDOFF_EVIDENCE
        or screening_plan.get("evidence_category")
        != FUNNEL_SCREENING_REFINEMENT_EVIDENCE
        or screening_plan.get("production_qualified") is not False
    ):
        raise ValidationError("Stage 3B source is not a Stage 3A screening run")
    analysis_path = screening_run_dir / "refinement_analysis" / "analysis_manifest.json"
    analysis = read_document(analysis_path, name="Stage 3A analysis manifest")
    if (
        analysis.get("schema") != "vela.validation-refinement-analysis-manifest/4"
        or analysis.get("status") != "completed"
        or analysis.get("source_evidence_category") != FUNNEL_SCREENING_HANDOFF_EVIDENCE
        or analysis.get("evidence_category") != FUNNEL_SCREENING_REFINEMENT_EVIDENCE
        or analysis.get("known_site_information_used") is not False
        or analysis.get("production_qualified") is not False
        or not is_vela_software_identity(analysis.get("sampling_software"))
        or not is_vela_software_identity(analysis.get("analysis_software"))
    ):
        raise ValidationError("Stage 3A analysis identity is invalid")
    try:
        cluster_record = object_mapping(
            analysis.get("refined_clusters"), name="Stage 3A cluster record"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3A cluster record is invalid") from exc
    cluster_path, _ = validate_record(
        root=analysis_path.parent,
        raw=cluster_record,
        name="Stage 3A refined clusters",
    )
    with cluster_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = [
            "cluster_id",
            "candidate_id",
            "receptor_id",
            "decoy_count",
            "refinement_seeds",
            "start_ids",
            "source_seeds",
            "task_cells",
            "evidence_status",
            "representative_decoy_id",
            "representative_path",
            "supported",
        ]
        if reader.fieldnames != expected_fields:
            raise ValidationError("Stage 3A refined cluster columns are invalid")
        rows = tuple(dict(row) for row in reader)
    if len(rows) != cluster_record.get("count"):
        raise ValidationError("Stage 3A refined cluster count is inconsistent")
    hit_clusters: dict[str, list[str]] = {}
    for row in rows:
        if row["evidence_status"] != "cross_source_screening_hit":
            continue
        candidate_id = safe_identifier(
            row["candidate_id"], name="Stage 3A hit candidate ID"
        )
        cluster_id = safe_identifier(row["cluster_id"], name="Stage 3A hit cluster ID")
        hit_clusters.setdefault(candidate_id, []).append(cluster_id)
    try:
        screening_selection = object_mapping(
            screening_plan.get("selection"), name="Stage 3A selection"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3A selection is invalid") from exc
    screening_candidates = _candidate_ids(
        screening_selection.get("selected_candidate_ids"),
        name="Stage 3A selected candidate IDs",
    )
    selected = tuple(
        candidate_id
        for candidate_id in screening_candidates
        if candidate_id in hit_clusters
    )[: config.validation.funnel.screening_promotion_budget]
    if not selected:
        raise ValidationError("Stage 3A has no cross-source screening hit")
    unknown_hits = set(hit_clusters) - set(screening_candidates)
    if unknown_hits:
        raise ValidationError("Stage 3A hit is outside the frozen screening candidates")
    starts = tuple(
        {
            task.start.start_id: task.start
            for task in screening_tasks
            if task.start.candidate_id in selected
        }.values()
    )
    expected_starts = config.validation.funnel.screening_starts_per_receptor_site
    counts: dict[tuple[str, str], int] = {}
    for start in starts:
        key = (start.candidate_id, start.receptor_site_id)
        counts[key] = counts.get(key, 0) + 1
    if not starts or any(count != expected_starts for count in counts.values()):
        raise ValidationError(
            "Stage 3B source does not preserve the frozen start budget"
        )
    try:
        inputs = object_mapping(screening_plan.get("inputs"), name="Stage 3A inputs")
        source_record = object_mapping(
            inputs.get("source_manifest"), name="Stage 3A source manifest"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3A source manifest record is invalid") from exc
    source_manifest = _source_manifest(
        config=config,
        record=source_record,
        category="handoffs",
        name="Stage 3A source manifest",
    )
    original_source, _ = read_refinement_source(
        config=config,
        source_run_dir=source_manifest.parent,
    )
    source = RefinementSource(
        kind=original_source.kind,
        category=original_source.category,
        manifest_path=original_source.manifest_path,
        source_evidence_category=original_source.source_evidence_category,
        evidence_category=FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE,
        known_site_information_used=False,
        production_qualified=False,
        selected_candidate_ids=selected,
        selection={
            "mode": "stage3b_source_seed_confirmation",
            "selected_candidate_ids": list(selected),
            "screening_hit_cluster_ids": {
                candidate_id: hit_clusters[candidate_id] for candidate_id in selected
            },
            "promotion_contract": funnel_screening_contract(config),
        },
    )
    return source, starts, analysis_path


def write_funnel_confirmation_plan(
    *, config: AppConfig, screening_run_dir: Path, run_id: str
) -> RefinementPlan:
    """只为Stage 3A命中起点冻结第二条Rosetta随机流。"""
    source, starts, analysis_path = _funnel_confirmation_source(
        config=config,
        screening_run_dir=screening_run_dir,
    )
    first_batch, second_batch, _ = config.validation.refinement.seed_batch_sizes
    seeds = config.validation.seeds[first_batch : first_batch + second_batch]
    return _write_refinement_plan(
        config=config,
        source=source,
        starts=starts,
        seeds=seeds,
        run_id=run_id,
        task_prefix="confirm",
        additional_inputs={
            "screening_analysis": {
                "path": _relative_output(
                    analysis_path,
                    config=config,
                    name="Stage 3A analysis manifest",
                ),
                "sha256": sha256_file(analysis_path),
            }
        },
    )


def _funnel_deep_source(
    *, config: AppConfig, confirmation_run_dir: Path
) -> tuple[RefinementSource, tuple[RefinementStart, ...], Path]:
    """从Stage 3B确认报告确定性读取Stage 3C候选和原起点。"""
    confirmation_plan, confirmation_tasks = verify_refinement_plan(
        config=config,
        run_dir=confirmation_run_dir,
        require_current_software=False,
    )
    if (
        confirmation_plan.get("source_evidence_category")
        != FUNNEL_SCREENING_HANDOFF_EVIDENCE
        or confirmation_plan.get("evidence_category")
        != FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE
        or confirmation_plan.get("production_qualified") is not False
    ):
        raise ValidationError("Stage 3C source is not a Stage 3B confirmation run")
    analysis_path = (
        confirmation_run_dir / "funnel_confirmation_analysis" / "analysis_manifest.json"
    )
    analysis = read_document(analysis_path, name="Stage 3B analysis manifest")
    if (
        analysis.get("schema") != "vela.validation-funnel-confirmation-analysis/1"
        or analysis.get("stage") != "stage3b_source_seed_confirmation"
        or analysis.get("status") != "completed"
        or analysis.get("evidence_category") != FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE
        or analysis.get("production_qualified") is not False
        or not is_vela_software_identity(analysis.get("analysis_software"))
    ):
        raise ValidationError("Stage 3B analysis identity is invalid")
    selected = _candidate_ids(
        analysis.get("confirmed_candidate_ids"),
        name="Stage 3B confirmed candidate IDs",
    )
    try:
        cluster_record = object_mapping(
            analysis.get("combined_refined_clusters"),
            name="Stage 3B combined cluster record",
        )
    except TypeError as exc:
        raise ValidationError("Stage 3B combined cluster record is invalid") from exc
    cluster_path, _ = validate_record(
        root=analysis_path.parent,
        raw=cluster_record,
        name="Stage 3B combined refined clusters",
    )
    with cluster_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        expected_fields = [
            "cluster_id",
            "candidate_id",
            "receptor_id",
            "decoy_count",
            "refinement_seeds",
            "start_ids",
            "source_seeds",
            "task_cells",
            "evidence_status",
            "representative_decoy_id",
            "representative_path",
        ]
        if reader.fieldnames != expected_fields:
            raise ValidationError("Stage 3B combined cluster columns are invalid")
        rows = tuple(dict(row) for row in reader)
    if len(rows) != cluster_record.get("count"):
        raise ValidationError("Stage 3B combined cluster count is inconsistent")
    confirmed_clusters: dict[str, list[str]] = {}
    for row in rows:
        if row["evidence_status"] != "source_seed_confirmation_hit":
            continue
        candidate_id = safe_identifier(
            row["candidate_id"], name="Stage 3B confirmed candidate ID"
        )
        cluster_id = safe_identifier(
            row["cluster_id"], name="Stage 3B confirmed cluster ID"
        )
        confirmed_clusters.setdefault(candidate_id, []).append(cluster_id)
    if set(confirmed_clusters) != set(selected):
        raise ValidationError("Stage 3B confirmed candidates and clusters differ")
    starts = tuple(
        {
            task.start.start_id: task.start
            for task in confirmation_tasks
            if task.start.candidate_id in set(selected)
        }.values()
    )
    if not starts:
        raise ValidationError("Stage 3B contains no starts for Stage 3C")
    try:
        inputs = object_mapping(confirmation_plan.get("inputs"), name="Stage 3B inputs")
        source_record = object_mapping(
            inputs.get("source_manifest"), name="Stage 3B source manifest"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3B source manifest record is invalid") from exc
    source_manifest = _source_manifest(
        config=config,
        record=source_record,
        category="handoffs",
        name="Stage 3B source manifest",
    )
    original_source, _ = read_refinement_source(
        config=config,
        source_run_dir=source_manifest.parent,
    )
    source = RefinementSource(
        kind=original_source.kind,
        category=original_source.category,
        manifest_path=original_source.manifest_path,
        source_evidence_category=original_source.source_evidence_category,
        evidence_category=FUNNEL_DEEP_REFINEMENT_EVIDENCE,
        known_site_information_used=False,
        production_qualified=False,
        selected_candidate_ids=selected,
        selection={
            "mode": "stage3c_deep_confirmation",
            "selected_candidate_ids": list(selected),
            "confirmation_hit_cluster_ids": {
                candidate_id: confirmed_clusters[candidate_id]
                for candidate_id in selected
            },
            "promotion_contract": funnel_screening_contract(config),
        },
    )
    return source, starts, analysis_path


def write_funnel_deep_plan(
    *, config: AppConfig, confirmation_run_dir: Path, run_id: str
) -> RefinementPlan:
    """为Stage 3B确认候选冻结第3和第4条Rosetta随机流。"""
    source, starts, analysis_path = _funnel_deep_source(
        config=config,
        confirmation_run_dir=confirmation_run_dir,
    )
    first_batch, second_batch, deep_batch = (
        config.validation.refinement.seed_batch_sizes
    )
    offset = first_batch + second_batch
    seeds = config.validation.seeds[offset : offset + deep_batch]
    return _write_refinement_plan(
        config=config,
        source=source,
        starts=starts,
        seeds=seeds,
        run_id=run_id,
        task_prefix="deep",
        additional_inputs={
            "confirmation_analysis": {
                "path": _relative_output(
                    analysis_path,
                    config=config,
                    name="Stage 3B analysis manifest",
                ),
                "sha256": sha256_file(analysis_path),
            }
        },
    )


def verify_refinement_plan(
    *, config: AppConfig, run_dir: Path, require_current_software: bool = True
) -> tuple[dict[str, object], tuple[RefinementTask, ...]]:
    """复核计划快照、资格报告、起点来源、任务和工具身份。"""
    if require_current_software:
        readiness = assess_validation_readiness(config)
        if not readiness.production_ready:
            raise ValidationError("Stage 3 candidate refinement is no longer ready")
    plan = read_document(run_dir / REFINEMENT_PLAN_NAME, name="refinement plan")
    software = plan.get("software")
    software_valid = (
        is_current_vela_software(software)
        if require_current_software
        else is_vela_software_identity(software)
    )
    if (
        plan.get("schema") != REFINEMENT_PLAN_SCHEMA
        or plan.get("stage") != "validation_local_refinement"
        or plan.get("status") != "planned"
        or not software_valid
        or plan.get("method_id") != config.validation.method_id
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
    ):
        raise ValidationError("refinement plan identity is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="refinement inputs")
        rows = object_list(plan.get("tasks"), name="refinement tasks")
    except TypeError as exc:
        raise ValidationError("refinement plan structure is invalid") from exc
    snapshot_path, _ = validate_record(
        root=run_dir, raw=inputs.get("config_snapshot"), name="config snapshot"
    )
    if sha256_file(snapshot_path) != config.source_snapshot_sha256:
        raise ValidationError("current project config differs from refinement plan")
    try:
        source_record = object_mapping(
            inputs.get("source_manifest"), name="refinement source manifest"
        )
        report_record = object_mapping(
            inputs.get("qualification_report"), name="qualification report"
        )
    except TypeError as exc:
        raise ValidationError("refinement source records are invalid") from exc
    source_kind = source_record.get("kind")
    source_category = source_record.get("category")
    resolved_category: str
    if source_kind == "blind_handoff" and source_category == "handoffs":
        resolved_category = "handoffs"
        source_evidence = plan.get("source_evidence_category")
        if source_evidence == MAIN_HANDOFF_EVIDENCE:
            expected_evidence = BLIND_REFINEMENT_EVIDENCE
            expected_production = True
        elif source_evidence == EXPLORATORY_HANDOFF_EVIDENCE:
            expected_evidence = EXPLORATORY_REFINEMENT_EVIDENCE
            expected_production = False
        elif source_evidence == SOURCE_SEED_CONFIRMATION_HANDOFF_EVIDENCE:
            expected_evidence = SOURCE_SEED_CONFIRMATION_REFINEMENT_EVIDENCE
            expected_production = False
        elif source_evidence == FUNNEL_SCREENING_HANDOFF_EVIDENCE:
            recorded_evidence = plan.get("evidence_category")
            if recorded_evidence not in {
                FUNNEL_SCREENING_REFINEMENT_EVIDENCE,
                FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE,
                FUNNEL_DEEP_REFINEMENT_EVIDENCE,
            }:
                raise ValidationError("refinement funnel evidence identity is invalid")
            expected_evidence = recorded_evidence
            expected_production = False
        else:
            raise ValidationError("refinement source evidence identity is invalid")
        expected_known = False
    elif source_kind == "guided_handoff" and source_category == "guided":
        resolved_category = "guided"
        source_evidence = GUIDED_EVIDENCE
        expected_evidence = GUIDED_REFINEMENT_EVIDENCE
        expected_known = True
        expected_production = False
    else:
        raise ValidationError("refinement source identity is invalid")
    if (
        plan.get("source_evidence_category") != source_evidence
        or plan.get("evidence_category") != expected_evidence
        or plan.get("known_site_information_used") is not expected_known
        or plan.get("production_qualified") is not expected_production
    ):
        raise ValidationError("refinement evidence identity is invalid")
    source_manifest = _source_manifest(
        config=config,
        record=source_record,
        category=resolved_category,
        name="refinement source manifest",
    )
    _source_manifest(
        config=config,
        record=report_record,
        category="controls",
        name="qualification report",
    )
    source, starts = read_refinement_source(
        config=config, source_run_dir=source_manifest.parent
    )
    task_prefix = "refine"
    if expected_evidence == FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE:
        try:
            screening_record = object_mapping(
                inputs.get("screening_analysis"),
                name="Stage 3A analysis manifest",
            )
        except TypeError as exc:
            raise ValidationError("Stage 3A analysis record is invalid") from exc
        screening_analysis = _source_manifest(
            config=config,
            record=screening_record,
            category="refinements",
            name="Stage 3A analysis manifest",
        )
        source, starts, expected_analysis = _funnel_confirmation_source(
            config=config,
            screening_run_dir=screening_analysis.parent.parent,
        )
        if screening_analysis != expected_analysis:
            raise ValidationError("Stage 3A analysis path is invalid")
        first_batch, second_batch, _ = config.validation.refinement.seed_batch_sizes
        seeds = config.validation.seeds[first_batch : first_batch + second_batch]
        task_prefix = "confirm"
    elif expected_evidence == FUNNEL_DEEP_REFINEMENT_EVIDENCE:
        try:
            confirmation_record = object_mapping(
                inputs.get("confirmation_analysis"),
                name="Stage 3B analysis manifest",
            )
        except TypeError as exc:
            raise ValidationError("Stage 3B analysis record is invalid") from exc
        confirmation_analysis = _source_manifest(
            config=config,
            record=confirmation_record,
            category="refinements",
            name="Stage 3B analysis manifest",
        )
        source, starts, expected_analysis = _funnel_deep_source(
            config=config,
            confirmation_run_dir=confirmation_analysis.parent.parent,
        )
        if confirmation_analysis != expected_analysis:
            raise ValidationError("Stage 3B analysis path is invalid")
        first_batch, second_batch, deep_batch = (
            config.validation.refinement.seed_batch_sizes
        )
        offset = first_batch + second_batch
        seeds = config.validation.seeds[offset : offset + deep_batch]
        task_prefix = "deep"
    else:
        seeds = _refinement_seeds(config=config, source=source)
    if (
        source.kind != source_kind
        or source.category != resolved_category
        or source.manifest_path != source_manifest
        or source.source_evidence_category != source_evidence
        or source.evidence_category != expected_evidence
        or source.known_site_information_used is not expected_known
        or source.production_qualified is not expected_production
    ):
        raise ValidationError("current refinement source differs from the frozen plan")
    tasks = _tasks_from_starts(
        starts=starts,
        seeds=seeds,
        task_prefix=task_prefix,
    )
    recorded: list[tuple[str, str, int]] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="refinement task")
        except TypeError as exc:
            raise ValidationError("refinement task is invalid") from exc
        if row.get("status") != "planned":
            raise ValidationError("refinement task status is invalid")
        recorded.append(
            (
                safe_identifier(row.get("task_id"), name="refinement task ID"),
                safe_identifier(row.get("start_id"), name="refinement start ID"),
                nonnegative_integer(row.get("refinement_seed"), name="refinement seed"),
            )
        )
    expected = tuple((task.task_id, task.start.start_id, task.seed) for task in tasks)
    expected_budget = {
        "candidate_count": len(source.selected_candidate_ids),
        "start_count": len({task.start.start_id for task in tasks}),
        "seed_count": len(seeds),
        "task_count": len(tasks),
        "decoys_per_task": config.validation.rosetta.decoys_per_seed,
        "total_decoy_count": len(tasks) * config.validation.rosetta.decoys_per_seed,
    }
    if (
        tuple(recorded) != expected
        or plan.get("parameters") != refinement_parameters(config, seeds=seeds)
        or plan.get("selection") != source.selection
        or plan.get("budget") != expected_budget
    ):
        raise ValidationError("current refinement tasks differ from frozen plan")
    if require_current_software:
        verify_flexpepdock_tool(config.validation.rosetta)
        verify_rosetta_scripts_tool(config.validation.rosetta)
    return plan, tasks
