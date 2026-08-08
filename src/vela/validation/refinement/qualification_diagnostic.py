"""不同受体构象之间的 native-aware cross-docking 开发诊断。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.errors import VelaError
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    is_current_vela_software,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.run_identity import validate_run_id
from vela.core.typed_data import object_list, object_mapping
from vela.preparation.chemistry import (
    ChemistryDefinition,
    chemistry_identity_record,
)
from vela.validation.bound_states.assets import standard_peptide_chemistry
from vela.validation.models import LocalRecoveryControl, ValidationError
from vela.validation.records import (
    file_record,
    read_document,
    resume_completed_result,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.receptor_flexibility import (
    ReceptorBackboneMode,
    resolve_receptor_backbone_mode,
    select_local_receptor_backbone,
    write_local_receptor_movemap,
)
from vela.validation.refinement.reconstruction import (
    validate_flexpepdock_input,
    write_chemistry_prepack_protocol,
    write_chemistry_production_refine_protocol,
    write_disulfide_indices,
)
from vela.validation.rosetta import (
    build_chemistry_flexpepdock_command,
    rosetta_crash_log_dir,
    run_rosetta_command,
    single_rosetta_pdb_output,
    verify_rosetta_scripts_tool,
)
from vela.validation.scores import (
    index_rosetta_pdb_outputs,
    read_rosetta_scorefile,
    write_refinement_decoy_manifest,
)

PLAN_NAME = "qualification_refinement_plan.json"
MANIFEST_NAME = "qualification_refinement_manifest.json"
PLAN_SCHEMA = "vela.validation-qualification-refinement-plan/4"
MANIFEST_SCHEMA = "vela.validation-qualification-refinement-manifest/3"
TASK_SCHEMA = "vela.validation-qualification-refinement-task-result/2"
PREPACK_SCHEMA = "vela.validation-qualification-refinement-prepack-result/2"
EVIDENCE_CATEGORY = "cross_receptor_pose_robustness_diagnostic"


@dataclass(frozen=True, slots=True)
class DiagnosticStart:
    """一个因事后 native 信息而被选中的全原子开发起点。"""

    start_id: str
    receptor_site_id: str
    receptor_id: str
    target: str
    source_seed: int
    input_path: Path
    input_sha256: str
    receptor_residue_count: int
    fixed_histidine_pose_indices: tuple[int, ...]
    direct_contact_receptor_pose_indices: tuple[int, ...]
    flexible_receptor_pose_indices: tuple[int, ...]


@dataclass(frozen=True, slots=True)
class DiagnosticTask:
    """一个开发起点与正式阶段三随机流的组合。"""

    task_id: str
    start: DiagnosticStart
    refinement_seed: int


@dataclass(frozen=True, slots=True)
class DiagnosticPlan:
    """已经冻结的 native-aware 恢复诊断。"""

    run_dir: Path
    tasks: tuple[DiagnosticTask, ...]


@dataclass(frozen=True, slots=True)
class DiagnosticOutcome:
    """完成诊断采样后的任务清单。"""

    manifest_path: Path
    task_count: int


@dataclass(frozen=True, slots=True)
class PreparedDiagnosticStart:
    """预打包输入、二硫键和可选局部受体 MoveMap。"""

    prepacked_path: Path
    disulfide_path: Path
    movemap_path: Path | None


def _control_definition(
    config: AppConfig,
) -> tuple[LocalRecoveryControl, ChemistryDefinition]:
    if len(config.validation.local_controls) != 1:
        raise ValidationError("qualification diagnostic requires one local control")
    control = config.validation.local_controls[0]
    states = {item.state_id: item for item in config.validation.bound_states}
    definition = states.get(control.bound_state_id)
    if definition is None or definition.local_control_kind != "standard_cyclic_peptide":
        raise ValidationError(
            "qualification diagnostic control is not a standard peptide"
        )
    return control, standard_peptide_chemistry(definition)


def _output_record(path: Path, *, config: AppConfig, name: str) -> dict[str, str]:
    resolved = path.resolve()
    try:
        relative = resolved.relative_to(config.paths.outputs_dir.resolve())
    except ValueError as exc:
        raise ValidationError(
            f"{name} is outside the configured outputs directory"
        ) from exc
    return {"path": relative.as_posix(), "sha256": sha256_file(resolved)}


def _resolve_output_record(
    raw: object, *, config: AppConfig, root: Path, name: str
) -> Path:
    try:
        record = object_mapping(raw, name=name)
    except TypeError as exc:
        raise ValidationError(f"{name} is invalid") from exc
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError(f"{name} is invalid")
    path = (config.paths.outputs_dir / relative).resolve()
    if not path.is_relative_to(root.resolve()):
        raise ValidationError(f"{name} is outside its allowed run root")
    if not path.is_file() or sha256_file(path) != expected_hash:
        raise ValidationError(f"{name} is missing or has changed")
    return path


def _native_site_ids(report: dict[str, object], *, site_budget: int) -> tuple[str, ...]:
    try:
        recovery = object_mapping(
            report.get("control_recovery"), name="control recovery"
        )
        delivery = object_mapping(
            recovery.get("budgeted_site_delivery"), name="budgeted site delivery"
        )
        rows = object_list(delivery.get("native_sites"), name="native sites")
    except TypeError as exc:
        raise ValidationError(
            "discovery qualification native-site evidence is invalid"
        ) from exc
    selected: list[tuple[int, str]] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="native site")
        except TypeError as exc:
            raise ValidationError("native site record is invalid") from exc
        rank = row.get("rank")
        site_id = row.get("site_id")
        if (
            not isinstance(rank, int)
            or isinstance(rank, bool)
            or rank < 1
            or not isinstance(site_id, str)
        ):
            raise ValidationError("native site rank or ID is invalid")
        if rank <= site_budget:
            selected.append((rank, safe_identifier(site_id, name="native site ID")))
    if not selected:
        raise ValidationError("frozen handoff budget contains no native-recovered site")
    return tuple(site_id for _, site_id in sorted(selected))


def validated_handoff_task_counts(
    manifest: dict[str, object],
) -> tuple[int, int, int]:
    """核对交接任务状态与清单汇总; 区分科学淘汰和技术无效。"""
    try:
        rows = object_list(manifest.get("tasks"), name="qualification handoff tasks")
    except TypeError as exc:
        raise ValidationError("qualification handoff tasks are invalid") from exc
    passed = failed = invalid = 0
    for raw in rows:
        try:
            row = object_mapping(raw, name="qualification handoff task")
        except TypeError as exc:
            raise ValidationError("qualification handoff task is invalid") from exc
        status = (row.get("execution_status"), row.get("reconstruction_status"))
        if status == ("completed", "passed"):
            passed += 1
        elif status == ("completed", "failed"):
            failed += 1
        elif status == ("invalid", "not_assessed"):
            invalid += 1
        else:
            raise ValidationError("qualification handoff task status is invalid")
    counts = (passed, failed, invalid)
    recorded = (
        manifest.get("passed_task_count"),
        manifest.get("failed_task_count"),
        manifest.get("invalid_task_count"),
    )
    if recorded != counts or sum(counts) != len(rows):
        raise ValidationError("qualification handoff task counts are inconsistent")
    return counts


def _handoff_qc_record(manifest: dict[str, object]) -> dict[str, int | float]:
    passed, failed, invalid = validated_handoff_task_counts(manifest)
    total = passed + failed + invalid
    return {
        "task_count": total,
        "passed_task_count": passed,
        "failed_task_count": failed,
        "invalid_task_count": invalid,
        "passed_fraction": round(passed / total, 6),
    }


def _validated_sources(
    *, config: AppConfig, handoff_run_dir: Path, control_run_dir: Path
) -> tuple[
    dict[str, object],
    dict[str, object],
    dict[str, object],
    Path,
    Path,
    tuple[str, ...],
    ChemistryDefinition,
]:
    handoff_manifest_path = handoff_run_dir / "qualification_handoff_manifest.json"
    handoff_manifest = read_document(
        handoff_manifest_path, name="qualification handoff manifest"
    )
    if (
        handoff_manifest.get("schema")
        != "vela.validation-qualification-handoff-manifest/5"
        or handoff_manifest.get("status") != "completed"
        or handoff_manifest.get("development_only") is not True
        or handoff_manifest.get("formal_qualification_gate") is not False
        or handoff_manifest.get("known_site_information_used_for_selection")
        is not False
    ):
        raise ValidationError("qualification handoff manifest is not eligible")
    _, _, invalid_count = validated_handoff_task_counts(handoff_manifest)
    if invalid_count:
        raise ValidationError(
            "qualification handoff contains technically invalid tasks"
        )
    handoff_plan_path, _ = validate_record(
        root=handoff_run_dir,
        raw=handoff_manifest.get("qualification_handoff_plan"),
        name="qualification handoff plan",
    )
    handoff_plan = read_document(handoff_plan_path, name="qualification handoff plan")
    try:
        selection = object_mapping(
            handoff_plan.get("selection"), name="handoff selection"
        )
        inputs = object_mapping(handoff_plan.get("inputs"), name="handoff inputs")
        qualification_run = object_mapping(
            inputs.get("qualification_run"), name="qualification run"
        )
    except TypeError as exc:
        raise ValidationError(
            "qualification handoff plan structure is invalid"
        ) from exc
    site_budget = selection.get("site_budget")
    qualification_relative = qualification_run.get("path")
    if (
        not isinstance(site_budget, int)
        or isinstance(site_budget, bool)
        or site_budget < 1
        or not isinstance(qualification_relative, str)
    ):
        raise ValidationError("qualification handoff budget or source is invalid")
    qualification_dir = (config.paths.outputs_dir / qualification_relative).resolve()
    allowed_qualification_root = (
        config.paths.outputs_dir / "discovery" / "qualifications"
    ).resolve()
    if not qualification_dir.is_relative_to(allowed_qualification_root):
        raise ValidationError("discovery qualification source is outside its run root")
    try:
        qualification_files = object_list(
            qualification_run.get("files"), name="qualification source files"
        )
    except TypeError as exc:
        raise ValidationError("qualification source file records are invalid") from exc
    recorded_files: dict[str, str] = {}
    for raw in qualification_files:
        try:
            record = object_mapping(raw, name="qualification source file")
        except TypeError as exc:
            raise ValidationError(
                "qualification source file record is invalid"
            ) from exc
        relative = record.get("path")
        expected_hash = record.get("sha256")
        if not isinstance(relative, str) or not isinstance(expected_hash, str):
            raise ValidationError("qualification source file record is invalid")
        path = (qualification_dir / relative).resolve()
        if not path.is_relative_to(qualification_dir) or not path.is_file():
            raise ValidationError(
                "qualification source file is missing or outside its run"
            )
        if sha256_file(path) != expected_hash:
            raise ValidationError("qualification source file hash mismatch")
        recorded_files[relative] = expected_hash
    discovery_report_path = qualification_dir / "qualification_report.json"
    if discovery_report_path.name not in recorded_files:
        raise ValidationError("qualification report is absent from handoff provenance")
    discovery_report = read_document(
        discovery_report_path, name="discovery qualification report"
    )
    if (
        discovery_report.get("schema") != "vela.discovery-qualification-report/8"
        or discovery_report.get("status") != "unqualified"
    ):
        raise ValidationError(
            "development diagnostic requires the frozen failed holdout"
        )
    native_sites = _native_site_ids(discovery_report, site_budget=site_budget)

    control_report_path = control_run_dir / "qualification_report.json"
    control_report = read_document(
        control_report_path, name="method qualification report"
    )
    if (
        control_report.get("schema") != "vela.validation-qualification-report/3"
        or control_report.get("status") != "qualified"
        or control_report.get("method_id") != config.validation.method_id
    ):
        raise ValidationError("chemistry-aware FlexPepDock method is not qualified")
    _, chemistry = _control_definition(config)
    if handoff_plan.get("chemistry") != chemistry_identity_record(chemistry):
        raise ValidationError("handoff chemistry differs from the local control")
    native_reference = (
        config.paths.data_dir
        / "validation"
        / "bound_states"
        / config.validation.local_controls[0].bound_state_id
        / "pair_reference.cif"
    )
    if not native_reference.is_file():
        raise ValidationError("qualification native reference is missing")
    return (
        handoff_manifest,
        handoff_plan,
        discovery_report,
        discovery_report_path,
        native_reference,
        native_sites,
        chemistry,
    )


def _starts(
    *,
    config: AppConfig,
    handoff_run_dir: Path,
    handoff_manifest: dict[str, object],
    native_sites: tuple[str, ...],
    chemistry: ChemistryDefinition,
    requested_start_ids: tuple[str, ...],
    receptor_backbone_mode: ReceptorBackboneMode,
) -> tuple[DiagnosticStart, ...]:
    requested_start_ids = tuple(
        safe_identifier(value, name="diagnostic start ID")
        for value in requested_start_ids
    )
    if not requested_start_ids or len(requested_start_ids) != len(
        set(requested_start_ids)
    ):
        raise ValidationError("diagnostic start IDs must be non-empty and unique")
    try:
        rows = object_list(handoff_manifest.get("tasks"), name="handoff tasks")
    except TypeError as exc:
        raise ValidationError("qualification handoff tasks are invalid") from exc
    requested = frozenset(requested_start_ids)
    selected: dict[str, DiagnosticStart] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="handoff task")
        except TypeError as exc:
            raise ValidationError("qualification handoff task is invalid") from exc
        task_id = row.get("task_id")
        if task_id not in requested:
            continue
        start_id = safe_identifier(task_id, name="diagnostic start ID")
        site_id = row.get("receptor_site_id")
        if site_id not in native_sites:
            raise ValidationError(
                f"diagnostic start is outside the frozen native sites: {start_id}"
            )
        if (
            row.get("execution_status") != "completed"
            or row.get("reconstruction_status") != "passed"
        ):
            raise ValidationError(
                f"selected qualification handoff task did not pass: {start_id}"
            )
        input_path, input_hash = validate_record(
            root=handoff_run_dir,
            raw=row.get("flexpepdock_input"),
            name=f"{start_id} FlexPepDock input",
        )
        receptor_count, histidines = validate_flexpepdock_input(
            path=input_path,
            chemistry=chemistry,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        )
        if receptor_backbone_mode == "local_constrained":
            backbone_selection = select_local_receptor_backbone(
                path=input_path,
                contact_A=(config.validation.refinement.receptor_backbone_contact_A),
                sequence_padding=(
                    config.validation.refinement.receptor_backbone_sequence_padding
                ),
            )
            direct_contacts = backbone_selection.direct_contact_pose_indices
            flexible_receptor = backbone_selection.flexible_pose_indices
        else:
            direct_contacts = ()
            flexible_receptor = ()
        seed = row.get("source_seed")
        if not isinstance(seed, int) or isinstance(seed, bool) or seed < 0:
            raise ValidationError("diagnostic source seed is invalid")
        selected[start_id] = DiagnosticStart(
            start_id=start_id,
            receptor_site_id=safe_identifier(site_id, name="receptor site ID"),
            receptor_id=safe_identifier(row.get("receptor_id"), name="receptor ID"),
            target=safe_identifier(row.get("target"), name="target ID"),
            source_seed=seed,
            input_path=input_path,
            input_sha256=input_hash,
            receptor_residue_count=receptor_count,
            fixed_histidine_pose_indices=histidines,
            direct_contact_receptor_pose_indices=direct_contacts,
            flexible_receptor_pose_indices=flexible_receptor,
        )
    missing = tuple(value for value in requested_start_ids if value not in selected)
    if missing:
        raise ValidationError("unknown diagnostic start IDs: " + ", ".join(missing))
    return tuple(selected[value] for value in requested_start_ids)


def _backbone_selection_records(
    starts: tuple[DiagnosticStart, ...],
) -> list[dict[str, JsonValue]]:
    """冻结每个起点独立、且不读取 native 的局部受体自由度。"""
    return [
        {
            "start_id": start.start_id,
            "direct_contact_receptor_pose_indices": list(
                start.direct_contact_receptor_pose_indices
            ),
            "flexible_receptor_pose_indices": list(
                start.flexible_receptor_pose_indices
            ),
        }
        for start in starts
    ]


def build_diagnostic_tasks(
    starts: tuple[DiagnosticStart, ...], seeds: tuple[int, ...]
) -> tuple[DiagnosticTask, ...]:
    tasks: list[DiagnosticTask] = []
    for start_index, start in enumerate(starts, 1):
        for seed_index, seed in enumerate(seeds, 1):
            tasks.append(
                DiagnosticTask(
                    task_id=f"diagnostic_{start_index:02d}_{seed_index:02d}",
                    start=start,
                    refinement_seed=seed,
                )
            )
    if not starts or not seeds:
        raise ValidationError("qualification diagnostic requires starts and seeds")
    return tuple(tasks)


def _task_records(tasks: tuple[DiagnosticTask, ...]) -> list[dict[str, JsonValue]]:
    """生成计划写入和执行复核共用的诊断任务合同。"""
    return [
        {
            "task_id": task.task_id,
            "start_id": task.start.start_id,
            "receptor_site_id": task.start.receptor_site_id,
            "receptor_id": task.start.receptor_id,
            "target": task.start.target,
            "source_seed": task.start.source_seed,
            "refinement_seed": task.refinement_seed,
            "input_sha256": task.start.input_sha256,
            "status": "planned",
        }
        for task in tasks
    ]


def write_qualification_refinement_plan(
    *,
    config: AppConfig,
    handoff_run_dir: Path,
    control_run_dir: Path,
    run_id: str,
    start_ids: tuple[str, ...],
    receptor_backbone_mode: str,
) -> DiagnosticPlan:
    """冻结只用于定位阶段二/三边界瓶颈的 native-aware 诊断。"""
    try:
        validate_run_id(run_id)
    except VelaError as exc:
        raise ValidationError(str(exc)) from exc
    resolved_backbone_mode = resolve_receptor_backbone_mode(receptor_backbone_mode)
    (
        handoff_manifest,
        _,
        _,
        discovery_report_path,
        native_reference,
        native_sites,
        chemistry,
    ) = _validated_sources(
        config=config,
        handoff_run_dir=handoff_run_dir,
        control_run_dir=control_run_dir,
    )
    starts = _starts(
        config=config,
        handoff_run_dir=handoff_run_dir,
        handoff_manifest=handoff_manifest,
        native_sites=native_sites,
        chemistry=chemistry,
        requested_start_ids=start_ids,
        receptor_backbone_mode=resolved_backbone_mode,
    )
    tasks = build_diagnostic_tasks(starts, config.validation.seeds)
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    control, _ = _control_definition(config)
    run_dir = (
        config.paths.outputs_dir / "validation" / "qualification_refinements" / run_id
    )
    if run_dir.exists():
        raise ValidationError(f"qualification refinement run already exists: {run_dir}")
    snapshot = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot, config.source_snapshot_text)
    handoff_manifest_path = handoff_run_dir / "qualification_handoff_manifest.json"
    control_report_path = control_run_dir / "qualification_report.json"
    atomic_write_json(
        run_dir / PLAN_NAME,
        {
            "schema": PLAN_SCHEMA,
            "stage": "validation_qualification_refinement_diagnostic",
            "status": "planned",
            "run_id": run_id,
            "planned_at": utc_now(),
            "development_only": True,
            "formal_qualification_gate": False,
            "native_information_used_for_task_selection": True,
            "evidence_category": EVIDENCE_CATEGORY,
            "receptor_relation": {
                "kind": "cross_receptor_cross_docking",
                "source_receptor_ids": sorted({start.receptor_id for start in starts}),
                "native_bound_state_id": control.bound_state_id,
            },
            "method_id": config.validation.method_id,
            "chemistry": chemistry_identity_record(chemistry),
            "software": {
                **vela_software_identity(),
                "rosetta_version": tool.version,
                "rosetta_scripts_sha256": tool.executable_sha256,
            },
            "inputs": {
                "config_snapshot": file_record(snapshot, root=run_dir),
                "handoff_manifest": _output_record(
                    handoff_manifest_path, config=config, name="handoff manifest"
                ),
                "method_qualification_report": _output_record(
                    control_report_path,
                    config=config,
                    name="method qualification report",
                ),
                "discovery_qualification_report": _output_record(
                    discovery_report_path,
                    config=config,
                    name="discovery qualification report",
                ),
                "native_reference": {
                    "path": native_reference.relative_to(
                        config.paths.data_dir
                    ).as_posix(),
                    "sha256": sha256_file(native_reference),
                },
            },
            "selection": {
                "selected_native_site_ids": list(native_sites),
                "requested_start_ids": list(start_ids),
                "selected_start_count": len(starts),
                "source_start_selection_was_native_free": True,
                "diagnostic_subset_selection_used_native_information": True,
                "source_holdout_status": "unqualified",
                "receptor_backbone_selections": _backbone_selection_records(starts),
            },
            "source_handoff_qc": _handoff_qc_record(handoff_manifest),
            "parameters": {
                "seeds": list(config.validation.seeds),
                "decoys_per_seed": config.validation.rosetta.decoys_per_seed,
                "parallel_tasks": config.validation.rosetta.parallel_tasks,
                "random_translation_A": config.validation.refinement.random_translation_A,
                "random_rotation_degrees": config.validation.refinement.random_rotation_degrees,
                "ranking_score": config.validation.refinement.ranking_score,
                "max_native_backbone_rmsd_A": control.max_recovery_rmsd_A,
                "max_cluster_backbone_rmsd_A": control.max_cluster_backbone_rmsd_A,
                "min_refinement_seed_support": 2,
                "min_source_start_support": min(2, len(starts)),
                "receptor_backbone": {
                    "mode": resolved_backbone_mode,
                    "selection": (
                        "start_pose_heavy_atom_contact_shell_with_sequence_padding"
                    ),
                    "contact_A": (
                        config.validation.refinement.receptor_backbone_contact_A
                    ),
                    "sequence_padding": (
                        config.validation.refinement.receptor_backbone_sequence_padding
                    ),
                    "movemap_policy": (
                        "all_receptor_chi; selected_receptor_bb_chi; "
                        "all_peptide_bb_chi; docking_jump_only"
                    ),
                    "coordinate_constraint": (
                        "FlexPepDock min_receptor_bb receptor-CA harmonic; "
                        "mean=0_A; sd=1_A"
                    ),
                },
            },
            "tasks": _task_records(tasks),
        },
    )
    return DiagnosticPlan(run_dir, tasks)


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], tuple[DiagnosticTask, ...], ChemistryDefinition, Path]:
    plan = read_document(run_dir / PLAN_NAME, name="qualification refinement plan")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or plan.get("stage") != "validation_qualification_refinement_diagnostic"
        or plan.get("status") != "planned"
        or plan.get("development_only") is not True
        or plan.get("formal_qualification_gate") is not False
        or plan.get("native_information_used_for_task_selection") is not True
        or not is_current_vela_software(plan.get("software"))
        or plan.get("method_id") != config.validation.method_id
    ):
        raise ValidationError("qualification refinement plan identity is invalid")
    try:
        inputs = object_mapping(plan.get("inputs"), name="diagnostic inputs")
        parameters = object_mapping(
            plan.get("parameters"), name="diagnostic parameters"
        )
        backbone = object_mapping(
            parameters.get("receptor_backbone"),
            name="diagnostic receptor backbone protocol",
        )
    except TypeError as exc:
        raise ValidationError("qualification refinement plan is invalid") from exc
    receptor_backbone_mode = resolve_receptor_backbone_mode(backbone.get("mode"))
    if (
        backbone.get("contact_A")
        != config.validation.refinement.receptor_backbone_contact_A
        or backbone.get("sequence_padding")
        != config.validation.refinement.receptor_backbone_sequence_padding
    ):
        raise ValidationError(
            "diagnostic receptor backbone parameters differ from current config"
        )
    snapshot, _ = validate_record(
        root=run_dir, raw=inputs.get("config_snapshot"), name="config snapshot"
    )
    if sha256_file(snapshot) != config.source_snapshot_sha256:
        raise ValidationError("current config differs from the diagnostic plan")
    handoff_root = config.paths.outputs_dir / "validation" / "qualification_handoffs"
    handoff_manifest = _resolve_output_record(
        inputs.get("handoff_manifest"),
        config=config,
        root=handoff_root,
        name="handoff manifest",
    )
    control_root = config.paths.outputs_dir / "validation" / "controls"
    control_report = _resolve_output_record(
        inputs.get("method_qualification_report"),
        config=config,
        root=control_root,
        name="method qualification report",
    )
    (
        source_manifest,
        _,
        _,
        discovery_report_path,
        native_reference,
        native_sites,
        chemistry,
    ) = _validated_sources(
        config=config,
        handoff_run_dir=handoff_manifest.parent,
        control_run_dir=control_report.parent,
    )
    discovery_root = config.paths.outputs_dir / "discovery" / "qualifications"
    frozen_discovery_report = _resolve_output_record(
        inputs.get("discovery_qualification_report"),
        config=config,
        root=discovery_root,
        name="discovery qualification report",
    )
    if frozen_discovery_report != discovery_report_path:
        raise ValidationError(
            "diagnostic discovery report differs from the frozen plan"
        )
    try:
        native_record = object_mapping(
            inputs.get("native_reference"), name="native reference"
        )
    except TypeError as exc:
        raise ValidationError("diagnostic native reference is invalid") from exc
    native_relative = native_record.get("path")
    native_hash = native_record.get("sha256")
    if not isinstance(native_relative, str) or not isinstance(native_hash, str):
        raise ValidationError("diagnostic native reference is invalid")
    frozen_native = (config.paths.data_dir / native_relative).resolve()
    if (
        not frozen_native.is_relative_to(config.paths.data_dir.resolve())
        or frozen_native != native_reference.resolve()
        or not frozen_native.is_file()
        or sha256_file(frozen_native) != native_hash
    ):
        raise ValidationError(
            "diagnostic native reference differs from the frozen plan"
        )
    if plan.get("chemistry") != chemistry_identity_record(chemistry):
        raise ValidationError("diagnostic chemistry differs from the frozen plan")
    if plan.get("source_handoff_qc") != _handoff_qc_record(source_manifest):
        raise ValidationError("source handoff QC differs from the frozen plan")
    try:
        selection = object_mapping(plan.get("selection"), name="diagnostic selection")
        requested_rows = object_list(
            selection.get("requested_start_ids"), name="diagnostic start IDs"
        )
    except TypeError as exc:
        raise ValidationError("diagnostic selection is invalid") from exc
    if any(not isinstance(value, str) for value in requested_rows):
        raise ValidationError("diagnostic start IDs are invalid")
    requested_start_ids = tuple(
        value for value in requested_rows if isinstance(value, str)
    )
    starts = _starts(
        config=config,
        handoff_run_dir=handoff_manifest.parent,
        handoff_manifest=source_manifest,
        native_sites=native_sites,
        chemistry=chemistry,
        requested_start_ids=requested_start_ids,
        receptor_backbone_mode=receptor_backbone_mode,
    )
    if selection.get("selected_start_count") != len(starts):
        raise ValidationError("diagnostic selected-start count changed")
    if selection.get("receptor_backbone_selections") != _backbone_selection_records(
        starts
    ):
        raise ValidationError("diagnostic receptor backbone selection changed")
    tasks = build_diagnostic_tasks(starts, config.validation.seeds)
    try:
        rows = object_list(plan.get("tasks"), name="diagnostic tasks")
    except TypeError as exc:
        raise ValidationError("diagnostic tasks are invalid") from exc
    recorded = [object_mapping(raw, name="diagnostic task") for raw in rows]
    if recorded != _task_records(tasks):
        raise ValidationError("diagnostic tasks differ from the frozen plan")
    tool = verify_rosetta_scripts_tool(config.validation.rosetta)
    try:
        software = object_mapping(plan.get("software"), name="diagnostic software")
    except TypeError as exc:
        raise ValidationError("diagnostic software identity is invalid") from exc
    if (
        software.get("rosetta_version") != tool.version
        or software.get("rosetta_scripts_sha256") != tool.executable_sha256
    ):
        raise ValidationError("current RosettaScripts differs from the diagnostic plan")
    return plan, tasks, chemistry, native_reference


def _prepare_start(
    *,
    config: AppConfig,
    chemistry: ChemistryDefinition,
    start: DiagnosticStart,
    run_dir: Path,
    plan_hash: str,
) -> PreparedDiagnosticStart:
    start_dir = run_dir / "starts" / start.start_id
    local_backbone = bool(start.flexible_receptor_pose_indices)
    records = {
        "output": "prepack output",
        "fix_disulfide": "prepack disulfide",
        "protocol": "prepack protocol",
        "scorefile": "prepack scorefile",
        "log": "prepack log",
    }
    if local_backbone:
        records["movemap"] = "local receptor MoveMap"
    resumed = resume_completed_result(
        directory=start_dir,
        filename="prepack_result.json",
        document_name="diagnostic prepack result",
        schema=PREPACK_SCHEMA,
        identity={"start_id": start.start_id},
        plan_hash_key="diagnostic_plan_sha256",
        plan_hash=plan_hash,
        records=records,
        stale_label="diagnostic prepack",
    )
    if resumed is not None:
        return PreparedDiagnosticStart(
            resumed.files["output"],
            resumed.files["fix_disulfide"],
            resumed.files.get("movemap"),
        )
    if start_dir.exists():
        raise ValidationError(
            f"incomplete diagnostic prepack requires review: {start_dir}"
        )
    start_dir.mkdir(parents=True)
    disulfide = start_dir / "fix_disulfide.txt"
    protocol = start_dir / "prepack.xml"
    movemap = start_dir / "local_receptor.movemap" if local_backbone else None
    write_disulfide_indices(
        destination=disulfide,
        receptor_residue_count=start.receptor_residue_count,
        chemistry=chemistry,
    )
    write_chemistry_prepack_protocol(
        destination=protocol,
        receptor_residue_count=start.receptor_residue_count,
        chemistry=chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    if movemap is not None:
        write_local_receptor_movemap(
            destination=movemap,
            receptor_residue_count=start.receptor_residue_count,
            peptide_residue_count=len(chemistry.sequence),
            flexible_receptor_pose_indices=(start.flexible_receptor_pose_indices),
        )
    command = build_chemistry_flexpepdock_command(
        settings=config.validation.rosetta,
        input_path=start.input_path,
        protocol_path=protocol,
        disulfide_path=disulfide,
        output_dir=start_dir,
        seed=config.validation.refinement.prepack_seed,
        fixed_histidine_pose_indices=start.fixed_histidine_pose_indices,
        nstruct=1,
        scorefile_name="prepack.sc",
        native_path=None,
        movemap_path=None,
    )
    log = start_dir / "prepack.log"
    run_rosetta_command(
        command=command,
        log_path=log,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        thread_count=1,
    )
    output = single_rosetta_pdb_output(start_dir)
    scorefile = start_dir / "prepack.sc"
    read_rosetta_scorefile(scorefile)
    validate_flexpepdock_input(
        path=output,
        chemistry=chemistry,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    result: dict[str, JsonValue] = {
        "schema": PREPACK_SCHEMA,
        "status": "completed",
        "start_id": start.start_id,
        "diagnostic_plan_sha256": plan_hash,
        "command": list(command),
        "output": file_record(output, root=start_dir),
        "fix_disulfide": file_record(disulfide, root=start_dir),
        "protocol": file_record(protocol, root=start_dir),
        "scorefile": file_record(scorefile, root=start_dir),
        "log": file_record(log, root=start_dir),
    }
    if movemap is not None:
        result["movemap"] = file_record(movemap, root=start_dir)
    atomic_write_json(
        start_dir / "prepack_result.json",
        result,
    )
    return PreparedDiagnosticStart(output, disulfide, movemap)


def _run_task(
    *,
    config: AppConfig,
    chemistry: ChemistryDefinition,
    native_reference: Path,
    task: DiagnosticTask,
    prepared: PreparedDiagnosticStart,
    run_dir: Path,
    plan_hash: str,
) -> Path:
    task_dir = run_dir / "tasks" / task.task_id
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="diagnostic task result",
        schema=TASK_SCHEMA,
        identity={"task_id": task.task_id},
        plan_hash_key="diagnostic_plan_sha256",
        plan_hash=plan_hash,
        records={
            "scorefile": "task scorefile",
            "decoy_manifest": "task decoy manifest",
            "protocol": "refinement protocol",
            "log": "task log",
        },
        stale_label="diagnostic task",
    )
    if resumed is not None:
        return resumed.path
    if task_dir.exists():
        raise ValidationError(f"incomplete diagnostic task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    protocol = task_dir / "refine.xml"
    write_chemistry_production_refine_protocol(
        destination=protocol,
        receptor_residue_count=task.start.receptor_residue_count,
        chemistry=chemistry,
        score_function=config.validation.rosetta.score_function,
        random_translation_A=config.validation.refinement.random_translation_A,
        random_rotation_degrees=config.validation.refinement.random_rotation_degrees,
        lowres_preoptimize=config.validation.rosetta.lowres_preoptimize,
        min_receptor_backbone=prepared.movemap_path is not None,
    )
    command = build_chemistry_flexpepdock_command(
        settings=config.validation.rosetta,
        input_path=prepared.prepacked_path,
        protocol_path=protocol,
        disulfide_path=prepared.disulfide_path,
        output_dir=task_dir,
        seed=task.refinement_seed,
        fixed_histidine_pose_indices=task.start.fixed_histidine_pose_indices,
        nstruct=config.validation.rosetta.decoys_per_seed,
        scorefile_name="refine.sc",
        native_path=native_reference,
        movemap_path=prepared.movemap_path,
    )
    log = task_dir / "refine.log"
    started_at = utc_now()
    run_rosetta_command(
        command=command,
        log_path=log,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
    )
    scorefile = task_dir / "refine.sc"
    rows = read_rosetta_scorefile(scorefile)
    if len(rows) != config.validation.rosetta.decoys_per_seed:
        raise ValidationError(
            f"{task.task_id} produced {len(rows)} decoys; expected "
            f"{config.validation.rosetta.decoys_per_seed}"
        )
    decoy_paths = index_rosetta_pdb_outputs(task_dir)
    chemistry_failures: list[dict[str, str]] = []
    for description, path in decoy_paths.items():
        try:
            validate_flexpepdock_input(
                path=path,
                chemistry=chemistry,
                min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
                max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
            )
        except ValidationError as exc:
            chemistry_failures.append({"description": description, "failure": str(exc)})
    decoy_manifest = write_refinement_decoy_manifest(
        rows=rows,
        ranking_score=config.validation.refinement.ranking_score,
        decoy_paths=decoy_paths,
        task_dir=task_dir,
    )
    result = task_dir / "task_result.json"
    atomic_write_json(
        result,
        {
            "schema": TASK_SCHEMA,
            "status": "completed",
            "task_id": task.task_id,
            "start_id": task.start.start_id,
            "receptor_site_id": task.start.receptor_site_id,
            "source_seed": task.start.source_seed,
            "refinement_seed": task.refinement_seed,
            "receptor_backbone_mode": (
                "local_constrained" if prepared.movemap_path is not None else "fixed"
            ),
            "diagnostic_plan_sha256": plan_hash,
            "started_at": started_at,
            "completed_at": utc_now(),
            "chemistry_qc": {
                "status": "passed" if not chemistry_failures else "failed",
                "passed_decoy_count": len(decoy_paths) - len(chemistry_failures),
                "failed_decoy_count": len(chemistry_failures),
                "failures": chemistry_failures,
            },
            "command": list(command),
            "scorefile": file_record(scorefile, root=task_dir),
            "decoy_manifest": file_record(decoy_manifest, root=task_dir),
            "protocol": file_record(protocol, root=task_dir),
            "log": file_record(log, root=task_dir),
        },
    )
    return result


def run_qualification_refinement(
    *, config: AppConfig, run_dir: Path
) -> DiagnosticOutcome:
    """执行冻结的显式起点开发诊断; 不将结果晋升为方法资格。"""
    manifest_path = run_dir / MANIFEST_NAME
    if manifest_path.exists():
        raise ValidationError(f"diagnostic manifest already exists: {manifest_path}")
    plan, tasks, chemistry, native_reference = _verify_plan(
        config=config, run_dir=run_dir
    )
    try:
        parameters = object_mapping(
            plan.get("parameters"), name="diagnostic parameters"
        )
        backbone = object_mapping(
            parameters.get("receptor_backbone"),
            name="diagnostic receptor backbone protocol",
        )
    except TypeError as exc:
        raise ValidationError(
            "diagnostic receptor backbone protocol is invalid"
        ) from exc
    receptor_backbone_mode = resolve_receptor_backbone_mode(backbone.get("mode"))
    plan_path = run_dir / PLAN_NAME
    plan_hash = sha256_file(plan_path)
    starts = {task.start.start_id: task.start for task in tasks}
    prepared = {
        start_id: _prepare_start(
            config=config,
            chemistry=chemistry,
            start=start,
            run_dir=run_dir,
            plan_hash=plan_hash,
        )
        for start_id, start in starts.items()
    }

    def execute(task: DiagnosticTask) -> Path:
        return _run_task(
            config=config,
            chemistry=chemistry,
            native_reference=native_reference,
            task=task,
            prepared=prepared[task.start.start_id],
            run_dir=run_dir,
            plan_hash=plan_hash,
        )

    with ThreadPoolExecutor(
        max_workers=config.validation.rosetta.parallel_tasks
    ) as executor:
        results = tuple(executor.map(execute, tasks))
    atomic_write_json(
        manifest_path,
        {
            "schema": MANIFEST_SCHEMA,
            "stage": "validation_qualification_refinement_diagnostic",
            "status": "completed",
            "completed_at": utc_now(),
            "development_only": True,
            "formal_qualification_gate": False,
            "native_information_used_for_task_selection": True,
            "evidence_category": EVIDENCE_CATEGORY,
            "method_id": config.validation.method_id,
            "receptor_backbone_mode": receptor_backbone_mode,
            "chemistry": chemistry_identity_record(chemistry),
            "qualification_refinement_plan": file_record(plan_path, root=run_dir),
            "task_count": len(tasks),
            "decoy_count": len(tasks) * config.validation.rosetta.decoys_per_seed,
            "tasks": [
                {
                    "task_id": task.task_id,
                    "start_id": task.start.start_id,
                    "receptor_site_id": task.start.receptor_site_id,
                    "source_seed": task.start.source_seed,
                    "refinement_seed": task.refinement_seed,
                    "task_result": file_record(result, root=run_dir),
                }
                for task, result in zip(tasks, results, strict=True)
            ],
            "interpretation": (
                "development sampling only; native-aware task selection prevents "
                "formal qualification or biological promotion"
            ),
        },
    )
    return DiagnosticOutcome(manifest_path, len(tasks))
