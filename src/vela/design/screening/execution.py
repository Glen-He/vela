"""阶段四固定骨架初筛的命令、并行执行和恢复。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import atomic_write_json, sha256_file, utc_now
from vela.design.models import DesignError, ScreenTask
from vela.design.scores import SCREEN_SCORE_COLUMNS, screen_metrics
from vela.design.screening.records import SCREEN_PLAN_NAME, ScreenPlan, read_screen_plan
from vela.preparation.chemistry import ChemistryDefinition, HistidineState
from vela.validation.records import file_record, resume_completed_result
from vela.validation.refinement.reconstruction import (
    validate_flexpepdock_input,
    write_disulfide_indices,
)
from vela.validation.rosetta import (
    rosetta_crash_log_dir,
    run_rosetta_command,
    single_rosetta_pdb_output,
)
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class ScreenOutcome:
    """全部成对筛查任务完成后的结果索引。"""

    manifest_path: Path
    task_count: int


def candidate_chemistry(*, config: AppConfig, task: ScreenTask) -> ChemistryDefinition:
    """为 WT 或候选序列建立端基和二硫键不变的显式化学定义。"""
    if task.state == "wt":
        return config.chemistry
    sequence = task.candidate.sequence
    histidines = tuple(
        HistidineState(position, config.design.sequence.candidate_histidine_state)
        for position, residue in enumerate(sequence, 1)
        if residue == "H"
    )
    return ChemistryDefinition(
        ligand_id=config.chemistry.ligand_id,
        chemistry_id=(
            f"{config.chemistry.chemistry_id}__{task.candidate.candidate_id}"
        ),
        sequence=sequence,
        chirality=config.chemistry.chirality,
        disulfide_bonds=config.chemistry.disulfide_bonds,
        n_terminus=config.chemistry.n_terminus,
        c_terminus=config.chemistry.c_terminus,
        target_ph=config.chemistry.target_ph,
        net_charge=task.candidate.net_charge,
        histidines=histidines,
        other_modifications_status=config.chemistry.other_modifications_status,
        other_modifications=config.chemistry.other_modifications,
        decision_sources=config.chemistry.decision_sources,
    )


def task_histidine_pose_indices(
    *, config: AppConfig, task: ScreenTask
) -> tuple[int, ...]:
    """返回当前 WT 或候选结构中必须保持 tautomer 的全部 His pose 编号。"""
    sequence = (
        task.candidate.sequence if task.state == "mutant" else config.chemistry.sequence
    )
    receptor_indices = tuple(
        index
        for index in task.template.fixed_histidine_pose_indices
        if index <= task.template.receptor_residue_count
    )
    peptide_indices = tuple(
        task.template.receptor_residue_count + position
        for position, residue in enumerate(sequence, 1)
        if residue == "H"
    )
    return tuple(sorted((*receptor_indices, *peptide_indices)))


def build_screen_command(
    *,
    config: AppConfig,
    plan: ScreenPlan,
    task: ScreenTask,
    output_dir: Path,
    disulfide_path: Path,
) -> tuple[str, ...]:
    """构造成对局部重排和 InterfaceAnalyzer 的单结构命令。"""
    settings = config.validation.rosetta
    command = [
        str(settings.scripts_executable),
        "-database",
        str(settings.database),
        "-s",
        str(task.template.path),
        "-parser:protocol",
        str(plan.protocol_path),
        "-parser:script_vars",
        f"score_function={config.design.screen.score_function}",
        f"resfile_path={task.resfile_path}",
        ("pack_separated=" + str(config.design.screen.pack_separated).lower()),
        "-in:fix_disulf",
        str(disulfide_path),
        "-constant_seed",
        "-jran",
        str(task.seed),
        "-nstruct",
        "1",
        "-ex1",
        "-ex2aro",
        "-use_input_sc",
        "-packing:no_optH",
        "true",
    ]
    histidines = task_histidine_pose_indices(config=config, task=task)
    if histidines:
        command.extend(
            ["-packing:fix_his_tautomer", *(str(index) for index in histidines)]
        )
    command.extend(
        [
            "-out:path:pdb",
            str(output_dir),
            "-out:path:score",
            str(output_dir),
            "-out:file:scorefile",
            "screen.sc",
            "-overwrite",
        ]
    )
    return tuple(command)


def _run_task(
    *,
    config: AppConfig,
    plan: ScreenPlan,
    task: ScreenTask,
    plan_sha256: str,
) -> Path:
    task_dir = plan.run_dir / "tasks" / task.task_id
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="design screen task result",
        schema="vela.design-screen-task-result/2",
        identity={"task_id": task.task_id, "pair_id": task.pair_id},
        plan_hash_key="screen_plan_sha256",
        plan_hash=plan_sha256,
        records={
            "output": "screen task output",
            "scorefile": "screen task scorefile",
            "fix_disulfide": "screen task fix_disulfide",
            "log": "screen task log",
        },
        stale_label="design screen task",
        stale_error=DesignError,
    )
    if resumed is not None:
        return resumed.path
    if task_dir.exists():
        raise DesignError(f"incomplete design screen task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    if sha256_file(task.resfile_path) != task.resfile_sha256:
        raise DesignError(f"design resfile changed: {task.resfile_path}")
    disulfide_path = task_dir / "fix_disulfide.txt"
    chemistry = candidate_chemistry(config=config, task=task)
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=task.template.receptor_residue_count,
        chemistry=chemistry,
    )
    command = build_screen_command(
        config=config,
        plan=plan,
        task=task,
        output_dir=task_dir,
        disulfide_path=disulfide_path,
    )
    log_path = task_dir / "screen.log"
    started_at = utc_now()
    run_rosetta_command(
        command=command,
        log_path=log_path,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        thread_count=1,
    )
    score_path = task_dir / "screen.sc"
    rows = read_rosetta_scorefile(score_path)
    if len(rows) != 1:
        raise DesignError(f"design screen task must produce one score: {task.task_id}")
    metrics = screen_metrics(rows[0])
    if config.design.screen.ranking_score != "dG_separated":
        raise DesignError("design screen ranking_score must be dG_separated")
    output_path = single_rosetta_pdb_output(task_dir)
    receptor_count, _ = validate_flexpepdock_input(
        path=output_path,
        chemistry=chemistry,
        min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
        max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
    )
    if receptor_count != task.template.receptor_residue_count:
        raise DesignError("design screen changed the receptor residue count")
    result_path = task_dir / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": "vela.design-screen-task-result/2",
            "status": "completed",
            "task_id": task.task_id,
            "pair_id": task.pair_id,
            "wt_context_id": task.wt_context_id,
            "state": task.state,
            "candidate_id": task.candidate.candidate_id,
            "template_id": task.template.template_id,
            "evidence_role": task.template.evidence_role,
            "target": task.template.target,
            "seed": task.seed,
            "screen_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "score_units": "Rosetta score units (REU) unless stated otherwise",
            "score_columns": SCREEN_SCORE_COLUMNS,
            "interface_metrics": metrics.as_dict(),
            "command": list(command),
            "output": file_record(output_path, root=task_dir),
            "scorefile": file_record(score_path, root=task_dir),
            "fix_disulfide": file_record(disulfide_path, root=task_dir),
            "log": file_record(log_path, root=task_dir),
        },
    )
    return result_path


def run_screen(*, config: AppConfig, run_dir: Path) -> ScreenOutcome:
    """并行执行或恢复全部成对 WT/候选界面筛查任务。"""
    manifest_path = run_dir / "screen_manifest.json"
    if manifest_path.exists():
        raise DesignError(f"design screen manifest already exists: {manifest_path}")
    plan = read_screen_plan(config=config, run_dir=run_dir)
    plan_path = run_dir / SCREEN_PLAN_NAME
    plan_sha256 = sha256_file(plan_path)

    execution_tasks: list[ScreenTask] = []
    canonical_wt_tasks: dict[str, ScreenTask] = {}
    for task in plan.tasks:
        if task.wt_context_id is None:
            execution_tasks.append(task)
            continue
        if task.wt_context_id not in canonical_wt_tasks:
            canonical_wt_tasks[task.wt_context_id] = task
            execution_tasks.append(task)

    def execute(task: ScreenTask) -> tuple[str, Path]:
        result = _run_task(config=config, plan=plan, task=task, plan_sha256=plan_sha256)
        key = task.wt_context_id or task.task_id
        return key, result

    with ThreadPoolExecutor(
        max_workers=config.design.screen.parallel_tasks
    ) as executor:
        executed = dict(executor.map(execute, execution_tasks))
    results = tuple(executed[task.wt_context_id or task.task_id] for task in plan.tasks)
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.design-screen-manifest/2",
            "stage": "design_interface_screen",
            "status": "completed",
            "completed_at": utc_now(),
            "design_round": plan.design_round,
            "method_id": config.design.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "logical_pair_count": len(plan.tasks) // 2,
            "executed_task_count": len(execution_tasks),
            "shared_wt_context_count": len(canonical_wt_tasks),
            "screen_plan": file_record(plan_path, root=run_dir),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "pair_id": task.pair_id,
                    "state": task.state,
                    "candidate_id": task.candidate.candidate_id,
                    "template_id": task.template.template_id,
                    "seed": task.seed,
                    "wt_context_id": task.wt_context_id,
                    "task_result": file_record(result, root=run_dir),
                }
                for task, result in zip(plan.tasks, results, strict=True)
            ],
        },
    )
    return ScreenOutcome(manifest_path, len(execution_tasks))
