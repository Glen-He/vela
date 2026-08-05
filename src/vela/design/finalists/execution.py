"""阶段四候选 FlexPepDock 柔性复核执行与恢复。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import atomic_write_json, sha256_file, utc_now
from vela.design.finalists.records import (
    FINALIST_PLAN_NAME,
    FinalistPlan,
    read_finalist_plan,
)
from vela.design.models import DesignError, FinalistStart, FinalistTask
from vela.design.scores import FINALIST_SCORE_COLUMNS, write_finalist_decoy_manifest
from vela.preparation.chemistry import ChemistryDefinition, HistidineState
from vela.validation.records import file_record, resume_completed_result
from vela.validation.refinement.reconstruction import write_disulfide_indices
from vela.validation.rosetta import (
    build_refine_command,
    rosetta_crash_log_dir,
    run_rosetta_command,
)
from vela.validation.scores import (
    index_rosetta_pdb_outputs,
    read_rosetta_scorefile,
)


@dataclass(frozen=True, slots=True)
class FinalistOutcome:
    """全部候选柔性复核任务完成后的结果索引。"""

    manifest_path: Path
    task_count: int


def finalist_chemistry(
    *, config: AppConfig, start: FinalistStart
) -> ChemistryDefinition:
    """为一份 WT 或候选柔性起点建立显式配体化学定义。"""
    if start.state == "wt":
        return config.chemistry
    sequence = start.candidate.sequence
    histidines = tuple(
        HistidineState(position, config.design.sequence.candidate_histidine_state)
        for position, residue in enumerate(sequence, 1)
        if residue == "H"
    )
    return ChemistryDefinition(
        ligand_id=config.chemistry.ligand_id,
        chemistry_id=(
            f"{config.chemistry.chemistry_id}__{start.candidate.candidate_id}"
        ),
        sequence=sequence,
        chirality=config.chemistry.chirality,
        disulfide_bonds=config.chemistry.disulfide_bonds,
        n_terminus=config.chemistry.n_terminus,
        c_terminus=config.chemistry.c_terminus,
        target_ph=config.chemistry.target_ph,
        net_charge=start.candidate.net_charge,
        histidines=histidines,
        other_modifications_status=config.chemistry.other_modifications_status,
        other_modifications=config.chemistry.other_modifications,
        decision_sources=config.chemistry.decision_sources,
    )


def _run_task(
    *,
    config: AppConfig,
    plan: FinalistPlan,
    task: FinalistTask,
    plan_sha256: str,
) -> Path:
    task_dir = plan.run_dir / "tasks" / task.task_id
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="finalist task result",
        schema="vela.design-finalist-task-result/2",
        identity={"task_id": task.task_id, "pair_id": task.pair_id},
        plan_hash_key="finalist_plan_sha256",
        plan_hash=plan_sha256,
        records={
            "scorefile": "finalist scorefile",
            "decoy_manifest": "finalist decoy_manifest",
            "fix_disulfide": "finalist fix_disulfide",
            "log": "finalist log",
        },
        stale_label="finalist task",
        stale_error=DesignError,
    )
    if resumed is not None:
        return resumed.path
    if task_dir.exists():
        raise DesignError(f"incomplete finalist task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    chemistry = finalist_chemistry(config=config, start=task.start)
    disulfide_path = task_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=task.start.template.receptor_residue_count,
        chemistry=chemistry,
    )
    command = build_refine_command(
        settings=config.validation.rosetta,
        input_path=task.start.path,
        disulfide_path=disulfide_path,
        output_dir=task_dir,
        seed=task.seed,
        native_path=task.start.path,
        random_translation_A=config.validation.refinement.random_translation_A,
        random_rotation_degrees=config.validation.refinement.random_rotation_degrees,
        fixed_histidine_pose_indices=task.start.histidine_pose_indices,
    )
    log_path = task_dir / "refine.log"
    started_at = utc_now()
    run_rosetta_command(
        command=command,
        log_path=log_path,
        crash_dir=rosetta_crash_log_dir(outputs_dir=config.paths.outputs_dir),
        thread_count=1,
    )
    score_path = task_dir / "refine.sc"
    rows = read_rosetta_scorefile(score_path)
    expected = config.validation.rosetta.decoys_per_seed
    if len(rows) != expected:
        raise DesignError(
            f"{task.task_id} produced {len(rows)} decoys; expected {expected}"
        )
    decoy_manifest = write_finalist_decoy_manifest(
        rows=rows,
        decoy_paths=index_rosetta_pdb_outputs(task_dir),
        task_dir=task_dir,
    )
    result_path = task_dir / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": "vela.design-finalist-task-result/2",
            "status": "completed",
            "task_id": task.task_id,
            "pair_id": task.pair_id,
            "state": task.start.state,
            "candidate_id": task.start.candidate.candidate_id,
            "template_id": task.start.template.template_id,
            "target": task.start.template.target,
            "evidence_role": task.start.template.evidence_role,
            "seed": task.seed,
            "finalist_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "score_units": "Rosetta score units (REU) unless stated otherwise",
            "score_columns": FINALIST_SCORE_COLUMNS,
            "command": list(command),
            "scorefile": file_record(score_path, root=task_dir),
            "decoy_manifest": file_record(decoy_manifest, root=task_dir),
            "fix_disulfide": file_record(disulfide_path, root=task_dir),
            "log": file_record(log_path, root=task_dir),
        },
    )
    return result_path


def run_finalists(*, config: AppConfig, run_dir: Path) -> FinalistOutcome:
    """并行执行或恢复全部 WT/候选柔性复核任务。"""
    manifest_path = run_dir / "finalist_manifest.json"
    if manifest_path.exists():
        raise DesignError(f"finalist manifest already exists: {manifest_path}")
    plan = read_finalist_plan(config=config, run_dir=run_dir)
    plan_path = run_dir / FINALIST_PLAN_NAME
    plan_sha256 = sha256_file(plan_path)

    def execute(task: FinalistTask) -> Path:
        return _run_task(config=config, plan=plan, task=task, plan_sha256=plan_sha256)

    with ThreadPoolExecutor(
        max_workers=config.design.finalists.parallel_tasks
    ) as executor:
        results = tuple(executor.map(execute, plan.tasks))
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.design-finalist-manifest/1",
            "stage": "design_flexible_verification",
            "status": "completed",
            "completed_at": utc_now(),
            "method_id": config.design.method_id,
            "flexpepdock_method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "objective": config.design.objective,
            "finalist_plan": file_record(plan_path, root=run_dir),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "pair_id": task.pair_id,
                    "state": task.start.state,
                    "candidate_id": task.start.candidate.candidate_id,
                    "template_id": task.start.template.template_id,
                    "seed": task.seed,
                    "task_result": file_record(result, root=run_dir),
                }
                for task, result in zip(plan.tasks, results, strict=True)
            ],
        },
    )
    return FinalistOutcome(manifest_path, len(plan.tasks))
