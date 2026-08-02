"""阶段三已冻结局部精修任务的执行和逐任务恢复。"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import atomic_write_json, sha256_file, utc_now
from vela.validation.models import ValidationError
from vela.validation.records import file_record, resume_completed_result
from vela.validation.refinement.planning import (
    REFINEMENT_PLAN_NAME,
    RefinementStart,
    RefinementTask,
    refinement_identity,
    verify_refinement_plan,
)
from vela.validation.refinement.reconstruction import write_disulfide_indices
from vela.validation.rosetta import (
    build_prepack_command,
    build_refine_command,
    run_rosetta_command,
)
from vela.validation.scores import (
    index_rosetta_pdb_outputs,
    read_rosetta_scorefile,
    write_refinement_decoy_manifest,
)


@dataclass(frozen=True, slots=True)
class PreparedStart:
    """一个可供全部局部精修 seed 共用的 prepack 起点。"""

    start: RefinementStart
    prepacked_path: Path
    disulfide_path: Path


@dataclass(frozen=True, slots=True)
class RefinementOutcome:
    """局部精修运行完成后的结果索引。"""

    manifest_path: Path
    task_count: int


def _single_output(directory: Path, *, score_name: str) -> tuple[Path, Path]:
    score_path = directory / score_name
    rows = read_rosetta_scorefile(score_path)
    outputs = tuple(directory.glob("*.pdb"))
    if len(rows) != 1 or len(outputs) != 1:
        raise ValidationError(
            f"Rosetta prepack must produce one score and one PDB: {directory}"
        )
    return outputs[0], score_path


def _resume_prepack(
    *, start: RefinementStart, start_dir: Path, plan_hash: str
) -> PreparedStart | None:
    result = resume_completed_result(
        directory=start_dir,
        filename="prepack_result.json",
        document_name="refinement prepack result",
        schema="vela.validation-refinement-prepack-result/1",
        identity={"start_id": start.start_id},
        plan_hash_key="refinement_plan_sha256",
        plan_hash=plan_hash,
        records={
            "output": "prepack output",
            "fix_disulfide": "prepack disulfide",
            "scorefile": "prepack scorefile",
            "log": "prepack log",
        },
        stale_label="refinement prepack",
    )
    if result is None:
        return None
    return PreparedStart(
        start,
        result.files["output"],
        result.files["fix_disulfide"],
    )


def _prepare_start(
    *, config: AppConfig, start: RefinementStart, run_dir: Path, plan_hash: str
) -> PreparedStart:
    start_dir = run_dir / "starts" / start.start_id
    resumed = _resume_prepack(start=start, start_dir=start_dir, plan_hash=plan_hash)
    if resumed is not None:
        return resumed
    if start_dir.exists():
        raise ValidationError(f"incomplete prepack requires review: {start_dir}")
    start_dir.mkdir(parents=True)
    disulfide_path = start_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=start.receptor_residue_count,
        chemistry=config.chemistry,
    )
    command = build_prepack_command(
        settings=config.validation.rosetta,
        input_path=start.input_path,
        disulfide_path=disulfide_path,
        output_dir=start_dir,
        seed=config.validation.refinement.prepack_seed,
        fixed_histidine_pose_indices=start.fixed_histidine_pose_indices,
    )
    log_path = start_dir / "prepack.log"
    run_rosetta_command(command=command, log_path=log_path, thread_count=1)
    output, score_path = _single_output(start_dir, score_name="prepack.sc")
    atomic_write_json(
        start_dir / "prepack_result.json",
        {
            "schema": "vela.validation-refinement-prepack-result/1",
            "status": "completed",
            "start_id": start.start_id,
            "refinement_plan_sha256": plan_hash,
            "command": list(command),
            "output": file_record(output, root=start_dir),
            "scorefile": file_record(score_path, root=start_dir),
            "fix_disulfide": file_record(disulfide_path, root=start_dir),
            "log": file_record(log_path, root=start_dir),
        },
    )
    return PreparedStart(start, output, disulfide_path)


def _run_task(
    *,
    config: AppConfig,
    task: RefinementTask,
    prepared: PreparedStart,
    run_dir: Path,
    plan_hash: str,
) -> Path:
    task_dir = run_dir / "tasks" / task.task_id
    resumed = resume_completed_result(
        directory=task_dir,
        filename="task_result.json",
        document_name="refinement task result",
        schema="vela.validation-refinement-task-result/1",
        identity={"task_id": task.task_id},
        plan_hash_key="refinement_plan_sha256",
        plan_hash=plan_hash,
        records={
            "scorefile": "task scorefile",
            "decoy_manifest": "task decoy_manifest",
            "log": "task log",
        },
        stale_label="refinement task",
    )
    if resumed is not None:
        return resumed.path
    if task_dir.exists():
        raise ValidationError(f"incomplete refinement task requires review: {task_dir}")
    task_dir.mkdir(parents=True)
    command = build_refine_command(
        settings=config.validation.rosetta,
        input_path=prepared.prepacked_path,
        disulfide_path=prepared.disulfide_path,
        output_dir=task_dir,
        seed=task.seed,
        native_path=task.start.input_path,
        random_translation_A=config.validation.refinement.random_translation_A,
        random_rotation_degrees=config.validation.refinement.random_rotation_degrees,
        fixed_histidine_pose_indices=task.start.fixed_histidine_pose_indices,
    )
    log_path = task_dir / "refine.log"
    started_at = utc_now()
    run_rosetta_command(command=command, log_path=log_path)
    score_path = task_dir / "refine.sc"
    rows = read_rosetta_scorefile(score_path)
    if len(rows) != config.validation.rosetta.decoys_per_seed:
        raise ValidationError(
            f"{task.task_id} produced {len(rows)} decoys; expected "
            f"{config.validation.rosetta.decoys_per_seed}"
        )
    decoy_manifest = write_refinement_decoy_manifest(
        rows=rows,
        ranking_score=config.validation.refinement.ranking_score,
        decoy_paths=index_rosetta_pdb_outputs(task_dir),
        task_dir=task_dir,
    )
    result_path = task_dir / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": "vela.validation-refinement-task-result/1",
            "status": "completed",
            "task_id": task.task_id,
            "start_id": task.start.start_id,
            "candidate_id": task.start.candidate_id,
            "receptor_site_id": task.start.receptor_site_id,
            "pose_id": task.start.pose_id,
            "receptor_id": task.start.receptor_id,
            "target": task.start.target,
            "source_seed": task.start.source_seed,
            "refinement_seed": task.seed,
            "refinement_plan_sha256": plan_hash,
            "started_at": started_at,
            "completed_at": utc_now(),
            "command": list(command),
            "scorefile": file_record(score_path, root=task_dir),
            "decoy_manifest": file_record(decoy_manifest, root=task_dir),
            "log": file_record(log_path, root=task_dir),
        },
    )
    return result_path


def run_refinement(*, config: AppConfig, run_dir: Path) -> RefinementOutcome:
    """执行或恢复冻结的全部 candidate prepack 和局部精修任务。"""
    manifest_path = run_dir / "refinement_manifest.json"
    if manifest_path.exists():
        raise ValidationError(f"refinement manifest already exists: {manifest_path}")
    plan, tasks = verify_refinement_plan(config=config, run_dir=run_dir)
    evidence_category, known_site_information_used = refinement_identity(plan)
    plan_path = run_dir / REFINEMENT_PLAN_NAME
    plan_hash = sha256_file(plan_path)
    starts = {task.start.start_id: task.start for task in tasks}
    prepared = {
        start_id: _prepare_start(
            config=config, start=start, run_dir=run_dir, plan_hash=plan_hash
        )
        for start_id, start in starts.items()
    }

    def execute(task: RefinementTask) -> Path:
        return _run_task(
            config=config,
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
            "schema": "vela.validation-refinement-manifest/1",
            "stage": "validation_local_refinement",
            "status": "completed",
            "completed_at": utc_now(),
            "method_id": config.validation.method_id,
            "chemistry_id": config.chemistry.chemistry_id,
            "evidence_category": evidence_category,
            "known_site_information_used": known_site_information_used,
            "refinement_plan": file_record(plan_path, root=run_dir),
            "tasks": [
                {
                    "task_id": task.task_id,
                    "start_id": task.start.start_id,
                    "candidate_id": task.start.candidate_id,
                    "receptor_id": task.start.receptor_id,
                    "refinement_seed": task.seed,
                    "task_result": file_record(result, root=run_dir),
                }
                for task, result in zip(tasks, results, strict=True)
            ],
        },
    )
    return RefinementOutcome(manifest_path, len(tasks))
