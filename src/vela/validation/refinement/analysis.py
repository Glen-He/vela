"""阶段三局部精修运行的完整性检查、分析和报告。"""

from __future__ import annotations

import csv
import io
import math
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
from vela.validation.models import ValidationError
from vela.validation.records import read_document, validate_record
from vela.validation.refinement.geometry import (
    ComplexGeometry,
    RefinedDecoy,
    ResolvedAnalysisSettings,
    assess_refined_decoy,
    cluster_refined_decoys,
    read_complex_geometry,
    resolve_analysis_settings,
)
from vela.validation.refinement.planning import (
    RefinementTask,
    refinement_identity,
    verify_refinement_plan,
)
from vela.validation.scores import read_rosetta_scorefile


@dataclass(frozen=True, slots=True)
class RefinementAnalysisOutcome:
    """局部精修分析报告的路径和规模。"""

    manifest_path: Path
    decoy_count: int
    cluster_count: int


def _decoy_rows(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != [
            "description",
            "ranking_score",
            "path",
            "sha256",
        ]:
            raise ValidationError(f"invalid refinement decoy columns: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows:
        raise ValidationError(f"refinement decoy manifest contains no rows: {path}")
    return rows


def _manifest_tasks(
    *,
    run_dir: Path,
    tasks: tuple[RefinementTask, ...],
    plan: dict[str, object],
) -> tuple[tuple[RefinementTask, Path, Path], ...]:
    manifest = read_document(
        run_dir / "refinement_manifest.json", name="refinement manifest"
    )
    if (
        manifest.get("schema") != "vela.validation-refinement-manifest/2"
        or manifest.get("stage") != "validation_local_refinement"
        or manifest.get("status") != "completed"
        or manifest.get("evidence_category") != plan.get("evidence_category")
        or manifest.get("known_site_information_used")
        is not plan.get("known_site_information_used")
    ):
        raise ValidationError("refinement manifest identity is invalid")
    validate_record(
        root=run_dir,
        raw=manifest.get("refinement_plan"),
        name="refinement plan",
    )
    try:
        rows = object_list(manifest.get("tasks"), name="refinement manifest tasks")
    except TypeError as exc:
        raise ValidationError("refinement manifest tasks are invalid") from exc
    by_id = {task.task_id: task for task in tasks}
    results: list[tuple[RefinementTask, Path, Path]] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="refinement manifest task")
        except TypeError as exc:
            raise ValidationError("refinement manifest task is invalid") from exc
        task_id = row.get("task_id")
        task = by_id.get(task_id) if isinstance(task_id, str) else None
        if task is None:
            raise ValidationError("refinement manifest contains an unknown task")
        result_path, _ = validate_record(
            root=run_dir, raw=row.get("task_result"), name=f"{task_id} task result"
        )
        result = read_document(result_path, name="refinement task result")
        task_dir = result_path.parent
        decoy_manifest, _ = validate_record(
            root=task_dir,
            raw=result.get("decoy_manifest"),
            name=f"{task_id} decoy manifest",
        )
        scorefile, _ = validate_record(
            root=task_dir,
            raw=result.get("scorefile"),
            name=f"{task_id} scorefile",
        )
        results.append((task, decoy_manifest, scorefile))
    if len(results) != len(tasks) or len({item[0].task_id for item in results}) != len(
        tasks
    ):
        raise ValidationError("refinement manifest task coverage is incomplete")
    return tuple(results)


def _start_geometries(
    *, config: AppConfig, tasks: tuple[tuple[RefinementTask, Path, Path], ...]
) -> tuple[dict[str, ComplexGeometry], dict[tuple[str, str], ComplexGeometry]]:
    starts = {
        task.start.start_id: read_complex_geometry(
            path=task.start.input_path,
            interface_contact_A=config.validation.interface_contact_A,
        )
        for task, _, _ in tasks
    }
    references: dict[tuple[str, str], ComplexGeometry] = {}
    for task, _, _ in tasks:
        key = (task.start.candidate_id, task.start.receptor_id)
        references.setdefault(key, starts[task.start.start_id])
    return starts, references


def _analyzed_decoys(
    *,
    config: AppConfig,
    tasks: tuple[tuple[RefinementTask, Path, Path], ...],
    settings: ResolvedAnalysisSettings,
) -> tuple[RefinedDecoy, ...]:
    starts, references = _start_geometries(config=config, tasks=tasks)
    results: list[RefinedDecoy] = []
    seen_ids: set[str] = set()
    for task, manifest_path, score_path in tasks:
        score_by_id = {
            row.description: row.score(config.validation.refinement.ranking_score)
            for row in read_rosetta_scorefile(score_path)
        }
        for row in _decoy_rows(manifest_path):
            description = row["description"]
            expected_hash = row["sha256"]
            path = (manifest_path.parent / row["path"]).resolve()
            try:
                path.relative_to(manifest_path.parent.resolve())
                score = float(row["ranking_score"])
            except (ValueError, OverflowError) as exc:
                raise ValidationError(
                    f"invalid refinement decoy record: {description}"
                ) from exc
            reference_score = score_by_id.get(description)
            if (
                reference_score is None
                or not math.isfinite(score)
                or row["ranking_score"] != f"{reference_score:.6f}"
                or not path.is_file()
                or sha256_file(path) != expected_hash
            ):
                raise ValidationError(
                    f"refinement decoy record mismatch: {description}"
                )
            decoy_id = f"{task.task_id}__{description}"
            if decoy_id in seen_ids:
                raise ValidationError(f"duplicate refined decoy ID: {decoy_id}")
            seen_ids.add(decoy_id)
            results.append(
                assess_refined_decoy(
                    task=task,
                    decoy_id=decoy_id,
                    path=path,
                    path_sha256=expected_hash,
                    ranking_score=score,
                    start=starts[task.start.start_id],
                    cluster_reference=references[
                        (task.start.candidate_id, task.start.receptor_id)
                    ],
                    config=config,
                    settings=settings,
                )
            )
    return tuple(results)


def _tsv(fields: tuple[str, ...], rows: list[dict[str, str]]) -> str:
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def _write_decoys(
    *, path: Path, run_dir: Path, decoys: tuple[RefinedDecoy, ...]
) -> None:
    fields = (
        "decoy_id",
        "task_id",
        "start_id",
        "candidate_id",
        "receptor_id",
        "source_seed",
        "refinement_seed",
        "ranking_score",
        "interface_contact_pairs",
        "interface_receptor_residues",
        "minimum_interface_distance_A",
        "receptor_ca_rmsd_A",
        "start_contact_overlap",
        "start_site_displacement_A",
        "qc_status",
        "path",
        "sha256",
    )
    rows = [
        {
            "decoy_id": item.decoy_id,
            "task_id": item.task_id,
            "start_id": item.start_id,
            "candidate_id": item.candidate_id,
            "receptor_id": item.receptor_id,
            "source_seed": "" if item.source_seed is None else str(item.source_seed),
            "refinement_seed": str(item.refinement_seed),
            "ranking_score": f"{item.ranking_score:.6f}",
            "interface_contact_pairs": str(item.interface_contact_pairs),
            "interface_receptor_residues": str(item.interface_receptor_residues),
            "minimum_interface_distance_A": (
                f"{item.minimum_interface_distance_A:.6f}"
            ),
            "receptor_ca_rmsd_A": f"{item.receptor_ca_rmsd_A:.6f}",
            "start_contact_overlap": f"{item.start_contact_overlap:.6f}",
            "start_site_displacement_A": f"{item.start_site_displacement_A:.6f}",
            "qc_status": item.qc_status,
            "path": item.path.relative_to(run_dir).as_posix(),
            "sha256": item.sha256,
        }
        for item in decoys
    ]
    atomic_write_text(path, _tsv(fields, rows))


def _parameters(
    *, config: AppConfig, settings: ResolvedAnalysisSettings
) -> dict[str, JsonValue]:
    return {
        "interface_contact_A": config.validation.interface_contact_A,
        "min_interface_contact_pairs": settings.min_interface_contact_pairs,
        "min_interface_receptor_residues": (settings.min_interface_receptor_residues),
        "max_receptor_ca_rmsd_A": settings.max_receptor_ca_rmsd_A,
        "min_start_contact_overlap": settings.min_start_contact_overlap,
        "max_start_site_displacement_A": settings.max_start_site_displacement_A,
        "max_cluster_backbone_rmsd_A": settings.max_cluster_backbone_rmsd_A,
        "min_heavy_atom_distance_A": settings.min_heavy_atom_distance_A,
        "min_refinement_seed_support": settings.min_refinement_seed_support,
        "min_refinement_start_support": settings.min_refinement_start_support,
    }


def analyze_refinement_run(
    *, config: AppConfig, run_dir: Path
) -> RefinementAnalysisOutcome:
    """验证完整运行并写出不可覆盖的 decoy QC 与构象簇报告。"""
    output_dir = run_dir / "refinement_analysis"
    if output_dir.exists():
        raise ValidationError(f"refinement analysis already exists: {output_dir}")
    plan, planned_tasks = verify_refinement_plan(config=config, run_dir=run_dir)
    evidence_category, known_site_information_used = refinement_identity(plan)
    settings = resolve_analysis_settings(config.validation.analysis)
    task_results = _manifest_tasks(run_dir=run_dir, tasks=planned_tasks, plan=plan)
    decoys = _analyzed_decoys(config=config, tasks=task_results, settings=settings)
    clusters = cluster_refined_decoys(decoys=decoys, settings=settings)
    decoy_path = output_dir / "refined_decoys.tsv"
    cluster_path = output_dir / "refined_clusters.tsv"
    manifest_path = output_dir / "analysis_manifest.json"
    _write_decoys(path=decoy_path, run_dir=run_dir, decoys=decoys)
    decoy_by_id = {item.decoy_id: item for item in decoys}
    cluster_fields = (
        "cluster_id",
        "candidate_id",
        "receptor_id",
        "decoy_count",
        "refinement_seeds",
        "start_ids",
        "representative_decoy_id",
        "representative_path",
        "supported",
    )
    cluster_rows = [
        {
            "cluster_id": item.cluster_id,
            "candidate_id": item.candidate_id,
            "receptor_id": item.receptor_id,
            "decoy_count": str(len(item.decoy_ids)),
            "refinement_seeds": ";".join(map(str, item.refinement_seeds)),
            "start_ids": ";".join(item.start_ids),
            "representative_decoy_id": item.representative_decoy_id,
            "representative_path": decoy_by_id[item.representative_decoy_id]
            .path.relative_to(run_dir)
            .as_posix(),
            "supported": str(item.supported).lower(),
        }
        for item in clusters
    ]
    atomic_write_text(cluster_path, _tsv(cluster_fields, cluster_rows))
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-refinement-analysis-manifest/1",
            "stage": "validation_candidate_refinement_analysis",
            "status": "completed",
            "generated_at": utc_now(),
            "evidence_category": evidence_category,
            "known_site_information_used": known_site_information_used,
            "parameters": _parameters(config=config, settings=settings),
            "refined_decoys": {
                "path": decoy_path.name,
                "sha256": sha256_file(decoy_path),
                "count": len(decoys),
                "passed_count": sum(item.qc_status == "passed" for item in decoys),
            },
            "refined_clusters": {
                "path": cluster_path.name,
                "sha256": sha256_file(cluster_path),
                "count": len(clusters),
                "supported_count": sum(item.supported for item in clusters),
            },
        },
    )
    return RefinementAnalysisOutcome(manifest_path, len(decoys), len(clusters))
