"""阶段三局部精修运行的完整性检查、分析和报告。"""

from __future__ import annotations

import csv
import io
import math
from collections import Counter
from dataclasses import dataclass
from pathlib import Path

import gemmi

from vela.config import AppConfig
from vela.core.provenance import (
    JsonValue,
    atomic_write_json,
    atomic_write_text,
    sha256_file,
    utc_now,
    vela_software_identity,
)
from vela.core.typed_data import json_value, object_list, object_mapping
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.evidence import (
    align_receptor,
    read_structure,
    required_atom,
    split_model,
)
from vela.validation.models import ValidationError
from vela.validation.records import read_document, validate_record
from vela.validation.refinement.geometry import (
    ComplexGeometry,
    RefinedCluster,
    RefinedDecoy,
    ResolvedAnalysisSettings,
    assess_refined_decoy,
    cluster_refined_decoys,
    read_complex_geometry,
    resolve_analysis_settings,
)
from vela.validation.refinement.planning import (
    FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE,
    FUNNEL_DEEP_REFINEMENT_EVIDENCE,
    FUNNEL_SCREENING_REFINEMENT_EVIDENCE,
    RefinementTask,
    refinement_authorization,
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


@dataclass(frozen=True, slots=True)
class FunnelConfirmationOutcome:
    """Stage 3A与3B合并分析后的候选晋级结果。"""

    manifest_path: Path
    decoy_count: int
    cluster_count: int
    confirmed_candidate_ids: tuple[str, ...]


@dataclass(frozen=True, slots=True)
class FunnelDeepOutcome:
    """Stage 3A至3C合并分析后的最终原子姿态假设。"""

    manifest_path: Path
    decoy_count: int
    cluster_count: int
    final_hypothesis_cluster_ids: tuple[str, ...]


def is_funnel_confirmation_hit(
    *,
    cluster: RefinedCluster,
    settings: ResolvedAnalysisSettings,
    minimum_task_cells: int,
) -> bool:
    """应用Stage 3B冻结的跨source、跨seed和任务单元联合门槛。"""
    return (
        len(cluster.start_ids) >= settings.min_refinement_start_support
        and len(cluster.source_seeds) >= settings.min_refinement_source_seed_support
        and len(cluster.refinement_seeds) >= settings.min_refinement_seed_support
        and len(cluster.task_cells) >= minimum_task_cells
    )


def funnel_deep_source_seed_counts(cluster: RefinedCluster) -> dict[str, int]:
    """按全原子起点统计一个簇覆盖的独立Rosetta随机流数。"""
    seeds_by_start: dict[str, set[int]] = {}
    for cell in cluster.task_cells:
        start_id, separator, raw_seed = cell.rpartition(":")
        if not separator or not start_id:
            raise ValidationError("refined cluster task cell is invalid")
        try:
            seed = int(raw_seed)
        except ValueError as exc:
            raise ValidationError("refined cluster task cell seed is invalid") from exc
        seeds_by_start.setdefault(start_id, set()).add(seed)
    return {start_id: len(seeds) for start_id, seeds in sorted(seeds_by_start.items())}


def is_funnel_deep_hit(
    *,
    cluster: RefinedCluster,
    minimum_source_starts: int,
    minimum_rosetta_seeds_per_source: int,
) -> bool:
    """应用Stage 3C冻结的逐source独立随机流门槛。"""
    counts = funnel_deep_source_seed_counts(cluster)
    qualifying_starts = sum(
        count >= minimum_rosetta_seeds_per_source for count in counts.values()
    )
    return (
        len(cluster.source_seeds) >= minimum_source_starts
        and qualifying_starts >= minimum_source_starts
    )


def _sampling_software(plan: dict[str, object]) -> dict[str, JsonValue]:
    """保留历史精修计划中实际执行采样的软件身份。"""
    try:
        raw = object_mapping(plan.get("software"), name="refinement sampling software")
    except TypeError as exc:
        raise ValidationError("refinement sampling software is invalid") from exc
    result: dict[str, JsonValue] = {}
    for key in (
        "vela_version",
        "vela_source_sha256",
        "rosetta_version",
        "flexpepdock_sha256",
        "rosetta_scripts_sha256",
    ):
        value = raw.get(key)
        if not isinstance(value, str) or not value:
            raise ValidationError("refinement sampling software is invalid")
        result[key] = value
    return result


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
    allow_manifest_superset: bool = False,
) -> tuple[tuple[RefinementTask, Path, Path], ...]:
    manifest = read_document(
        run_dir / "refinement_manifest.json", name="refinement manifest"
    )
    if (
        manifest.get("schema") != "vela.validation-refinement-manifest/3"
        or manifest.get("stage") != "validation_local_refinement"
        or manifest.get("status") != "completed"
        or manifest.get("evidence_category") != plan.get("evidence_category")
        or manifest.get("known_site_information_used")
        is not plan.get("known_site_information_used")
        or manifest.get("source_evidence_category")
        != plan.get("source_evidence_category")
        or manifest.get("production_qualified") is not plan.get("production_qualified")
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
            if allow_manifest_superset:
                continue
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
        "qc_failures",
        "chemistry_failure",
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
            "qc_failures": ";".join(item.qc_failures),
            "chemistry_failure": item.chemistry_failure or "",
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
        "min_refinement_source_seed_support": (
            settings.min_refinement_source_seed_support
        ),
    }


def analyze_refinement_run(
    *, config: AppConfig, run_dir: Path
) -> RefinementAnalysisOutcome:
    """验证完整运行并写出不可覆盖的 decoy QC 与构象簇报告。"""
    output_dir = run_dir / "refinement_analysis"
    if output_dir.exists():
        raise ValidationError(f"refinement analysis already exists: {output_dir}")
    plan, planned_tasks = verify_refinement_plan(
        config=config,
        run_dir=run_dir,
        require_current_software=False,
    )
    evidence_category, known_site_information_used = refinement_identity(plan)
    source_evidence_category, production_qualified = refinement_authorization(plan)
    sampling_software = _sampling_software(plan)
    settings = resolve_analysis_settings(config.validation.analysis)
    task_results = _manifest_tasks(run_dir=run_dir, tasks=planned_tasks, plan=plan)
    decoys = _analyzed_decoys(config=config, tasks=task_results, settings=settings)
    clusters = cluster_refined_decoys(
        decoys=decoys,
        settings=settings,
        require_source_seed_support=not known_site_information_used,
    )
    failure_counts = Counter(
        failure for decoy in decoys for failure in decoy.qc_failures
    )
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
        "source_seeds",
        "task_cells",
        "evidence_status",
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
            "source_seeds": ";".join(map(str, item.source_seeds)),
            "task_cells": ";".join(item.task_cells),
            "evidence_status": (
                "cross_source_screening_hit"
                if evidence_category == FUNNEL_SCREENING_REFINEMENT_EVIDENCE
                and len(item.start_ids) >= 2
                and len(item.source_seeds) >= 2
                else ("supported" if item.supported else "not_supported")
            ),
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
            "schema": "vela.validation-refinement-analysis-manifest/4",
            "stage": "validation_candidate_refinement_analysis",
            "status": "completed",
            "generated_at": utc_now(),
            "source_evidence_category": source_evidence_category,
            "evidence_category": evidence_category,
            "known_site_information_used": known_site_information_used,
            "production_qualified": production_qualified,
            "sampling_software": sampling_software,
            "analysis_software": vela_software_identity(),
            "parameters": _parameters(config=config, settings=settings),
            "refined_decoys": {
                "path": decoy_path.name,
                "sha256": sha256_file(decoy_path),
                "count": len(decoys),
                "passed_count": sum(item.qc_status == "passed" for item in decoys),
                "failed_count": sum(item.qc_status == "failed" for item in decoys),
                "chemistry_invalid_count": sum(
                    item.chemistry_failure is not None for item in decoys
                ),
                "failure_counts": dict(sorted(failure_counts.items())),
            },
            "refined_clusters": {
                "path": cluster_path.name,
                "sha256": sha256_file(cluster_path),
                "count": len(clusters),
                "supported_count": sum(item.supported for item in clusters),
                "cross_source_screening_hit_count": sum(
                    evidence_category == FUNNEL_SCREENING_REFINEMENT_EVIDENCE
                    and len(item.start_ids) >= 2
                    and len(item.source_seeds) >= 2
                    for item in clusters
                ),
            },
        },
    )
    return RefinementAnalysisOutcome(manifest_path, len(decoys), len(clusters))


def _screening_analysis_from_confirmation_plan(
    *, config: AppConfig, plan: dict[str, object]
) -> Path:
    try:
        inputs = object_mapping(plan.get("inputs"), name="Stage 3B inputs")
        record = object_mapping(
            inputs.get("screening_analysis"), name="Stage 3A analysis record"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3A analysis record is invalid") from exc
    path, _ = validate_record(
        root=config.paths.outputs_dir,
        raw=record,
        name="Stage 3A analysis manifest",
    )
    root = (config.paths.outputs_dir / "validation" / "refinements").resolve()
    if (
        not path.is_relative_to(root)
        or path.name != "analysis_manifest.json"
        or path.parent.name != "refinement_analysis"
    ):
        raise ValidationError("Stage 3A analysis path is invalid")
    return path


def _confirmation_analysis_from_deep_plan(
    *, config: AppConfig, plan: dict[str, object]
) -> Path:
    try:
        inputs = object_mapping(plan.get("inputs"), name="Stage 3C inputs")
        record = object_mapping(
            inputs.get("confirmation_analysis"), name="Stage 3B analysis record"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3B analysis record is invalid") from exc
    path, _ = validate_record(
        root=config.paths.outputs_dir,
        raw=record,
        name="Stage 3B analysis manifest",
    )
    root = (config.paths.outputs_dir / "validation" / "refinements").resolve()
    if (
        not path.is_relative_to(root)
        or path.name != "analysis_manifest.json"
        or path.parent.name != "funnel_confirmation_analysis"
    ):
        raise ValidationError("Stage 3B analysis path is invalid")
    return path


def _cross_receptor_atomic_metrics(
    *, config: AppConfig, first: RefinedDecoy, second: RefinedDecoy
) -> tuple[float, float]:
    """在共同受体C-alpha坐标系中比较两个全原子代表姿态。"""
    try:
        first_structure = read_structure(first.path)
        second_structure = read_structure(second.path)
        if len(first_structure) != 1 or len(second_structure) != 1:
            raise ValidationError("refined representative must contain one model")
        first_receptor, first_peptide = split_model(
            first_structure[0], peptide_sequence=config.chemistry.sequence
        )
        second_receptor, second_peptide = split_model(
            second_structure[0], peptide_sequence=config.chemistry.sequence
        )
        alignment = align_receptor(
            receptor=second_receptor,
            reference_chains=(first_receptor,),
        )
    except DiscoveryError as exc:
        raise ValidationError("cross-receptor representative is invalid") from exc
    first_backbone = tuple(
        required_atom(residue, atom_name).pos
        for residue in first_peptide
        for atom_name in ("N", "CA", "C")
    )
    second_backbone = tuple(
        gemmi.Position(alignment.transform.apply(required_atom(residue, atom_name).pos))
        for residue in second_peptide
        for atom_name in ("N", "CA", "C")
    )
    if len(first_backbone) != len(second_backbone) or not first_backbone:
        raise ValidationError("cross-receptor peptide backbones do not correspond")
    backbone_rmsd = math.sqrt(
        sum(
            left.dist(right) ** 2
            for left, right in zip(first_backbone, second_backbone, strict=True)
        )
        / len(first_backbone)
    )
    first_contacts = read_complex_geometry(
        path=first.path,
        interface_contact_A=config.validation.interface_contact_A,
    ).receptor_contacts
    second_contacts = read_complex_geometry(
        path=second.path,
        interface_contact_A=config.validation.interface_contact_A,
    ).receptor_contacts
    union = first_contacts | second_contacts
    contact_jaccard = (
        len(first_contacts & second_contacts) / len(union) if union else 0.0
    )
    return backbone_rmsd, contact_jaccard


def analyze_funnel_confirmation(
    *, config: AppConfig, run_dir: Path
) -> FunnelConfirmationOutcome:
    """合并Stage 3A与3B decoy并应用冻结的3/4任务单元门槛。"""
    output_dir = run_dir / "funnel_confirmation_analysis"
    if output_dir.exists():
        raise ValidationError(
            f"funnel confirmation analysis already exists: {output_dir}"
        )
    confirmation_plan, confirmation_tasks = verify_refinement_plan(
        config=config,
        run_dir=run_dir,
        require_current_software=False,
    )
    if (
        confirmation_plan.get("source_evidence_category")
        != "exploratory_funnel_screening_handoff"
        or confirmation_plan.get("evidence_category")
        != FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE
        or confirmation_plan.get("known_site_information_used") is not False
        or confirmation_plan.get("production_qualified") is not False
    ):
        raise ValidationError("Stage 3B refinement identity is invalid")
    screening_analysis_path = _screening_analysis_from_confirmation_plan(
        config=config,
        plan=confirmation_plan,
    )
    screening_run_dir = screening_analysis_path.parent.parent
    screening_plan, screening_tasks = verify_refinement_plan(
        config=config,
        run_dir=screening_run_dir,
        require_current_software=False,
    )
    screening_analysis = read_document(
        screening_analysis_path,
        name="Stage 3A analysis manifest",
    )
    if (
        screening_plan.get("evidence_category") != FUNNEL_SCREENING_REFINEMENT_EVIDENCE
        or screening_analysis.get("evidence_category")
        != FUNNEL_SCREENING_REFINEMENT_EVIDENCE
        or screening_analysis.get("status") != "completed"
    ):
        raise ValidationError("Stage 3A evidence identity is invalid")
    try:
        selection = object_mapping(
            confirmation_plan.get("selection"), name="Stage 3B selection"
        )
        selected_candidates = tuple(
            value
            for value in object_list(
                selection.get("selected_candidate_ids"),
                name="Stage 3B selected candidates",
            )
            if isinstance(value, str)
        )
        promotion_contract = json_value(
            selection.get("promotion_contract"), name="Stage 3B promotion contract"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3B selection is invalid") from exc
    if (
        not selected_candidates
        or len(selected_candidates)
        != len(object_list(selection.get("selected_candidate_ids"), name="candidates"))
        or len(set(selected_candidates)) != len(selected_candidates)
    ):
        raise ValidationError("Stage 3B selected candidate IDs are invalid")
    selected_set = set(selected_candidates)
    screening_tasks = tuple(
        task for task in screening_tasks if task.start.candidate_id in selected_set
    )
    if not screening_tasks:
        raise ValidationError("Stage 3A contains no tasks for Stage 3B candidates")
    settings = resolve_analysis_settings(config.validation.analysis)
    screening_results = _manifest_tasks(
        run_dir=screening_run_dir,
        tasks=screening_tasks,
        plan=screening_plan,
        allow_manifest_superset=True,
    )
    confirmation_results = _manifest_tasks(
        run_dir=run_dir,
        tasks=confirmation_tasks,
        plan=confirmation_plan,
    )
    decoys = _analyzed_decoys(
        config=config,
        tasks=screening_results,
        settings=settings,
    ) + _analyzed_decoys(
        config=config,
        tasks=confirmation_results,
        settings=settings,
    )
    if len({item.decoy_id for item in decoys}) != len(decoys):
        raise ValidationError("Stage 3A and Stage 3B decoy identities overlap")
    clusters = cluster_refined_decoys(
        decoys=decoys,
        settings=settings,
        require_source_seed_support=True,
    )
    minimum_cells = config.validation.funnel.confirmation_min_task_cells
    hit_clusters = tuple(
        cluster
        for cluster in clusters
        if is_funnel_confirmation_hit(
            cluster=cluster,
            settings=settings,
            minimum_task_cells=minimum_cells,
        )
    )
    hit_candidates = {cluster.candidate_id for cluster in hit_clusters}
    confirmed_candidates = tuple(
        candidate_id
        for candidate_id in selected_candidates
        if candidate_id in hit_candidates
    )[: config.validation.funnel.confirmation_promotion_budget]
    failure_counts = Counter(
        failure for decoy in decoys for failure in decoy.qc_failures
    )
    decoy_path = output_dir / "combined_refined_decoys.tsv"
    cluster_path = output_dir / "combined_refined_clusters.tsv"
    manifest_path = output_dir / "analysis_manifest.json"
    _write_decoys(path=decoy_path, run_dir=config.paths.outputs_dir, decoys=decoys)
    decoy_by_id = {item.decoy_id: item for item in decoys}
    cluster_fields = (
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
    )
    hit_ids = {cluster.cluster_id for cluster in hit_clusters}
    cluster_rows = [
        {
            "cluster_id": cluster.cluster_id,
            "candidate_id": cluster.candidate_id,
            "receptor_id": cluster.receptor_id,
            "decoy_count": str(len(cluster.decoy_ids)),
            "refinement_seeds": ";".join(map(str, cluster.refinement_seeds)),
            "start_ids": ";".join(cluster.start_ids),
            "source_seeds": ";".join(map(str, cluster.source_seeds)),
            "task_cells": ";".join(cluster.task_cells),
            "evidence_status": (
                "source_seed_confirmation_hit"
                if cluster.cluster_id in hit_ids
                else "not_confirmed"
            ),
            "representative_decoy_id": cluster.representative_decoy_id,
            "representative_path": decoy_by_id[cluster.representative_decoy_id]
            .path.relative_to(config.paths.outputs_dir)
            .as_posix(),
        }
        for cluster in clusters
    ]
    atomic_write_text(cluster_path, _tsv(cluster_fields, cluster_rows))
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-funnel-confirmation-analysis/1",
            "stage": "stage3b_source_seed_confirmation",
            "status": "completed",
            "generated_at": utc_now(),
            "evidence_category": FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE,
            "production_qualified": False,
            "screening_analysis": {
                "path": screening_analysis_path.relative_to(
                    config.paths.outputs_dir
                ).as_posix(),
                "sha256": sha256_file(screening_analysis_path),
            },
            "confirmation_refinement_manifest": {
                "path": (run_dir / "refinement_manifest.json")
                .relative_to(config.paths.outputs_dir)
                .as_posix(),
                "sha256": sha256_file(run_dir / "refinement_manifest.json"),
            },
            "sampling_software": {
                "screening": _sampling_software(screening_plan),
                "confirmation": _sampling_software(confirmation_plan),
            },
            "analysis_software": vela_software_identity(),
            "parameters": {
                **_parameters(config=config, settings=settings),
                "minimum_source_seed_task_cells": minimum_cells,
                "promotion_budget": config.validation.funnel.confirmation_promotion_budget,
            },
            "promotion_contract": promotion_contract,
            "selected_candidate_ids": list(selected_candidates),
            "confirmed_candidate_ids": list(confirmed_candidates),
            "combined_refined_decoys": {
                "path": decoy_path.name,
                "sha256": sha256_file(decoy_path),
                "count": len(decoys),
                "passed_count": sum(item.qc_status == "passed" for item in decoys),
                "failed_count": sum(item.qc_status == "failed" for item in decoys),
                "failure_counts": dict(sorted(failure_counts.items())),
            },
            "combined_refined_clusters": {
                "path": cluster_path.name,
                "sha256": sha256_file(cluster_path),
                "count": len(clusters),
                "confirmation_hit_count": len(hit_clusters),
            },
        },
    )
    return FunnelConfirmationOutcome(
        manifest_path,
        len(decoys),
        len(clusters),
        confirmed_candidates,
    )


def analyze_funnel_deep_confirmation(
    *, config: AppConfig, run_dir: Path
) -> FunnelDeepOutcome:
    """合并Stage 3A至3C并生成最多两个最终原子姿态假设。"""
    output_dir = run_dir / "funnel_deep_analysis"
    if output_dir.exists():
        raise ValidationError(f"funnel deep analysis already exists: {output_dir}")
    deep_plan, deep_tasks = verify_refinement_plan(
        config=config,
        run_dir=run_dir,
        require_current_software=False,
    )
    if (
        deep_plan.get("source_evidence_category")
        != "exploratory_funnel_screening_handoff"
        or deep_plan.get("evidence_category") != FUNNEL_DEEP_REFINEMENT_EVIDENCE
        or deep_plan.get("known_site_information_used") is not False
        or deep_plan.get("production_qualified") is not False
    ):
        raise ValidationError("Stage 3C refinement identity is invalid")
    confirmation_analysis_path = _confirmation_analysis_from_deep_plan(
        config=config,
        plan=deep_plan,
    )
    confirmation_run_dir = confirmation_analysis_path.parent.parent
    confirmation_plan, confirmation_tasks = verify_refinement_plan(
        config=config,
        run_dir=confirmation_run_dir,
        require_current_software=False,
    )
    confirmation_analysis = read_document(
        confirmation_analysis_path,
        name="Stage 3B analysis manifest",
    )
    if (
        confirmation_plan.get("evidence_category")
        != FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE
        or confirmation_analysis.get("schema")
        != "vela.validation-funnel-confirmation-analysis/1"
        or confirmation_analysis.get("status") != "completed"
        or confirmation_analysis.get("evidence_category")
        != FUNNEL_CONFIRMATION_REFINEMENT_EVIDENCE
    ):
        raise ValidationError("Stage 3B evidence identity is invalid")
    screening_analysis_path = _screening_analysis_from_confirmation_plan(
        config=config,
        plan=confirmation_plan,
    )
    screening_run_dir = screening_analysis_path.parent.parent
    screening_plan, screening_tasks = verify_refinement_plan(
        config=config,
        run_dir=screening_run_dir,
        require_current_software=False,
    )
    if screening_plan.get("evidence_category") != FUNNEL_SCREENING_REFINEMENT_EVIDENCE:
        raise ValidationError("Stage 3A evidence identity is invalid")
    try:
        selection = object_mapping(
            deep_plan.get("selection"), name="Stage 3C selection"
        )
        selected_raw = object_list(
            selection.get("selected_candidate_ids"),
            name="Stage 3C selected candidates",
        )
        selected_candidates = tuple(
            value for value in selected_raw if isinstance(value, str)
        )
        promotion_contract = json_value(
            selection.get("promotion_contract"), name="Stage 3C promotion contract"
        )
    except TypeError as exc:
        raise ValidationError("Stage 3C selection is invalid") from exc
    if (
        not selected_candidates
        or len(selected_candidates) != len(selected_raw)
        or len(set(selected_candidates)) != len(selected_candidates)
    ):
        raise ValidationError("Stage 3C selected candidate IDs are invalid")
    selected_set = set(selected_candidates)
    screening_tasks = tuple(
        task for task in screening_tasks if task.start.candidate_id in selected_set
    )
    confirmation_tasks = tuple(
        task for task in confirmation_tasks if task.start.candidate_id in selected_set
    )
    if not screening_tasks or not confirmation_tasks or not deep_tasks:
        raise ValidationError("Stage 3C evidence chain contains no selected tasks")
    settings = resolve_analysis_settings(config.validation.analysis)
    screening_results = _manifest_tasks(
        run_dir=screening_run_dir,
        tasks=screening_tasks,
        plan=screening_plan,
        allow_manifest_superset=True,
    )
    confirmation_results = _manifest_tasks(
        run_dir=confirmation_run_dir,
        tasks=confirmation_tasks,
        plan=confirmation_plan,
    )
    deep_results = _manifest_tasks(
        run_dir=run_dir,
        tasks=deep_tasks,
        plan=deep_plan,
    )
    decoys = (
        _analyzed_decoys(config=config, tasks=screening_results, settings=settings)
        + _analyzed_decoys(
            config=config,
            tasks=confirmation_results,
            settings=settings,
        )
        + _analyzed_decoys(config=config, tasks=deep_results, settings=settings)
    )
    if len({item.decoy_id for item in decoys}) != len(decoys):
        raise ValidationError("Stage 3A, 3B, and 3C decoy identities overlap")
    clusters = cluster_refined_decoys(
        decoys=decoys,
        settings=settings,
        require_source_seed_support=True,
    )
    try:
        contract = object_mapping(promotion_contract, name="Stage 3 promotion contract")
        deep_contract = object_mapping(
            contract.get("stage3c_deep_confirmation"),
            name="Stage 3C promotion contract",
        )
    except TypeError as exc:
        raise ValidationError("Stage 3C promotion contract is invalid") from exc
    minimum_sources = deep_contract.get("minimum_source_starts")
    minimum_seeds = deep_contract.get("minimum_rosetta_seeds_per_source")
    if (
        not isinstance(minimum_sources, int)
        or isinstance(minimum_sources, bool)
        or minimum_sources < 1
        or not isinstance(minimum_seeds, int)
        or isinstance(minimum_seeds, bool)
        or minimum_seeds < 1
    ):
        raise ValidationError("Stage 3C source and seed thresholds are invalid")
    hit_clusters = tuple(
        cluster
        for cluster in clusters
        if is_funnel_deep_hit(
            cluster=cluster,
            minimum_source_starts=minimum_sources,
            minimum_rosetta_seeds_per_source=minimum_seeds,
        )
    )
    decoy_by_id = {item.decoy_id: item for item in decoys}
    max_backbone_rmsd = config.validation.funnel.cross_receptor_max_backbone_rmsd_A
    min_contact_jaccard = config.validation.funnel.cross_receptor_min_contact_jaccard
    compatibility_rows: list[dict[str, str]] = []
    compatible_cluster_ids: set[str] = set()
    for index, first in enumerate(hit_clusters):
        for second in hit_clusters[index + 1 :]:
            if (
                first.candidate_id != second.candidate_id
                or first.receptor_id == second.receptor_id
            ):
                continue
            backbone_rmsd, contact_jaccard = _cross_receptor_atomic_metrics(
                config=config,
                first=decoy_by_id[first.representative_decoy_id],
                second=decoy_by_id[second.representative_decoy_id],
            )
            compatible = (
                backbone_rmsd <= max_backbone_rmsd
                and contact_jaccard >= min_contact_jaccard
            )
            if compatible:
                compatible_cluster_ids.update((first.cluster_id, second.cluster_id))
            compatibility_rows.append(
                {
                    "candidate_id": first.candidate_id,
                    "first_cluster_id": first.cluster_id,
                    "second_cluster_id": second.cluster_id,
                    "peptide_backbone_rmsd_A": f"{backbone_rmsd:.6f}",
                    "receptor_contact_jaccard": f"{contact_jaccard:.6f}",
                    "compatible": str(compatible).lower(),
                }
            )
    hypothesis_budget = config.validation.funnel.final_hypothesis_budget
    final_clusters = hit_clusters[:hypothesis_budget]
    final_ids = tuple(cluster.cluster_id for cluster in final_clusters)
    failure_counts = Counter(
        failure for decoy in decoys for failure in decoy.qc_failures
    )
    decoy_path = output_dir / "combined_refined_decoys.tsv"
    cluster_path = output_dir / "combined_refined_clusters.tsv"
    compatibility_path = output_dir / "cross_receptor_compatibility.tsv"
    manifest_path = output_dir / "analysis_manifest.json"
    _write_decoys(path=decoy_path, run_dir=config.paths.outputs_dir, decoys=decoys)
    hit_ids = {cluster.cluster_id for cluster in hit_clusters}
    cluster_fields = (
        "cluster_id",
        "candidate_id",
        "receptor_id",
        "decoy_count",
        "refinement_seeds",
        "start_ids",
        "source_seeds",
        "task_cells",
        "rosetta_seed_counts_by_start",
        "evidence_status",
        "representative_decoy_id",
        "representative_path",
    )
    cluster_rows = [
        {
            "cluster_id": cluster.cluster_id,
            "candidate_id": cluster.candidate_id,
            "receptor_id": cluster.receptor_id,
            "decoy_count": str(len(cluster.decoy_ids)),
            "refinement_seeds": ";".join(map(str, cluster.refinement_seeds)),
            "start_ids": ";".join(cluster.start_ids),
            "source_seeds": ";".join(map(str, cluster.source_seeds)),
            "task_cells": ";".join(cluster.task_cells),
            "rosetta_seed_counts_by_start": ";".join(
                f"{start_id}:{count}"
                for start_id, count in funnel_deep_source_seed_counts(cluster).items()
            ),
            "evidence_status": (
                "deep_confirmation_hit"
                if cluster.cluster_id in hit_ids
                else "not_confirmed"
            ),
            "representative_decoy_id": cluster.representative_decoy_id,
            "representative_path": decoy_by_id[cluster.representative_decoy_id]
            .path.relative_to(config.paths.outputs_dir)
            .as_posix(),
        }
        for cluster in clusters
    ]
    atomic_write_text(cluster_path, _tsv(cluster_fields, cluster_rows))
    compatibility_fields = (
        "candidate_id",
        "first_cluster_id",
        "second_cluster_id",
        "peptide_backbone_rmsd_A",
        "receptor_contact_jaccard",
        "compatible",
    )
    atomic_write_text(
        compatibility_path,
        _tsv(compatibility_fields, compatibility_rows),
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.validation-funnel-deep-analysis/1",
            "stage": "stage3c_deep_confirmation",
            "status": "completed",
            "generated_at": utc_now(),
            "evidence_category": FUNNEL_DEEP_REFINEMENT_EVIDENCE,
            "production_qualified": False,
            "confirmation_analysis": {
                "path": confirmation_analysis_path.relative_to(
                    config.paths.outputs_dir
                ).as_posix(),
                "sha256": sha256_file(confirmation_analysis_path),
            },
            "deep_refinement_manifest": {
                "path": (run_dir / "refinement_manifest.json")
                .relative_to(config.paths.outputs_dir)
                .as_posix(),
                "sha256": sha256_file(run_dir / "refinement_manifest.json"),
            },
            "sampling_software": {
                "screening": _sampling_software(screening_plan),
                "confirmation": _sampling_software(confirmation_plan),
                "deep": _sampling_software(deep_plan),
            },
            "analysis_software": vela_software_identity(),
            "parameters": {
                **_parameters(config=config, settings=settings),
                "minimum_source_starts": minimum_sources,
                "minimum_rosetta_seeds_per_source": minimum_seeds,
                "cross_receptor_maximum_peptide_backbone_rmsd_A": (max_backbone_rmsd),
                "cross_receptor_minimum_receptor_contact_jaccard": (
                    min_contact_jaccard
                ),
                "final_hypothesis_budget": hypothesis_budget,
            },
            "promotion_contract": promotion_contract,
            "selected_candidate_ids": list(selected_candidates),
            "final_hypothesis_cluster_ids": list(final_ids),
            "final_hypothesis_candidate_ids": list(
                dict.fromkeys(cluster.candidate_id for cluster in final_clusters)
            ),
            "combined_refined_decoys": {
                "path": decoy_path.name,
                "sha256": sha256_file(decoy_path),
                "count": len(decoys),
                "passed_count": sum(item.qc_status == "passed" for item in decoys),
                "failed_count": sum(item.qc_status == "failed" for item in decoys),
                "failure_counts": dict(sorted(failure_counts.items())),
            },
            "combined_refined_clusters": {
                "path": cluster_path.name,
                "sha256": sha256_file(cluster_path),
                "count": len(clusters),
                "deep_confirmation_hit_count": len(hit_clusters),
            },
            "cross_receptor_compatibility": {
                "path": compatibility_path.name,
                "sha256": sha256_file(compatibility_path),
                "count": len(compatibility_rows),
                "compatible_pair_count": sum(
                    row["compatible"] == "true" for row in compatibility_rows
                ),
                "compatible_cluster_ids": sorted(compatible_cluster_ids),
            },
        },
    )
    return FunnelDeepOutcome(manifest_path, len(decoys), len(clusters), final_ids)
