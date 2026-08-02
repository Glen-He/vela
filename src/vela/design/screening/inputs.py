"""阶段三精修证据到阶段四初筛模板的严格边界。"""

from __future__ import annotations

import csv
from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.design.models import DesignError, DesignTemplate
from vela.validation.records import read_document, safe_identifier, validate_record
from vela.validation.refinement.planning import BLIND_REFINEMENT_EVIDENCE
from vela.validation.refinement.reconstruction import validate_flexpepdock_input


def _table(path: Path, *, required: frozenset[str]) -> tuple[dict[str, str], ...]:
    if not path.is_file():
        raise DesignError(f"design source table does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise DesignError(f"design source table columns are invalid: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows:
        raise DesignError(f"design source table contains no rows: {path}")
    return rows


def _source_identity(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, str], dict[tuple[str, str], tuple[str, str]]]:
    plan_path = run_dir / "refinement_plan.json"
    plan = read_document(plan_path, name="design source refinement plan")
    if (
        plan.get("schema") != "vela.validation-refinement-plan/1"
        or plan.get("stage") != "validation_local_refinement"
        or plan.get("status") != "planned"
        or plan.get("chemistry_id") != config.chemistry.chemistry_id
        or plan.get("evidence_category") != BLIND_REFINEMENT_EVIDENCE
        or plan.get("known_site_information_used") is not False
    ):
        raise DesignError("Stage 4 requires a blind Stage 3 refinement plan")
    try:
        task_rows = object_list(plan.get("tasks"), name="refinement plan tasks")
    except TypeError as exc:
        raise DesignError("design source refinement tasks are invalid") from exc
    target_by_key: dict[tuple[str, str], tuple[str, str]] = {}
    for raw in task_rows:
        try:
            task = object_mapping(raw, name="refinement plan task")
        except TypeError as exc:
            raise DesignError("design source refinement task is invalid") from exc
        candidate_id = safe_identifier(task.get("candidate_id"), name="candidate ID")
        receptor_id = safe_identifier(task.get("receptor_id"), name="receptor ID")
        target = safe_identifier(task.get("target"), name="target ID")
        key = (candidate_id, receptor_id)
        previous = target_by_key.setdefault(key, (target, receptor_id))
        if previous[0] != target:
            raise DesignError("refinement plan target identity is inconsistent")

    manifest = read_document(
        run_dir / "refinement_manifest.json", name="design source refinement manifest"
    )
    if (
        manifest.get("schema") != "vela.validation-refinement-manifest/1"
        or manifest.get("stage") != "validation_local_refinement"
        or manifest.get("status") != "completed"
        or manifest.get("chemistry_id") != config.chemistry.chemistry_id
        or manifest.get("evidence_category") != BLIND_REFINEMENT_EVIDENCE
        or manifest.get("known_site_information_used") is not False
    ):
        raise DesignError("Stage 3 refinement manifest is not eligible for design")
    recorded_plan, recorded_hash = validate_record(
        root=run_dir,
        raw=manifest.get("refinement_plan"),
        name="design source refinement plan",
    )
    if recorded_plan != plan_path or recorded_hash != sha256_file(plan_path):
        raise DesignError("Stage 3 refinement plan record is inconsistent")

    review_root = run_dir / "candidate_review"
    review = read_document(
        review_root / "review_manifest.json", name="candidate review manifest"
    )
    if (
        review.get("schema") != "vela.candidate-review-manifest/1"
        or review.get("status") != "completed"
        or review.get("known_site_information_used_for_discovery") is not False
        or review.get("classification_applied") is not False
    ):
        raise DesignError("Stage 3 candidate review identity is invalid")
    evidence_path, _ = validate_record(
        root=review_root,
        raw=review.get("candidate_evidence"),
        name="candidate evidence table",
    )
    evidence_rows = _table(
        evidence_path,
        required=frozenset({"candidate_id", "target", "evidence_summary_status"}),
    )
    evidence = {
        row["candidate_id"]: row["evidence_summary_status"] for row in evidence_rows
    }
    return evidence, target_by_key


def selected_templates(
    *,
    config: AppConfig,
    refinement_run_dir: Path,
    target_cluster_ids: tuple[str, ...],
) -> tuple[DesignTemplate, ...]:
    """核对单亚型的显式所选 cluster, 并冻结其代表结构。"""
    if not target_cluster_ids:
        raise DesignError("at least one target Stage 3 cluster is required")
    selected = target_cluster_ids
    if len(selected) != len(set(selected)):
        raise DesignError("design cluster selections must be unique")
    evidence, target_by_key = _source_identity(
        config=config, run_dir=refinement_run_dir
    )

    analysis_root = refinement_run_dir / "refinement_analysis"
    analysis = read_document(
        analysis_root / "analysis_manifest.json", name="refinement analysis manifest"
    )
    if (
        analysis.get("schema") != "vela.validation-refinement-analysis-manifest/1"
        or analysis.get("status") != "completed"
        or analysis.get("evidence_category") != BLIND_REFINEMENT_EVIDENCE
        or analysis.get("known_site_information_used") is not False
    ):
        raise DesignError("Stage 3 refinement analysis is not blind completed evidence")
    cluster_path, _ = validate_record(
        root=analysis_root,
        raw=analysis.get("refined_clusters"),
        name="refined clusters",
    )
    decoy_path, _ = validate_record(
        root=analysis_root,
        raw=analysis.get("refined_decoys"),
        name="refined decoys",
    )
    clusters = _table(
        cluster_path,
        required=frozenset(
            {
                "cluster_id",
                "candidate_id",
                "receptor_id",
                "representative_decoy_id",
                "representative_path",
                "supported",
            }
        ),
    )
    decoys = _table(
        decoy_path,
        required=frozenset({"decoy_id", "path", "sha256", "qc_status"}),
    )
    decoy_by_id = {row["decoy_id"]: row for row in decoys}
    if len(decoy_by_id) != len(decoys):
        raise DesignError("refinement decoy identities are duplicated")
    cluster_by_id = {row["cluster_id"]: row for row in clusters}
    if len(cluster_by_id) != len(clusters):
        raise DesignError("refinement cluster identities are duplicated")
    missing = tuple(
        cluster_id for cluster_id in selected if cluster_id not in cluster_by_id
    )
    if missing:
        raise DesignError(
            "selected Stage 3 clusters do not exist: " + ", ".join(missing)
        )

    templates: list[DesignTemplate] = []
    for index, cluster_id in enumerate(selected, 1):
        row = cluster_by_id[cluster_id]
        if row["supported"] != "true":
            raise DesignError(
                f"selected Stage 3 cluster is not supported: {cluster_id}"
            )
        candidate_id = safe_identifier(row["candidate_id"], name="candidate ID")
        receptor_id = safe_identifier(row["receptor_id"], name="receptor ID")
        if not evidence.get(candidate_id, "").startswith("locally_supported_"):
            raise DesignError(
                f"selected candidate lacks a locally supported review state: {candidate_id}"
            )
        target_record = target_by_key.get((candidate_id, receptor_id))
        if target_record is None:
            raise DesignError(
                "selected cluster is absent from the refinement task matrix"
            )
        representative = decoy_by_id.get(row["representative_decoy_id"])
        if representative is None or representative["qc_status"] != "passed":
            raise DesignError("selected cluster representative failed structural QC")
        if representative["path"] != row["representative_path"]:
            raise DesignError("cluster and decoy representative paths disagree")
        path = (refinement_run_dir / representative["path"]).resolve()
        try:
            path.relative_to(refinement_run_dir.resolve())
        except ValueError as exc:
            raise DesignError("selected template path escapes its Stage 3 run") from exc
        if not path.is_file() or sha256_file(path) != representative["sha256"]:
            raise DesignError(f"selected template hash mismatch: {path}")
        receptor_count, histidines = validate_flexpepdock_input(
            path=path,
            chemistry=config.chemistry,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
        )
        templates.append(
            DesignTemplate(
                template_id=f"template_{index:03d}",
                evidence_role="positive",
                cluster_id=safe_identifier(cluster_id, name="cluster ID"),
                candidate_id=candidate_id,
                receptor_id=receptor_id,
                target=target_record[0],
                path=path,
                sha256=representative["sha256"],
                receptor_residue_count=receptor_count,
                fixed_histidine_pose_indices=histidines,
            )
        )
    validate_objective(config=config, templates=tuple(templates))
    return tuple(templates)


def validate_objective(
    *, config: AppConfig, templates: tuple[DesignTemplate, ...]
) -> None:
    """要求所有模板只属于一个获得支持的目标亚型。"""
    if config.design.objective is None:
        raise DesignError("design objective is unresolved")
    targets = {item.target for item in templates}
    if len(targets) != 1 or not targets <= {"ck2_alpha", "ck2_alpha_prime"}:
        raise DesignError(
            "single_supported_target requires templates from exactly one CK2 subtype"
        )
