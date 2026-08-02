"""阶段三 blind candidate 的跨证据综合复核。"""

from __future__ import annotations

import csv
import io
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
from vela.discovery.analysis.reports import read_site_analysis_report
from vela.discovery.models import DiscoveryError
from vela.discovery.sampling.planning import MAIN_DISCOVERY_EVIDENCE
from vela.validation.models import ValidationError
from vela.validation.records import file_record, read_document, validate_record
from vela.validation.refinement.planning import (
    BLIND_REFINEMENT_EVIDENCE,
    refinement_identity,
    verify_refinement_plan,
)


@dataclass(frozen=True, slots=True)
class CandidateReviewOutcome:
    """候选事实卡清单及其规模。"""

    manifest_path: Path
    candidate_count: int


def _table(path: Path, *, required: frozenset[str]) -> tuple[dict[str, str], ...]:
    if not path.is_file():
        raise ValidationError(f"candidate review table does not exist: {path}")
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValidationError(f"candidate review table columns are invalid: {path}")
        return tuple(dict(row) for row in reader)


def _replication(
    *, discovery_manifest: Path, replication_run_dir: Path
) -> tuple[Path, dict[str, tuple[str, ...]]]:
    root = replication_run_dir / "comparison"
    manifest_path = root / "comparison_manifest.json"
    manifest = read_document(manifest_path, name="replication comparison manifest")
    if (
        manifest.get("schema") != "vela.validation-replication-comparison-manifest/1"
        or manifest.get("status") != "completed"
        or manifest.get("known_site_information_used") is not False
        or manifest.get("classification_applied") is not False
    ):
        raise ValidationError("replication comparison identity is invalid")
    try:
        inputs = object_mapping(manifest.get("inputs"), name="comparison inputs")
        main = object_mapping(
            inputs.get("main_discovery_analysis"), name="main discovery analysis"
        )
    except TypeError as exc:
        raise ValidationError("replication comparison inputs are invalid") from exc
    if main.get("sha256") != sha256_file(discovery_manifest):
        raise ValidationError(
            "replication comparison uses a different main discovery analysis"
        )
    table_path, _ = validate_record(
        root=root,
        raw=manifest.get("candidate_replication"),
        name="candidate replication",
    )
    rows = _table(
        table_path,
        required=frozenset({"candidate_id", "replication_state_ids", "matched"}),
    )
    result: dict[str, tuple[str, ...]] = {}
    for row in rows:
        if row["matched"] not in {"true", "false"}:
            raise ValidationError("candidate replication matched value is invalid")
        states = tuple(item for item in row["replication_state_ids"].split(";") if item)
        if (row["matched"] == "true") != bool(states):
            raise ValidationError("candidate replication support is inconsistent")
        result[row["candidate_id"]] = states
    return manifest_path, result


def _refinement(
    *, config: AppConfig, run_dir: Path
) -> tuple[
    dict[str, object],
    Path,
    set[str],
    dict[str, tuple[str, ...]],
]:
    plan, _ = verify_refinement_plan(config=config, run_dir=run_dir)
    evidence, known = refinement_identity(plan)
    if evidence != BLIND_REFINEMENT_EVIDENCE or known:
        raise ValidationError("candidate review requires a blind refinement run")
    try:
        tasks = object_list(plan.get("tasks"), name="refinement tasks")
    except TypeError as exc:
        raise ValidationError("refinement tasks are invalid") from exc
    evaluated: set[str] = set()
    for raw in tasks:
        try:
            task = object_mapping(raw, name="refinement task")
        except TypeError as exc:
            raise ValidationError("refinement task is invalid") from exc
        candidate_id = task.get("candidate_id")
        if not isinstance(candidate_id, str):
            raise ValidationError("refinement candidate ID is invalid")
        evaluated.add(candidate_id)
    analysis_root = run_dir / "refinement_analysis"
    manifest_path = analysis_root / "analysis_manifest.json"
    manifest = read_document(manifest_path, name="refinement analysis manifest")
    if (
        manifest.get("schema") != "vela.validation-refinement-analysis-manifest/1"
        or manifest.get("status") != "completed"
        or manifest.get("evidence_category") != evidence
        or manifest.get("known_site_information_used") is not False
    ):
        raise ValidationError("blind refinement analysis identity is invalid")
    cluster_path, _ = validate_record(
        root=analysis_root,
        raw=manifest.get("refined_clusters"),
        name="refined clusters",
    )
    rows = _table(
        cluster_path,
        required=frozenset({"candidate_id", "receptor_id", "supported"}),
    )
    supported: dict[str, set[str]] = defaultdict(set)
    for row in rows:
        if row["supported"] == "true":
            supported[row["candidate_id"]].add(row["receptor_id"])
        elif row["supported"] != "false":
            raise ValidationError("refined cluster supported value is invalid")
    return (
        plan,
        manifest_path,
        evaluated,
        {key: tuple(sorted(values)) for key, values in supported.items()},
    )


def _environment(
    *, run_dir: Path, expected_source_evidence: str
) -> tuple[Path | None, dict[str, tuple[int, float, tuple[str, ...]]]]:
    root = run_dir / "environment_mapping"
    manifest_path = root / "mapping_manifest.json"
    if not manifest_path.is_file():
        return None, {}
    manifest = read_document(manifest_path, name="environment mapping manifest")
    if (
        manifest.get("schema") != "vela.environment-mapping-manifest/1"
        or manifest.get("status") != "completed"
        or manifest.get("source_evidence_category") != expected_source_evidence
        or manifest.get("classification_applied") is not False
    ):
        raise ValidationError("environment mapping identity is invalid")
    table_path, _ = validate_record(
        root=root,
        raw=manifest.get("mappings"),
        name="environment mappings",
    )
    rows = _table(
        table_path,
        required=frozenset(
            {
                "candidate_id",
                "layout_kind",
                "beta_contact_pairs",
                "minimum_beta_distance_A",
            }
        ),
    )
    grouped: dict[str, list[tuple[int, float, str]]] = defaultdict(list)
    for row in rows:
        try:
            pairs = int(row["beta_contact_pairs"])
            minimum = float(row["minimum_beta_distance_A"])
        except (ValueError, OverflowError) as exc:
            raise ValidationError("environment mapping metrics are invalid") from exc
        grouped[row["candidate_id"]].append((pairs, minimum, row["layout_kind"]))
    return manifest_path, {
        candidate_id: (
            max(item[0] for item in values),
            min(item[1] for item in values),
            tuple(sorted({item[2] for item in values})),
        )
        for candidate_id, values in grouped.items()
    }


def _status(*, evaluated: bool, locally_supported: bool, replicated: bool) -> str:
    if not evaluated:
        return "not_selected_for_local_refinement"
    if not locally_supported:
        return "locally_evaluated_without_supported_cluster"
    if replicated:
        return "locally_supported_with_bound_state_match"
    return "locally_supported_without_bound_state_match"


def _tsv(rows: list[dict[str, str]]) -> str:
    fields: list[str] = [
        "candidate_id",
        "target",
        "main_receptor_support",
        "main_receptor_ids",
        "replication_state_support",
        "replication_state_ids",
        "local_refinement_evaluated",
        "supported_refinement_receptor_count",
        "supported_refinement_receptor_ids",
        "environment_mapping_available",
        "maximum_beta_contact_pairs",
        "minimum_beta_distance_A",
        "environment_layout_kinds",
        "evidence_summary_status",
    ]
    output = io.StringIO()
    writer = csv.DictWriter(
        output, fieldnames=fields, delimiter="\t", lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(rows)
    return output.getvalue()


def write_candidate_review(
    *,
    config: AppConfig,
    discovery_run_dir: Path,
    replication_run_dir: Path,
    refinement_run_dir: Path,
) -> CandidateReviewOutcome:
    """汇总盲发现、结合态复现、局部精修和全酶映射事实。"""
    output_dir = refinement_run_dir / "candidate_review"
    if output_dir.exists():
        raise ValidationError(f"candidate review already exists: {output_dir}")
    try:
        discovery = read_site_analysis_report(
            run_dir=discovery_run_dir,
            expected_evidence_category=MAIN_DISCOVERY_EVIDENCE,
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
    replication_manifest, replication = _replication(
        discovery_manifest=discovery.manifest_path,
        replication_run_dir=replication_run_dir,
    )
    _, refinement_manifest, evaluated, refined = _refinement(
        config=config, run_dir=refinement_run_dir
    )
    environment_manifest, environment = _environment(
        run_dir=refinement_run_dir,
        expected_source_evidence=BLIND_REFINEMENT_EVIDENCE,
    )
    rows: list[dict[str, str]] = []
    for candidate in sorted(
        discovery.candidate_sites.values(), key=lambda item: item.candidate_id
    ):
        if not candidate.supported:
            continue
        replication_states = replication.get(candidate.candidate_id, ())
        refined_receptors = refined.get(candidate.candidate_id, ())
        environment_metrics = environment.get(candidate.candidate_id)
        rows.append(
            {
                "candidate_id": candidate.candidate_id,
                "target": candidate.target,
                "main_receptor_support": str(candidate.receptor_support),
                "main_receptor_ids": ";".join(candidate.receptor_ids),
                "replication_state_support": str(len(replication_states)),
                "replication_state_ids": ";".join(replication_states),
                "local_refinement_evaluated": str(
                    candidate.candidate_id in evaluated
                ).lower(),
                "supported_refinement_receptor_count": str(len(refined_receptors)),
                "supported_refinement_receptor_ids": ";".join(refined_receptors),
                "environment_mapping_available": str(
                    environment_metrics is not None
                ).lower(),
                "maximum_beta_contact_pairs": (
                    "" if environment_metrics is None else str(environment_metrics[0])
                ),
                "minimum_beta_distance_A": (
                    ""
                    if environment_metrics is None
                    else f"{environment_metrics[1]:.6f}"
                ),
                "environment_layout_kinds": (
                    ""
                    if environment_metrics is None
                    else ";".join(environment_metrics[2])
                ),
                "evidence_summary_status": _status(
                    evaluated=candidate.candidate_id in evaluated,
                    locally_supported=bool(refined_receptors),
                    replicated=bool(replication_states),
                ),
            }
        )
    if not rows:
        raise ValidationError("main discovery analysis has no supported candidates")
    table_path = output_dir / "candidate_evidence.tsv"
    atomic_write_text(table_path, _tsv(rows))
    inputs: dict[str, JsonValue] = {
        "main_discovery_analysis": {
            "path": discovery.manifest_path.resolve()
            .relative_to(config.paths.outputs_dir.resolve())
            .as_posix(),
            "sha256": sha256_file(discovery.manifest_path),
        },
        "replication_comparison": file_record(
            replication_manifest, root=config.paths.outputs_dir
        ),
        "refinement_analysis": file_record(
            refinement_manifest, root=config.paths.outputs_dir
        ),
        "environment_mapping": (
            None
            if environment_manifest is None
            else file_record(environment_manifest, root=config.paths.outputs_dir)
        ),
    }
    manifest_path = output_dir / "review_manifest.json"
    atomic_write_json(
        manifest_path,
        {
            "schema": "vela.candidate-review-manifest/1",
            "stage": "validation_candidate_review",
            "status": "completed",
            "generated_at": utc_now(),
            "evidence_category": "blind_candidate_cross_evidence_review",
            "known_site_information_used_for_discovery": False,
            "classification_applied": False,
            "interpretation": (
                "Evidence summary statuses report completed checks only; they are not "
                "affinity, inhibition, antitumor, or final candidate grades."
            ),
            "inputs": inputs,
            "candidate_evidence": {
                "path": table_path.name,
                "sha256": sha256_file(table_path),
                "count": len(rows),
            },
        },
    )
    return CandidateReviewOutcome(manifest_path, len(rows))
