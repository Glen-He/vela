"""Cross-receptor native-aware 开发诊断的事后分析。"""

from __future__ import annotations

import csv
import io
import math
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from statistics import median

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
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.sampling.evidence import (
    align_receptor,
    read_structure,
    required_atom,
    split_model,
)
from vela.validation.models import ValidationError
from vela.validation.records import (
    read_document,
    safe_identifier,
    validate_record,
)
from vela.validation.refinement.qualification_diagnostic import (
    EVIDENCE_CATEGORY,
    MANIFEST_NAME,
    MANIFEST_SCHEMA,
    PLAN_NAME,
    PLAN_SCHEMA,
    PREPACK_SCHEMA,
    TASK_SCHEMA,
)
from vela.validation.scores import read_rosetta_scorefile

ANALYSIS_DIRECTORY = "qualification_refinement_analysis"
REPORT_NAME = "analysis_report.json"
REPORT_SCHEMA = "vela.validation-qualification-refinement-analysis-report/3"


@dataclass(frozen=True, slots=True)
class AnalysisTask:
    """一个冻结诊断任务的分析身份。"""

    task_id: str
    start_id: str
    receptor_site_id: str
    source_seed: int
    refinement_seed: int


@dataclass(frozen=True, slots=True)
class StartAssessment:
    """实际送入 FlexPepDock 的预打包起点相对 native 的几何。"""

    start_id: str
    receptor_site_id: str
    source_seed: int
    path: Path
    sha256: str
    receptor_alignment_rmsd_A: float
    native_backbone_rmsd_A: float


@dataclass(frozen=True, slots=True)
class DiagnosticDecoy:
    """带化学有效性、native 几何和来源身份的精修 decoy。"""

    decoy_id: str
    description: str
    task_id: str
    start_id: str
    receptor_site_id: str
    source_seed: int
    refinement_seed: int
    path: Path
    sha256: str
    ranking_score: float
    chemistry_valid: bool
    chemistry_failure: str | None
    receptor_alignment_rmsd_A: float
    native_backbone_rmsd_A: float
    aligned_backbone: tuple[gemmi.Position, ...]


@dataclass(frozen=True, slots=True)
class DiagnosticCluster:
    """同一受体 site 内按肽骨架聚合的分数有序构象簇。"""

    cluster_id: str
    receptor_site_id: str
    rank: int
    members: tuple[DiagnosticDecoy, ...]
    energy_representative_id: str
    geometric_medoid_id: str
    refinement_seeds: tuple[int, ...]
    start_ids: tuple[str, ...]
    recovered_seeds: tuple[int, ...]
    recovered_start_ids: tuple[str, ...]
    minimum_native_backbone_rmsd_A: float
    supported: bool


@dataclass(frozen=True, slots=True)
class QualificationRefinementAnalysisOutcome:
    """开发诊断分析报告及有效模型规模。"""

    report_path: Path
    valid_decoy_count: int
    cluster_count: int
    recovery_supported: bool


def _integer(value: object, *, name: str, minimum: int = 0) -> int:
    if not isinstance(value, int) or isinstance(value, bool) or value < minimum:
        raise ValidationError(f"{name} must be an integer >= {minimum}")
    return value


def _number(value: object, *, name: str, positive: bool = False) -> float:
    if isinstance(value, bool) or not isinstance(value, int | float):
        raise ValidationError(f"{name} must be numeric")
    result = float(value)
    if not math.isfinite(result) or (positive and result <= 0):
        raise ValidationError(f"{name} is outside its allowed range")
    return result


def _text(value: object, *, name: str) -> str:
    if not isinstance(value, str) or not value:
        raise ValidationError(f"{name} must be non-empty text")
    return value


def _parameters(plan: dict[str, object]) -> dict[str, object]:
    try:
        return object_mapping(plan.get("parameters"), name="diagnostic parameters")
    except TypeError as exc:
        raise ValidationError("diagnostic parameters are invalid") from exc


def _receptor_backbone_parameters(
    parameters: dict[str, object],
) -> dict[str, JsonValue]:
    try:
        raw = object_mapping(
            parameters.get("receptor_backbone"),
            name="diagnostic receptor backbone protocol",
        )
    except TypeError as exc:
        raise ValidationError(
            "diagnostic receptor backbone protocol is invalid"
        ) from exc
    mode = _text(raw.get("mode"), name="receptor backbone mode")
    if mode not in {"fixed", "local_constrained"}:
        raise ValidationError("diagnostic receptor backbone mode is invalid")
    return {
        "mode": mode,
        "selection": _text(raw.get("selection"), name="receptor backbone selection"),
        "contact_A": _number(
            raw.get("contact_A"), name="receptor backbone contact", positive=True
        ),
        "sequence_padding": _integer(
            raw.get("sequence_padding"), name="receptor backbone sequence padding"
        ),
        "movemap_policy": _text(
            raw.get("movemap_policy"), name="receptor backbone MoveMap policy"
        ),
        "coordinate_constraint": _text(
            raw.get("coordinate_constraint"),
            name="receptor backbone coordinate constraint",
        ),
    }


def _plan_tasks(plan: dict[str, object]) -> tuple[AnalysisTask, ...]:
    try:
        rows = object_list(plan.get("tasks"), name="diagnostic plan tasks")
    except TypeError as exc:
        raise ValidationError("diagnostic plan tasks are invalid") from exc
    tasks: list[AnalysisTask] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="diagnostic plan task")
        except TypeError as exc:
            raise ValidationError("diagnostic plan task is invalid") from exc
        if row.get("status") != "planned":
            raise ValidationError("diagnostic plan task is not frozen as planned")
        tasks.append(
            AnalysisTask(
                safe_identifier(row.get("task_id"), name="diagnostic task ID"),
                safe_identifier(row.get("start_id"), name="diagnostic start ID"),
                safe_identifier(
                    row.get("receptor_site_id"), name="diagnostic receptor site ID"
                ),
                _integer(row.get("source_seed"), name="source seed"),
                _integer(row.get("refinement_seed"), name="refinement seed"),
            )
        )
    if not tasks or len({task.task_id for task in tasks}) != len(tasks):
        raise ValidationError("diagnostic plan task identities are incomplete")
    return tuple(tasks)


def _validated_run(
    *, config: AppConfig, run_dir: Path
) -> tuple[
    dict[str, object],
    dict[str, object],
    tuple[AnalysisTask, ...],
    Path,
    str,
]:
    plan_path = run_dir / PLAN_NAME
    plan = read_document(plan_path, name="qualification refinement plan")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or plan.get("stage") != "validation_qualification_refinement_diagnostic"
        or plan.get("status") != "planned"
        or plan.get("development_only") is not True
        or plan.get("formal_qualification_gate") is not False
        or plan.get("native_information_used_for_task_selection") is not True
    ):
        raise ValidationError("qualification refinement plan identity is invalid")
    manifest = read_document(
        run_dir / MANIFEST_NAME, name="qualification refinement manifest"
    )
    if (
        manifest.get("schema") != MANIFEST_SCHEMA
        or manifest.get("stage") != "validation_qualification_refinement_diagnostic"
        or manifest.get("status") != "completed"
        or manifest.get("development_only") is not True
        or manifest.get("formal_qualification_gate") is not False
        or manifest.get("native_information_used_for_task_selection") is not True
        or manifest.get("method_id") != plan.get("method_id")
        or manifest.get("chemistry") != plan.get("chemistry")
    ):
        raise ValidationError("qualification refinement manifest identity is invalid")
    parameters = _parameters(plan)
    backbone = _receptor_backbone_parameters(parameters)
    if manifest.get("receptor_backbone_mode") != backbone.get("mode"):
        raise ValidationError("diagnostic receptor backbone identity is inconsistent")
    frozen_plan, plan_hash = validate_record(
        root=run_dir,
        raw=manifest.get("qualification_refinement_plan"),
        name="qualification refinement plan",
    )
    if frozen_plan != plan_path.resolve():
        raise ValidationError("diagnostic manifest references a different plan")
    tasks = _plan_tasks(plan)
    if manifest.get("task_count") != len(tasks):
        raise ValidationError("diagnostic manifest task count is inconsistent")
    expected_decoys = len(tasks) * _integer(
        parameters.get("decoys_per_seed"),
        name="decoys per seed",
        minimum=1,
    )
    if manifest.get("decoy_count") != expected_decoys:
        raise ValidationError("diagnostic manifest decoy count is inconsistent")
    try:
        inputs = object_mapping(plan.get("inputs"), name="diagnostic inputs")
        native_record = object_mapping(
            inputs.get("native_reference"), name="native reference"
        )
    except TypeError as exc:
        raise ValidationError("diagnostic native reference is invalid") from exc
    relative = native_record.get("path")
    expected_hash = native_record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError("diagnostic native reference is invalid")
    data_root = config.paths.data_dir.resolve()
    native_path = (data_root / relative).resolve()
    if (
        not native_path.is_relative_to(data_root)
        or not native_path.is_file()
        or sha256_file(native_path) != expected_hash
    ):
        raise ValidationError("diagnostic native reference is missing or changed")
    return plan, manifest, tasks, native_path, plan_hash


def _native_coordinates(
    *, path: Path, peptide_sequence: str
) -> tuple[tuple[gemmi.Residue, ...], tuple[gemmi.Position, ...]]:
    structure = read_structure(path)
    if len(structure) != 1:
        raise ValidationError("diagnostic native reference must contain one model")
    receptor, peptide = split_model(structure[0], peptide_sequence=peptide_sequence)
    backbone = tuple(
        required_atom(residue, atom_name).pos
        for residue in peptide
        for atom_name in ("N", "CA", "C")
    )
    if not backbone:
        raise ValidationError("diagnostic native peptide backbone is empty")
    return receptor, backbone


def _aligned_backbone(
    *,
    path: Path,
    peptide_sequence: str,
    native_receptor: tuple[gemmi.Residue, ...],
) -> tuple[float, tuple[gemmi.Position, ...]]:
    structure = read_structure(path)
    if len(structure) != 1:
        raise ValidationError(f"diagnostic structure must contain one model: {path}")
    receptor, peptide = split_model(structure[0], peptide_sequence=peptide_sequence)
    alignment = align_receptor(receptor=receptor, reference_chains=(native_receptor,))
    backbone = tuple(
        gemmi.Position(alignment.transform.apply(required_atom(residue, atom_name).pos))
        for residue in peptide
        for atom_name in ("N", "CA", "C")
    )
    return alignment.rmsd, backbone


def _rmsd(
    first: tuple[gemmi.Position, ...], second: tuple[gemmi.Position, ...]
) -> float:
    if len(first) != len(second) or not first:
        raise ValidationError("diagnostic peptide backbone correspondence is invalid")
    return math.sqrt(
        sum(left.dist(right) ** 2 for left, right in zip(first, second, strict=True))
        / len(first)
    )


def _starts(
    *,
    run_dir: Path,
    tasks: tuple[AnalysisTask, ...],
    plan_hash: str,
    peptide_sequence: str,
    native_receptor: tuple[gemmi.Residue, ...],
    native_backbone: tuple[gemmi.Position, ...],
) -> tuple[StartAssessment, ...]:
    by_start: dict[str, AnalysisTask] = {}
    for task in tasks:
        existing = by_start.setdefault(task.start_id, task)
        if (
            existing.receptor_site_id != task.receptor_site_id
            or existing.source_seed != task.source_seed
        ):
            raise ValidationError("diagnostic start identity is inconsistent")
    results: list[StartAssessment] = []
    for start_id, task in sorted(by_start.items()):
        start_dir = run_dir / "starts" / start_id
        document = read_document(
            start_dir / "prepack_result.json", name="diagnostic prepack result"
        )
        if (
            document.get("schema") != PREPACK_SCHEMA
            or document.get("status") != "completed"
            or document.get("start_id") != start_id
            or document.get("diagnostic_plan_sha256") != plan_hash
        ):
            raise ValidationError(f"diagnostic prepack identity is invalid: {start_id}")
        output, output_hash = validate_record(
            root=start_dir,
            raw=document.get("output"),
            name=f"{start_id} prepack output",
        )
        receptor_rmsd, backbone = _aligned_backbone(
            path=output,
            peptide_sequence=peptide_sequence,
            native_receptor=native_receptor,
        )
        results.append(
            StartAssessment(
                start_id,
                task.receptor_site_id,
                task.source_seed,
                output,
                output_hash,
                receptor_rmsd,
                _rmsd(backbone, native_backbone),
            )
        )
    return tuple(results)


def _decoy_manifest(path: Path) -> tuple[dict[str, str], ...]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames != [
            "description",
            "ranking_score",
            "path",
            "sha256",
        ]:
            raise ValidationError(f"invalid diagnostic decoy columns: {path}")
        rows = tuple(dict(row) for row in reader)
    if not rows:
        raise ValidationError(f"diagnostic decoy manifest contains no rows: {path}")
    return rows


def _chemistry_failures(
    document: dict[str, object], *, expected_count: int
) -> dict[str, str]:
    try:
        qc = object_mapping(document.get("chemistry_qc"), name="chemistry QC")
        rows = object_list(qc.get("failures"), name="chemistry QC failures")
    except TypeError as exc:
        raise ValidationError("diagnostic chemistry QC is invalid") from exc
    passed = _integer(qc.get("passed_decoy_count"), name="passed decoy count")
    failed = _integer(qc.get("failed_decoy_count"), name="failed decoy count")
    failures: dict[str, str] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="chemistry QC failure")
        except TypeError as exc:
            raise ValidationError("diagnostic chemistry QC failure is invalid") from exc
        description = safe_identifier(
            row.get("description"), name="failed decoy description"
        )
        reason = _text(row.get("failure"), name="chemistry failure reason")
        if description in failures:
            raise ValidationError("duplicate diagnostic chemistry QC failure")
        failures[description] = reason
    expected_status = "passed" if not failures else "failed"
    if (
        passed + failed != expected_count
        or failed != len(failures)
        or qc.get("status") != expected_status
    ):
        raise ValidationError("diagnostic chemistry QC counts are inconsistent")
    return failures


def _manifest_task_results(
    *,
    run_dir: Path,
    manifest: dict[str, object],
    tasks: tuple[AnalysisTask, ...],
) -> dict[str, Path]:
    try:
        rows = object_list(manifest.get("tasks"), name="diagnostic manifest tasks")
    except TypeError as exc:
        raise ValidationError("diagnostic manifest tasks are invalid") from exc
    by_id = {task.task_id: task for task in tasks}
    results: dict[str, Path] = {}
    for raw in rows:
        try:
            row = object_mapping(raw, name="diagnostic manifest task")
        except TypeError as exc:
            raise ValidationError("diagnostic manifest task is invalid") from exc
        task_id = row.get("task_id")
        task = by_id.get(task_id) if isinstance(task_id, str) else None
        if task is None or any(
            row.get(key) != expected
            for key, expected in (
                ("start_id", task.start_id),
                ("receptor_site_id", task.receptor_site_id),
                ("source_seed", task.source_seed),
                ("refinement_seed", task.refinement_seed),
            )
        ):
            raise ValidationError("diagnostic manifest task identity is inconsistent")
        result, _ = validate_record(
            root=run_dir,
            raw=row.get("task_result"),
            name=f"{task.task_id} task result",
        )
        if task.task_id in results:
            raise ValidationError("duplicate diagnostic manifest task")
        results[task.task_id] = result
    if set(results) != set(by_id):
        raise ValidationError("diagnostic manifest task coverage is incomplete")
    return results


def _decoys(
    *,
    run_dir: Path,
    manifest: dict[str, object],
    tasks: tuple[AnalysisTask, ...],
    plan_hash: str,
    expected_per_task: int,
    ranking_score_name: str,
    peptide_sequence: str,
    native_receptor: tuple[gemmi.Residue, ...],
    native_backbone: tuple[gemmi.Position, ...],
) -> tuple[DiagnosticDecoy, ...]:
    result_paths = _manifest_task_results(
        run_dir=run_dir, manifest=manifest, tasks=tasks
    )
    decoys: list[DiagnosticDecoy] = []
    seen_ids: set[str] = set()
    for task in tasks:
        result_path = result_paths[task.task_id]
        task_dir = result_path.parent
        document = read_document(result_path, name="diagnostic task result")
        if (
            document.get("schema") != TASK_SCHEMA
            or document.get("status") != "completed"
            or document.get("task_id") != task.task_id
            or document.get("start_id") != task.start_id
            or document.get("receptor_site_id") != task.receptor_site_id
            or document.get("source_seed") != task.source_seed
            or document.get("refinement_seed") != task.refinement_seed
            or document.get("diagnostic_plan_sha256") != plan_hash
        ):
            raise ValidationError(
                f"diagnostic task identity is invalid: {task.task_id}"
            )
        scorefile, _ = validate_record(
            root=task_dir,
            raw=document.get("scorefile"),
            name=f"{task.task_id} scorefile",
        )
        manifest_path, _ = validate_record(
            root=task_dir,
            raw=document.get("decoy_manifest"),
            name=f"{task.task_id} decoy manifest",
        )
        rows = _decoy_manifest(manifest_path)
        if len(rows) != expected_per_task:
            raise ValidationError(
                f"{task.task_id} has {len(rows)} decoys; expected {expected_per_task}"
            )
        scores = {
            row.description: row.score(ranking_score_name)
            for row in read_rosetta_scorefile(scorefile)
        }
        if len(scores) != expected_per_task:
            raise ValidationError(f"{task.task_id} score count is inconsistent")
        failures = _chemistry_failures(document, expected_count=expected_per_task)
        descriptions = {row["description"] for row in rows}
        if set(scores) != descriptions or not set(failures).issubset(descriptions):
            raise ValidationError(f"{task.task_id} decoy identities are inconsistent")
        for row in rows:
            description = safe_identifier(
                row["description"], name="diagnostic decoy description"
            )
            decoy_id = f"{task.task_id}__{description}"
            if decoy_id in seen_ids:
                raise ValidationError(f"duplicate diagnostic decoy ID: {decoy_id}")
            seen_ids.add(decoy_id)
            expected_hash = row["sha256"]
            path = (task_dir / row["path"]).resolve()
            try:
                ranking_score = float(row["ranking_score"])
            except (ValueError, OverflowError) as exc:
                raise ValidationError(
                    f"invalid diagnostic ranking score: {decoy_id}"
                ) from exc
            reference_score = scores.get(description)
            if (
                not path.is_relative_to(task_dir.resolve())
                or not path.is_file()
                or sha256_file(path) != expected_hash
                or reference_score is None
                or not math.isfinite(ranking_score)
                or row["ranking_score"] != f"{reference_score:.6f}"
            ):
                raise ValidationError(f"diagnostic decoy record mismatch: {decoy_id}")
            receptor_rmsd, backbone = _aligned_backbone(
                path=path,
                peptide_sequence=peptide_sequence,
                native_receptor=native_receptor,
            )
            decoys.append(
                DiagnosticDecoy(
                    decoy_id,
                    description,
                    task.task_id,
                    task.start_id,
                    task.receptor_site_id,
                    task.source_seed,
                    task.refinement_seed,
                    path,
                    expected_hash,
                    ranking_score,
                    description not in failures,
                    failures.get(description),
                    receptor_rmsd,
                    _rmsd(backbone, native_backbone),
                    backbone,
                )
            )
    return tuple(decoys)


def cluster_diagnostic_decoys(
    *,
    decoys: tuple[DiagnosticDecoy, ...],
    cluster_threshold_A: float,
    recovery_threshold_A: float,
    min_seed_support: int,
    min_start_support: int,
) -> tuple[DiagnosticCluster, ...]:
    """按分数顺序执行确定性的 complete-link 阈值聚类。"""
    grouped: dict[str, list[DiagnosticDecoy]] = defaultdict(list)
    for decoy in decoys:
        if decoy.chemistry_valid:
            grouped[decoy.receptor_site_id].append(decoy)
    results: list[DiagnosticCluster] = []
    for site_id, members in sorted(grouped.items()):
        groups: list[list[DiagnosticDecoy]] = []
        for decoy in sorted(
            members, key=lambda item: (item.ranking_score, item.decoy_id)
        ):
            for group in groups:
                if all(
                    _rmsd(decoy.aligned_backbone, member.aligned_backbone)
                    <= cluster_threshold_A
                    for member in group
                ):
                    group.append(decoy)
                    break
            else:
                groups.append([decoy])
        for rank, group in enumerate(groups, 1):
            cluster_members = tuple(group)
            recovered = tuple(
                item
                for item in cluster_members
                if item.native_backbone_rmsd_A <= recovery_threshold_A
            )
            medoid = min(
                cluster_members,
                key=lambda candidate: (
                    sum(
                        _rmsd(candidate.aligned_backbone, other.aligned_backbone)
                        for other in cluster_members
                    ),
                    candidate.decoy_id,
                ),
            )
            recovered_seeds = tuple(
                sorted({item.refinement_seed for item in recovered})
            )
            recovered_starts = tuple(sorted({item.start_id for item in recovered}))
            results.append(
                DiagnosticCluster(
                    f"{site_id}__R{rank:03d}",
                    site_id,
                    rank,
                    cluster_members,
                    cluster_members[0].decoy_id,
                    medoid.decoy_id,
                    tuple(sorted({item.refinement_seed for item in cluster_members})),
                    tuple(sorted({item.start_id for item in cluster_members})),
                    recovered_seeds,
                    recovered_starts,
                    min(item.native_backbone_rmsd_A for item in cluster_members),
                    len(recovered_seeds) >= min_seed_support
                    and len(recovered_starts) >= min_start_support,
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


def _write_starts(
    *, path: Path, run_dir: Path, starts: tuple[StartAssessment, ...]
) -> None:
    fields = (
        "start_id",
        "receptor_site_id",
        "source_seed",
        "receptor_alignment_rmsd_A",
        "native_backbone_rmsd_A",
        "path",
        "sha256",
    )
    atomic_write_text(
        path,
        _tsv(
            fields,
            [
                {
                    "start_id": item.start_id,
                    "receptor_site_id": item.receptor_site_id,
                    "source_seed": str(item.source_seed),
                    "receptor_alignment_rmsd_A": (
                        f"{item.receptor_alignment_rmsd_A:.6f}"
                    ),
                    "native_backbone_rmsd_A": (f"{item.native_backbone_rmsd_A:.6f}"),
                    "path": item.path.relative_to(run_dir).as_posix(),
                    "sha256": item.sha256,
                }
                for item in starts
            ],
        ),
    )


def _write_decoys(
    *,
    path: Path,
    run_dir: Path,
    decoys: tuple[DiagnosticDecoy, ...],
    clusters: tuple[DiagnosticCluster, ...],
    recovery_threshold_A: float,
) -> None:
    cluster_rank = {
        member.decoy_id: cluster.rank
        for cluster in clusters
        for member in cluster.members
    }
    fields = (
        "decoy_id",
        "task_id",
        "start_id",
        "receptor_site_id",
        "source_seed",
        "refinement_seed",
        "ranking_score",
        "chemistry_valid",
        "chemistry_failure",
        "receptor_alignment_rmsd_A",
        "native_backbone_rmsd_A",
        "recovered",
        "cluster_rank",
        "path",
        "sha256",
    )
    atomic_write_text(
        path,
        _tsv(
            fields,
            [
                {
                    "decoy_id": item.decoy_id,
                    "task_id": item.task_id,
                    "start_id": item.start_id,
                    "receptor_site_id": item.receptor_site_id,
                    "source_seed": str(item.source_seed),
                    "refinement_seed": str(item.refinement_seed),
                    "ranking_score": f"{item.ranking_score:.6f}",
                    "chemistry_valid": str(item.chemistry_valid).lower(),
                    "chemistry_failure": item.chemistry_failure or "",
                    "receptor_alignment_rmsd_A": (
                        f"{item.receptor_alignment_rmsd_A:.6f}"
                    ),
                    "native_backbone_rmsd_A": (f"{item.native_backbone_rmsd_A:.6f}"),
                    "recovered": str(
                        item.chemistry_valid
                        and item.native_backbone_rmsd_A <= recovery_threshold_A
                    ).lower(),
                    "cluster_rank": str(cluster_rank.get(item.decoy_id, "")),
                    "path": item.path.relative_to(run_dir).as_posix(),
                    "sha256": item.sha256,
                }
                for item in decoys
            ],
        ),
    )


def _write_clusters(
    *, path: Path, clusters: tuple[DiagnosticCluster, ...], threshold_A: float
) -> None:
    fields = (
        "cluster_id",
        "receptor_site_id",
        "cluster_rank",
        "decoy_count",
        "refinement_seeds",
        "start_ids",
        "recovered_refinement_seeds",
        "recovered_start_ids",
        "best_ranking_score",
        "minimum_native_backbone_rmsd_A",
        "near_native",
        "supported",
        "energy_representative_decoy_id",
        "geometric_medoid_decoy_id",
    )
    atomic_write_text(
        path,
        _tsv(
            fields,
            [
                {
                    "cluster_id": item.cluster_id,
                    "receptor_site_id": item.receptor_site_id,
                    "cluster_rank": str(item.rank),
                    "decoy_count": str(len(item.members)),
                    "refinement_seeds": ";".join(map(str, item.refinement_seeds)),
                    "start_ids": ";".join(item.start_ids),
                    "recovered_refinement_seeds": ";".join(
                        map(str, item.recovered_seeds)
                    ),
                    "recovered_start_ids": ";".join(item.recovered_start_ids),
                    "best_ranking_score": (f"{item.members[0].ranking_score:.6f}"),
                    "minimum_native_backbone_rmsd_A": (
                        f"{item.minimum_native_backbone_rmsd_A:.6f}"
                    ),
                    "near_native": str(
                        item.minimum_native_backbone_rmsd_A <= threshold_A
                    ).lower(),
                    "supported": str(item.supported).lower(),
                    "energy_representative_decoy_id": (item.energy_representative_id),
                    "geometric_medoid_decoy_id": item.geometric_medoid_id,
                }
                for item in clusters
            ],
        ),
    )


def _source_exact_pose_seeds(
    *, config: AppConfig, plan: dict[str, object]
) -> tuple[int, ...]:
    try:
        inputs = object_mapping(plan.get("inputs"), name="diagnostic inputs")
    except TypeError as exc:
        raise ValidationError("diagnostic inputs are invalid") from exc
    outputs_root = config.paths.outputs_dir.resolve()
    report = inputs.get("discovery_qualification_report")
    try:
        record = object_mapping(report, name="discovery qualification report")
    except TypeError as exc:
        raise ValidationError(
            "discovery qualification report record is invalid"
        ) from exc
    relative = record.get("path")
    expected_hash = record.get("sha256")
    if not isinstance(relative, str) or not isinstance(expected_hash, str):
        raise ValidationError("discovery qualification report record is invalid")
    path = (outputs_root / relative).resolve()
    if (
        not path.is_relative_to(
            (outputs_root / "discovery" / "qualifications").resolve()
        )
        or not path.is_file()
        or sha256_file(path) != expected_hash
    ):
        raise ValidationError("discovery qualification report is missing or changed")
    document = read_document(path, name="discovery qualification report")
    try:
        recovery = object_mapping(
            document.get("control_recovery"), name="control recovery"
        )
        selection = object_mapping(
            recovery.get("candidate_selection"), name="candidate selection"
        )
        seeds = object_list(
            selection.get("exact_pose_successful_seeds"),
            name="exact-pose successful seeds",
        )
    except TypeError as exc:
        raise ValidationError("discovery exact-pose seed evidence is invalid") from exc
    return tuple(_integer(seed, name="exact-pose successful seed") for seed in seeds)


def _site_summaries(
    *,
    starts: tuple[StartAssessment, ...],
    decoys: tuple[DiagnosticDecoy, ...],
    clusters: tuple[DiagnosticCluster, ...],
    recovery_threshold_A: float,
) -> list[dict[str, JsonValue]]:
    site_ids = sorted({item.receptor_site_id for item in starts})
    summaries: list[dict[str, JsonValue]] = []
    for site_id in site_ids:
        site_starts = tuple(item for item in starts if item.receptor_site_id == site_id)
        site_decoys = tuple(item for item in decoys if item.receptor_site_id == site_id)
        valid = tuple(item for item in site_decoys if item.chemistry_valid)
        recovered = tuple(
            item
            for item in valid
            if item.native_backbone_rmsd_A <= recovery_threshold_A
        )
        site_clusters = tuple(
            item for item in clusters if item.receptor_site_id == site_id
        )
        ranked = sorted(valid, key=lambda item: (item.ranking_score, item.decoy_id))
        first_recovered_rank = next(
            (
                rank
                for rank, item in enumerate(ranked, 1)
                if item.native_backbone_rmsd_A <= recovery_threshold_A
            ),
            None,
        )
        supported_cluster = next(
            (item for item in site_clusters if item.supported), None
        )
        best_start = min(item.native_backbone_rmsd_A for item in site_starts)
        best_decoy = min(item.native_backbone_rmsd_A for item in valid)
        summaries.append(
            {
                "receptor_site_id": site_id,
                "start_count": len(site_starts),
                "source_seeds": sorted({item.source_seed for item in site_starts}),
                "best_start_native_backbone_rmsd_A": round(best_start, 6),
                "median_start_native_backbone_rmsd_A": round(
                    median(item.native_backbone_rmsd_A for item in site_starts), 6
                ),
                "valid_decoy_count": len(valid),
                "chemistry_invalid_decoy_count": len(site_decoys) - len(valid),
                "best_decoy_native_backbone_rmsd_A": round(best_decoy, 6),
                "median_decoy_native_backbone_rmsd_A": round(
                    median(item.native_backbone_rmsd_A for item in valid), 6
                ),
                "best_rmsd_improvement_A": round(best_start - best_decoy, 6),
                "recovered_decoy_count": len(recovered),
                "recovered_refinement_seeds": sorted(
                    {item.refinement_seed for item in recovered}
                ),
                "recovered_start_ids": sorted({item.start_id for item in recovered}),
                "cluster_count": len(site_clusters),
                "supported_near_native_cluster_count": sum(
                    item.supported for item in site_clusters
                ),
                "first_recovered_score_rank": first_recovered_rank,
                "first_supported_cluster_rank": (
                    supported_cluster.rank if supported_cluster is not None else None
                ),
                "sampling_success": bool(recovered),
                "repeatable_recovery_success": supported_cluster is not None,
            }
        )
    return summaries


def analyze_qualification_refinement(
    *, config: AppConfig, run_dir: Path
) -> QualificationRefinementAnalysisOutcome:
    """验证完整诊断并写出不覆盖原始运行的 native recovery 报告。"""
    output_dir = run_dir / ANALYSIS_DIRECTORY
    if output_dir.exists():
        raise ValidationError(
            f"qualification refinement analysis already exists: {output_dir}"
        )
    plan, manifest, tasks, native_path, plan_hash = _validated_run(
        config=config, run_dir=run_dir
    )
    parameters = _parameters(plan)
    receptor_backbone = _receptor_backbone_parameters(parameters)
    expected_per_task = _integer(
        parameters.get("decoys_per_seed"), name="decoys per seed", minimum=1
    )
    ranking_score_name = _text(parameters.get("ranking_score"), name="ranking score")
    recovery_threshold = _number(
        parameters.get("max_native_backbone_rmsd_A"),
        name="native recovery RMSD",
        positive=True,
    )
    cluster_threshold = _number(
        parameters.get("max_cluster_backbone_rmsd_A"),
        name="cluster RMSD",
        positive=True,
    )
    min_seed_support = _integer(
        parameters.get("min_refinement_seed_support"),
        name="minimum refinement seed support",
        minimum=1,
    )
    min_start_support = _integer(
        parameters.get("min_source_start_support"),
        name="minimum source start support",
        minimum=1,
    )
    try:
        chemistry = object_mapping(plan.get("chemistry"), name="diagnostic chemistry")
    except TypeError as exc:
        raise ValidationError("diagnostic chemistry is invalid") from exc
    peptide_sequence = _text(chemistry.get("sequence"), name="peptide sequence")
    native_receptor, native_backbone = _native_coordinates(
        path=native_path, peptide_sequence=peptide_sequence
    )
    starts = _starts(
        run_dir=run_dir,
        tasks=tasks,
        plan_hash=plan_hash,
        peptide_sequence=peptide_sequence,
        native_receptor=native_receptor,
        native_backbone=native_backbone,
    )
    decoys = _decoys(
        run_dir=run_dir,
        manifest=manifest,
        tasks=tasks,
        plan_hash=plan_hash,
        expected_per_task=expected_per_task,
        ranking_score_name=ranking_score_name,
        peptide_sequence=peptide_sequence,
        native_receptor=native_receptor,
        native_backbone=native_backbone,
    )
    clusters = cluster_diagnostic_decoys(
        decoys=decoys,
        cluster_threshold_A=cluster_threshold,
        recovery_threshold_A=recovery_threshold,
        min_seed_support=min_seed_support,
        min_start_support=min_start_support,
    )
    starts_path = output_dir / "starts.tsv"
    decoys_path = output_dir / "decoys.tsv"
    clusters_path = output_dir / "clusters.tsv"
    report_path = output_dir / REPORT_NAME
    _write_starts(path=starts_path, run_dir=run_dir, starts=starts)
    _write_decoys(
        path=decoys_path,
        run_dir=run_dir,
        decoys=decoys,
        clusters=clusters,
        recovery_threshold_A=recovery_threshold,
    )
    _write_clusters(
        path=clusters_path,
        clusters=clusters,
        threshold_A=recovery_threshold,
    )
    site_summaries = _site_summaries(
        starts=starts,
        decoys=decoys,
        clusters=clusters,
        recovery_threshold_A=recovery_threshold,
    )
    valid = tuple(item for item in decoys if item.chemistry_valid)
    supported = tuple(item for item in clusters if item.supported)
    source_exact_seeds = _source_exact_pose_seeds(config=config, plan=plan)
    delivered_source_seeds = tuple(sorted({item.source_seed for item in starts}))
    delivered_exact_seeds = tuple(
        sorted(set(source_exact_seeds) & set(delivered_source_seeds))
    )
    try:
        source_handoff_qc = object_mapping(
            plan.get("source_handoff_qc"), name="source handoff QC"
        )
    except TypeError as exc:
        raise ValidationError("source handoff QC is invalid") from exc
    source_task_count = _integer(
        source_handoff_qc.get("task_count"), name="source handoff task count", minimum=1
    )
    source_passed_count = _integer(
        source_handoff_qc.get("passed_task_count"),
        name="source handoff passed task count",
    )
    source_failed_count = _integer(
        source_handoff_qc.get("failed_task_count"),
        name="source handoff failed task count",
    )
    source_invalid_count = _integer(
        source_handoff_qc.get("invalid_task_count"),
        name="source handoff invalid task count",
    )
    source_passed_fraction = _number(
        source_handoff_qc.get("passed_fraction"),
        name="source handoff passed fraction",
    )
    if (
        source_passed_count + source_failed_count + source_invalid_count
        != source_task_count
        or source_invalid_count
        or not 0.0 <= source_passed_fraction <= 1.0
        or not math.isclose(
            source_passed_fraction,
            round(source_passed_count / source_task_count, 6),
        )
    ):
        raise ValidationError("source handoff QC is inconsistent")
    source_handoff_qc_record: dict[str, JsonValue] = {
        "task_count": source_task_count,
        "passed_task_count": source_passed_count,
        "failed_task_count": source_failed_count,
        "invalid_task_count": source_invalid_count,
        "passed_fraction": source_passed_fraction,
        "scientific_qc_failures_excluded_before_refinement": True,
    }
    recovery_supported = bool(supported)
    failure_reasons: list[str] = []
    if not any(item.native_backbone_rmsd_A <= recovery_threshold for item in valid):
        failure_reasons.append("no_near_native_refined_decoy")
    if not supported:
        failure_reasons.append("no_repeatable_near_native_cluster")
    if source_exact_seeds and not delivered_exact_seeds:
        failure_reasons.append("source_exact_pose_seed_not_delivered")
    atomic_write_json(
        report_path,
        {
            "schema": REPORT_SCHEMA,
            "stage": "validation_qualification_refinement_diagnostic_analysis",
            "status": "completed",
            "generated_at": utc_now(),
            "analysis_software": vela_software_identity(),
            "development_only": True,
            "formal_qualification_gate": False,
            "native_information_used": True,
            "evidence_category": EVIDENCE_CATEGORY,
            "scientific_result": (
                "recovery_supported" if recovery_supported else "recovery_not_supported"
            ),
            "failure_reasons": failure_reasons,
            "parameters": {
                "ranking_score": ranking_score_name,
                "max_native_backbone_rmsd_A": recovery_threshold,
                "max_cluster_backbone_rmsd_A": cluster_threshold,
                "min_refinement_seed_support": min_seed_support,
                "min_source_start_support": min_start_support,
                "cluster_ranking": (
                    "best_ranking_score_asc_then_decoy_id; no frozen Top-K gate"
                ),
                "receptor_backbone": receptor_backbone,
            },
            "source": {
                "qualification_refinement_plan": {
                    "path": f"../{PLAN_NAME}",
                    "sha256": plan_hash,
                },
                "qualification_refinement_manifest": {
                    "path": f"../{MANIFEST_NAME}",
                    "sha256": sha256_file(run_dir / MANIFEST_NAME),
                },
                "native_reference": {
                    "path": native_path.relative_to(config.paths.data_dir).as_posix(),
                    "sha256": sha256_file(native_path),
                },
            },
            "data_quality": {
                "source_handoff": source_handoff_qc_record,
                "task_count": len(tasks),
                "expected_decoy_count": len(tasks) * expected_per_task,
                "observed_decoy_count": len(decoys),
                "chemistry_valid_decoy_count": len(valid),
                "chemistry_invalid_decoy_count": len(decoys) - len(valid),
                "chemistry_valid_fraction": round(len(valid) / len(decoys), 6),
                "invalid_decoys_excluded_from_ranking_and_clustering": True,
            },
            "handoff_diagnostic": {
                "source_exact_pose_successful_seeds": list(source_exact_seeds),
                "delivered_source_seeds": list(delivered_source_seeds),
                "delivered_exact_pose_successful_seeds": list(delivered_exact_seeds),
                "source_exact_pose_seed_delivered": bool(delivered_exact_seeds),
                "interpretation": (
                    "The frozen native-free handoff delivered "
                    "at least one start from a strict CABS recovery seed."
                    if delivered_exact_seeds
                    else "The frozen native-free handoff did not "
                    "deliver a start from any strict CABS recovery seed."
                ),
            },
            "site_assessments": site_summaries,
            "overall": {
                "near_native_decoy_count": sum(
                    item.native_backbone_rmsd_A <= recovery_threshold for item in valid
                ),
                "supported_near_native_cluster_count": len(supported),
                "recovery_supported": recovery_supported,
                "ranking_gate_evaluable": False,
                "ranking_gate_limitation": (
                    "The frozen development plan defines no final Top-K ranking "
                    "window; this report records ranks descriptively."
                ),
                "recommended_next_action": (
                    "freeze_joint_contract_then_run_unseen_holdout"
                    if recovery_supported
                    else "no_required_action_cross_receptor_robustness_optional"
                ),
            },
            "artifacts": {
                "starts": {
                    "path": starts_path.name,
                    "sha256": sha256_file(starts_path),
                    "count": len(starts),
                },
                "decoys": {
                    "path": decoys_path.name,
                    "sha256": sha256_file(decoys_path),
                    "count": len(decoys),
                },
                "clusters": {
                    "path": clusters_path.name,
                    "sha256": sha256_file(clusters_path),
                    "count": len(clusters),
                },
            },
            "interpretation": (
                "This completed native-aware development diagnostic does not alter "
                "the frozen Stage 2 holdout decision and cannot qualify the joint "
                "workflow."
            ),
        },
    )
    return QualificationRefinementAnalysisOutcome(
        report_path,
        len(valid),
        len(clusters),
        recovery_supported,
    )
