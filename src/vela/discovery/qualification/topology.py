"""用开发 seed 独立校准 CABS 二硫环拓扑的全原子可恢复性。"""

from __future__ import annotations

import json
import math
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

from vela.config import AppConfig
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
from vela.core.typed_data import json_value, object_list, object_mapping
from vela.discovery.analysis.cluster_engine import (
    bounded_leader_clusters,
    normalized_site_distance,
)
from vela.discovery.models import DiscoveryError
from vela.discovery.qualification.control import control_bound_state, control_chemistry
from vela.discovery.qualification.planning import CONTROL_RECOVERY
from vela.discovery.sampling.cabsdock import cabsdock_archive_path
from vela.discovery.sampling.evidence import (
    align_receptor,
    ca_contact_residues,
    centroid,
    disulfide_ca_distance,
    disulfide_cabs_sc_distance,
    read_reference_chains,
    read_structure,
    required_atom,
    split_model,
)
from vela.discovery.sampling.trajectory import audit_cabs_trajectory
from vela.validation.models import ValidationError
from vela.validation.records import file_record, validate_record
from vela.validation.refinement.reconstruction import (
    assess_topology_reconstruction,
    build_cg2all_command,
    verify_cg2all_tool,
    write_cg2all_input,
    write_disulfide_indices,
    write_reference_receptor_complex,
    write_topology_rebuild_protocol,
)
from vela.validation.rosetta import (
    build_chemistry_command,
    build_prepack_command,
    build_topology_refine_command,
    run_rosetta_command,
    single_rosetta_pdb_output,
    verify_flexpepdock_tool,
    verify_rosetta_scripts_tool,
)
from vela.validation.scores import read_rosetta_scorefile

PLAN_NAME = "topology_calibration_plan.json"
MANIFEST_NAME = "topology_calibration_manifest.json"
REPORT_NAME = "topology_calibration_report.json"
PLAN_SCHEMA = "vela.disulfide-topology-calibration-plan/2"
TASK_SCHEMA = "vela.disulfide-topology-calibration-task-result/2"
MANIFEST_SCHEMA = "vela.disulfide-topology-calibration-manifest/2"
REPORT_SCHEMA = "vela.disulfide-topology-calibration-report/2"


@dataclass(frozen=True, slots=True)
class _Candidate:
    """一个带真实 CABS 能量和位点特征的 Top-1000 模型。"""

    task_id: str
    seed: int
    model_path: Path
    model_sha256: str
    model_index: int
    archive_path: Path
    archive_sha256: str
    ca_distance_A: float
    cabs_sc_distance_A: float
    interaction_energy: float
    contact_residues: frozenset[str]
    local_position: tuple[float, float, float]
    stratum_index: int

    @property
    def identity(self) -> str:
        return f"{self.task_id}__model_{self.model_index:04d}"


@dataclass(frozen=True, slots=True)
class TopologyCalibrationPlan:
    """已经写入磁盘的拓扑校准计划。"""

    run_dir: Path
    task_count: int


@dataclass(frozen=True, slots=True)
class _ControlSource:
    """一个开发控制任务及其不含配体的坐标参考。"""

    task_id: str
    seed: int
    reference_path: Path
    reference_sha256: str


def topology_calibration_contract(config: AppConfig) -> dict[str, JsonValue]:
    """返回资格计划和校准报告必须完全一致的方法合同。"""
    calibration = config.discovery.qualification.topology_calibration
    return {
        "coarse_grained_representation": "CABS_CA_SC",
        "topology_feasibility": {
            "metric": "maximum_disulfide_endpoint_CA_distance",
            "candidate_thresholds_A": list(calibration.candidate_ca_thresholds_A),
            "above_threshold_comparator_upper_bound_A": (
                calibration.above_threshold_comparator_upper_bound_A
            ),
            "threshold_selection": "highest_contiguous_passing_candidate",
            "chemical_bond_claimed": False,
        },
        "stratified_sampling": {
            "ca_strata_upper_bounds_A": list(calibration.strata_upper_bounds_A),
            "models_per_stratum": calibration.models_per_stratum,
            "coverage_dimensions": [
                "seed",
                "binding_site",
                "cabsdock_interaction_energy",
                "cabs_SC_pseudocenter_distance",
            ],
            "native_information_used": False,
        },
        "all_atom_reconstruction": {
            "pipeline": (
                "cg2all_peptide_then_aligned_experimental_receptor_graft_then_"
                "RosettaScripts_ForceDisulfides_repack_then_FlexPepDock_prepack_"
                "then_single_local_refine"
            ),
            "cg2all_representation": config.validation.cg2all.representation,
            "cg2all_checkpoint_sha256": (config.validation.cg2all.checkpoint_sha256),
            "max_receptor_CA_RMSD_A": config.validation.cg2all.max_ca_rmsd_A,
            "max_peptide_internal_CA_RMSD_A": (
                calibration.max_peptide_internal_ca_rmsd_A
            ),
            "max_ligand_centroid_displacement_A": (
                calibration.max_ligand_centroid_displacement_A
            ),
            "receptor_contact_representation": "receptor_CA_to_peptide_CA",
            "receptor_contact_distance_A": calibration.contact_ca_threshold_A,
            "min_receptor_contact_retention_fraction": (
                calibration.min_receptor_contact_retention_fraction
            ),
            "min_disulfide_SG_distance_A": (config.validation.min_disulfide_sg_A),
            "max_disulfide_SG_distance_A": (config.validation.max_disulfide_sg_A),
            "min_interchain_heavy_atom_distance_A": (
                calibration.min_interchain_heavy_atom_distance_A
            ),
            "max_refine_weighted_fa_rep_per_residue": (
                calibration.max_refine_fa_rep_per_residue
            ),
            "max_refine_weighted_backbone_strain_per_residue": (
                calibration.max_refine_backbone_strain_per_residue
            ),
            "terminal_chemistry_assessed": False,
        },
        "qualification_rule": {
            "eligible_threshold": (
                "all_strata_up_to_candidate_must_pass; comparator_is_descriptive"
            ),
            "min_success_fraction_per_stratum": (
                calibration.min_success_fraction_per_stratum
            ),
            "min_successful_seeds_per_stratum": (
                calibration.min_successful_seeds_per_stratum
            ),
            "above_threshold_comparator_is_qualification_gate": False,
        },
    }


def _document(path: Path, *, name: str) -> dict[str, object]:
    if not path.is_file():
        raise DiscoveryError(f"{name} does not exist: {path}")
    try:
        value: object = json.loads(path.read_text(encoding="utf-8"))
        return object_mapping(value, name=name)
    except (UnicodeDecodeError, json.JSONDecodeError, TypeError) as exc:
        raise DiscoveryError(f"invalid {name}: {path}") from exc


def _required_integer(source: dict[str, object], key: str, *, name: str) -> int:
    value = source.get(key)
    if not isinstance(value, int) or isinstance(value, bool):
        raise DiscoveryError(f"{name} must be an integer")
    return value


def _required_text(source: dict[str, object], key: str, *, name: str) -> str:
    value = source.get(key)
    if not isinstance(value, str) or not value:
        raise DiscoveryError(f"{name} must be a non-empty string")
    return value


def _calibration_root(config: AppConfig) -> Path:
    return (config.paths.outputs_dir / "discovery" / "topology_calibrations").resolve()


def resolve_topology_calibration_run(*, config: AppConfig, run_dir: Path) -> Path:
    """把运行目录限制在拓扑校准输出根。"""
    resolved = run_dir.expanduser().resolve()
    if not resolved.is_relative_to(_calibration_root(config)):
        raise DiscoveryError(
            f"topology calibration is outside its output root: {resolved}"
        )
    return resolved


def _source_run(*, config: AppConfig, source_run: Path) -> Path:
    resolved = source_run.expanduser().resolve()
    root = (config.paths.outputs_dir / "discovery" / "qualifications").resolve()
    if not resolved.is_relative_to(root):
        raise DiscoveryError(
            f"topology calibration source is outside qualifications: {resolved}"
        )
    return resolved


def _control_rows(
    *, source_plan: dict[str, object], state_id: str, data_dir: Path
) -> tuple[_ControlSource, ...]:
    if (
        source_plan.get("stage") != "discovery_qualification"
        or source_plan.get("schema") != "vela.discovery-qualification-plan/4"
    ):
        raise DiscoveryError(
            "topology calibration source is not a recognized development qualification"
        )
    try:
        rows = object_list(source_plan.get("tasks"), name="source qualification tasks")
    except TypeError as exc:
        raise DiscoveryError("source qualification tasks are invalid") from exc
    controls: list[_ControlSource] = []
    for raw in rows:
        try:
            row = object_mapping(raw, name="source qualification task")
        except TypeError as exc:
            raise DiscoveryError("source qualification task is invalid") from exc
        if row.get("case") != CONTROL_RECOVERY:
            continue
        task_id = row.get("task_id")
        seed = row.get("seed")
        if (
            not isinstance(task_id, str)
            or not isinstance(seed, int)
            or isinstance(seed, bool)
            or row.get("receptor_id") != state_id
        ):
            raise DiscoveryError("source control task identity is invalid")
        try:
            reference_path, reference_sha256 = validate_record(
                root=data_dir,
                raw=row.get("reference"),
                name=f"{task_id} receptor-only reference",
            )
        except ValidationError as exc:
            raise DiscoveryError(str(exc)) from exc
        controls.append(_ControlSource(task_id, seed, reference_path, reference_sha256))
    if not controls or len({item.seed for item in controls}) != len(controls):
        raise DiscoveryError(
            "source control tasks must contain unique development seeds"
        )
    reference_hashes = {item.reference_sha256 for item in controls}
    if len(reference_hashes) != 1:
        raise DiscoveryError(
            "source control tasks must use one receptor-only coordinate reference"
        )
    return tuple(sorted(controls, key=lambda item: item.seed))


def _stratum_index(distance_A: float, bounds: tuple[float, ...]) -> int | None:
    return next(
        (index for index, upper in enumerate(bounds) if distance_A <= upper), None
    )


def _filtered_candidates(
    *,
    config: AppConfig,
    source_run: Path,
    task_id: str,
    seed: int,
    reference_path: Path,
) -> tuple[_Candidate, ...]:
    state = control_bound_state(config)
    chemistry = control_chemistry(state)
    task_dir = source_run / "tasks" / task_id
    model_path = task_dir / "output_pdbs" / "top1000.pdb"
    archive_path = cabsdock_archive_path(task_dir)
    structure = read_structure(model_path)
    settings = config.discovery.cabsdock
    if len(structure) != settings.filtering_count:
        raise DiscoveryError(f"calibration source model count is invalid: {task_id}")
    reference_chains = read_reference_chains(reference_path)
    ca_distances: list[float] = []
    sc_distances: list[float] = []
    topology_by_replica = [0] * settings.replicas
    per_replica = settings.filtering_count // settings.replicas
    for model_index, model in enumerate(structure, 1):
        _, peptide = split_model(model, peptide_sequence=chemistry.sequence)
        ca_distance = disulfide_ca_distance(peptide=peptide, chemistry=chemistry)
        ca_distances.append(ca_distance)
        sc_distances.append(
            disulfide_cabs_sc_distance(peptide=peptide, chemistry=chemistry)
        )
        if ca_distance <= settings.max_disulfide_ca_distance_A:
            topology_by_replica[(model_index - 1) // per_replica] += 1
    try:
        audit = audit_cabs_trajectory(
            archive_path=archive_path,
            chemistry=chemistry,
            replicas=settings.replicas,
            filtering_count=settings.filtering_count,
            max_disulfide_ca_distance_A=settings.max_disulfide_ca_distance_A,
            filtered_topology_feasible_by_replica=tuple(topology_by_replica),
        )
    except DiscoveryError as exc:
        raise DiscoveryError(f"{task_id}: {exc}") from exc
    model_hash = sha256_file(model_path)
    archive_hash = sha256_file(archive_path)
    bounds = config.discovery.qualification.topology_calibration.strata_upper_bounds_A
    candidates: list[_Candidate] = []
    for model_index, model in enumerate(structure, 1):
        stratum = _stratum_index(ca_distances[model_index - 1], bounds)
        if stratum is None:
            continue
        receptor, peptide = split_model(model, peptide_sequence=chemistry.sequence)
        contacts = ca_contact_residues(
            receptor=receptor,
            peptide=peptide,
            threshold_A=(
                config.discovery.qualification.topology_calibration.contact_ca_threshold_A
            ),
        )
        if not contacts:
            continue
        alignment = align_receptor(
            receptor=receptor,
            reference_chains=reference_chains,
        )
        peptide_positions = tuple(
            alignment.transform.apply(required_atom(residue, "CA").pos)
            for residue in peptide
        )
        local_position = centroid(
            tuple(
                (position.x, position.y, position.z) for position in peptide_positions
            )
        )
        candidates.append(
            _Candidate(
                task_id=task_id,
                seed=seed,
                model_path=model_path,
                model_sha256=model_hash,
                model_index=model_index,
                archive_path=archive_path,
                archive_sha256=archive_hash,
                ca_distance_A=ca_distances[model_index - 1],
                cabs_sc_distance_A=sc_distances[model_index - 1],
                interaction_energy=audit.filtered_interaction_energies[model_index - 1],
                contact_residues=contacts,
                local_position=local_position,
                stratum_index=stratum,
            )
        )
    return tuple(candidates)


def _percentile_ranks(
    candidates: tuple[_Candidate, ...], *, attribute: str
) -> dict[str, float]:
    ordered = sorted(
        candidates,
        key=lambda candidate: (
            float(getattr(candidate, attribute)),
            candidate.identity,
        ),
    )
    denominator = max(1, len(ordered) - 1)
    return {
        candidate.identity: index / denominator
        for index, candidate in enumerate(ordered)
    }


def _site_distance(first: _Candidate, second: _Candidate, config: AppConfig) -> float:
    selection = config.discovery.cabsdock
    return normalized_site_distance(
        first_contacts=first.contact_residues,
        first_position=first.local_position,
        second_contacts=second.contact_residues,
        second_position=second.local_position,
        contact_limit=selection.selection_contact_jaccard_distance,
        position_limit=selection.selection_position_distance_A,
    )


def _select_stratum(
    *, candidates: tuple[_Candidate, ...], config: AppConfig
) -> tuple[tuple[_Candidate, int], ...]:
    budget = config.discovery.qualification.topology_calibration.models_per_stratum
    if len(candidates) < budget:
        raise DiscoveryError("topology calibration stratum lacks enough models")
    energy_ranks = _percentile_ranks(candidates, attribute="interaction_energy")
    sc_ranks = _percentile_ranks(candidates, attribute="cabs_sc_distance_A")
    selected: list[_Candidate] = []
    remaining = set(candidates)
    while len(selected) < budget:
        used_seeds = frozenset(candidate.seed for candidate in selected)

        def key(
            candidate: _Candidate,
            *,
            used_seed_ids: frozenset[int] = used_seeds,
        ) -> tuple[float, float, float, float, float, int, str]:
            site_distance = min(
                (_site_distance(candidate, existing, config) for existing in selected),
                default=math.inf,
            )
            feature_distance = min(
                (
                    math.dist(
                        (
                            energy_ranks[candidate.identity],
                            sc_ranks[candidate.identity],
                        ),
                        (
                            energy_ranks[existing.identity],
                            sc_ranks[existing.identity],
                        ),
                    )
                    for existing in selected
                ),
                default=math.inf,
            )
            return (
                float(candidate.seed not in used_seed_ids),
                float(site_distance > 1.0),
                min(site_distance, 2.0),
                min(feature_distance, 2.0),
                -candidate.ca_distance_A,
                -candidate.seed,
                candidate.identity,
            )

        chosen = max(remaining, key=key)
        selected.append(chosen)
        remaining.remove(chosen)
    clusters = bounded_leader_clusters(
        tuple(selected),
        distance=lambda first, second: _site_distance(first, second, config),
        identity=lambda candidate: candidate.identity,
        maximum_distance=1.0,
    )
    site_by_identity = {
        candidate.identity: index
        for index, cluster in enumerate(clusters, 1)
        for candidate in cluster
    }
    return tuple(
        (candidate, site_by_identity[candidate.identity]) for candidate in selected
    )


def write_topology_calibration_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> TopologyCalibrationPlan:
    """从开发 seed 冻结分层、native-free 的全原子重建样本。"""
    validate_run_id(run_id)
    source = _source_run(config=config, source_run=source_run)
    source_plan_path = source / "qualification_plan.json"
    source_plan = _document(source_plan_path, name="source qualification plan")
    source_schema = source_plan.get("schema")
    if not isinstance(source_schema, str):
        raise DiscoveryError("source qualification plan schema is invalid")
    state = control_bound_state(config)
    controls = _control_rows(
        source_plan=source_plan,
        state_id=state.state_id,
        data_dir=config.paths.data_dir,
    )
    calibration = config.discovery.qualification.topology_calibration
    if len(controls) < calibration.min_successful_seeds_per_stratum:
        raise DiscoveryError(
            "topology calibration source has too few development seeds"
        )
    candidates = tuple(
        candidate
        for control in controls
        for candidate in _filtered_candidates(
            config=config,
            source_run=source,
            task_id=control.task_id,
            seed=control.seed,
            reference_path=control.reference_path,
        )
    )
    selected: list[tuple[_Candidate, int]] = []
    for stratum_index in range(len(calibration.strata_upper_bounds_A)):
        pool = tuple(
            candidate
            for candidate in candidates
            if candidate.stratum_index == stratum_index
        )
        selected.extend(_select_stratum(candidates=pool, config=config))
    run_dir = _calibration_root(config) / run_id
    if run_dir.exists():
        raise DiscoveryError(
            f"topology calibration run directory already exists: {run_dir}"
        )
    snapshot_path = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot_path, config.source_snapshot_text)
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    rosetta_scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    bounds = calibration.strata_upper_bounds_A
    task_rows: list[dict[str, JsonValue]] = []
    for ordinal, (candidate, site_group) in enumerate(selected):
        lower = bounds[candidate.stratum_index - 1] if candidate.stratum_index else None
        upper = bounds[candidate.stratum_index]
        task_rows.append(
            {
                "task_id": f"stratum_{candidate.stratum_index + 1:02d}__{ordinal + 1:03d}",
                "source_task_id": candidate.task_id,
                "source_seed": candidate.seed,
                "source_model": {
                    "path": candidate.model_path.relative_to(source).as_posix(),
                    "sha256": candidate.model_sha256,
                    "model_index": candidate.model_index,
                },
                "source_archive": {
                    "path": candidate.archive_path.relative_to(source).as_posix(),
                    "sha256": candidate.archive_sha256,
                },
                "stratum": {
                    "index": candidate.stratum_index + 1,
                    "lower_exclusive_A": lower,
                    "upper_inclusive_A": upper,
                    "threshold_candidate": (
                        upper in calibration.candidate_ca_thresholds_A
                    ),
                    "above_threshold_comparator": (
                        upper == calibration.above_threshold_comparator_upper_bound_A
                    ),
                },
                "selection_features": {
                    "disulfide_ca_distance_A": round(candidate.ca_distance_A, 6),
                    "disulfide_cabs_sc_distance_A": round(
                        candidate.cabs_sc_distance_A, 6
                    ),
                    "cabsdock_interaction_energy": candidate.interaction_energy,
                    "site_group": site_group,
                    "contact_residues": sorted(candidate.contact_residues),
                    "aligned_peptide_ca_centroid": [
                        round(value, 6) for value in candidate.local_position
                    ],
                },
                "rosetta_seed": config.validation.handoff.chemistry_seed + ordinal,
            }
        )
    source_record: dict[str, JsonValue] = {
        "run_path": source.relative_to(config.paths.outputs_dir).as_posix(),
        "qualification_plan": {
            "path": source_plan_path.name,
            "sha256": sha256_file(source_plan_path),
            "schema": source_schema,
        },
        "control_bound_state_id": state.state_id,
        "development_seeds": [control.seed for control in controls],
        "receptor_only_reference": {
            "path": controls[0]
            .reference_path.relative_to(config.paths.data_dir)
            .as_posix(),
            "sha256": controls[0].reference_sha256,
        },
    }
    tools_record: dict[str, JsonValue] = {
        "cg2all": {
            "version": cg2all.version,
            "executable_sha256": cg2all.executable_sha256,
            "checkpoint_sha256": cg2all.checkpoint_sha256,
        },
        "flexpepdock": {
            "version": flexpepdock.version,
            "executable_sha256": flexpepdock.executable_sha256,
        },
        "rosetta_scripts": {
            "version": rosetta_scripts.version,
            "executable_sha256": rosetta_scripts.executable_sha256,
        },
    }
    plan_document: dict[str, JsonValue] = {
        "schema": PLAN_SCHEMA,
        "stage": "disulfide_topology_calibration",
        "status": "planned",
        "run_id": run_id,
        "planned_at": utc_now(),
        "software": vela_software_identity(),
        "evidence_scope": "development_calibration_only",
        "native_information_used": False,
        "config_snapshot": file_record(snapshot_path, root=run_dir),
        "source": source_record,
        "contract": topology_calibration_contract(config),
        "tools": tools_record,
        "task_count": len(task_rows),
        "tasks": task_rows,
    }
    atomic_write_json(run_dir / PLAN_NAME, plan_document)
    return TopologyCalibrationPlan(run_dir, len(task_rows))


def _result_artifacts(task_dir: Path, result: dict[str, object]) -> None:
    try:
        artifacts = object_mapping(result.get("artifacts"), name="task artifacts")
    except TypeError as exc:
        raise DiscoveryError("topology calibration task artifacts are invalid") from exc
    for name, record in artifacts.items():
        try:
            validate_record(root=task_dir, raw=record, name=name)
        except ValidationError as exc:
            raise DiscoveryError(str(exc)) from exc


def _resume_task(*, task_dir: Path, task_id: str, plan_sha256: str) -> Path | None:
    path = task_dir / "task_result.json"
    if not path.is_file():
        return None
    result = _document(path, name="topology calibration task result")
    if (
        result.get("schema") != TASK_SCHEMA
        or result.get("status") != "completed"
        or result.get("task_id") != task_id
        or result.get("topology_calibration_plan_sha256") != plan_sha256
    ):
        raise DiscoveryError(f"stale topology calibration task result: {task_id}")
    _result_artifacts(task_dir, result)
    return path


def _write_task_result(
    *,
    task_dir: Path,
    task_id: str,
    plan_sha256: str,
    started_at: str,
    reconstruction_status: str,
    failure_reasons: tuple[str, ...],
    failure_detail: str | None,
    metrics: dict[str, JsonValue] | None,
    commands: dict[str, JsonValue],
    artifacts: dict[str, JsonValue],
) -> Path:
    result_path = task_dir / "task_result.json"
    atomic_write_json(
        result_path,
        {
            "schema": TASK_SCHEMA,
            "status": "completed",
            "task_id": task_id,
            "topology_calibration_plan_sha256": plan_sha256,
            "started_at": started_at,
            "completed_at": utc_now(),
            "reconstruction_status": reconstruction_status,
            "failure_reasons": list(failure_reasons),
            "failure_detail": failure_detail,
            "metrics": metrics,
            "commands": commands,
            "artifacts": artifacts,
        },
    )
    return result_path


def _run_task(
    *,
    config: AppConfig,
    run_dir: Path,
    source_run: Path,
    reference_receptor_path: Path,
    raw_task: object,
    plan_sha256: str,
) -> Path:
    try:
        task = object_mapping(raw_task, name="topology calibration task")
        source_model = object_mapping(
            task.get("source_model"), name="topology calibration source model"
        )
    except TypeError as exc:
        raise DiscoveryError("topology calibration task is invalid") from exc
    task_id = task.get("task_id")
    model_index = source_model.get("model_index")
    rosetta_seed = task.get("rosetta_seed")
    if (
        not isinstance(task_id, str)
        or not isinstance(model_index, int)
        or isinstance(model_index, bool)
        or not isinstance(rosetta_seed, int)
        or isinstance(rosetta_seed, bool)
    ):
        raise DiscoveryError("topology calibration task identity is invalid")
    task_dir = run_dir / "tasks" / task_id
    resumed = _resume_task(task_dir=task_dir, task_id=task_id, plan_sha256=plan_sha256)
    if resumed is not None:
        return resumed
    if task_dir.exists():
        raise DiscoveryError(
            f"incomplete topology calibration task requires review: {task_dir}"
        )
    task_dir.mkdir(parents=True)
    source_path, _ = validate_record(
        root=source_run, raw=source_model, name="topology calibration source model"
    )
    state = control_bound_state(config)
    chemistry = control_chemistry(state)
    started_at = utc_now()
    artifacts: dict[str, JsonValue] = {}
    commands: dict[str, JsonValue] = {}
    cg_input_path = task_dir / "cg2all_input.pdb"
    cg_input = write_cg2all_input(
        source_path=source_path,
        model_index=model_index,
        destination=cg_input_path,
        chemistry=chemistry,
        settings=config.validation.cg2all,
    )
    artifacts["cg2all_input"] = file_record(cg_input_path, root=task_dir)
    cg_output_path = task_dir / "cg2all_raw.pdb"
    cg_command = build_cg2all_command(
        settings=config.validation.cg2all,
        input_path=cg_input_path,
        output_path=cg_output_path,
    )
    commands["cg2all"] = list(cg_command)
    cg_log_path = task_dir / "cg2all.log"
    try:
        run_rosetta_command(
            command=cg_command,
            log_path=cg_log_path,
            thread_count=config.validation.cg2all.processes,
        )
    except ValidationError as exc:
        artifacts["cg2all_log"] = file_record(cg_log_path, root=task_dir)
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("cg2all_reconstruction_failed",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    artifacts["cg2all_log"] = file_record(cg_log_path, root=task_dir)
    artifacts["cg2all_raw"] = file_record(cg_output_path, root=task_dir)
    grafted_path = task_dir / "all_atom_reference_receptor.pdb"
    try:
        write_reference_receptor_complex(
            coarse_pose_path=cg_input_path,
            reconstructed_path=cg_output_path,
            reference_receptor_path=reference_receptor_path,
            destination=grafted_path,
            chemistry=chemistry,
        )
    except ValidationError as exc:
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("cg2all_output_invalid",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    artifacts["all_atom_reference_receptor"] = file_record(grafted_path, root=task_dir)
    disulfide_path = task_dir / "fix_disulfide.txt"
    write_disulfide_indices(
        destination=disulfide_path,
        receptor_residue_count=cg_input.receptor_residue_count,
        chemistry=chemistry,
    )
    artifacts["fix_disulfide"] = file_record(disulfide_path, root=task_dir)
    topology_protocol_path = task_dir / "rebuild_topology.xml"
    write_topology_rebuild_protocol(
        destination=topology_protocol_path,
        receptor_residue_count=cg_input.receptor_residue_count,
        chemistry=chemistry,
        score_function=config.validation.rosetta.score_function,
    )
    artifacts["topology_rebuild_protocol"] = file_record(
        topology_protocol_path, root=task_dir
    )
    rebuild_dir = task_dir / "topology_rebuild"
    rebuild_dir.mkdir()
    rebuild_command = build_chemistry_command(
        settings=config.validation.rosetta,
        input_path=grafted_path,
        protocol_path=topology_protocol_path,
        disulfide_path=disulfide_path,
        output_dir=rebuild_dir,
        seed=rosetta_seed,
    )
    commands["disulfide_rebuild"] = list(rebuild_command)
    rebuild_log_path = task_dir / "disulfide_rebuild.log"
    try:
        run_rosetta_command(command=rebuild_command, log_path=rebuild_log_path)
        rebuild_scores = read_rosetta_scorefile(rebuild_dir / "chemistry.sc")
        if len(rebuild_scores) != 1:
            raise ValidationError(
                "disulfide rebuild must produce exactly one score row"
            )
        rebuilt_output = single_rosetta_pdb_output(rebuild_dir)
    except ValidationError as exc:
        artifacts["disulfide_rebuild_log"] = file_record(
            rebuild_log_path, root=task_dir
        )
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("rosetta_disulfide_rebuild_failed",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    artifacts["disulfide_rebuild_log"] = file_record(rebuild_log_path, root=task_dir)
    rebuilt_path = task_dir / "all_atom_disulfide_rebuilt.pdb"
    atomic_write_text(rebuilt_path, rebuilt_output.read_text(encoding="utf-8"))
    artifacts["all_atom_disulfide_rebuilt"] = file_record(rebuilt_path, root=task_dir)
    prepack_dir = task_dir / "prepack"
    prepack_dir.mkdir()
    prepack_command = build_prepack_command(
        settings=config.validation.rosetta,
        input_path=rebuilt_path,
        disulfide_path=disulfide_path,
        output_dir=prepack_dir,
        seed=rosetta_seed,
        fixed_histidine_pose_indices=cg_input.fixed_histidine_pose_indices,
    )
    commands["disulfide_prepack"] = list(prepack_command)
    prepack_log_path = task_dir / "disulfide_prepack.log"
    try:
        run_rosetta_command(command=prepack_command, log_path=prepack_log_path)
        prepack_scores = read_rosetta_scorefile(prepack_dir / "prepack.sc")
        if len(prepack_scores) != 1:
            raise ValidationError(
                "disulfide prepack must produce exactly one score row"
            )
        rosetta_output = single_rosetta_pdb_output(prepack_dir)
    except ValidationError as exc:
        artifacts["disulfide_prepack_log"] = file_record(
            prepack_log_path, root=task_dir
        )
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("rosetta_disulfide_prepack_failed",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    artifacts["disulfide_prepack_log"] = file_record(prepack_log_path, root=task_dir)
    prepacked_path = task_dir / "all_atom_prepacked.pdb"
    atomic_write_text(prepacked_path, rosetta_output.read_text(encoding="utf-8"))
    artifacts["all_atom_prepacked"] = file_record(prepacked_path, root=task_dir)
    refine_dir = task_dir / "topology_refine"
    refine_dir.mkdir()
    refine_command = build_topology_refine_command(
        settings=config.validation.rosetta,
        input_path=prepacked_path,
        disulfide_path=disulfide_path,
        output_dir=refine_dir,
        seed=rosetta_seed,
        fixed_histidine_pose_indices=cg_input.fixed_histidine_pose_indices,
    )
    commands["topology_refine"] = list(refine_command)
    refine_log_path = task_dir / "topology_refine.log"
    try:
        run_rosetta_command(command=refine_command, log_path=refine_log_path)
        refine_scores = read_rosetta_scorefile(refine_dir / "refine.sc")
        if len(refine_scores) != 1:
            raise ValidationError("topology refine must produce exactly one score row")
        refined_output = single_rosetta_pdb_output(refine_dir)
    except ValidationError as exc:
        artifacts["topology_refine_log"] = file_record(refine_log_path, root=task_dir)
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("flexpepdock_topology_refine_failed",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    artifacts["topology_refine_log"] = file_record(refine_log_path, root=task_dir)
    final_path = task_dir / "all_atom_topology.pdb"
    atomic_write_text(final_path, refined_output.read_text(encoding="utf-8"))
    artifacts["all_atom_topology"] = file_record(final_path, root=task_dir)
    try:
        assessment = assess_topology_reconstruction(
            input_path=cg_input_path,
            output_path=final_path,
            chemistry=chemistry,
            settings=config.validation.cg2all,
            min_disulfide_sg_A=config.validation.min_disulfide_sg_A,
            max_disulfide_sg_A=config.validation.max_disulfide_sg_A,
            min_interchain_heavy_atom_distance_A=(
                config.discovery.qualification.topology_calibration.min_interchain_heavy_atom_distance_A
            ),
            contact_ca_threshold_A=(
                config.discovery.qualification.topology_calibration.contact_ca_threshold_A
            ),
            max_peptide_internal_ca_rmsd_A=(
                config.discovery.qualification.topology_calibration.max_peptide_internal_ca_rmsd_A
            ),
            max_ligand_centroid_displacement_A=(
                config.discovery.qualification.topology_calibration.max_ligand_centroid_displacement_A
            ),
            min_receptor_contact_retention_fraction=(
                config.discovery.qualification.topology_calibration.min_receptor_contact_retention_fraction
            ),
        )
    except ValidationError as exc:
        return _write_task_result(
            task_dir=task_dir,
            task_id=task_id,
            plan_sha256=plan_sha256,
            started_at=started_at,
            reconstruction_status="failed",
            failure_reasons=("all_atom_output_invalid",),
            failure_detail=str(exc),
            metrics=None,
            commands=commands,
            artifacts=artifacts,
        )
    residue_count = cg_input.receptor_residue_count + cg_input.peptide_residue_count
    refine_fa_rep = refine_scores[0].score("fa_rep")
    refine_backbone_strain = max(refine_scores[0].score("omega"), 0.0) + max(
        refine_scores[0].score("rama_prepro"), 0.0
    )
    refine_fa_rep_per_residue = refine_fa_rep / residue_count
    refine_backbone_strain_per_residue = refine_backbone_strain / residue_count
    calibration = config.discovery.qualification.topology_calibration
    score_failures: list[str] = []
    if refine_fa_rep_per_residue > calibration.max_refine_fa_rep_per_residue:
        score_failures.append("refined_fa_rep_exceeds_normalized_limit")
    if (
        refine_backbone_strain_per_residue
        > calibration.max_refine_backbone_strain_per_residue
    ):
        score_failures.append("refined_backbone_strain_exceeds_normalized_limit")
    failures = tuple((*assessment.failures, *score_failures))
    metrics: dict[str, JsonValue] = {
        "receptor_ca_rmsd_A": round(assessment.receptor_ca_rmsd_A, 6),
        "peptide_pose_ca_rmsd_A": round(assessment.peptide_pose_ca_rmsd_A, 6),
        "peptide_internal_ca_rmsd_A": round(assessment.peptide_internal_ca_rmsd_A, 6),
        "ligand_centroid_displacement_A": round(
            assessment.ligand_centroid_displacement_A, 6
        ),
        "receptor_contact_retention_fraction": round(
            assessment.receptor_contact_retention_fraction, 6
        ),
        "disulfide_sg_distances_A": [
            round(value, 6) for value in assessment.disulfide_sg_distances_A
        ],
        "min_interchain_heavy_atom_distance_A": round(
            assessment.min_interchain_heavy_atom_distance_A, 6
        ),
        "rosetta_scores": {
            "rebuild_total_score": rebuild_scores[0].score("total_score"),
            "rebuild_dslf_fa13": rebuild_scores[0].score("dslf_fa13"),
            "prepack_total_score": prepack_scores[0].score("total_score"),
            "prepack_dslf_fa13": prepack_scores[0].score("dslf_fa13"),
            "prepack_fa_rep": prepack_scores[0].score("fa_rep"),
            "refine_total_score": refine_scores[0].score("total_score"),
            "refine_dslf_fa13": refine_scores[0].score("dslf_fa13"),
            "refine_fa_rep": refine_fa_rep,
            "refine_omega": refine_scores[0].score("omega"),
            "refine_rama_prepro": refine_scores[0].score("rama_prepro"),
            "refine_fa_rep_per_residue": round(refine_fa_rep_per_residue, 6),
            "refine_backbone_strain_per_residue": round(
                refine_backbone_strain_per_residue, 6
            ),
        },
    }
    return _write_task_result(
        task_dir=task_dir,
        task_id=task_id,
        plan_sha256=plan_sha256,
        started_at=started_at,
        reconstruction_status=("passed" if not failures else "failed"),
        failure_reasons=failures,
        failure_detail=None,
        metrics=metrics,
        commands=commands,
        artifacts=artifacts,
    )


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], Path, Path, tuple[object, ...]]:
    plan_path = run_dir / PLAN_NAME
    plan = _document(plan_path, name="topology calibration plan")
    if (
        plan.get("schema") != PLAN_SCHEMA
        or plan.get("stage") != "disulfide_topology_calibration"
        or plan.get("status") != "planned"
        or not is_current_vela_software(plan.get("software"))
        or plan.get("contract") != topology_calibration_contract(config)
    ):
        raise DiscoveryError("topology calibration plan identity is invalid")
    try:
        snapshot = object_mapping(plan.get("config_snapshot"), name="config snapshot")
        source = object_mapping(plan.get("source"), name="calibration source")
        tasks = tuple(object_list(plan.get("tasks"), name="calibration tasks"))
    except TypeError as exc:
        raise DiscoveryError("topology calibration plan structure is invalid") from exc
    snapshot_path, _ = validate_record(
        root=run_dir, raw=snapshot, name="config snapshot"
    )
    if sha256_file(snapshot_path) != config.source_snapshot_sha256:
        raise DiscoveryError("current config differs from topology calibration plan")
    relative_source = source.get("run_path")
    if not isinstance(relative_source, str):
        raise DiscoveryError("topology calibration source path is invalid")
    source_run = _source_run(
        config=config, source_run=config.paths.outputs_dir / relative_source
    )
    source_plan = source_run / "qualification_plan.json"
    plan_record = object_mapping(
        source.get("qualification_plan"), name="source qualification plan"
    )
    if plan_record.get("path") != source_plan.name or plan_record.get(
        "sha256"
    ) != sha256_file(source_plan):
        raise DiscoveryError("source qualification plan has changed")
    try:
        reference_receptor_path, _ = validate_record(
            root=config.paths.data_dir,
            raw=source.get("receptor_only_reference"),
            name="topology calibration receptor-only reference",
        )
    except ValidationError as exc:
        raise DiscoveryError(str(exc)) from exc
    if plan.get("task_count") != len(tasks):
        raise DiscoveryError("topology calibration task count is invalid")
    return plan, source_run, reference_receptor_path, tasks


def run_topology_calibration(*, config: AppConfig, run_dir: Path) -> Path:
    """执行或恢复冻结的全原子拓扑校准任务。"""
    plan, source_run, reference_receptor_path, tasks = _verify_plan(
        config=config, run_dir=run_dir
    )
    manifest_path = run_dir / MANIFEST_NAME
    if manifest_path.exists():
        raise DiscoveryError(
            f"topology calibration manifest already exists: {manifest_path}"
        )
    tools = object_mapping(plan.get("tools"), name="calibration tools")
    planned_cg2all = object_mapping(tools.get("cg2all"), name="planned cg2all")
    planned_rosetta = object_mapping(
        tools.get("flexpepdock"), name="planned FlexPepDock"
    )
    planned_rosetta_scripts = object_mapping(
        tools.get("rosetta_scripts"), name="planned RosettaScripts"
    )
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    rosetta_scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    if (
        planned_cg2all
        != {
            "version": cg2all.version,
            "executable_sha256": cg2all.executable_sha256,
            "checkpoint_sha256": cg2all.checkpoint_sha256,
        }
        or planned_rosetta
        != {
            "version": flexpepdock.version,
            "executable_sha256": flexpepdock.executable_sha256,
        }
        or planned_rosetta_scripts
        != {
            "version": rosetta_scripts.version,
            "executable_sha256": rosetta_scripts.executable_sha256,
        }
    ):
        raise DiscoveryError("topology calibration tool identity has changed")
    plan_sha256 = sha256_file(run_dir / PLAN_NAME)
    worker_count = min(config.validation.rosetta.parallel_tasks, len(tasks))
    with ThreadPoolExecutor(max_workers=worker_count) as executor:
        futures = tuple(
            executor.submit(
                _run_task,
                config=config,
                run_dir=run_dir,
                source_run=source_run,
                reference_receptor_path=reference_receptor_path,
                raw_task=task,
                plan_sha256=plan_sha256,
            )
            for task in tasks
        )
        results = tuple(future.result() for future in futures)
    atomic_write_json(
        manifest_path,
        {
            "schema": MANIFEST_SCHEMA,
            "stage": "disulfide_topology_calibration",
            "status": "completed",
            "completed_at": utc_now(),
            "topology_calibration_plan": {
                "path": PLAN_NAME,
                "sha256": plan_sha256,
            },
            "execution": {
                "parallel_tasks": worker_count,
                "task_count": len(results),
            },
            "task_results": [file_record(path, root=run_dir) for path in results],
        },
    )
    return manifest_path


def analyze_topology_calibration(*, config: AppConfig, run_dir: Path) -> Path:
    """按预注册分层规则汇总全原子可恢复性并写出资格证据。"""
    plan, _, _, tasks = _verify_plan(config=config, run_dir=run_dir)
    manifest_path = run_dir / MANIFEST_NAME
    manifest = _document(manifest_path, name="topology calibration manifest")
    if (
        manifest.get("schema") != MANIFEST_SCHEMA
        or manifest.get("status") != "completed"
        or manifest.get("topology_calibration_plan")
        != {"path": PLAN_NAME, "sha256": sha256_file(run_dir / PLAN_NAME)}
    ):
        raise DiscoveryError("topology calibration manifest identity is invalid")
    try:
        result_records = object_list(
            manifest.get("task_results"), name="topology calibration task results"
        )
    except TypeError as exc:
        raise DiscoveryError("topology calibration result records are invalid") from exc
    results: dict[str, dict[str, object]] = {}
    for raw in result_records:
        path, _ = validate_record(
            root=run_dir, raw=raw, name="topology calibration task result"
        )
        result = _document(path, name="topology calibration task result")
        task_id = result.get("task_id")
        if (
            result.get("schema") != TASK_SCHEMA
            or result.get("status") != "completed"
            or not isinstance(task_id, str)
            or task_id in results
        ):
            raise DiscoveryError("topology calibration task result is invalid")
        _result_artifacts(path.parent, result)
        results[task_id] = result
    task_rows: list[dict[str, object]] = []
    for raw in tasks:
        task = object_mapping(raw, name="topology calibration task")
        task_id = task.get("task_id")
        if not isinstance(task_id, str) or task_id not in results:
            raise DiscoveryError("topology calibration task results are incomplete")
        task_rows.append(task)
    calibration = config.discovery.qualification.topology_calibration
    bounds = calibration.strata_upper_bounds_A
    strata: list[dict[str, JsonValue]] = []
    threshold_candidates: list[dict[str, JsonValue]] = []
    contiguous_passed = True
    calibrated_threshold: float | None = None
    for index, upper in enumerate(bounds, 1):
        planned = [
            row
            for row in task_rows
            if object_mapping(row.get("stratum"), name="task stratum").get("index")
            == index
        ]
        passed = [
            row
            for row in planned
            if results[str(row["task_id"])].get("reconstruction_status") == "passed"
        ]
        successful_seeds = sorted(
            {
                _required_integer(
                    row, "source_seed", name="topology calibration source seed"
                )
                for row in passed
            }
        )
        fraction = len(passed) / len(planned) if planned else 0.0
        is_threshold_candidate = upper in calibration.candidate_ca_thresholds_A
        stratum_passed = (
            len(planned) == calibration.models_per_stratum
            and fraction >= calibration.min_success_fraction_per_stratum
            and len(successful_seeds) >= calibration.min_successful_seeds_per_stratum
        )
        if is_threshold_candidate:
            contiguous_passed = contiguous_passed and stratum_passed
            threshold_candidates.append(
                {
                    "maximum_disulfide_ca_distance_A": upper,
                    "all_included_strata_passed": contiguous_passed,
                }
            )
            if contiguous_passed:
                calibrated_threshold = upper
        lower = bounds[index - 2] if index > 1 else None
        strata.append(
            {
                "index": index,
                "lower_exclusive_A": lower,
                "upper_inclusive_A": upper,
                "threshold_candidate": is_threshold_candidate,
                "above_threshold_comparator": not is_threshold_candidate,
                "contributes_to_threshold_selection": is_threshold_candidate,
                "model_count": len(planned),
                "successful_model_count": len(passed),
                "success_fraction": round(fraction, 6),
                "successful_seeds": successful_seeds,
                "successful_seed_count": len(successful_seeds),
                "passed": stratum_passed,
            }
        )
    task_outcomes: list[dict[str, JsonValue]] = []
    for row in task_rows:
        task_id = _required_text(row, "task_id", name="topology calibration task ID")
        stratum = object_mapping(row.get("stratum"), name="task stratum")
        stratum_index = _required_integer(
            stratum, "index", name="topology calibration stratum index"
        )
        result = results[task_id]
        task_outcomes.append(
            {
                "task_id": task_id,
                "source_seed": _required_integer(
                    row, "source_seed", name="topology calibration source seed"
                ),
                "stratum_index": stratum_index,
                "reconstruction_status": json_value(
                    result.get("reconstruction_status"),
                    name="reconstruction status",
                ),
                "failure_reasons": json_value(
                    result.get("failure_reasons"), name="failure reasons"
                ),
                "metrics": json_value(result.get("metrics"), name="task metrics"),
            }
        )
    report_path = run_dir / REPORT_NAME
    qualified = calibrated_threshold is not None
    report_document: dict[str, JsonValue] = {
        "schema": REPORT_SCHEMA,
        "stage": "disulfide_topology_calibration",
        "status": "qualified" if qualified else "unqualified",
        "generated_at": utc_now(),
        "evidence_scope": "development_calibration_only",
        "native_information_used": False,
        "contract": topology_calibration_contract(config),
        "topology_calibration_plan": {
            "path": PLAN_NAME,
            "sha256": sha256_file(run_dir / PLAN_NAME),
        },
        "topology_calibration_manifest": {
            "path": MANIFEST_NAME,
            "sha256": sha256_file(manifest_path),
        },
        "source": json_value(plan.get("source"), name="calibration source"),
        "tools": json_value(plan.get("tools"), name="calibration tools"),
        "strata": strata,
        "threshold_candidates": threshold_candidates,
        "calibrated_max_disulfide_ca_distance_A": calibrated_threshold,
        "all_atom_disulfide_rebuild_success_definition": [
            "cg2all_command_completed",
            "RosettaScripts_ForceDisulfides_repack_completed",
            "FlexPepDock_prepack_completed",
            "single_native_free_FlexPepDock_local_refine_completed",
            "Rosetta_scores_are_finite_and_include_dslf_fa13",
            "refined_weighted_fa_rep_per_residue_within_contract",
            "refined_weighted_backbone_strain_per_residue_within_contract",
            "receptor_CA_pose_RMSD_within_contract",
            "peptide_internal_CA_RMSD_within_contract",
            "ligand_centroid_displacement_within_contract",
            "receptor_contact_retention_within_contract",
            "all_atom_backbones_complete",
            "all_disulfide_SG_distances_within_contract",
            "no_interchain_heavy_atom_distance_below_contract",
        ],
        "limitations": [
            "calibrates standard-amino-acid disulfide topology only",
            "does not validate peptide terminal chemistry or affinity ranking",
            "does not qualify transfer to a non-disulfide macrocycle",
        ],
        "task_outcomes": task_outcomes,
    }
    atomic_write_json(report_path, report_document)
    return report_path
