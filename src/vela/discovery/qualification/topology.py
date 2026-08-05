"""用开发 seed 独立校准 CABS 二硫环拓扑的全原子可恢复性。"""

from __future__ import annotations

import json
import math
from collections.abc import Mapping
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import gemmi

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
from vela.discovery.qualification.schemas import (
    PLAN_SCHEMA as QUALIFICATION_PLAN_SCHEMA,
)
from vela.discovery.qualification.schemas import (
    SAMPLING_SCHEMA as QUALIFICATION_SAMPLING_SCHEMA,
)
from vela.discovery.sampling.cabsdock import (
    CABS_TASK_RESULT_SCHEMA,
    cabsdock_archive_path,
    verify_cabsdock_tool,
)
from vela.discovery.sampling.evidence import (
    align_trajectory_receptor,
    centroid,
    disulfide_ca_distance,
    disulfide_cabs_sc_distance,
    read_reference_chains,
    read_structure,
    split_model,
    trajectory_ca_contact_residues,
)
from vela.discovery.sampling.materialization import (
    CabsFrameIdentity,
    materialize_cabs_frames,
)
from vela.discovery.sampling.trajectory import (
    iter_cabs_trajectory,
    read_cabs_sequence,
    trajectory_disulfide_ca_distance,
)
from vela.validation.models import ValidationError
from vela.validation.records import file_record, validate_record
from vela.validation.refinement.reconstruction import verify_cg2all_tool
from vela.validation.refinement.topology_recovery import recover_topology
from vela.validation.rosetta import (
    verify_flexpepdock_tool,
    verify_rosetta_scripts_tool,
)

PLAN_NAME = "topology_calibration_plan.json"
MANIFEST_NAME = "topology_calibration_manifest.json"
REPORT_NAME = "topology_calibration_report.json"
PLAN_SCHEMA = "vela.disulfide-topology-calibration-plan/6"
TASK_SCHEMA = "vela.disulfide-topology-calibration-task-result/6"
MANIFEST_SCHEMA = "vela.disulfide-topology-calibration-manifest/6"
REPORT_SCHEMA = "vela.disulfide-topology-calibration-report/6"


@dataclass(frozen=True, slots=True)
class _Candidate:
    """一个来自完整 TRAF 的 native-free 拓扑校准候选。"""

    task_id: str
    seed: int
    frame_identity: CabsFrameIdentity
    archive_path: Path
    archive_sha256: str
    ca_distance_A: float
    interaction_energy: float
    contact_residues: frozenset[str]
    local_position: tuple[float, float, float]
    stratum_index: int

    @property
    def identity(self) -> str:
        return (
            f"{self.task_id}__traf_r{self.frame_identity.replica:02d}_"
            f"m{self.frame_identity.model:04d}"
        )


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


@dataclass(frozen=True, slots=True)
class _SourceEvidence:
    """经完整哈希链验证的资格采样来源。"""

    sampling_path: Path
    task_result_paths: tuple[Path, ...]


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
            "source_pool": "complete_CABS_TRAF_before_energy_filtering",
            "coverage_dimensions": [
                "seed",
                "binding_site",
                "cabsdock_interaction_energy",
                "disulfide_endpoint_CA_distance",
            ],
            "native_information_used": False,
        },
        "all_atom_reconstruction": {
            "pipeline": (
                "cg2all_peptide_then_aligned_experimental_receptor_graft_then_"
                "RosettaScripts_ForceDisulfides_repack_then_FlexPepDock_prepack_"
                "then_site_coordinate_constrained_single_local_refine"
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
            "min_nonlocal_peptide_heavy_atom_distance_A": (
                calibration.min_nonlocal_peptide_heavy_atom_distance_A
            ),
            "site_coordinate_constraints": {
                "atoms": "all_peptide_CA",
                "reference_frame": "fixed_receptor_first_CA",
                "function": "FLAT_HARMONIC",
                "flat_width_A": (calibration.site_coordinate_constraint_flat_width_A),
                "standard_deviation_A": (calibration.site_coordinate_constraint_sd_A),
                "score_weight": calibration.site_coordinate_constraint_weight,
            },
            "rosetta_score_terms_are_descriptive": True,
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


def validate_topology_source_evidence(
    *, source_run: Path, source_plan: Mapping[str, object]
) -> _SourceEvidence:
    """验证资格计划、采样清单、任务结果和 CABS 归档的完整证据链。"""
    plan_path = source_run / "qualification_plan.json"
    sampling_path = source_run / "qualification_sampling.json"
    sampling = _document(sampling_path, name="source qualification sampling")
    if (
        sampling.get("schema") != QUALIFICATION_SAMPLING_SCHEMA
        or sampling.get("stage") != "discovery_qualification"
        or sampling.get("status") != "sampling_completed"
        or sampling.get("target_id") != source_plan.get("target_id")
        or sampling.get("qualification_plan_sha256") != sha256_file(plan_path)
    ):
        raise DiscoveryError("source qualification sampling identity is invalid")
    try:
        planned_rows = object_list(
            source_plan.get("tasks"), name="source qualification plan tasks"
        )
        sampled_rows = object_list(
            sampling.get("tasks"), name="source qualification sampling tasks"
        )
    except TypeError as exc:
        raise DiscoveryError("source qualification task records are invalid") from exc
    planned_ids = tuple(
        _required_text(
            object_mapping(row, name="source planned task"),
            "task_id",
            name="source planned task_id",
        )
        for row in planned_rows
    )
    sampled_by_id: dict[str, dict[str, object]] = {}
    for raw_row in sampled_rows:
        try:
            row = object_mapping(raw_row, name="source sampled task")
        except TypeError as exc:
            raise DiscoveryError("source sampled task is invalid") from exc
        task_id = _required_text(row, "task_id", name="source sampled task_id")
        if task_id in sampled_by_id or row.get("execution_status") != "completed":
            raise DiscoveryError("source qualification sampling tasks are incomplete")
        sampled_by_id[task_id] = row
    if len(set(planned_ids)) != len(planned_ids) or set(sampled_by_id) != set(
        planned_ids
    ):
        raise DiscoveryError("source qualification tasks differ across evidence files")

    plan_sha256 = sha256_file(plan_path)
    result_paths: list[Path] = []
    for task_id in planned_ids:
        task_dir = source_run / "tasks" / task_id
        result_path = task_dir / "task_result.json"
        result = _document(result_path, name=f"source task result {task_id}")
        if (
            result.get("schema") != CABS_TASK_RESULT_SCHEMA
            or result.get("execution_status") != "completed"
            or result.get("task_id") != task_id
            or result.get("run_manifest_sha256") != plan_sha256
        ):
            raise DiscoveryError(f"source task result identity is invalid: {task_id}")
        try:
            outputs = object_mapping(
                result.get("outputs"), name=f"source task outputs {task_id}"
            )
        except TypeError as exc:
            raise DiscoveryError(f"source task outputs are invalid: {task_id}") from exc
        if "cabs_archive" not in outputs:
            raise DiscoveryError(f"source task lacks its CABS archive: {task_id}")
        archive_path: Path | None = None
        for output_name, record in outputs.items():
            try:
                output_path, _ = validate_record(
                    root=task_dir,
                    raw=record,
                    name=f"source task {task_id} output {output_name}",
                )
            except ValidationError as exc:
                raise DiscoveryError(str(exc)) from exc
            if output_name == "cabs_archive":
                archive_path = output_path
        if archive_path is None or archive_path != cabsdock_archive_path(task_dir):
            raise DiscoveryError(f"source CABS archive identity changed: {task_id}")
        result_paths.append(result_path)
    return _SourceEvidence(sampling_path, tuple(result_paths))


def _control_rows(
    *, source_plan: dict[str, object], receptor_id: str, data_dir: Path
) -> tuple[_ControlSource, ...]:
    if (
        source_plan.get("stage") != "discovery_qualification"
        or source_plan.get("schema") != QUALIFICATION_PLAN_SCHEMA
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
            or row.get("receptor_id") != receptor_id
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


def _trajectory_candidates(
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
    archive_path = cabsdock_archive_path(task_dir)
    chains = read_cabs_sequence(archive_path)
    if chains[-1].sequence != chemistry.sequence:
        raise DiscoveryError("topology calibration TRAF peptide is invalid")
    reference_chains = read_reference_chains(reference_path)
    archive_hash = sha256_file(archive_path)
    bounds = config.discovery.qualification.topology_calibration.strata_upper_bounds_A
    candidates: list[_Candidate] = []
    for frame in iter_cabs_trajectory(archive_path=archive_path, chains=chains):
        ca_distance = trajectory_disulfide_ca_distance(frame, chemistry=chemistry)
        stratum = _stratum_index(ca_distance, bounds)
        if stratum is None:
            continue
        receptor = frame.chain_ca[0]
        peptide = frame.chain_ca[-1]
        contacts = trajectory_ca_contact_residues(
            sequence=chains[0],
            receptor=receptor,
            peptide=peptide,
            threshold_A=(
                config.discovery.qualification.topology_calibration.contact_ca_threshold_A
            ),
        )
        if not contacts:
            continue
        alignment = align_trajectory_receptor(
            sequence=chains[0],
            positions=receptor,
            reference_chains=reference_chains,
        )
        peptide_positions = tuple(
            alignment.transform.apply(gemmi.Position(*position)) for position in peptide
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
                frame_identity=CabsFrameIdentity(frame.replica, frame.model),
                archive_path=archive_path,
                archive_sha256=archive_hash,
                ca_distance_A=ca_distance,
                interaction_energy=frame.interaction_energy,
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
    ca_ranks = _percentile_ranks(candidates, attribute="ca_distance_A")
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
                            ca_ranks[candidate.identity],
                        ),
                        (
                            energy_ranks[existing.identity],
                            ca_ranks[existing.identity],
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


def _materialize_calibration_sources(
    *,
    selected: tuple[tuple[_Candidate, int], ...],
    run_dir: Path,
    config: AppConfig,
) -> dict[str, tuple[Path, str, int, float]]:
    """按源任务批量物化少量校准帧, 并核对物化后的粗粒化几何。"""
    chemistry = control_chemistry(control_bound_state(config))
    by_task: dict[str, list[_Candidate]] = {}
    for candidate, _ in selected:
        by_task.setdefault(candidate.task_id, []).append(candidate)
    result: dict[str, tuple[Path, str, int, float]] = {}
    for task_id, candidates in sorted(by_task.items()):
        identities = tuple(candidate.frame_identity for candidate in candidates)
        task_dir = run_dir / "source_frames" / task_id
        _, model_path = materialize_cabs_frames(
            archive_path=candidates[0].archive_path,
            identities=identities,
            task_dir=task_dir,
            settings=config.discovery.cabsdock,
        )
        structure = read_structure(model_path)
        if len(structure) != len(candidates):
            raise DiscoveryError("materialized topology calibration count is invalid")
        model_sha256 = sha256_file(model_path)
        for model_index, (candidate, model) in enumerate(
            zip(candidates, structure, strict=True), 1
        ):
            _, peptide = split_model(model, peptide_sequence=chemistry.sequence)
            materialized_ca_distance = disulfide_ca_distance(
                peptide=peptide, chemistry=chemistry
            )
            if abs(materialized_ca_distance - candidate.ca_distance_A) > 0.002:
                raise DiscoveryError(
                    "materialized CABS frame differs from its TRAF CA geometry"
                )
            result[candidate.identity] = (
                model_path,
                model_sha256,
                model_index,
                disulfide_cabs_sc_distance(peptide=peptide, chemistry=chemistry),
            )
    if set(result) != {candidate.identity for candidate, _ in selected}:
        raise DiscoveryError(
            "topology calibration source materialization is incomplete"
        )
    return result


def write_topology_calibration_plan(
    *, config: AppConfig, source_run: Path, run_id: str
) -> TopologyCalibrationPlan:
    """从开发 seed 冻结分层、native-free 的全原子重建样本。"""
    validate_run_id(run_id)
    source = _source_run(config=config, source_run=source_run)
    source_plan_path = source / "qualification_plan.json"
    source_plan = _document(source_plan_path, name="source qualification plan")
    source_evidence = validate_topology_source_evidence(
        source_run=source, source_plan=source_plan
    )
    source_schema = source_plan.get("schema")
    if not isinstance(source_schema, str):
        raise DiscoveryError("source qualification plan schema is invalid")
    state = control_bound_state(config)
    try:
        control_scope = object_mapping(
            source_plan.get("control_scope"), name="source control scope"
        )
    except TypeError as exc:
        raise DiscoveryError("source qualification control scope is invalid") from exc
    if (
        control_scope.get("native_bound_state_id") != state.state_id
        or control_scope.get("control_receptor_id")
        != config.discovery.qualification.control_receptor_id
    ):
        raise DiscoveryError("source qualification uses a different control system")
    controls = _control_rows(
        source_plan=source_plan,
        receptor_id=config.discovery.qualification.control_receptor_id,
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
        for candidate in _trajectory_candidates(
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
    verify_cabsdock_tool(config.discovery.cabsdock)
    snapshot_path = run_dir / "config.snapshot.txt"
    atomic_write_text(snapshot_path, config.source_snapshot_text)
    materialized = _materialize_calibration_sources(
        selected=tuple(selected), run_dir=run_dir, config=config
    )
    cg2all = verify_cg2all_tool(config.validation.cg2all)
    flexpepdock = verify_flexpepdock_tool(config.validation.rosetta)
    rosetta_scripts = verify_rosetta_scripts_tool(config.validation.rosetta)
    bounds = calibration.strata_upper_bounds_A
    task_rows: list[dict[str, JsonValue]] = []
    for ordinal, (candidate, site_group) in enumerate(selected):
        model_path, model_sha256, model_index, cabs_sc_distance = materialized[
            candidate.identity
        ]
        lower = bounds[candidate.stratum_index - 1] if candidate.stratum_index else None
        upper = bounds[candidate.stratum_index]
        task_rows.append(
            {
                "task_id": f"stratum_{candidate.stratum_index + 1:02d}__{ordinal + 1:03d}",
                "source_task_id": candidate.task_id,
                "source_seed": candidate.seed,
                "source_model": {
                    "path": model_path.relative_to(run_dir).as_posix(),
                    "sha256": model_sha256,
                    "model_index": model_index,
                    "frame": {
                        "replica": candidate.frame_identity.replica,
                        "model": candidate.frame_identity.model,
                    },
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
                    "disulfide_cabs_sc_distance_A": round(cabs_sc_distance, 6),
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
        "qualification_sampling": file_record(
            source_evidence.sampling_path, root=source
        ),
        "task_results": [
            file_record(path, root=source) for path in source_evidence.task_result_paths
        ],
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
    execution_status: str,
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
            "execution_status": execution_status,
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
        root=run_dir, raw=source_model, name="topology calibration source model"
    )
    state = control_bound_state(config)
    chemistry = control_chemistry(state)
    started_at = utc_now()
    recovery = recover_topology(
        config=config,
        source_path=source_path,
        model_index=model_index,
        reference_receptor_path=reference_receptor_path,
        chemistry=chemistry,
        task_dir=task_dir,
        rosetta_seed=rosetta_seed,
    )
    return _write_task_result(
        task_dir=task_dir,
        task_id=task_id,
        plan_sha256=plan_sha256,
        started_at=started_at,
        execution_status=recovery.execution_status,
        reconstruction_status=recovery.reconstruction_status,
        failure_reasons=recovery.failure_reasons,
        failure_detail=recovery.failure_detail,
        metrics=recovery.metrics,
        commands=recovery.commands,
        artifacts=recovery.artifacts,
    )


def _verify_plan(
    *, config: AppConfig, run_dir: Path
) -> tuple[dict[str, object], Path, tuple[object, ...]]:
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
    return plan, reference_receptor_path, tasks


def run_topology_calibration(*, config: AppConfig, run_dir: Path) -> Path:
    """执行或恢复冻结的全原子拓扑校准任务。"""
    plan, reference_receptor_path, tasks = _verify_plan(config=config, run_dir=run_dir)
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
                reference_receptor_path=reference_receptor_path,
                raw_task=task,
                plan_sha256=plan_sha256,
            )
            for task in tasks
        )
        results = tuple(future.result() for future in futures)
    result_documents = tuple(
        _document(path, name="topology calibration task result") for path in results
    )
    atomic_write_json(
        manifest_path,
        {
            "schema": MANIFEST_SCHEMA,
            "stage": "disulfide_topology_calibration",
            "status": (
                "invalid"
                if any(
                    result.get("execution_status") == "invalid"
                    for result in result_documents
                )
                else "completed"
            ),
            "completed_at": utc_now(),
            "topology_calibration_plan": {
                "path": PLAN_NAME,
                "sha256": plan_sha256,
            },
            "execution": {
                "parallel_tasks": worker_count,
                "task_count": len(results),
                "invalid_task_count": sum(
                    result.get("execution_status") == "invalid"
                    for result in result_documents
                ),
            },
            "task_results": [file_record(path, root=run_dir) for path in results],
        },
    )
    return manifest_path


def analyze_topology_calibration(*, config: AppConfig, run_dir: Path) -> Path:
    """按预注册分层规则汇总全原子可恢复性并写出资格证据。"""
    plan, _, tasks = _verify_plan(config=config, run_dir=run_dir)
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
        "calibrated_max_reconstructable_disulfide_ca_distance_A": (
            calibrated_threshold
        ),
        "all_atom_disulfide_rebuild_success_definition": [
            "cg2all_command_completed",
            "RosettaScripts_ForceDisulfides_repack_completed",
            "FlexPepDock_prepack_completed",
            "single_native_free_FlexPepDock_local_refine_completed",
            "peptide_CA_site_coordinate_constraints_applied",
            "required_Rosetta_scores_are_finite_and_include_dslf_fa13",
            "receptor_CA_pose_RMSD_within_contract",
            "peptide_internal_CA_RMSD_within_contract",
            "ligand_centroid_displacement_within_contract",
            "receptor_contact_retention_within_contract",
            "all_atom_backbones_complete",
            "all_disulfide_SG_distances_within_contract",
            "no_interchain_heavy_atom_distance_below_contract",
            "no_nonlocal_peptide_heavy_atom_distance_below_contract",
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
