"""实验结合态去配体后的全表面复现计划。"""

from __future__ import annotations

from pathlib import Path

from vela.config import AppConfig
from vela.core.provenance import JsonValue, sha256_file
from vela.core.typed_data import object_list, object_mapping
from vela.discovery.models import DiscoveryError, DiscoveryTask
from vela.discovery.readiness import assess_discovery_readiness
from vela.discovery.sampling.planning import (
    DiscoveryPlan,
    production_authorization,
    write_sampling_plan,
)
from vela.validation.models import BoundStateDefinition, ValidationError
from vela.validation.records import read_document, safe_identifier, validate_record

REPLICATION_ROLE = "bound_state_blind_replication"
REPLICATION_EVIDENCE = "bound_state_blind_replication"
PREPARATION_MANIFEST = "preparation_manifest.json"


def _require_global_method(config: AppConfig, *, target_id: str) -> tuple[str, str]:
    """复用阶段二已资格化的全表面方法, 不引入另一套隐藏参数。"""
    readiness = assess_discovery_readiness(
        target_id=target_id,
        chemistry=config.chemistry,
        settings=config.discovery,
        receptors=config.receptors,
        audit=config.audit,
        preparation=config.preparation,
        data_dir=config.paths.data_dir,
    )
    if not readiness.ready:
        raise ValidationError(
            "bound-state blind replication is not ready: "
            + "; ".join(f"{item.code}: {item.message}" for item in readiness.issues)
        )
    method_id = config.discovery.method_id
    adapter_id = config.discovery.adapter_id
    if method_id is None or adapter_id is None:
        raise ValidationError("qualified global method identity is unresolved")
    return method_id, adapter_id


def _manifest_entries(
    *, config: AppConfig
) -> tuple[Path, dict[str, dict[str, object]]]:
    manifest_path = (
        config.paths.data_dir / "validation" / "bound_states" / PREPARATION_MANIFEST
    )
    manifest = read_document(manifest_path, name="bound-state preparation manifest")
    if manifest.get("schema") != "vela.bound-state-preparation-manifest/1":
        raise ValidationError("bound-state preparation manifest schema is invalid")
    try:
        raw_entries = object_list(
            manifest.get("entries"), name="bound-state preparation entries"
        )
    except TypeError as exc:
        raise ValidationError("bound-state preparation entries are invalid") from exc
    entries: dict[str, dict[str, object]] = {}
    for raw in raw_entries:
        try:
            entry = object_mapping(raw, name="bound-state preparation entry")
        except TypeError as exc:
            raise ValidationError("bound-state preparation entry is invalid") from exc
        state_id = safe_identifier(entry.get("state_id"), name="bound-state ID")
        if state_id in entries:
            raise ValidationError(f"duplicate prepared bound state: {state_id}")
        entries[state_id] = entry
    configured_ids = {state.state_id for state in config.validation.bound_states}
    if set(entries) != configured_ids:
        raise ValidationError(
            "prepared bound-state identities differ from the current config"
        )
    return manifest_path, entries


def _prepared_receptor(
    *,
    config: AppConfig,
    state: BoundStateDefinition,
    entry: dict[str, object],
) -> tuple[Path, str]:
    expected_fields = {
        "receptor_id": state.receptor_id,
        "ligand_id": state.ligand_id,
        "source_ligand_chain": state.ligand_author_chain_id,
        "local_control_kind": state.local_control_kind,
    }
    if any(entry.get(name) != value for name, value in expected_fields.items()):
        raise ValidationError(
            f"prepared bound-state metadata is stale: {state.state_id}"
        )
    try:
        outputs = object_mapping(
            entry.get("outputs"), name=f"{state.state_id} prepared outputs"
        )
    except TypeError as exc:
        raise ValidationError(
            f"prepared outputs are invalid: {state.state_id}"
        ) from exc
    return validate_record(
        root=config.paths.data_dir,
        raw=outputs.get("receptor_only"),
        name=f"{state.state_id} stripped receptor",
    )


def build_replication_tasks(
    config: AppConfig, *, target_id: str
) -> tuple[DiscoveryTask, ...]:
    """按 seed 展开一个目标的去配体结合态受体。"""
    method_id, adapter_id = _require_global_method(config, target_id=target_id)
    _, entries = _manifest_entries(config=config)
    receptors = {receptor.receptor_id: receptor for receptor in config.receptors}
    tasks: list[DiscoveryTask] = []
    selected_states: list[tuple[BoundStateDefinition, str]] = []
    for state in config.validation.bound_states:
        receptor = receptors.get(state.receptor_id)
        if receptor is None:
            raise ValidationError(f"unknown bound-state receptor: {state.receptor_id}")
        if REPLICATION_ROLE in receptor.roles and receptor.target == target_id:
            selected_states.append((state, receptor.target))
    for seed in config.discovery.seeds:
        for state, target in selected_states:
            state_id = safe_identifier(state.state_id, name="replication state ID")
            receptor_path, receptor_hash = _prepared_receptor(
                config=config,
                state=state,
                entry=entries[state.state_id],
            )
            task_id = safe_identifier(
                f"{state_id}__seed_{seed}", name="replication task ID"
            )
            tasks.append(
                DiscoveryTask(
                    task_id=task_id,
                    receptor_id=state_id,
                    target=target,
                    receptor_path=receptor_path,
                    receptor_sha256=receptor_hash,
                    chemistry_id=config.chemistry.chemistry_id,
                    method_id=method_id,
                    adapter_id=adapter_id,
                    seed=seed,
                    evidence_category=REPLICATION_EVIDENCE,
                )
            )
    if not tasks:
        raise ValidationError(
            f"no {REPLICATION_ROLE} receptor is configured for target {target_id}"
        )
    return tuple(tasks)


def write_replication_plan(
    *, config: AppConfig, run_id: str, target_id: str
) -> DiscoveryPlan:
    """冻结去配体结合态复现; 不读取原配体位置作为当前配体约束。"""
    tasks = build_replication_tasks(config, target_id=target_id)
    manifest_path, _ = _manifest_entries(config=config)
    receptor_by_id = {receptor.receptor_id: receptor for receptor in config.receptors}
    states: list[dict[str, JsonValue]] = []
    for state in config.validation.bound_states:
        receptor = receptor_by_id[state.receptor_id]
        if REPLICATION_ROLE in receptor.roles and receptor.target == target_id:
            states.append(
                {
                    "state_id": state.state_id,
                    "source_receptor_id": state.receptor_id,
                    "target": receptor.target,
                }
            )
    try:
        return write_sampling_plan(
            config=config,
            run_id=run_id,
            target_id=target_id,
            run_dir=(config.paths.outputs_dir / "validation" / "replications" / run_id),
            tasks=tasks,
            evidence_category=REPLICATION_EVIDENCE,
            receptor_selection={
                "mode": "configured_bound_state_role",
                "role": REPLICATION_ROLE,
                "target_id": target_id,
                "states": states,
                "reference_receptor": (
                    config.discovery.target(target_id).reference_receptor
                ),
            },
            method_authorization=production_authorization(
                config=config, target_id=target_id
            ),
            additional_inputs={
                "bound_state_preparation_manifest": {
                    "path": manifest_path.relative_to(config.paths.data_dir).as_posix(),
                    "sha256": sha256_file(manifest_path),
                }
            },
        )
    except DiscoveryError as exc:
        raise ValidationError(str(exc)) from exc
