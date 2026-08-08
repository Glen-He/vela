"""Vela 的唯一进程入口, 参数和业务处理位于 commands 包。"""

from __future__ import annotations

import sys
from collections.abc import Sequence
from pathlib import Path

from vela.commands.handlers import execute
from vela.commands.parser import build_parser
from vela.config import load_config
from vela.core.errors import VelaError
from vela.core.typed_data import object_list


def _unique_strings(raw: object, *, name: str) -> tuple[str, ...]:
    """把 argparse 的可重复字符串参数收窄为有序唯一元组。"""
    try:
        values = object_list(raw, name=name)
    except TypeError as exc:
        raise TypeError(f"parsed {name} are invalid") from exc
    if any(not isinstance(value, str) for value in values):
        raise TypeError(f"parsed {name} are invalid")
    result = tuple(value for value in values if isinstance(value, str))
    if len(result) != len(set(result)):
        raise TypeError(f"parsed {name} are invalid")
    return result


def _unique_paths(raw: object, *, name: str) -> tuple[Path, ...]:
    """把argparse的可重复路径参数收窄为有序唯一元组。"""
    try:
        values = object_list(raw, name=name)
    except TypeError as exc:
        raise TypeError(f"parsed {name} are invalid") from exc
    if any(not isinstance(value, Path) for value in values):
        raise TypeError(f"parsed {name} are invalid")
    result = tuple(value for value in values if isinstance(value, Path))
    if len(result) != len(set(result)):
        raise TypeError(f"parsed {name} are invalid")
    return result


def run(argv: Sequence[str] | None = None) -> int:
    """解析参数并执行一个命令, 便于测试复用。"""
    arguments = build_parser().parse_args(argv)
    path = arguments.config
    group = arguments.group
    command = arguments.command
    run_id: object = getattr(arguments, "run_id", None)
    run_dir: object = getattr(arguments, "run_dir", None)
    discovery_run: object = getattr(arguments, "discovery_run", None)
    control_run: object = getattr(arguments, "control_run", None)
    replication_run: object = getattr(arguments, "replication_run", None)
    refinement_source: object = getattr(arguments, "source_run", None)
    topology_source: object = getattr(arguments, "topology_source_run", None)
    exploration_basis: object = getattr(arguments, "exploration_basis_run", None)
    qualification_source: object = getattr(arguments, "qualification_run", None)
    site_budget: object = getattr(arguments, "site_budget", None)
    design_source: object = getattr(arguments, "design_source_run", None)
    candidate_ids: object = getattr(arguments, "candidate_id", [])
    blind_candidate_ids: object = getattr(arguments, "blind_candidate_id", [])
    functional_candidate_ids: object = getattr(arguments, "functional_candidate_id", [])
    negative_refinement_runs: object = getattr(arguments, "negative_refinement_run", [])
    start_ids: object = getattr(arguments, "start_id", [])
    receptor_backbone_mode: object = getattr(arguments, "receptor_backbone_mode", None)
    target_clusters: object = getattr(arguments, "target_cluster", [])
    target_id: object = getattr(arguments, "target", None)
    if not isinstance(path, Path):
        raise TypeError("parsed config path is not a Path")
    if not isinstance(group, str) or not isinstance(command, str):
        raise TypeError("parsed command identity is invalid")
    if run_id is not None and not isinstance(run_id, str):
        raise TypeError("parsed run_id is invalid")
    if run_dir is not None and not isinstance(run_dir, Path):
        raise TypeError("parsed run_dir is invalid")
    if discovery_run is not None and not isinstance(discovery_run, Path):
        raise TypeError("parsed discovery_run is invalid")
    if control_run is not None and not isinstance(control_run, Path):
        raise TypeError("parsed control_run is invalid")
    if replication_run is not None and not isinstance(replication_run, Path):
        raise TypeError("parsed replication_run is invalid")
    if refinement_source is not None and not isinstance(refinement_source, Path):
        raise TypeError("parsed source_run is invalid")
    if topology_source is not None and not isinstance(topology_source, Path):
        raise TypeError("parsed topology source_run is invalid")
    if exploration_basis is not None and not isinstance(exploration_basis, Path):
        raise TypeError("parsed exploration basis_run is invalid")
    if qualification_source is not None and not isinstance(qualification_source, Path):
        raise TypeError("parsed qualification_run is invalid")
    if site_budget is not None and (
        not isinstance(site_budget, int) or isinstance(site_budget, bool)
    ):
        raise TypeError("parsed site_budget is invalid")
    if design_source is not None and not isinstance(design_source, Path):
        raise TypeError("parsed design source_run is invalid")
    if target_id is not None and not isinstance(target_id, str):
        raise TypeError("parsed target is invalid")
    if receptor_backbone_mode is not None and not isinstance(
        receptor_backbone_mode, str
    ):
        raise TypeError("parsed receptor_backbone_mode is invalid")
    parsed_candidate_ids = _unique_strings(candidate_ids, name="candidate IDs")
    parsed_blind_candidate_ids = _unique_strings(
        blind_candidate_ids, name="blind candidate IDs"
    )
    parsed_functional_candidate_ids = _unique_strings(
        functional_candidate_ids, name="functional candidate IDs"
    )
    parsed_start_ids = _unique_strings(start_ids, name="start IDs")
    parsed_negative_refinement_runs = _unique_paths(
        negative_refinement_runs, name="negative refinement runs"
    )
    parsed_clusters = _unique_strings(target_clusters, name="target cluster IDs")
    return execute(
        group=group,
        command=command,
        config=load_config(path),
        run_id=run_id,
        run_dir=run_dir,
        source_run=discovery_run,
        control_run=control_run,
        replication_run=replication_run,
        refinement_source=refinement_source,
        topology_source=topology_source,
        exploration_basis=exploration_basis,
        qualification_source=qualification_source,
        site_budget=site_budget,
        candidate_ids=parsed_candidate_ids,
        blind_candidate_ids=parsed_blind_candidate_ids,
        functional_candidate_ids=parsed_functional_candidate_ids,
        negative_refinement_runs=parsed_negative_refinement_runs,
        start_ids=parsed_start_ids,
        receptor_backbone_mode=receptor_backbone_mode,
        design_source=design_source,
        target_cluster_ids=parsed_clusters,
        target_id=target_id,
    )


def main() -> None:
    """控制进程退出状态并保持终端错误为英文。"""
    try:
        exit_code = run()
    except (VelaError, OSError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc
    raise SystemExit(exit_code)
