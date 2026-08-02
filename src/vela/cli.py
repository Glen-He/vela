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
    design_source: object = getattr(arguments, "design_source_run", None)
    candidate_ids: object = getattr(arguments, "candidate_id", [])
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
    if design_source is not None and not isinstance(design_source, Path):
        raise TypeError("parsed design source_run is invalid")
    if target_id is not None and not isinstance(target_id, str):
        raise TypeError("parsed target is invalid")
    try:
        raw_candidate_ids = object_list(candidate_ids, name="candidate IDs")
    except TypeError as exc:
        raise TypeError("parsed candidate IDs are invalid") from exc
    parsed_candidate_ids: list[str] = []
    for value in raw_candidate_ids:
        if not isinstance(value, str):
            raise TypeError("parsed candidate IDs are invalid")
        parsed_candidate_ids.append(value)
    if len(parsed_candidate_ids) != len(set(parsed_candidate_ids)):
        raise TypeError("parsed candidate IDs are invalid")
    try:
        raw_clusters = object_list(target_clusters, name="target cluster IDs")
    except TypeError as exc:
        raise TypeError("parsed target cluster IDs are invalid") from exc
    parsed_clusters: list[str] = []
    for value in raw_clusters:
        if not isinstance(value, str):
            raise TypeError("parsed target cluster IDs are invalid")
        parsed_clusters.append(value)
    if len(parsed_clusters) != len(set(parsed_clusters)):
        raise TypeError("parsed target cluster IDs are invalid")
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
        candidate_ids=tuple(parsed_candidate_ids),
        design_source=design_source,
        target_cluster_ids=tuple(parsed_clusters),
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
