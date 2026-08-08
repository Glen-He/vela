"""阶段二分析共用的 site 距离和确定性凝聚聚类。"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from heapq import heappop, heappush

from vela.discovery.analysis.evidence import Point3D


def contact_jaccard_distance(first: frozenset[str], second: frozenset[str]) -> float:
    """返回接触残基集合的 Jaccard 距离。"""
    union = first | second
    if not union:
        return 1.0
    return 1.0 - len(first & second) / len(union)


def normalized_site_distance(
    *,
    first_contacts: frozenset[str],
    first_position: Point3D,
    second_contacts: frozenset[str],
    second_position: Point3D,
    contact_limit: float,
    position_limit: float,
) -> float:
    """以各自门槛归一化接触和空间距离; 使用较差者。"""
    return max(
        contact_jaccard_distance(first_contacts, second_contacts) / contact_limit,
        math.dist(first_position, second_position) / position_limit,
    )


def complete_linkage[T](
    items: Sequence[T],
    *,
    distance: Callable[[T, T], float],
    identity: Callable[[T], str],
    can_merge: Callable[[tuple[T, ...], tuple[T, ...]], bool] | None = None,
) -> list[tuple[T, ...]]:
    """以归一化距离 1.0 为门槛执行确定性 complete-linkage。"""
    ordered = tuple(sorted(items, key=identity))
    names = tuple(identity(item) for item in ordered)
    if len(set(names)) != len(names):
        raise ValueError("complete-linkage item identities must be unique")
    active: dict[int, tuple[T, ...]] = {
        index: (item,) for index, item in enumerate(ordered)
    }
    labels = {index: index for index in range(len(ordered))}
    distances: dict[tuple[int, int], float] = {}
    candidates: list[tuple[float, int, int, int, int]] = []

    def pair_key(first_id: int, second_id: int) -> tuple[int, int]:
        return (first_id, second_id) if first_id < second_id else (second_id, first_id)

    def ordered_pair(first_id: int, second_id: int) -> tuple[int, int]:
        return (
            (first_id, second_id)
            if labels[first_id] < labels[second_id]
            else (second_id, first_id)
        )

    def add_candidate(first_id: int, second_id: int, maximum: float) -> None:
        first_id, second_id = ordered_pair(first_id, second_id)
        first = active[first_id]
        second = active[second_id]
        if maximum <= 1.0 and (can_merge is None or can_merge(first, second)):
            heappush(
                candidates,
                (
                    maximum,
                    labels[first_id],
                    labels[second_id],
                    first_id,
                    second_id,
                ),
            )

    for first_id, first in enumerate(ordered):
        for second_id in range(first_id + 1, len(ordered)):
            maximum = distance(first, ordered[second_id])
            distances[(first_id, second_id)] = maximum
            add_candidate(first_id, second_id, maximum)

    next_id = len(ordered)
    while candidates:
        _, _, _, first_id, second_id = heappop(candidates)
        if first_id not in active or second_id not in active:
            continue
        first = active[first_id]
        second = active[second_id]
        if can_merge is not None and not can_merge(first, second):
            continue
        merged = tuple(sorted(first + second, key=identity))
        other_ids = tuple(
            cluster_id
            for cluster_id in active
            if cluster_id not in {first_id, second_id}
        )
        distances.pop(pair_key(first_id, second_id))
        new_distances: dict[int, float] = {}
        for other_id in other_ids:
            first_distance = distances.pop(pair_key(first_id, other_id))
            second_distance = distances.pop(pair_key(second_id, other_id))
            new_distances[other_id] = max(first_distance, second_distance)
        del active[first_id]
        del active[second_id]
        active[next_id] = merged
        labels[next_id] = min(labels[first_id], labels[second_id])
        for other_id, maximum in new_distances.items():
            distances[pair_key(next_id, other_id)] = maximum
            add_candidate(next_id, other_id, maximum)
        next_id += 1

    return [
        active[cluster_id]
        for cluster_id in sorted(active, key=lambda item: labels[item])
    ]


def bounded_leader_clusters[T](
    items: Sequence[T],
    *,
    distance: Callable[[T, T], float],
    identity: Callable[[T], str],
    maximum_distance: float,
) -> list[tuple[T, ...]]:
    """确定性地建立直径受限簇, 供较大候选池保留多样性后压缩。"""
    if maximum_distance <= 0:
        raise ValueError("maximum_distance must be positive")
    clusters: list[list[T]] = []
    for item in sorted(items, key=identity):
        candidates: list[tuple[float, str, int]] = []
        for index, cluster in enumerate(clusters):
            maximum = max(distance(item, member) for member in cluster)
            if maximum <= maximum_distance:
                candidates.append((maximum, identity(cluster[0]), index))
        if not candidates:
            clusters.append([item])
            continue
        _, _, index = min(candidates)
        clusters[index].append(item)
    return [tuple(cluster) for cluster in clusters]
