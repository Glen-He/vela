"""阶段二分析共用的 site 距离和确定性凝聚聚类。"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence

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
    clusters: list[tuple[T, ...]] = [(item,) for item in sorted(items, key=identity)]
    while True:
        candidates: list[tuple[float, str, str, int, int]] = []
        for first_index, first in enumerate(clusters):
            for second_index in range(first_index + 1, len(clusters)):
                second = clusters[second_index]
                if can_merge is not None and not can_merge(first, second):
                    continue
                maximum = max(
                    distance(left, right) for left in first for right in second
                )
                if maximum <= 1.0:
                    candidates.append(
                        (
                            maximum,
                            min(identity(item) for item in first),
                            min(identity(item) for item in second),
                            first_index,
                            second_index,
                        )
                    )
        if not candidates:
            return clusters
        _, _, _, first_index, second_index = min(candidates)
        merged = tuple(
            sorted(clusters[first_index] + clusters[second_index], key=identity)
        )
        clusters = [
            cluster
            for index, cluster in enumerate(clusters)
            if index not in {first_index, second_index}
        ]
        clusters.append(merged)
        clusters.sort(key=lambda cluster: min(identity(item) for item in cluster))


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
