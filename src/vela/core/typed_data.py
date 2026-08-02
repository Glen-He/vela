"""把外部解析器返回的未类型化容器收窄为 object 边界。"""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from typing import Protocol, runtime_checkable

from vela.core.provenance import JsonValue


@runtime_checkable
class _ObjectMapping(Protocol):
    def items(self) -> Iterable[tuple[object, object]]:
        """迭代键和值。"""
        ...


@runtime_checkable
class _ObjectIterable(Protocol):
    def __iter__(self) -> Iterator[object]:
        """迭代值。"""
        ...


def object_mapping(value: object, *, name: str) -> dict[str, object]:
    """校验字符串键映射并复制到严格类型容器。"""
    if not isinstance(value, _ObjectMapping):
        raise TypeError(f"{name} must be a mapping")
    result: dict[str, object] = {}
    for key, item in value.items():
        if not isinstance(key, str):
            raise TypeError(f"{name} contains a non-string key")
        result[key] = item
    return result


def object_list(value: object, *, name: str) -> list[object]:
    """校验列表并复制到严格类型容器。"""
    if not isinstance(value, _ObjectIterable) or value.__class__ is not list:
        raise TypeError(f"{name} must be a list")
    return list(value)


def json_value(value: object, *, name: str) -> JsonValue:
    """递归校验外部对象能够无损写入项目 JSON 合同。"""
    if value is None or isinstance(value, (str, bool, int, float)):
        return value
    if isinstance(value, _ObjectMapping):
        mapping = object_mapping(value, name=name)
        return {
            key: json_value(item, name=f"{name}.{key}") for key, item in mapping.items()
        }
    if isinstance(value, _ObjectIterable) and value.__class__ is list:
        return [
            json_value(item, name=f"{name}[{index}]")
            for index, item in enumerate(value)
        ]
    raise TypeError(f"{name} is not a JSON value")
