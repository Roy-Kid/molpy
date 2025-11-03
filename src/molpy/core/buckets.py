from __future__ import annotations

from typing import Any, Iterator, TypeVar, Generic

T = TypeVar("T")


class TypeBucket(Generic[T]):
    """A simple set-backed bucket for instances of a registered type."""

    def __init__(self) -> None:
        self.items: set[T] = set()

    def add(self, *objs: T) -> None:
        for o in objs:
            self.items.add(o)

    def discard(self, *objs: T) -> None:
        for o in objs:
            self.items.discard(o)

    def __iter__(self) -> Iterator[T]:
        return iter(self.items)


class BucketRegistry:
    """Registry mapping classes to their instance buckets.

    Notes
    -----
    - Registration is by class; instances are added based on the nearest
      registered superclass in MRO.
    - Discard removes the object from all buckets that contain it.
    """

    def __init__(self) -> None:
        self._buckets: dict[type[Any], TypeBucket[Any]] = {}
        self._order: list[type[Any]] = []

    def register(self, cls: type[Any]) -> None:
        if cls not in self._buckets:
            self._buckets[cls] = TypeBucket()
            self._order.append(cls)

    def bucket(self, cls: type[T]) -> TypeBucket[T]:
        b = self._buckets.get(cls)
        if b is None:
            raise KeyError(f"Class not registered: {cls!r}")
        return b  # type: ignore[return-value]

    def classes(self) -> list[type[Any]]:
        return list(self._order)

    def _find_registered_for_instance(self, obj: Any) -> type[Any] | None:
        for c in obj.__class__.mro():
            if c in self._buckets:
                return c
        return None

    def add_by_instance(self, obj: Any) -> None:
        cls = self._find_registered_for_instance(obj)
        if cls is None:
            raise KeyError(f"No registered class for instance: {obj!r}")
        self._buckets[cls].add(obj)

    def discard_any(self, obj: Any) -> None:
        for b in self._buckets.values():
            b.discard(obj)
