from __future__ import annotations

from collections import UserDict
from copy import deepcopy
from typing import Generic, Iterable, Iterator, TypeVar


class Entity(UserDict):
    """Dictionary-like base object for all structure elements.

    Minimal by design; no IDs, no persistence, no global context.
    """

    # Keep identity-based hashing/equality from UserDict's object semantics
    def __hash__(self) -> int:  # pragma: no cover - trivial identity
        return id(self)


class Struct(UserDict):
    """Lightweight dict-like container with cloning helpers.

    Behaves like a dict, plus:
    - repr shows the "name" field when present
    - clone(**updates) returns a deep-copied Struct with updates applied
    - calling the instance acts like clone: s(name="foo")
    - to_dict() returns a plain dict snapshot
    """

    def __repr__(self) -> str:  # pragma: no cover - representation helper
        name = self.data.get("name")
        return f"<Struct: {name}>" if name else "<Struct>"

    def clone(self, /, **updates) -> "Struct":
        cloned = Struct(**deepcopy(self.data))
        for k, v in updates.items():
            cloned[k] = v
        return cloned

    def __call__(self, /, **updates) -> "Struct":
        return self.clone(**updates)

    def to_dict(self) -> dict:
        return deepcopy(self.data)


T = TypeVar("T")


class Entities(list, Generic[T]):
    """Simple typed container for structure element collections.

    Provides list semantics plus an .add(item) convenience that returns the item.
    Subscription like Entities[Atom] is supported at runtime via __class_getitem__.
    """

    def __init__(self, items: Iterable[T] | None = None):
        super().__init__(items or [])

    def add(self, item: T) -> T:
        self.append(item)
        return item

    # Explicit iterator return type for typing clarity
    def __iter__(self) -> Iterator[T]:  # type: ignore[override]
        return super().__iter__()
