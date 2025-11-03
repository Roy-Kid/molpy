from __future__ import annotations

from collections import UserDict
from typing import Any, Iterable

from .entity import Entity


class Link(UserDict):
    """Connectivity object holding direct references to endpoint entities.

    Attributes
    ----------
    endpoints: list[Entity]
        The ordered list of endpoint entity references.
    """

    endpoints: list[Entity]

    def __init__(self, endpoints: Iterable[Entity], /, **attrs: Any):
        super().__init__()
        self.endpoints = list(endpoints)
        # store remaining attributes in the mapping
        for k, v in attrs.items():
            self.data[k] = v

    def replace_endpoint(self, old: Entity, new: Entity) -> None:
        """Replace one endpoint reference with another in-place."""
        self.endpoints = [new if e is old else e for e in self.endpoints]

    def __hash__(self) -> int:  # pragma: no cover - trivial identity
        return id(self)
