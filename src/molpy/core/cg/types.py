from __future__ import annotations

from typing import Any

from ..entity import Entity
from ..link import Link


class Bead(Entity):
    """Coarse-grain bead; may include {"type": "X", "pos": [...], "members": list[Entity] | None}."""


class Bond(Link):
    def __init__(self, a: Entity, b: Entity, /, **attrs: Any):
        super().__init__([a, b], **attrs)
