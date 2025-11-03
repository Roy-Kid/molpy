from __future__ import annotations

from typing import Any

from ..entity import Entity
from ..link import Link


class Atom(Entity):
    """Atom entity (expects optional keys like {"type": "C", "pos": [...]})"""


class Bond(Link):
    def __init__(self, a: Entity, b: Entity, /, **attrs: Any):
        super().__init__([a, b], **attrs)
