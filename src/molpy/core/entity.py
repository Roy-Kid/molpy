from typing import Any

"""Identity exports for the molrs-owned live graph view layer.

Molrs exposes NodeRef / RelationRef / Refs only. MolPy keeps Entity / Link /
Entities as **local** domain names (single definition site here).
"""

from molrs.views import GraphViews, NodeRef, Refs, RelationRef, _GraphViews

# Domain vocabulary — not dual APIs on molrs.
Entity = NodeRef
Link = RelationRef
Entities = Refs

__all__ = [
    "Entities",
    "Entity",
    "GraphViews",
    "Link",
    "NodeRef",
    "Refs",
    "RelationRef",
    "_GraphViews",
]


class NotPublic:
    """Class attribute that hides a native molrs builder behind its ``def_*`` door.

    Reading the attribute raises ``AttributeError`` naming the public method;
    every other attribute lookup on the class is untouched.
    """

    __slots__ = ("_door", "_name")

    def __init__(self, door: str) -> None:
        self._door = door
        self._name = ""

    def __set_name__(self, owner: type, name: str) -> None:
        self._name = name

    def __get__(self, obj: Any, objtype: type | None = None) -> Any:
        raise AttributeError(f"{self._name} is not public; use {self._door} instead")
