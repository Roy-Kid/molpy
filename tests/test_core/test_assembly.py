from __future__ import annotations

from typing import Any

from molpy.core.assembly import Assembly
from molpy.core.entity import Entity
from molpy.core.link import Link


class Node(Entity):
    ...


class Edge(Link):
    def __init__(self, a: Entity, b: Entity, /, **attrs: Any):
        super().__init__([a, b], **attrs)


class TestAssembly:
    def test_register_and_bucket(self) -> None:
        asm = Assembly()
        asm.register_entity(Node)
        asm.register_link(Edge)
        assert asm.bucket(Node) is asm.entities.bucket(Node)
        assert asm.bucket(Edge) is asm.links.bucket(Edge)

    def test_copy_returns_deep_copy(self) -> None:
        asm = Assembly()
        asm.register_entity(Node)
        asm.register_link(Edge)
        a, b = Node({"pos": [0.0, 0.0, 0.0]}), Node({"pos": [1.0, 0.0, 0.0]})
        asm.entities.add_by_instance(a)
        asm.entities.add_by_instance(b)
        e = Edge(a, b, order=1)
        asm.links.add_by_instance(e)
        cpy = asm.copy()
        # entities and links are different objects
        assert next(iter(cpy.entities.bucket(Node))).get("pos") == [0.0, 0.0, 0.0]
        assert a not in set(cpy.entities.bucket(Node))

    def test_merge_returns_mapping(self) -> None:
        left = Assembly()
        left.register_entity(Node)
        left.register_link(Edge)

        right = Assembly()
        right.register_entity(Node)
        right.register_link(Edge)
        a, b = Node({"x": 1}), Node({"x": 2})
        right.entities.add_by_instance(a)
        right.entities.add_by_instance(b)
        right.links.add_by_instance(Edge(a, b))

        emap = left.merge(right)
        assert a in emap and b in emap
        assert emap[a] in set(left.entities.bucket(Node))

    def test_attach_creates_link(self) -> None:
        asm = Assembly()
        asm.register_entity(Node)
        asm.register_link(Edge)
        a, b = Node(), Node()
        link = asm.attach(a, b, Edge, order=2)
        assert link in set(asm.links.bucket(Edge))
