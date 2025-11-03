from __future__ import annotations

from collections import UserDict

from molpy.core.entity import Entity


class TestEntity:
    def test_identity_hash(self) -> None:
        e1 = Entity({"x": 1})
        e2 = Entity({"x": 1})
        assert isinstance(e1, UserDict)
        assert e1 is not e2
        assert hash(e1) != hash(e2)
