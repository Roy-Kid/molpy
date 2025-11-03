from __future__ import annotations

from molpy.core.aa.assembly import Atomistic
from molpy.core.aa.types import Atom, Bond


class TestConnectivity:
    def test_neighbors(self) -> None:
        asm = Atomistic()
        a, b, c = Atom(), Atom(), Atom()
        asm.add_link(Bond(a, b))
        asm.add_link(Bond(a, c))
        neigh = set(asm.neighbors(a))
        assert neigh == {b, c}
