from __future__ import annotations

from molpy.core.cg.assembly import CoarseGrain
from molpy.core.cg.types import Bead, Bond


class TestCoarseGrain:
    def test_basic_buckets_and_attach(self) -> None:
        cg = CoarseGrain()
        a, b = Bead(), Bead()
        l = cg.attach(a, b, Bond)
        assert a in set(cg.beads) and b in set(cg.beads)
        assert l in set(cg.bonds)
