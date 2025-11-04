from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom, Bond


class TestAtomistic:
    def test_basic_buckets_and_attach(self) -> None:
        aa = Atomistic()
        a, b = Atom(), Atom()
        l = Bond(a, b)
        aa.add_link(l)
        assert a in set(aa.atoms) and b in set(aa.atoms)
        assert l in set(aa.bonds)

    def test_copy_merge(self) -> None:
        aa = Atomistic()
        a, b = Atom({"pos": [0.0, 0.0, 0.0]}), Atom({"pos": [1.0, 0.0, 0.0]})
        aa.add_link(Bond(a, b))
        cpy = aa.copy()
        m = aa.merge(cpy)
        # mapping is from cpy (source) -> aa (dest) new entries
        srcs = set(m.keys())
        assert len(srcs) == 2
        # keys should come from the copied assembly (not original a, b)
        assert all(s is not a and s is not b for s in srcs)
        # values should be present in destination assembly
        assert all(t in set(aa.atoms) for t in m.values())
            
