from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom, Bond
from molpy.core.wrappers.monomer import Monomer


class TestMonomer:
    def test_ports_and_copy_and_resolve(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        asm.add_link(Bond(a, b))
        mono = Monomer(asm)
        mono.set_port("A", a)
        mono.set_port("B", b)
        names = mono.port_names()
        assert set(names) == {"A", "B"}

        cloned = mono.copy_graph()
        emap = cloned.merge(asm)  # map from asm ents to new in cloned
        resolved = mono.resolve_ports(emap)
        assert set(resolved.keys()) == {"A", "B"}
        assert resolved["A"] is emap[a] and resolved["B"] is emap[b]
