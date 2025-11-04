from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom, Bond


class TestConnectivity:
    def test_neighbors(self) -> None:
        asm = Atomistic()
        a, b, c = Atom(), Atom(), Atom()
        asm.add_link(Bond(a, b))
        asm.add_link(Bond(a, c))
        neigh = set(asm.get_neighbors(a))
        assert neigh == {b, c}
