from molpy.core.aa.base import Atomistic
from molpy.core.aa import Atom, Bond


class TestMembershipMixin:
    def test_add_entity(self) -> None:
        asm = Atomistic()
        a = Atom({"id": 1})
        asm.add_entity(a)
        assert a in set(asm.atoms)

    def test_add_link_include_endpoints(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        asm.add_link(Bond(a, b))
        assert a in set(asm.atoms) and b in set(asm.atoms)
        assert any(True for _ in asm.bonds)

    def test_remove_entity_drop_incident_links(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        bond = Bond(a, b)
        asm.add_link(bond)
        asm.remove_entity(a, drop_incident_links=True)
        assert bond not in set(asm.bonds)
        assert a not in set(asm.atoms)

    def test_remove_link(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        bond = Bond(a, b)
        asm.add_link(bond)
        asm.remove_link(bond)
        assert bond not in set(asm.bonds)

    def test_normalize_prune_missing_endpoints(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        bond = Bond(a, b)
        # add link without endpoints
        asm.add_link(bond, include_endpoints=False)
        asm.normalize(include_missing_endpoints=False)
        assert bond not in set(asm.bonds)

    def test_normalize_include_missing_endpoints(self) -> None:
        asm = Atomistic()
        a, b = Atom(), Atom()
        bond = Bond(a, b)
        asm.add_link(bond, include_endpoints=False)
        asm.normalize(include_missing_endpoints=True)
        assert a in set(asm.atoms) and b in set(asm.atoms)
        assert bond in set(asm.bonds)
