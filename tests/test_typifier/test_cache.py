"""RetypeCache: one typing per distinct region structure, confirmed by isomorphism."""

import molpy as mp
from molpy.core.atomistic import Atom, Atomistic
from molpy.typifier.affected_region import AffectedRegion
from molpy.typifier.base import Match, Typifier
from molpy.typifier.cache import RetypeCache


class _CountingTypifier(Typifier[Atomistic]):
    """Types every atom by element and counts how often it was asked."""

    def __init__(self) -> None:
        self.calls = 0

    def match(self, graph: Atomistic) -> Match:
        self.calls += 1
        return Match(nodes=tuple({"type": f"t_{a['element']}"} for a in graph.atoms))


def _chain(length: int) -> tuple[Atomistic, list[Atom]]:
    mol = mp.Atomistic()
    carbons: list[Atom] = []
    prev: Atom | None = None
    for i in range(length):
        c = mol.def_atom(element="C", x=float(i), y=0.0, z=0.0)
        mol.def_bond(c, mol.def_atom(element="H", x=float(i), y=1.0, z=0.0))
        if prev is not None:
            mol.def_bond(prev, c, bond_type=1, bond_number=1)
        carbons.append(c)
        prev = c
    return mol, carbons


def _region(mol: Atomistic, centre: Atom) -> AffectedRegion:
    return AffectedRegion._from(mol, [centre], extract_radius=1, interior_reach=0)


def test_isomorphic_regions_are_typed_once():
    mol, carbons = _chain(7)
    typifier = _CountingTypifier()
    cache = RetypeCache(typifier)
    first = cache.retype(_region(mol, carbons[2]))
    second = cache.retype(_region(mol, carbons[4]))  # same local structure
    assert typifier.calls == 1
    assert second is first


def test_a_different_structure_is_a_miss():
    mol, carbons = _chain(7)
    typifier = _CountingTypifier()
    cache = RetypeCache(typifier)
    cache.retype(_region(mol, carbons[3]))  # interior carbon
    cache.retype(_region(mol, carbons[0]))  # chain end: different neighbourhood
    assert typifier.calls == 2


def test_retype_and_apply_writes_the_interior_types_onto_the_parent():
    mol, carbons = _chain(5)
    cache = RetypeCache(_CountingTypifier())
    cache.retype_and_apply(_region(mol, carbons[2]))
    assert carbons[2].get("type") == "t_C"
