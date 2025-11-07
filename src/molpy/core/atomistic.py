from .entity import ConnectivityMixin, Entity, Assembly, TypeBucket, Link, MembershipMixin, SpatialMixin
from typing import Any, override, cast


class Atom(Entity):
    """Atom entity (expects optional keys like {"type": "C", "pos": [...]})"""


class Bond(Link):
    def __init__(self, a: Atom, b: Atom, /, **attrs: Any):
        super().__init__([a, b], **attrs)

    @property
    def itom(self) -> Atom:
        return self.endpoints[0]

    @property
    def jtom(self) -> Atom:
        return self.endpoints[1]
    

class Atomistic(Assembly, MembershipMixin, SpatialMixin, ConnectivityMixin):

    def __init__(self) -> None:
        super().__init__()

    @property
    def atoms(self) -> list[Atom]:
        items = list(cast(list[Atom], self.entities.bucket(Atom)))
        return items

    @property
    def bonds(self) -> list[Bond]:
        return cast(list[Bond], self.links.bucket(Bond))

    @property
    def symbols(self) -> list[str]:
        atoms = list(self.atoms)
        return [str(a.get("symbol", "")) for a in atoms]
    
    def add_atom(self, **attrs: Any) -> Atom:
        atom = Atom(**attrs)
        self.entities.add(atom)
        return atom
    
    def add_bond(self, a: Atom, b: Atom, **attrs: Any) -> Bond:
        bond = Bond(a, b, **attrs)
        self.links.add(bond)
        return bond

    def move(self, delta: list[float], *, entity_type: type[Entity]=Atom) -> None:
        super().move(delta, entity_type=entity_type)

    def rotate(self, axis: list[float], angle: float, about: list[float] | None = None, *, entity_type: type[Entity]=Atom) -> None:
        super().rotate(axis, angle, about=about, entity_type=entity_type)

    def scale(self, factor: float, about: list[float] | None = None, *, entity_type: type[Entity]=Atom) -> None:
        super().scale(factor, about=about, entity_type=entity_type)

    def align(
        self,
        a: Entity,
        b: Entity,
        *,
        a_dir: list[float] | None = None,
        b_dir: list[float] | None = None,
        flip: bool = False,
        entity_type: type[Entity] = Atom,
    ) -> None:
        super().align(a, b, a_dir=a_dir, b_dir=b_dir, flip=flip, entity_type=entity_type)