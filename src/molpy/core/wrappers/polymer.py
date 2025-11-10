from typing import TypeVar, Generic, Self
from ..entity import Assembly, Entity
from .base import Wrapper

T = TypeVar("T", bound=Assembly)


class Polymer(Wrapper[T], Generic[T]):
    """Wrapper representing a polymer (an assembly with named ports).

    Polymer stores ports (remaining connection points) similarly to Monomer
    and wraps an Assembly containing the merged graph. This is intentionally
    minimal for topology-only operations.
    """

    def __init__(self, wrapped: T):
        super().__init__(wrapped)
        self._ports: dict[str, Entity] = {}
        # store port metadata dicts keyed by port name
        self._port_meta: dict[str, dict] = {}

    def set_port(
        self,
        name: str,
        target: Entity,
        *,
        role: str | None = None,
        bond_kind: str | None = None,
        compat: set | str | None = None,
        multiplicity: int = 1,
        priority: int = 0,
    ) -> None:
        """Create or update a port pointing to an entity in the wrapped assembly."""
        from .monomer import Port

        port = Port(name, target, role=role, bond_kind=bond_kind, compat=compat, multiplicity=multiplicity, priority=priority)
        self._ports[name] = port
        self._port_meta[name] = {
            "role": role,
            "bond_kind": bond_kind,
            "compat": compat if compat is not None else '*',
            "multiplicity": multiplicity,
            "priority": priority,
        }

    def get_port(self, name: str):
        return self._ports.get(name)

    @property
    def ports(self) -> dict[str, "Port"]:
        return dict(self._ports)

    def copy(self) -> Self:
        """Deep copy polymer: copy wrapped assembly and remap ports."""
        wrapped = self.unwrap()
        new_wrapped = wrapped.copy()
        new_poly = type(self)(new_wrapped)

        # Build entity map by assuming copy preserves iteration order
        emap: dict[Entity, Entity] = {}
        for et in wrapped.entities.classes():
            old_ents = wrapped.entities.bucket(et)
            new_ents = new_wrapped.entities.bucket(et)
            if len(old_ents) == len(new_ents):
                for a, b in zip(old_ents, new_ents):
                    emap[a] = b

        # Remap ports
        for name, port in self._ports.items():
            if port.target in emap:
                meta = self._port_meta.get(name, {})
                new_poly.set_port(
                    name,
                    emap[port.target],
                    role=meta.get('role'),
                    bond_kind=meta.get('bond_kind'),
                    compat=meta.get('compat'),
                    multiplicity=meta.get('multiplicity', 1),
                    priority=meta.get('priority', 0),
                )

        return new_poly

    def _repr_info(self) -> str:
        if self._ports:
            return f"{len(self._ports)} port(s)"
        return "no ports"
