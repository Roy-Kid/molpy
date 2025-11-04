from ..entity import Entity, Assembly
from .base import Wrapper
from typing import TypeVar
from ..entity import Assembly

T = TypeVar("T", bound=Assembly)

class Monomer[T: Assembly](Wrapper[T]):
    def __init__(self, core: T):
        super().__init__(core)
        self.ports: dict[str, Entity] = {}

    def set_port(self, name: str, ent: Entity) -> None:
        self.ports[name] = ent

    def port_names(self) -> list[str]:
        return list(self.ports.keys())

    def copy_graph(self) -> T:
        core = self.unwrap()
        return core.copy()

    def resolve_ports(self, emap: dict[Entity, Entity]) -> dict[str, Entity]:
        resolved: dict[str, Entity] = {}
        for name, ent in self.ports.items():
            if ent in emap:
                resolved[name] = emap[ent]
        return resolved
