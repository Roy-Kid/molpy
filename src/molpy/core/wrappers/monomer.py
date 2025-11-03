from __future__ import annotations

from typing import Any

from ..assembly import Assembly
from ..entity import Entity
from ..wrapper import Wrapper


class Monomer(Wrapper[Assembly]):
    def __init__(self, core: Assembly):
        super().__init__(core)
        self.ports: dict[str, Entity] = {}

    def set_port(self, name: str, ent: Entity) -> None:
        self.ports[name] = ent

    def port_names(self) -> list[str]:
        return list(self.ports.keys())

    def copy_graph(self) -> Assembly:
        core = self.unwrap()
        return core.copy()

    def resolve_ports(self, emap: dict[Entity, Entity]) -> dict[str, Entity]:
        resolved: dict[str, Entity] = {}
        for name, ent in self.ports.items():
            if ent in emap:
                resolved[name] = emap[ent]
        return resolved
