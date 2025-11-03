from __future__ import annotations

from typing import Any

from ..assembly import Assembly
from ..entity import Entity
from ..link import Link


class MembershipMixin:
    """CRUD operations for entities and links within an Assembly."""

    # Entities -------------------------------------------------------------
    def add_entity(self: Assembly, *ents: Entity) -> None:
        for e in ents:
            self.entities.add_by_instance(e)

    def remove_entity(self: Assembly, *ents: Entity, drop_incident_links: bool = True) -> None:
        to_remove = set(ents)
        # optionally drop incident links
        if drop_incident_links:
            for lcls in self.links.classes():
                bucket = self.links.bucket(lcls)
                doomed: list[Link] = []
                for l in bucket:
                    if any(ep in to_remove for ep in l.endpoints):
                        doomed.append(l)
                if doomed:
                    self.remove_link(*doomed)
        # finally discard entities
        for e in ents:
            self.entities.discard_any(e)

    # Links ----------------------------------------------------------------
    def add_link(self: Assembly, *links: Link, include_endpoints: bool = True) -> None:
        for l in links:
            self.links.add_by_instance(l)
            if include_endpoints:
                for ep in l.endpoints:
                    self.entities.add_by_instance(ep)

    def remove_link(self: Assembly, *links: Link) -> None:
        for l in links:
            self.links.discard_any(l)

    # Normalize ------------------------------------------------------------
    def normalize(self: Assembly, include_missing_endpoints: bool = False) -> None:
        present: set[Entity] = set()
        for ecls in self.entities.classes():
            present.update(self.entities.bucket(ecls))
        for lcls in self.links.classes():
            bucket = self.links.bucket(lcls)
            doomed: list[Link] = []
            for l in bucket:
                missing = [ep for ep in l.endpoints if ep not in present]
                if missing:
                    if include_missing_endpoints:
                        for ep in missing:
                            self.entities.add_by_instance(ep)
                            present.add(ep)
                    else:
                        doomed.append(l)
            if doomed:
                self.remove_link(*doomed)
