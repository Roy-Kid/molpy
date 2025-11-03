from __future__ import annotations

from copy import deepcopy
from typing import Any, Iterable, Iterator

from .buckets import BucketRegistry, TypeBucket
from .entity import Entity
from .link import Link


class Assembly:
    """Container holding entities and links via typed buckets."""

    def __init__(self) -> None:
        self.entities = BucketRegistry()
        self.links = BucketRegistry()

    # registration ---------------------------------------------------------
    def register_entity(self, cls: type[Entity]) -> None:
        self.entities.register(cls)

    def register_link(self, cls: type[Link]) -> None:
        self.links.register(cls)

    def bucket(self, cls: type[Any]) -> TypeBucket[Any]:
        try:
            return self.entities.bucket(cls)
        except KeyError:
            return self.links.bucket(cls)

    # connectivity ---------------------------------------------------------
    def neighbors(self, e: Entity) -> Iterator[Entity]:
        for lcls in self.links.classes():
            for link in self.links.bucket(lcls):
                if e in link.endpoints:
                    for other in link.endpoints:
                        if other is not e:
                            yield other

    # built-ins ------------------------------------------------------------
    def copy(self) -> "Assembly":
        new = type(self)()
        # replicate registrations
        for cls in self.entities.classes():
            new.register_entity(cls)
        for cls in self.links.classes():
            new.register_link(cls)
        # deep-copy entities
        emap: dict[Entity, Entity] = {}
        for ecls in self.entities.classes():
            for ent in self.entities.bucket(ecls):
                cloned = ent.__class__(deepcopy(ent.data))
                emap[ent] = cloned
                new.entities.add_by_instance(cloned)
        # deep-copy links with remapped endpoints
        for lcls in self.links.classes():
            for link in self.links.bucket(lcls):
                mapped_eps = [emap[ep] for ep in link.endpoints]
                attrs = deepcopy(link.data)
                try:
                    new_link = lcls(*mapped_eps, **attrs)  # type: ignore[misc]
                except TypeError:
                    new_link = lcls(mapped_eps, **attrs)  # type: ignore[misc]
                new.links.add_by_instance(new_link)
        return new

    def merge(self, sub: "Assembly") -> dict[Entity, Entity]:
        emap: dict[Entity, Entity] = {}
        # ensure registration compatibility
        for cls in sub.entities.classes():
            self.register_entity(cls)
        for cls in sub.links.classes():
            self.register_link(cls)
        # copy entities
        for ecls in sub.entities.classes():
            for ent in sub.entities.bucket(ecls):
                cloned = ent.__class__(deepcopy(ent.data))
                emap[ent] = cloned
                self.entities.add_by_instance(cloned)
        # copy links
        for lcls in sub.links.classes():
            for link in sub.links.bucket(lcls):
                mapped_eps = [emap[ep] for ep in link.endpoints]
                attrs = deepcopy(link.data)
                try:
                    new_link = lcls(*mapped_eps, **attrs)  # type: ignore[misc]
                except TypeError:
                    new_link = lcls(mapped_eps, **attrs)  # type: ignore[misc]
                self.links.add_by_instance(new_link)
        return emap

    def attach(self, a: Entity, b: Entity, link_cls: type[Link], /, **attrs: Any) -> Link:
        try:
            link = link_cls(a, b, **attrs)
        except TypeError:
            link = link_cls([a, b], **attrs)
        self.links.add_by_instance(link)
        # include endpoints by default when attaching a link
        self.entities.add_by_instance(a)
        self.entities.add_by_instance(b)
        return link
