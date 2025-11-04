from collections import UserDict
from copy import deepcopy
from typing import Any, Iterable, Protocol, Self, Iterator
from collections import defaultdict
from molpy.core.ops.geometry import (
    _dot,
    _norm,
    _unit,
    _cross,
    _rodrigues_rotate,
    _vec_add,
    _vec_scale,
    _vec_sub,
)


class EntityLike(Protocol):
    """Protocol for objects that can act as Entities (Entity or subclass)."""

    def __getitem__(self, key: str) -> Any: ...
    def __setitem__(self, key: str, value: Any) -> None: ...
    def __hash__(self) -> int: ...
    def get(self, key: str, default: Any = None) -> Any: ...


class Entity(UserDict):
    """Dictionary-like base object for all structure elements.

    Minimal by design; no IDs, no persistence, no global context.
    """

    # Keep identity-based hashing/equality from UserDict's object semantics
    def __hash__(self) -> int:  # pragma: no cover - trivial identity
        return id(self)


class LinkLike(Protocol):
    """Protocol for objects that can act as Links (Link or subclass)."""

    endpoints: tuple[EntityLike, ...]


class Link[T: Entity](UserDict):
    """Connectivity object holding direct references to endpoint entities.

    Attributes
    ----------
    endpoints: tuple[Entity]
        The ordered tuple of endpoint entity references.
    """

    endpoints: tuple[T, ...]

    def __init__(self, endpoints: Iterable[T], /, **attrs: Any):
        super().__init__()
        self.endpoints = tuple(endpoints)
        # store remaining attributes in the mapping
        for k, v in attrs.items():
            self.data[k] = v

    def replace_endpoint(self, old: T, new: T) -> None:
        """Replace one endpoint reference with another in-place."""
        self.endpoints = tuple(new if e is old else e for e in self.endpoints)

    def __hash__(self) -> int:  # pragma: no cover - trivial identity
        return id(self)


def get_nearest_type[U](item: U) -> type[U]:
    return type(item)


class TypeBucket[T]:
    def __init__(self) -> None:
        # Map concrete classes to a set of instances
        self._items: dict[type[Any], set[T]] = defaultdict(set)

    def add(self, item: T) -> None:
        cls = get_nearest_type(item)
        self._items[cls].add(item)

    def remove(self, item: T) -> None:
        cls = get_nearest_type(item)
        self._items[cls].discard(item)

    def bucket(self, cls: type[Any]) -> list[T]:
        # Return all items whose concrete class is the same as
        # or a subclass of the requested class.
        result: list[T] = []
        if cls in self._items:
            result.extend(self._items[cls])
        for k, items in self._items.items():
            if k is not cls and isinstance(k, type) and issubclass(k, cls):
                result.extend(items)
        return result

    def classes(self) -> Iterator[type[Any]]:
        return iter(self._items.keys())


class AssemblyLike(Protocol):
    """Protocol for objects that can act as Assemblies (Assembly or subclass)."""

    entities: TypeBucket[EntityLike]
    links: TypeBucket[LinkLike]


class Assembly:
    """Container holding entities and links via typed buckets."""

    def __init__(self) -> None:
        self.entities: TypeBucket[Entity] = TypeBucket()
        self.links: TypeBucket[Link] = TypeBucket()

    # ---------- helpers ----------
    def _iter_all_entities(self) -> Iterable[Entity]:
        for cls in self.entities.classes():
            yield from self.entities.bucket(cls)

    def _iter_all_links(self) -> Iterable[Link]:
        for cls in self.links.classes():
            yield from self.links.bucket(cls)

    # ---------- built-ins ----------
    def copy(self) -> Self:
        new = type(self)()
        # deep-copy entities
        emap: dict[Entity, Entity] = {}
        for ent in self._iter_all_entities():
            cloned = ent.__class__(deepcopy(getattr(ent, "data", None)))
            emap[ent] = cloned
            new.entities.add(cloned)

        # deep-copy links (remap endpoints)
        for link in self._iter_all_links():
            mapped_eps = [emap[ep] for ep in link.endpoints]
            attrs = deepcopy(getattr(link, "data", {}))
            lcls: type[Link] = type(link)
            try:
                new_link = lcls(*mapped_eps, **attrs)  # 两端点位置参数
            except TypeError:
                new_link = lcls(mapped_eps, **attrs)  # 或者列表形式
            new.links.add(new_link)

        return new

    def merge(self, sub: "Assembly") -> dict[Entity, Entity]:
        """Deep-copy `sub` into self, return entity mapping old->new."""
        emap: dict[Entity, Entity] = {}

        # copy entities
        for ent in sub._iter_all_entities():
            cloned = ent.__class__(deepcopy(getattr(ent, "data", None)))
            emap[ent] = cloned
            self.entities.add(cloned)

        # copy links with remapped endpoints
        for link in sub._iter_all_links():
            mapped_eps = [emap[ep] for ep in link.endpoints]
            attrs = deepcopy(getattr(link, "data", {}))
            lcls: type[Link] = type(link)
            try:
                new_link = lcls(*mapped_eps, **attrs)
            except TypeError:
                new_link = lcls(mapped_eps, **attrs)
            self.links.add(new_link)

        return emap

class SpatialMixin:
    """Geometry operations on entities with a "xyz" key only."""

    entities: TypeBucket[Entity]
    links: TypeBucket[Link]

    def move(self, delta: list[float], *, entity_type: type[Entity]) -> None:
        for e in self.entities.bucket(entity_type):
            e["xyz"] = _vec_add(e.get("xyz", [0, 0, 0]), delta)

    def rotate(
        self,
        axis: list[float],
        angle: float,
        about: list[float] | None = None,
        *,
        entity_type: type[Entity],
    ) -> None:
        k = _unit(axis)
        o = [0.0, 0.0, 0.0] if about is None else about
        for e in self.entities.bucket(entity_type):
            pos = e.get("xyz")
            if isinstance(pos, list) and len(pos) == 3:
                e["xyz"] = _rodrigues_rotate(pos, k, angle, o)

    def scale(
        self,
        factor: float,
        about: list[float] | None = None,
        *,
        entity_type: type[Entity],
    ) -> None:
        o = [0.0, 0.0, 0.0] if about is None else about
        for e in self.entities.bucket(entity_type):
            pos = e.get("xyz")
            if isinstance(pos, list) and len(pos) == 3:
                v = _vec_sub(pos, o)
                e["xyz"] = _vec_add(o, _vec_scale(v, factor))

    def align(
        self,
        a: Entity,
        b: Entity,
        *,
        a_dir: list[float] | None = None,
        b_dir: list[float] | None = None,
        flip: bool = False,
        entity_type: type[Entity],
    ) -> None:
        pa = a.get("xyz")
        pb = b.get("xyz")
        if not (
            isinstance(pa, list)
            and isinstance(pb, list)
            and len(pa) == 3
            and len(pb) == 3
        ):
            return  # silently skip if missing positions

        ents = self.entities.bucket(entity_type)

        # rotate if directions provided
        if a_dir is not None and b_dir is not None:
            va = _unit(a_dir)
            vb = _unit(b_dir)
            if flip:
                vb = _vec_scale(vb, -1.0)
            # axis = va x vb; angle = atan2(|axis|, dot)
            axis = _cross(va, vb)
            na = _norm(axis)
            if na > 0:
                # angle via sin/cos components
                from math import atan2

                angle = atan2(na, _dot(va, vb))
                for e in ents:
                    pos = e.get("xyz")
                    if isinstance(pos, list) and len(pos) == 3:
                        e["xyz"] = _rodrigues_rotate(
                            pos, _vec_scale(axis, 1.0 / na), angle, pa
                        )
        # translate so that a -> b
        new_pa = a.get("xyz")
        if isinstance(new_pa, list) and len(new_pa) == 3:
            delta = _vec_sub(pb, new_pa)
            self.move(delta, entity_type=entity_type)


class MembershipMixin:
    """CRUD operations for entities and links within an AssemblyLike."""

    entities: TypeBucket[Entity]
    links: TypeBucket[Link]

    # Entities -------------------------------------------------------------
    def add_entity(self, *ents: Entity) -> None:
        for e in ents:
            self.entities.add(e)

    def remove_entity(self, *ents: Entity, drop_incident_links: bool = True) -> None:
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
            self.entities.remove(e)

    # Links ----------------------------------------------------------------
    def add_link(self, *links: Link, include_endpoints: bool = True) -> None:
        for l in links:
            self.links.add(l)
            if include_endpoints:
                for ep in l.endpoints:
                    self.entities.add(ep)

    def remove_link(self, *links: Link) -> None:
        for l in links:
            self.links.remove(l)

    # Normalize ------------------------------------------------------------
    def normalize(self, include_missing_endpoints: bool = False) -> None:
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
                            self.entities.add(ep)
                            present.add(ep)
                    else:
                        doomed.append(l)
            if doomed:
                self.remove_link(*doomed)

class ConnectivityMixin:

    entities: TypeBucket[Entity]
    links: TypeBucket[Link]

    def get_neighbors(self, entity: Entity, link_type: type[Link] = Link) -> list[Entity]:
        neighbors: list[Entity] = []
        try:
            bucket = self.links.bucket(link_type)
        except KeyError:
            return neighbors
        for link in bucket:
            if entity in link.endpoints:
                for ep in link.endpoints:
                    if ep is not entity:
                        neighbors.append(ep)
        return neighbors