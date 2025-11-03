from __future__ import annotations

from ..assembly import Assembly
from ..buckets import TypeBucket
from ..mixins.membership import MembershipMixin
from ..mixins.spatial import SpatialMixin
from .types import Atom, Bond


class Atomistic(Assembly, MembershipMixin, SpatialMixin):
    def __init__(self) -> None:
        super().__init__()
        self.register_entity(Atom)
        self.register_link(Bond)

    @property
    def atoms(self) -> TypeBucket[Atom]:
        return self.entities.bucket(Atom)

    @property
    def bonds(self) -> TypeBucket[Bond]:
        return self.links.bucket(Bond)
