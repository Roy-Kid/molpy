from __future__ import annotations

from ..assembly import Assembly
from ..buckets import TypeBucket
from ..mixins.membership import MembershipMixin
from ..mixins.spatial import SpatialMixin
from .types import Bead, Bond


class CoarseGrain(Assembly, MembershipMixin, SpatialMixin):
    def __init__(self) -> None:
        super().__init__()
        self.register_entity(Bead)
        self.register_link(Bond)

    @property
    def beads(self) -> TypeBucket[Bead]:
        return self.entities.bucket(Bead)

    @property
    def bonds(self) -> TypeBucket[Bond]:
        return self.links.bucket(Bond)
