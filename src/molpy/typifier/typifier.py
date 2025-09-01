"""
Simple typifier implementations using mixin pattern.
"""

from molpy.core import ForceField
from molpy.typifier.base import (
    AngleTypifierMixin,
    AtomTypifierMixin,
    BaseTypifier,
    BondTypifierMixin,
    ForceFieldMatchingMixin,
    TopoTypifierMixin,
)


class NullAtomTypifier(BaseTypifier, AtomTypifierMixin):
    """Concrete typifier that doesn't change atom types - just uses existing ones."""

    def typify(self, struct):
        """Do nothing - assume atoms already have types."""
        return struct


class TopoTypifier(BaseTypifier, TopoTypifierMixin):
    """Concrete typifier that focuses on topology-based typing."""

    def typify(self, struct):
        """Apply topology-based typing to the structure."""
        struct = self.typify_bonds(struct)
        struct = self.typify_angles(struct)
        struct = self.typify_dihedrals(struct)
        return struct


class ForceFieldTypifier(BaseTypifier, ForceFieldMatchingMixin):
    """Typifier that matches topology to force field parameters."""

    def __init__(self, forcefield: ForceField):
        ForceFieldMatchingMixin.__init__(self, forcefield)

    def typify(self, struct):
        """Apply force field typing to the structure."""
        # Ensure atoms have types
        for atom in struct.atoms:
            if "type" not in atom:
                atom["type"] = atom.get("element", "unknown")

        # Match topology to force field
        return self.match_topology_to_forcefield(struct)


class SmartsTypifier(
    BaseTypifier, AtomTypifierMixin, BondTypifierMixin, AngleTypifierMixin
):
    """Typifier that uses SMARTS patterns for typing."""

    def __init__(self, forcefield: ForceField):
        self.forcefield = forcefield

    def typify(self, struct):
        """Apply SMARTS-based typing to the structure."""
        struct = self.typify_atoms(struct)
        struct = self.typify_bonds(struct)
        struct = self.typify_angles(struct)
        return struct


class TopologyForceFieldTypifier(
    BaseTypifier, TopoTypifierMixin, ForceFieldMatchingMixin
):
    """Typifier that combines topology typing with force field matching."""

    def __init__(self, forcefield: ForceField):
        ForceFieldMatchingMixin.__init__(self, forcefield)

    def typify(self, struct):
        """Apply topology typing then match to force field."""
        # Apply topology typing
        struct = self.typify_bonds(struct)
        struct = self.typify_angles(struct)
        struct = self.typify_dihedrals(struct)

        # Match to force field
        struct = self.match_topology_to_forcefield(struct)
        return struct
