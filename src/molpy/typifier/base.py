"""
Simple base typifier classes using mixin pattern.
"""

from abc import ABC, abstractmethod

from molpy.core._atomistic import Atomistic
from molpy.core.forcefield import ForceField


class BaseTypifier(ABC):
    """Base abstract class for all typifiers."""

    @abstractmethod
    def typify(self, struct: Atomistic) -> Atomistic:
        """Apply typing to a molecular structure."""
        pass


class AtomTypifierMixin:
    """Mixin for atom typing functionality."""

    def typify_atoms(self, struct: Atomistic) -> Atomistic:
        """Assign types to atoms in the structure."""
        return struct


class TopoTypifierMixin:
    """Mixin for topology-based typing - focuses on bonds, angles, etc."""

    def typify_bonds(self, struct: Atomistic) -> Atomistic:
        """Assign types to bonds based on topology."""
        return struct

    def typify_angles(self, struct: Atomistic) -> Atomistic:
        """Assign types to angles based on topology."""
        return struct

    def typify_dihedrals(self, struct: Atomistic) -> Atomistic:
        """Assign types to dihedrals based on topology."""
        return struct


class BondTypifierMixin:
    """Mixin for bond typing functionality."""

    def typify_bonds(self, struct: Atomistic) -> Atomistic:
        """Assign types to bonds in the structure."""
        return struct


class AngleTypifierMixin:
    """Mixin for angle typing functionality."""

    def typify_angles(self, struct: Atomistic) -> Atomistic:
        """Assign types to angles in the structure."""
        return struct


class ForceFieldMatchingMixin:
    """Mixin for force field parameter matching."""

    def __init__(self, forcefield: ForceField):
        self.forcefield = forcefield

    def match_bonds_to_forcefield(self, struct: Atomistic) -> Atomistic:
        """Match bonds to force field bond types based on atom types."""
        for bond in struct.bonds:
            i_type = bond.itom["type"]
            j_type = bond.jtom["type"]

            for bondstyle in self.forcefield.bondstyles:
                for bondtype in bondstyle.types:
                    if {bondtype.itype.name, bondtype.jtype.name} == {i_type, j_type}:
                        bond["type"] = bondtype.name
                        bond["parameters"] = bondtype.parms
                        break
        return struct

    def match_angles_to_forcefield(self, struct: Atomistic) -> Atomistic:
        """Match angles to force field angle types based on atom types."""
        for angle in struct.angles:
            i_type = angle.itom["type"]
            j_type = angle.jtom["type"]
            k_type = angle.ktom["type"]

            for anglestyle in self.forcefield.anglestyles:
                for angletype in anglestyle.types:
                    if {
                        angletype.itype.name,
                        angletype.jtype.name,
                        angletype.ktype.name,
                    } == {i_type, j_type, k_type}:
                        angle["type"] = angletype.name
                        angle["parameters"] = angletype.parms
                        break
        return struct

    def match_topology_to_forcefield(self, struct: Atomistic) -> Atomistic:
        """Match all topology elements to force field parameters."""
        struct = self.match_bonds_to_forcefield(struct)
        struct = self.match_angles_to_forcefield(struct)
        return struct
