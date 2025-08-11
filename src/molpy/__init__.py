"""
MolPy - Modular Python toolkit for molecular modeling, simulation setup, and structural analysis.

This package provides a unified framework for manipulating atomic structures,
generating force field inputs, and analyzing simulation data.
"""

__version__ = "0.1.0"

# Import all modules - they remain accessible via mp.module_name
from . import analysis, builder, core, io, op, pack, potential, typifier

# Only re-export core classes directly (mp.xxx)
from .core.atomistic import Angle, Atom, Atomistic, Bond, Dihedral, Improper, Topology
from .core.box import Box
from .core.element import Element
from .core.forcefield import ForceField
from .core.frame import Block, Frame
from .core.region import (
    AndRegion,
    BoxRegion,
    Cube,
    NotRegion,
    OrRegion,
    Region,
    SphereRegion,
)
from .core.selection import (
    AtomIndexSelection,
    AtomTypeSelection,
    CoordinateRangeSelection,
    DistanceSelection,
    ElementSelection,
    MaskPredicate,
)
from .core.system import FrameSystem, PeriodicSystem, StructSystem
from .core.trajectory import (
    FrameGenerator,
    FrameProvider,
    SplitStrategy,
    Trajectory,
    TrajectorySplitter,
)
from .core.units import Unit

# Re-export commonly used functions from wrapper
from .core.wrapper import (
    HierarchyWrapper,
    Spatial,
    Wrapper,
    is_wrapped,
    unwrap_all,
    wrap,
)

# Other modules remain accessible via mp.module_name (e.g., mp.io.read_pdb)
# No direct re-exports to avoid namespace pollution
