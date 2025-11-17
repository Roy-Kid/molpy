"""Lightweight molpy package initializer for refactored core tests.

Avoid eager imports to prevent legacy dependencies from loading during
unit tests that target the new core architecture.
"""

from .version import version  # noqa: F401

# Core atomistic classes
from .core.atomistic import Atomistic, Atom, Bond, Angle, Dihedral  # noqa: F401

# Core frame and box classes
from .core.frame import Frame, Block  # noqa: F401
from .core.box import Box  # noqa: F401

# Core forcefield classes
from .core.forcefield import (  # noqa: F401
    ForceField,
    AtomisticForcefield,
    AtomType,
    BondType,
    AngleType,
    DihedralType,
    ImproperType,
    PairType,
    AtomStyle,
    BondStyle,
    AngleStyle,
    DihedralStyle,
    ImproperStyle,
    PairStyle,
    Type,
    Style,
    Parameters,
    TypeBucket,
)

# Core topology class
from .core.topology import Topology  # noqa: F401

# Core script classes
from .core.script import Script, ScriptLanguage  # noqa: F401

# Submodules - Import these AFTER core classes to avoid circular imports
from . import typifier
from . import io
from . import data  # noqa: F401
