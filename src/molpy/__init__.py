"""Lightweight molpy package initializer for refactored core tests.

Avoid eager imports to prevent legacy dependencies from loading during
unit tests that target the new core architecture.
"""

from .version import version  # noqa: F401
from .core.atomistic import Atomistic, Atom, Bond, Angle, Dihedral  # noqa: F401
from . import typifier
from . import io