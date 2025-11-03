"""System classes for organizing molecular data.

This module defines three types of molecular systems:
- FrameSystem: Uses Frame/Block API for columnar data storage
- StructSystem: Uses Struct | Wrapper API for object graph storage
- PeriodicSystem: Wrapper for periodic systems supporting supercell operations
"""

from dataclasses import dataclass

import numpy as np

from .box import Box
from .forcefield import ForceField
from .frame import Frame
from .entity import Struct
from .wrapper import Wrapper


@dataclass
class FrameSystem:
    frame: Frame
    forcefield: ForceField | None = None
    box: Box | None = None
