"""Convenience exports for the builder subpackage.

The legacy ``PolymerBuilder`` class and bulk builders have been removed in
favour of the new declarative API documented in
``notebooks/reacter_polymerbuilder_integration.ipynb``.
"""

from .polymer import *  # re-export linear(), connectors, geometry utils
from .placer import (
    Placer,
    VdWSeparator,
    CovalentSeparator,
    LinearOrienter,
    create_vdw_linear_placer,
    create_covalent_linear_placer,
)
from .presets import *
