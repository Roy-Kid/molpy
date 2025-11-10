"""
Polymer assembly module.

Provides linear polymer assembly with topology-only connectors
and optional geometric placement via Placer strategies.
"""

from .linear import PolymerBuilder
from .connectors import (
    Connector,
    ConnectorContext,
    BondKind,
    AutoConnector,
    TableConnector,
    ChainConnector,
    CallbackConnector,
)
from .geom_utils import (
    Placer,
    NoOpPlacer,
    DockPlacer,
    rodrigues,
    get_vdw_radius,
)

__all__ = [
    # Builder
    "PolymerBuilder",
    # Connectors
    "Connector",
    "ConnectorContext",
    "BondKind",
    "AutoConnector",
    "TableConnector",
    "ChainConnector",
    "CallbackConnector",
    # Geometry
    "Placer",
    "NoOpPlacer",
    "DockPlacer",
    "rodrigues",
    "get_vdw_radius",
]
