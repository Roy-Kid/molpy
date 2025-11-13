"""
Polymer assembly module.

Provides linear polymer assembly with both topology-only and chemical reaction connectors,
plus optional geometric placement via Placer strategies.
"""

from .linear import linear
from .connectors import (
    Connector,
    ConnectorContext,
    BondKind,
    AutoConnector,
    TableConnector,
    ChainConnector,
    CallbackConnector,
    ReacterConnector,
    TopologyConnector,
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
    "linear",
    # Connectors
    "Connector",
    "ConnectorContext",
    "BondKind",
    "AutoConnector",
    "TableConnector",
    "ChainConnector",
    "CallbackConnector",
    "ReacterConnector",
    "TopologyConnector",
    # Geometry
    "Placer",
    "NoOpPlacer",
    "DockPlacer",
    "rodrigues",
    "get_vdw_radius",
]
