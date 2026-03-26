"""
Polymer assembly module.

Provides linear polymer assembly with both topology-only and chemical reaction connectors,
plus optional geometric placement via Placer strategies.
"""

from .connectors import (
    AutoConnector,
    BondKind,
    CallbackConnector,
    ChainConnector,
    Connector,
    ConnectorContext,
    ReacterConnector,
    TableConnector,
    TopologyConnector,
)
from .growth_kernel import GrowthKernel, ProbabilityTableKernel
from .polymer_builder import PolymerBuilder
from .sequence_generator import SequenceGenerator, WeightedSequenceGenerator
from .system import (
    Chain,
    FlorySchulzPolydisperse,
    PoissonPolydisperse,
    PolydisperseChainGenerator,
    SchulzZimmPolydisperse,
    SystemPlan,
    SystemPlanner,
    UniformPolydisperse,
)
from .types import (
    ConnectionMetadata,
    ConnectionResult,
    MonomerPlacement,
    MonomerTemplate,
    PolymerBuildResult,
    PortDescriptor,
    StochasticChain,
)

# Re-export tools from molpy.tool.polymer for backward compatibility
from molpy.tool.polymer import (
    BuildPolymer,
    BuildPolymerAmber,
    BuildSystem,
    PlanSystem,
    PrepareMonomer,
    polymer,
    polymer_system,
)

__all__ = [
    # Tools (from molpy.tool.polymer)
    "BuildPolymer",
    "BuildPolymerAmber",
    "BuildSystem",
    "PlanSystem",
    "PrepareMonomer",
    "polymer",
    "polymer_system",
    # Connectors
    "AutoConnector",
    "BondKind",
    "CallbackConnector",
    "ChainConnector",
    "Connector",
    "ConnectorContext",
    "ReacterConnector",
    "TableConnector",
    "TopologyConnector",
    # CGSmiles Builder
    "PolymerBuilder",
    # Sequence Generators
    "SequenceGenerator",
    "WeightedSequenceGenerator",
    # System-level (new three-layer architecture)
    "Chain",
    "FlorySchulzPolydisperse",
    "PoissonPolydisperse",
    "PolydisperseChainGenerator",
    "SchulzZimmPolydisperse",
    "SystemPlan",
    "SystemPlanner",
    "UniformPolydisperse",
    # Types
    "ConnectionMetadata",
    "ConnectionResult",
    "PolymerBuildResult",
    # G-BigSMILES Stochastic Growth Types
    "MonomerTemplate",
    "PortDescriptor",
    "MonomerPlacement",
    "StochasticChain",
    # G-BigSMILES Growth Kernel
    "GrowthKernel",
    "ProbabilityTableKernel",
]
