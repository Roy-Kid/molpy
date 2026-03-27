"""
Polymer assembly module.

Provides linear polymer assembly with both topology-only and chemical reaction connectors,
plus optional geometric placement via Placer strategies.
"""

from .connectors import (
    Connector,
    ConnectorContext,
)
from .growth_kernel import GrowthKernel, ProbabilityTableKernel
from .placer import CovalentSeparator, LinearOrienter, Placer, VdWSeparator
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

__all__ = [
    # Connector
    "Connector",
    "ConnectorContext",
    # Placer
    "CovalentSeparator",
    "LinearOrienter",
    "Placer",
    "VdWSeparator",
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
