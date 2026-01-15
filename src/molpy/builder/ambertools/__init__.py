"""
AmberTools polymer builder module.

This module provides tools for building polymers using the AmberTools suite
(antechamber, parmchk2, prepgen, tleap).

Public API:
    - AmberPolymerBuilder: Main entry point for polymer building
    - AmberPolymerBuilderConfig: Configuration for the builder
    - AmberBuildResult: Result dataclass containing Frame, ForceField, and paths

Example usage:
    >>> from molpy.builder.ambertools import AmberPolymerBuilder
    >>> from molpy.io import read_pdb
    >>> 
    >>> # Load monomer and mark ports
    >>> eo_monomer = read_pdb("PEO_initial.pdb")
    >>> eo_monomer.atoms[0]["port"] = "<"  # Head port
    >>> eo_monomer.atoms[6]["port"] = ">"  # Tail port
    >>> 
    >>> # Build polymer
    >>> builder = AmberPolymerBuilder(library={"EO": eo_monomer})
    >>> result = builder.build("{[#EO]|10}")
    >>> print(result.frame)
"""

from .amber_builder import AmberPolymerBuilder, AmberPolymerBuilderConfig
from .types import AmberBuildResult

__all__ = [
    "AmberPolymerBuilder",
    "AmberPolymerBuilderConfig",
    "AmberBuildResult",
]
