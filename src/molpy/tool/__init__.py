"""Tools and compute operations for molecular modeling.

Tools are packaged recipes that wrap multiple MolPy modules into
single-call operations. Compute classes are analysis operations
on trajectory data.

Tool examples::

    from molpy.tool import PrepareMonomer, BuildPolymer, polymer

    prep = PrepareMonomer()
    eo = prep.run("{[<]CCO[>]}")

    chain = polymer("{[<]CCO[>]}|10|")

Compute examples::

    from molpy.tool import MSD, DisplacementCorrelation

    msd = MSD(max_lag=3000)
    msd_values = msd(unwrapped_coords)           # -> NDArray (max_lag,)

    xdc = DisplacementCorrelation(max_lag=3000)
    corr = xdc(cation_coords, anion_coords)       # -> NDArray (max_lag,)
"""

from __future__ import annotations

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from molpy.core.atomistic import Atomistic

# Base classes
from .base import Compute, Tool, ToolRegistry

# Compute operations (analysis)
from .cross_correlation import DisplacementCorrelation, displacement_correlation
from .msd import MSD
from .msd import msd as compute_msd_shorthand
from .time_series import TimeAverage, TimeCache, compute_acf, compute_msd

# Tool operations (polymer building recipes)
from .polymer import (
    BuildPolymer,
    BuildPolymerAmber,
    BuildSystem,
    PlanSystem,
    PrepareMonomer,
    polymer,
    polymer_system,
)

# Optional RDKit compute nodes
try:  # pragma: no cover
    from .rdkit import Generate3D, OptimizeGeometry

    _HAS_RDKIT = True
except ModuleNotFoundError:  # rdkit missing
    _HAS_RDKIT = False
    Generate3D = None  # type: ignore[assignment]
    OptimizeGeometry = None  # type: ignore[assignment]


def generate_3d(
    mol: "Atomistic",
    add_hydrogens: bool = True,
    optimize: bool = True,
) -> "Atomistic":
    """Generate 3D coordinates for a molecular structure via RDKit.

    Wraps RDKitAdapter + Generate3D into a single convenience call.

    Args:
        mol: Atomistic structure (typically from parser.parse_molecule)
        add_hydrogens: Add implicit hydrogens before embedding
        optimize: Run force-field geometry optimization after embedding

    Returns:
        New Atomistic with 3D coordinates and (optionally) explicit hydrogens

    Raises:
        ImportError: if RDKit is not installed
    """
    from molpy.adapter import RDKitAdapter

    if RDKitAdapter is None or Generate3D is None:
        raise ImportError(
            "RDKit is required for 3D coordinate generation. "
            "Install with: pip install rdkit"
        )

    adapter = RDKitAdapter(internal=mol)
    compute = Generate3D(
        add_hydrogens=add_hydrogens,
        embed=True,
        optimize=optimize,
        update_internal=True,
    )
    adapter = compute(adapter)
    return adapter.get_internal()


__all__ = [
    # Base
    "Compute",
    "Tool",
    "ToolRegistry",
    # Compute operations
    "MSD",
    "DisplacementCorrelation",
    "TimeCache",
    "TimeAverage",
    "compute_msd",
    "compute_msd_shorthand",
    "compute_acf",
    "displacement_correlation",
    "generate_3d",
    # Tool operations (polymer building)
    "PrepareMonomer",
    "BuildPolymer",
    "PlanSystem",
    "BuildSystem",
    "BuildPolymerAmber",
    "polymer",
    "polymer_system",
]

if _HAS_RDKIT:
    __all__ += ["Generate3D", "OptimizeGeometry"]
