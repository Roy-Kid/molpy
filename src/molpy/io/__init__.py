"""
MolPy IO Module - Unified interface for molecular file I/O.

This module provides a clean, organized interface for reading and writing
various molecular file formats. It supports:
- Data files (PDB, XYZ, LAMMPS, GROMACS, AMBER, etc.)
- Force field files (LAMMPS, XML, AMBER prmtop, GROMACS top)
- Trajectory files (LAMMPS dump, XYZ)

Design Principles:
-----------------
1. Reader/Writer Pattern: Each format has dedicated Reader and Writer classes
2. Factory Functions: Convenient read_xxx/write_xxx functions for simple usage
3. Lazy Imports: Dependencies loaded only when needed
4. Unified Interface: All readers have read(), all writers have write()

Basic Usage:
-----------
    # Reading data files
    from molpy.io import read_pdb, read_lammps_data
    
    frame = read_pdb("structure.pdb")
    frame = read_lammps_data("data.lammps", atom_style="full")
    
    # Writing data files
    from molpy.io import write_pdb, write_lammps_data
    
    write_pdb("output.pdb", frame)
    write_lammps_data("output.data", frame, atom_style="full")
    
    # Reading force fields
    from molpy.io import read_xml_forcefield, read_lammps_forcefield
    
    ff = read_xml_forcefield("oplsaa.xml")
    ff = read_lammps_forcefield("forcefield.in")
    
    # Reading trajectories
    from molpy.io import read_lammps_trajectory, read_xyz_trajectory
    
    traj = read_lammps_trajectory("dump.lammpstrj")
    for frame in traj:
        process(frame)
"""

from pathlib import Path
from typing import Union
import numpy as np

# Type aliases
PathLike = Union[str, Path]

# =============================================================================
# Import organized reader/writer functions
# =============================================================================

# Data readers
from .readers import (
    read_lammps_data,
    read_lammps_molecule,
    read_pdb,
    read_amber_inpcrd,
    read_amber_ac,
    read_mol2,
    read_xsf,
    read_gro,
    read_xyz,
)

# Force field readers
from .readers import (
    read_lammps_forcefield,
    read_xml_forcefield,
    read_amber_prmtop as read_amber,  # Keep old name for compatibility
    read_top,
)

# Trajectory readers
from .readers import (
    read_lammps_trajectory,
    read_xyz_trajectory,
)

# Data writers
from .writers import (
    write_lammps_data,
    write_pdb,
    write_xsf,
    write_lammps_molecule,
)

# Force field writers
from .writers import (
    write_lammps_forcefield,
)

# Trajectory writers
from .writers import (
    write_lammps_trajectory,
    write_xyz_trajectory,
)

# System writers
from .writers import (
    write_lammps_system as write_lammps,  # Keep old name for compatibility
)

# =============================================================================
# Backward Compatibility: Complex system readers
# =============================================================================

def read_lammps(
    data: PathLike,
    scripts: PathLike | list[PathLike] | None = None,
    frame=None,
    atomstyle: str = "full",
):
    """
    Read LAMMPS data and optional force field files.
    
    Args:
        data: Path to LAMMPS data file
        scripts: Optional path(s) to LAMMPS force field scripts
        frame: Optional existing Frame to populate
        atomstyle: LAMMPS atom style (default: 'full')
        
    Returns:
        Frame object (force field is loaded but returned separately if needed)
    
    Note:
        For new code, prefer using read_lammps_data() and read_lammps_forcefield() separately.
    """
    if scripts is not None:
        # Load force field first (though return value not used in original)
        _ = read_lammps_forcefield(scripts)
    
    return read_lammps_data(data, atomstyle, frame)


def read_amber_system(
    prmtop: PathLike, 
    inpcrd: PathLike | None = None,
    system = None
):
    """
    Read AMBER prmtop and optional inpcrd files (legacy function).
    
    Args:
        prmtop: Path to AMBER prmtop file
        inpcrd: Optional path to AMBER inpcrd file
        system: Optional FrameSystem (unused, kept for compatibility)
        
    Returns:
        Tuple of (Frame, ForceField)
        
    Note:
        For new code, prefer using read_amber() directly.
    """
    frame, ff = read_amber(prmtop, inpcrd)
    
    # Original function returned FrameSystem namedtuple
    from collections import namedtuple
    FrameSystem = namedtuple("FrameSystem", ["frame", "forcefield", "box"])
    return FrameSystem(frame=frame, forcefield=ff, box=getattr(frame, 'box', None))


def read_gromacs_system(
    gro_file: PathLike,
    top_file: PathLike | None = None,
    system = None,
):
    """
    Read GROMACS structure and optional topology files.
    
    Args:
        gro_file: Path to GROMACS .gro file
        top_file: Optional path to GROMACS .top file
        system: Optional FrameSystem (unused, kept for compatibility)
        
    Returns:
        Frame if no topology, or tuple of (Frame, ForceField) if topology provided
        
    Note:
        For new code, prefer using read_gro() and read_top() separately.
    """
    frame = read_gro(gro_file)
    
    if top_file is not None:
        forcefield = read_top(top_file)
        from collections import namedtuple
        FrameSystem = namedtuple("FrameSystem", ["frame", "forcefield", "box"])
        return FrameSystem(frame=frame, forcefield=forcefield, box=getattr(frame, 'box', None))
    
    return frame


# =============================================================================
# Utility functions
# =============================================================================

# Numpy loadtxt shortcut
read_txt = np.loadtxt


# =============================================================================
# Module Exports
# =============================================================================

__all__ = [
    # Core types
    "PathLike",
    
    # Data readers
    "read_lammps_data",
    "read_lammps_molecule",
    "read_pdb",
    "read_amber_inpcrd",
    "read_amber_ac",
    "read_mol2",
    "read_xsf",
    "read_gro",
    "read_xyz",
    
    # Force field readers
    "read_lammps_forcefield",
    "read_xml_forcefield",
    "read_amber",
    "read_top",
    
    # System readers (backward compatibility)
    "read_lammps",
    "read_amber_system",
    "read_gromacs_system",
    
    # Data writers
    "write_lammps_data",
    "write_pdb",
    "write_xsf",
    "write_lammps_molecule",
    
    # Force field writers
    "write_lammps_forcefield",
    
    # System writers
    "write_lammps",
    
    # Trajectory functions
    "read_lammps_trajectory",
    "read_xyz_trajectory",
    "write_lammps_trajectory",
    "write_xyz_trajectory",
    
    # Utility functions
    "read_txt",
]
