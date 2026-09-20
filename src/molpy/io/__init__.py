"""
MolPy I/O — the **only** public file I/O surface (``mp.io.read_*`` / ``write_*``).

There is no package-root ``mp.read_*`` / ``mp.write_*``, and no ``MolStore`` /
Zarr layer. Kernels and formats the native core owns are reached through this module
(thin wrappers / readers that call into it); callers never import the core
extension themselves.

Supports:
- Data files (PDB, XYZ, LAMMPS, GROMACS, AMBER, …)
- Force field files (LAMMPS ``*.ff``, OpenMM/OPLS XML, AMBER prmtop, GROMACS top)
- Trajectory files (LAMMPS dump, XYZ, DCD/TRR/XTC where available)

Basic usage::

    import molpy as mp
    from molpy.data import get_forcefield_path

    frame = mp.io.read_pdb("structure.pdb")
    result = mp.io.read_lammps_data("data.lammps", atom_style="full")
    ff = mp.io.read_xml_forcefield(get_forcefield_path("oplsaa.xml"))
    traj = mp.io.read_lammps_trajectory("dump.lammpstrj")
"""

from pathlib import Path

import numpy as np

# Type aliases
PathLike = str | Path

# =============================================================================
# Import order: Deepest to shallowest to avoid circular dependencies
# =============================================================================

# 2. Data Readers and Writers
from .data.ac import AcReader
from .data.amber import AmberInpcrdReader

# 1. Deepest level: Base classes
from .data.base import DataReader, DataWriter
from .data.gro import GroReader, GroWriter
from .data.lammps import LammpsDataReader, LammpsDataResult, LammpsDataWriter
from .data.lammps_molecule import (
    LammpsMoleculeReader,
    LammpsMoleculeWriter,
)
from .data.mol2 import Mol2Reader, Mol2Writer
from .data.smiles import SmilesReader
from .data.pdb import PDBReader, PDBWriter
from .data.top import TopReader
from .data.xsf import XsfReader, XsfWriter
from .data.xyz import XYZReader

# ForceField Readers and Writers
from .forcefield.amber import AmberPrmtopReader
from .forcefield.base import ForceFieldReader, ForceFieldWriter
from .forcefield.lammps import LAMMPSForceFieldWriter
from .forcefield.moltemplate import MolTemplateReader
from .forcefield.top import GromacsTopReader
from .utils import ZipReader

# 5. Factory functions (use the classes above)
from .readers import (
    read_amber,
    read_amber_ac,
    read_amber_inpcrd,
    read_chgcar,
    read_cube,
    read_dcd_trajectory,
    read_gro,
    read_lammps_data,
    read_lammps_forcefield,
    read_lammps_molecule,
    read_lammps_trajectory,
    read_mol2,
    read_smiles,
    write_smarts,
    read_pdb,
    read_pdb_trajectory,
    read_top,
    read_trr_trajectory,
    read_xml_forcefield,
    read_xsf,
    read_xtc_trajectory,
    read_xyz,
    read_xyz_trajectory,
)
from .base import BaseReader
from .trajectory.base import (
    BaseTrajectoryReader,
    TrajectoryWriter,
)

# 3. Trajectory Readers and Writers
from .trajectory.lammps import (
    LammpsDumpLocalWriter,
    LammpsTrajectoryWriter,
)
from .trajectory.xyz import XYZTrajectoryWriter

# 4. Log Readers
from .log import (
    LammpsCpuUse,
    LammpsLoadBalance,
    LammpsLog,
    LammpsLogHeader,
    LammpsLoopTime,
    LammpsMemoryUsage,
    LammpsNeighborStatistics,
    LammpsPerformance,
    LammpsRun,
    LammpsThermo,
    LammpsTimingBreakdown,
    LammpsTimingRow,
    LammpsWarning,
    parse_lammps_log_text,
    read_lammps_log,
)
from .writers import (
    write_gro,
    write_lammps_data,
    write_lammps_data_coeffs,
    write_lammps_forcefield,
    write_lammps_molecule,
    write_bond_react_map,
    write_lammps_bond_react_system,
    write_lammps_system,
    write_lammps_trajectory,
    write_lammps_dump_local,
    write_mol2,
    write_pdb,
    write_top,
    write_trr,
    write_xsf,
    write_xtc,
    write_xyz,
    write_xyz_trajectory,
    write_dcd_trajectory,
    write_cube,
)

# 6. Utility functions (shallowest level)
read_txt = np.loadtxt

# 7. Scientific-record I/O (submodule; names stay off this package root)
from . import mrec

__all__ = [
    # Core types
    "PathLike",
    "mrec",
    # Factory functions - Readers
    "read_amber",
    "read_amber_ac",
    "read_amber_inpcrd",
    "read_gro",
    "read_lammps_log",
    "parse_lammps_log_text",
    "read_lammps_data",
    "read_lammps_forcefield",
    "read_lammps_molecule",
    "read_lammps_trajectory",
    "read_mol2",
    "read_smiles",
    "write_smarts",
    "read_pdb",
    "read_pdb_trajectory",
    "read_top",
    "read_xml_forcefield",
    "read_xsf",
    "read_xyz",
    "read_xyz_trajectory",
    "read_dcd_trajectory",
    "read_trr_trajectory",
    "read_xtc_trajectory",
    "read_cube",
    "read_chgcar",
    # Factory functions - Writers
    "write_gro",
    "write_lammps_data",
    "write_lammps_data_coeffs",
    "write_lammps_forcefield",
    "write_lammps_molecule",
    "write_bond_react_map",
    "write_lammps_bond_react_system",
    "write_lammps_system",
    "write_lammps_trajectory",
    "write_lammps_dump_local",
    "write_mol2",
    "write_pdb",
    "write_top",
    "write_xsf",
    "write_xyz",
    "write_xyz_trajectory",
    "write_trr",
    "write_xtc",
    "write_dcd_trajectory",
    "write_cube",
    # Utility functions
    "read_txt",
    # Data Readers
    "DataReader",
    "AcReader",
    "AmberInpcrdReader",
    "GroReader",
    "LammpsDataReader",
    "LammpsDataResult",
    "LammpsMoleculeReader",
    "Mol2Reader",
    "Mol2Writer",
    "SmilesReader",
    "PDBReader",
    "TopReader",
    "XsfReader",
    "XYZReader",
    # Data Writers
    "DataWriter",
    "GroWriter",
    "LammpsDataWriter",
    "LammpsMoleculeWriter",
    "PDBWriter",
    "XsfWriter",
    # ForceField Readers
    "ForceFieldReader",
    "AmberPrmtopReader",
    "GromacsTopReader",
    "MolTemplateReader",
    # ForceField Writers
    "ForceFieldWriter",
    "LAMMPSForceFieldWriter",
    # Trajectory Readers
    "BaseReader",
    "BaseTrajectoryReader",
    # Trajectory Writers
    "TrajectoryWriter",
    "LammpsDumpLocalWriter",
    "LammpsTrajectoryWriter",
    "XYZTrajectoryWriter",
    # Log Readers
    "LammpsCpuUse",
    "LammpsLoadBalance",
    "LammpsLog",
    "LammpsLogHeader",
    "LammpsLoopTime",
    "LammpsMemoryUsage",
    "LammpsNeighborStatistics",
    "LammpsPerformance",
    "LammpsRun",
    "LammpsThermo",
    "LammpsTimingBreakdown",
    "LammpsTimingRow",
    "LammpsWarning",
    # Utility Classes
    "ZipReader",
]
