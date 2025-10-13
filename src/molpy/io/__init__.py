"""
MolPy IO Module

This module provides a unified interface for reading and writing various molecular
file formats. It supports both data files and trajectory files with lazy loading
capabilities.
"""

from pathlib import Path
from typing import List, Union

import numpy as np

from molpy.core.forcefield import ForceField
from molpy.core.frame import Frame
from molpy.core.system import FrameSystem

from . import data, forcefield, log, trajectory

# Type aliases for convenience
PathLike = Union[str, Path]

# Utility function
read_txt = np.loadtxt

# =============================================================================
# Data File Readers
# =============================================================================


def read_lammps_data(
    file: PathLike, atom_style: str, frame: Frame | None = None
) -> Frame:
    """Read LAMMPS data file and return a molpy Frame object."""
    from .data.lammps import LammpsDataReader

    reader = LammpsDataReader(Path(file), atom_style)
    return reader.read(frame=frame)


def read_lammps_molecule(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read LAMMPS molecule file and return a molpy Frame object."""
    from .data.lammps_molecule import LammpsMoleculeReader

    reader = LammpsMoleculeReader(Path(file))
    return reader.read(frame=frame)


def read_pdb(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read a PDB file and return a molpy Frame object."""
    from .data.pdb import PDBReader

    reader = PDBReader(Path(file))
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_amber_inpcrd(inpcrd: PathLike, frame: Frame | None = None) -> Frame:
    """Read AMBER inpcrd file and return a molpy Frame object."""
    from .data.amber import AmberInpcrdReader

    reader = AmberInpcrdReader(Path(inpcrd))
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_amber_ac(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read an AC file and return a molpy Frame object."""
    from .data.ac import AcReader

    reader = AcReader(Path(file))
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_mol2(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read a mol2 file and return a molpy Frame object."""
    from .data.mol2 import Mol2Reader

    reader = Mol2Reader(Path(file))
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_xsf(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read an XSF file and return a molpy Frame object."""
    from .data.xsf import XsfReader

    reader = XsfReader(Path(file))
    return reader.read(frame)


def read_gro(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read a GROMACS gro file and return a Frame."""
    from .data.gro import GroReader

    if frame is None:
        frame = Frame()
    reader = GroReader(Path(file))
    frame = reader.read(frame)
    return frame


def read_xyz(file: PathLike, frame: Frame | None = None) -> Frame:
    """Read an XYZ file and return a molpy Frame object."""
    from .data.xyz import XYZReader

    reader = XYZReader(Path(file))
    if frame is None:
        frame = Frame()
    return reader.read(frame)


# =============================================================================
# Force Field Readers
# =============================================================================


def read_lammps_forcefield(
    scripts: PathLike | list[PathLike], forcefield: ForceField | None = None
) -> ForceField:
    """Read LAMMPS force field file and return a molpy ForceField object."""
    from .forcefield.lammps import LAMMPSForceFieldReader

    reader = LAMMPSForceFieldReader(scripts)
    return reader.read(forcefield=forcefield)


def read_amber(
    prmtop: PathLike, inpcrd: PathLike | None = None, frame: Frame | None = None
) -> tuple[Frame, ForceField]:
    """Read AMBER force field prmtop and return a molpy ForceField object."""
    from .forcefield.amber import AmberPrmtopReader

    prmtop = Path(prmtop)
    inpcrd = Path(inpcrd) if inpcrd is not None else None
    reader = AmberPrmtopReader(prmtop)
    if frame is None:
        frame = Frame()
    frame, ff = reader.read(frame)
    if inpcrd is not None:
        from .data.amber import AmberInpcrdReader

        reader = AmberInpcrdReader(inpcrd)
        frame = reader.read(frame)
    return frame, ff


def read_amber_system(
    prmtop: PathLike, inpcrd: PathLike | None = None, system: FrameSystem | None = None
) -> FrameSystem:
    """Read AMBER force field prmtop and coordinates and return a FrameSystem."""
    if system is None:
        system = FrameSystem()

    # Ensure we have a frame to work with
    frame = system._wrapped if system._wrapped is not None else Frame()
    frame, ff = read_amber(prmtop, inpcrd, frame)
    return FrameSystem(frame=frame, forcefield=ff, box=system.box)


def read_xml_forcefield(
    name_or_path: PathLike, system: FrameSystem | None = None
) -> FrameSystem:
    """Read an XML force field file and return a FrameSystem."""
    from .forcefield.xml import XMLForceFieldReader

    # Handle built-in force fields
    if isinstance(name_or_path, str) and not Path(name_or_path).exists():
        builtin_path = Path(__file__).parent / f"forcefield/xml/{name_or_path}.xml"
        if builtin_path.exists():
            xml_path = builtin_path
        else:
            # Try data directory
            import molpy as mp

            data_path = (
                Path(mp.__file__).parent / "data" / "forcefield" / f"{name_or_path}.xml"
            )
            if data_path.exists():
                xml_path = data_path
            else:
                raise FileNotFoundError(
                    f"Built-in force field '{name_or_path}' not found"
                )
    else:
        xml_path = Path(name_or_path)

    if system is None:
        system = FrameSystem()

    reader = XMLForceFieldReader(xml_path)
    reader.read(system)

    return system


def read_top(file: PathLike, forcefield: ForceField | None = None) -> ForceField:
    """Read a GROMACS top file and return a FrameSystem with force field data."""
    from .forcefield.top import GromacsTopReader

    if forcefield is None:
        forcefield = ForceField()

    reader = GromacsTopReader(Path(file))
    forcefield = reader.read(forcefield)
    # Create new FrameSystem with updated forcefield
    return forcefield


# =============================================================================
# System Readers (Combined Data + Force Field)
# =============================================================================


def read_lammps(
    data: PathLike,
    scripts: PathLike | list[PathLike] | None = None,
    frame: Frame | None = None,
    atomstyle: str = "full",
) -> Frame:
    """Read LAMMPS data and force field files and return a molpy Frame object. If data file is provided, only read model;
    If input file is provided, read force field.
    """
    if scripts is not None:  # read definition first
        forcefield = read_lammps_forcefield(scripts, frame)
        # Note: read_lammps_forcefield returns ForceField, not Frame
    frame = read_lammps_data(data, atomstyle, frame)
    return frame


def read_gromacs_system(
    gro_file: PathLike,
    top_file: PathLike | None = None,
    system: FrameSystem | None = None,
) -> FrameSystem:
    """Read GROMACS structure and topology files and return a complete FrameSystem.

    Args:
        gro_file: Path to the GROMACS .gro structure file
        top_file: Optional path to the GROMACS .top topology file
        system: Optional existing FrameSystem to populate

    Returns:
        FrameSystem with structure, box, and optional force field data
    """
    gro_frame = read_gro(gro_file, None)
    if system is None:
        system = FrameSystem(gro_frame)

    if top_file is not None:
        system = read_top(top_file, system)

    return system


# =============================================================================
# Data File Writers
# =============================================================================


def write_lammps_data(file: PathLike, frame: Frame) -> None:
    """Write a molpy Frame object to a LAMMPS data file."""
    from .data.lammps import LammpsDataWriter

    writer = LammpsDataWriter(Path(file))
    writer.write(frame)


def write_pdb(file: PathLike, frame: Frame) -> None:
    """Write a molpy Frame object to a PDB file."""
    from .data.pdb import PDBWriter

    writer = PDBWriter(Path(file))
    writer.write(frame)


def write_xsf(file: PathLike, frame: Frame) -> None:
    """Write a molpy Frame object to an XSF file."""
    from .data.xsf import XsfWriter

    writer = XsfWriter(Path(file))
    writer.write(frame)


def write_lammps_molecule(
    file: PathLike, frame: Frame, format_type: str = "native"
) -> None:
    """Write a molpy Frame object to a LAMMPS molecule file."""
    from .data.lammps_molecule import LammpsMoleculeWriter

    writer = LammpsMoleculeWriter(Path(file), format_type)
    writer.write(frame)


def write_lammps_forcefield(file: PathLike, forcefield: ForceField) -> None:
    """Write a molpy ForceField object to a LAMMPS force field file."""
    from .forcefield.lammps import LAMMPSForceFieldWriter

    writer = LAMMPSForceFieldWriter(
        Path(file),
    )
    writer.write(forcefield)


def write_lammps(workdir: PathLike, frame: Frame, forcefield: ForceField) -> None:
    """Write a molpy FrameSystem object to LAMMPS data and force field files."""
    if not Path(workdir).exists():
        Path(workdir).mkdir(parents=True, exist_ok=True)
    file_path = Path(workdir) / Path(workdir).stem
    write_lammps_data(file_path.with_suffix(".data"), frame)
    write_lammps_forcefield(file_path.with_suffix(".ff"), forcefield)


# =============================================================================
# Trajectory Readers and Writers
# =============================================================================


def read_lammps_trajectory(
    traj: PathLike, frame: Frame | None = None
) -> "trajectory.lammps.LammpsTrajectoryReader":
    """Read LAMMPS trajectory file and return a trajectory reader."""
    from .trajectory.lammps import LammpsTrajectoryReader

    return LammpsTrajectoryReader(Path(traj), frame)


def read_xyz_trajectory(
    file: PathLike,
) -> "trajectory.xyz.XYZTrajectoryReader":
    """Read XYZ trajectory file and return a trajectory reader."""
    from .trajectory.xyz import XYZTrajectoryReader

    return XYZTrajectoryReader(Path(file))


def write_lammps_trajectory(
    file: PathLike, frames: List[Frame], atom_style: str = "full"
) -> None:
    """Write frames to a LAMMPS trajectory file."""
    from .trajectory.lammps import LammpsTrajectoryWriter

    with LammpsTrajectoryWriter(Path(file), atom_style) as writer:
        for i, frame in enumerate(frames):
            timestep = getattr(frame, "step", i)
            writer.write_frame(frame, timestep)


def write_xyz_trajectory(file: PathLike, frames: List[Frame]) -> None:
    """Write frames to an XYZ trajectory file."""
    from .trajectory.xyz import XYZTrajectoryWriter

    with XYZTrajectoryWriter(file) as writer:
        for frame in frames:
            writer.write_frame(frame)


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
    "read_amber",
    "read_amber_system",
    "read_xml_forcefield",
    "read_top",
    # System readers
    "read_lammps",
    "read_gromacs_system",
    # Data writers
    "write_lammps_data",
    "write_pdb",
    "write_xsf",
    "write_lammps_molecule",
    "write_lammps_forcefield",
    "write_lammps",
    # Trajectory functions
    "read_lammps_trajectory",
    "read_xyz_trajectory",
    "write_lammps_trajectory",
    "write_xyz_trajectory",
    # Utility functions
    "read_txt",
]
