from pathlib import Path
from typing import List

import numpy as np

from molpy.core.forcefield import ForceField
from molpy.core.frame import Frame
from molpy.core.system import FrameSystem

from . import data, forcefield, log, trajectory

read_txt = np.loadtxt


def read_lammps_data(
    file: Path | str, atom_style: str, frame: Frame | None = None
) -> Frame:
    """Read LAMMPS data file and return a molpy System object."""
    from .data.lammps import LammpsDataReader

    reader = LammpsDataReader(file, atom_style)
    return reader.read(frame=frame)


def read_lammps_forcefield(
    scripts: Path | list[Path], frame: ForceField | None = None
) -> ForceField:
    """Read LAMMPS force field file and return a molpy ForceField object."""
    from .forcefield.lammps import LAMMPSForceFieldReader

    reader = LAMMPSForceFieldReader(scripts)
    return reader.read(frame=frame)


def read_lammps_molecule(file: Path | str, frame: Frame | None = None) -> Frame:
    """Read LAMMPS molecule file and return a molpy Frame object."""
    from .data.lammps_molecule import LammpsMoleculeReader

    reader = LammpsMoleculeReader(file)
    return reader.read(frame=frame)


def read_lammps(
    data: Path, scripts: Path | list[Path] | None = None, frame: Frame | None = None
) -> Frame:
    """Read LAMMPS data and force field files and return a molpy System object. If data file is provided, only read model;
    If input file is provided, read force field.
    """
    if scripts is not None:  # read defination first
        frame = read_lammps_forcefield(scripts, frame)
    frame = read_lammps_data(data, frame)
    return frame


def read_pdb(file: Path | str, frame: Frame | None = None) -> Frame:
    """Read a PDB file and return a molpy Frame object."""
    from .data.pdb import PDBReader

    reader = PDBReader(file)
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_amber(
    prmtop: Path, inpcrd: Path | None = None, frame: Frame | None = None
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


def read_amber_inpcrd(inpcrd: Path, frame: Frame | None = None) -> Frame:
    """Read AMBER inpcrd file and return a molpy Frame object."""
    from .data.amber import AmberInpcrdReader

    reader = AmberInpcrdReader(inpcrd)
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_amber_system(
    prmtop: Path, inpcrd: Path | None = None, system: FrameSystem | None = None
) -> FrameSystem:
    """Read AMBER force field prmtop and coordinates and return a FrameSystem."""
    if system is None:
        system = FrameSystem()

    # Ensure we have a frame to work with
    frame = system._wrapped if system._wrapped is not None else Frame()
    frame, ff = read_amber(prmtop, inpcrd, frame)
    return FrameSystem(frame=frame, forcefield=ff, box=system.box)


def read_amber_ac(file: Path | str, frame: Frame | None = None) -> Frame:
    """Read an AC file and return a molpy System object."""
    from .data.ac import AcReader

    reader = AcReader(file)
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_mol2(file: Path | str, frame: Frame | None = None) -> Frame:
    """Read a mol2 file and return a molpy System object."""
    from .data.mol2 import Mol2Reader

    reader = Mol2Reader(file)
    if frame is None:
        frame = Frame()
    return reader.read(frame)


def read_xsf(file: Path | str) -> FrameSystem:
    """Read an XSF file and return a molpy FrameSystem object."""
    from .data.xsf import XsfReader

    reader = XsfReader(file)
    return reader.read()


def read_xml_forcefield(
    name_or_path: Path | str, system: FrameSystem | None = None
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


def read_gro(file: Path | str, frame: Frame | None = None) -> Frame:
    """Read a GROMACS gro file and return a Frame."""
    from .data.gro import GroReader

    if frame is None:
        frame = Frame(frame=Frame())
    reader = GroReader(Path(file))
    frame = reader.read(frame)
    return frame


def read_top(file: Path | str, system: FrameSystem | None = None) -> FrameSystem:
    """Read a GROMACS top file and return a FrameSystem with force field data."""
    from .forcefield.top import GromacsTopReader

    if system is None:
        system = FrameSystem()

    reader = GromacsTopReader(Path(file))
    forcefield = reader.read(system.forcefield)
    # Create new FrameSystem with updated forcefield
    return FrameSystem(frame=system._wrapped, forcefield=forcefield, box=system.box)


def read_gromacs_system(
    gro_file: Path | str,
    top_file: Path | str | None = None,
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
    if system is None:
        system = FrameSystem()

    # Read structure file
    system = read_gro(gro_file, system)

    # Read topology file if provided
    if top_file is not None:
        system = read_top(top_file, system)

    return system


def write_lammps_data(file: Path, frame: Frame) -> None:
    """Write a molpy System object to a LAMMPS data file."""
    from .data.lammps import LammpsDataWriter

    writer = LammpsDataWriter(file)
    writer.write(frame)


def write_pdb(file: Path | str, frame: Frame) -> None:
    """Write a molpy System object to a PDB file."""
    from .data.pdb import PDBWriter

    writer = PDBWriter(Path(file))
    writer.write(frame)


def write_xsf(file: Path | str, system: FrameSystem) -> None:
    """Write a molpy FrameSystem object to an XSF file."""
    from .data.xsf import XsfWriter

    writer = XsfWriter(Path(file))
    writer.write(system)


def write_lammps_molecule(
    file: Path | str, frame: Frame, format_type: str = "native"
) -> None:
    """Write a molpy Frame object to a LAMMPS molecule file."""
    from .data.lammps_molecule import LammpsMoleculeWriter

    writer = LammpsMoleculeWriter(file, format_type)
    writer.write(frame)


def write_lammps_forcefield(file: Path | str, forcefield: ForceField) -> None:
    """Write a molpy System object to a LAMMPS force field file."""
    from .forcefield.lammps import LAMMPSForceFieldWriter

    writer = LAMMPSForceFieldWriter()
    writer.write(file, forcefield)


def write_lammps(workdir: Path, system: FrameSystem) -> None:
    """Write a molpy FrameSystem object to LAMMPS data and force field files."""
    if not workdir.exists():
        workdir.mkdir(parents=True, exist_ok=True)
    file_path = workdir / workdir.stem
    write_lammps_data(file_path.with_suffix(".data"), system._wrapped)
    write_lammps_forcefield(file_path.with_suffix(".ff"), system.forcefield)


def read_lammps_trajectory(
    file: Path | str,
) -> "trajectory.lammps.LammpsTrajectoryReader":
    """Read LAMMPS trajectory file and return a trajectory reader."""
    from .trajectory.lammps import LammpsTrajectoryReader

    return LammpsTrajectoryReader(Path(file))


def write_lammps_trajectory(
    file: Path | str, frames: List[Frame], atom_style: str = "full"
) -> None:
    """Write frames to a LAMMPS trajectory file."""
    from .trajectory.lammps import LammpsTrajectoryWriter

    with LammpsTrajectoryWriter(file, atom_style) as writer:
        for i, frame in enumerate(frames):
            timestep = getattr(frame, "timestep", i)
            writer.write_frame(frame, timestep)
