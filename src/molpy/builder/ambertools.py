"""
AmberTools-based molecular builder for molpy.
Provides clean, composable API for AmberTools workflows using molq for orchestration.
"""

import logging
from pathlib import Path
from typing import Any, Generator

from molq.submitor import LocalSubmitor

from molpy.core.atomistic import Angle, Atom, Atomistic, Bond, Dihedral
from molpy.core.element import Element
from molpy.core.frame import Frame
from molpy.core.protocol import Entities

from .polymer import Monomer
from .reacter_lammps import ReactantWrapper

logger = logging.getLogger(__name__)


class AmberToolsComponent:
    """Base class for all AmberTools components."""

    def __init__(self, conda_env: str = "AmberTools25", workdir: Path | None = None):
        """
        Initialize AmberTools component.

        Args:
            conda_env: Conda environment name for AmberTools
            workdir: Working directory for temporary files
        """
        self.conda_env = conda_env
        self.workdir = (
            workdir
            if workdir
            else Path.cwd() / f"{self.__class__.__name__.lower()}_work"
        )
        self.workdir.mkdir(parents=True, exist_ok=True)


class Antechamber(AmberToolsComponent):
    """
    AnteChamber component for atom type assignment and charge calculation.

    Usage:
        antech = Antechamber(conda_env="AmberTools25")
        result = antech(struct, net_charge=0.0, charge_type="bcc")
    """

    def __call__(
        self,
        struct: Atomistic | Monomer | ReactantWrapper,
        net_charge: float = 0.0,
        forcefield: str = "gaff",
        charge_type: str = "bcc",
        output_format: str = "ac",
        workdir: Path | None = None,
        **kwargs,
    ) -> tuple[Path, Frame]:
        """
        Run AnteChamber on a molecular structure.

        Args:
            struct: Molecular structure to process
            net_charge: Net charge of the system
            forcefield: Force field to use (e.g., "gaff")
            charge_type: Charge calculation method
            output_format: Output format for AnteChamber
            workdir: Optional working directory (overrides default)
            **kwargs: Additional arguments

        Returns:
            Tuple of (ac_file_path, updated_frame)
        """
        # Use provided workdir or default
        if workdir is not None:
            workdir = Path(workdir)
        else:
            workdir = self.workdir

        # Generate unique name for this structure
        struct_name = getattr(struct, "name", None) or f"struct_{id(struct)}"
        workdir = workdir / struct_name
        workdir.mkdir(parents=True, exist_ok=True)

        # Ensure workdir is Path type for the rest of the method
        workdir = Path(workdir)

        # Prepare input files
        pdb_name = f"{struct_name}.pdb"
        ac_name = f"{struct_name}.ac"
        pdb_path = workdir / pdb_name
        ac_path = workdir / ac_name

        # Convert structure to PDB format
        import molpy.io as mp_io

        if hasattr(struct, "to_frame"):
            frame = struct.to_frame()
        else:
            frame = struct

        # Ensure we have a proper Frame object
        if not isinstance(frame, Frame):
            raise TypeError(f"Expected Frame object, got {type(frame)}")

        mp_io.write_pdb(pdb_path, frame)

        # Run AnteChamber if output doesn't exist
        if not ac_path.exists():
            submitor = LocalSubmitor("local", {})
            submitor.local_submit(
                job_name="antechamber",
                cmd=f"antechamber -i {pdb_name} -fi pdb -o {ac_name} -fo {output_format} "
                f"-an y -at {forcefield} -c {charge_type} -nc {net_charge}",
                cwd=workdir,
                conda_env=self.conda_env,
                block=True,
            )

        # Read results and update structure
        result_frame = mp_io.read_amber_ac(ac_path, frame=None)

        # Update atom properties in the original structure
        if hasattr(struct, "atoms"):
            atom_names = result_frame["atoms"]["name"]
            atom_types = result_frame["atoms"]["type"]
            atom_charges = result_frame["atoms"]["q"]

            for i, (atom, typ, q, name) in enumerate(
                zip(struct.atoms, atom_types, atom_charges, atom_names)
            ):
                atom["name"] = name.item()
                atom["type"] = typ.item()
                atom["q"] = q.item()

        # Generate monomer connectivity file if applicable
        if isinstance(struct, Monomer):
            self._write_monomer_connectivity(struct, workdir, struct_name)

        return ac_path, result_frame

    def read_ac(self, ac_path: Path | str) -> Frame:
        """
        Manually read an existing .ac file.

        Args:
            ac_path: Path to the .ac file

        Returns:
            Frame object with atom types and charges
        """
        import molpy.io as mp_io

        return mp_io.read_amber_ac(ac_path, frame=None)

    def _write_monomer_connectivity(self, monomer: Monomer, workdir: Path, name: str):
        """Write monomer connectivity file for TLeap."""
        mc_path = workdir / f"{name}.mc"

        def get_atom_name(atom_ref) -> str:
            """Convert atom reference to name."""
            if isinstance(atom_ref, int):
                return monomer.atoms[atom_ref]["name"]
            return str(atom_ref)

        with open(mc_path, "w") as f:
            # Process head anchor
            if hasattr(monomer.anchors, "get") and monomer.anchors.get("head"):
                head_anchor = monomer.anchors["head"]
                head_name = get_atom_name(head_anchor.anchor)
                f.write(f"HEAD_NAME {head_name}\n")

                if hasattr(head_anchor, "deletes"):
                    for delete in head_anchor.deletes:
                        f.write(f"OMIT_NAME {get_atom_name(delete)}\n")

            # Process tail anchor
            if hasattr(monomer.anchors, "get") and monomer.anchors.get("tail"):
                tail_anchor = monomer.anchors["tail"]
                tail_name = get_atom_name(tail_anchor.anchor)
                f.write(f"TAIL_NAME {tail_name}\n")

                if hasattr(tail_anchor, "deletes"):
                    for delete in tail_anchor.deletes:
                        f.write(f"OMIT_NAME {get_atom_name(delete)}\n")


class Prepgen(AmberToolsComponent):
    """
    PrepGen component for generating prepi files.

    Usage:
        prepgen = Prepgen(conda_env="AmberTools25")
        prepi_path = prepgen(ac_file_path)
    """

    def __call__(self, ac_path: str | Path, **kwargs) -> Path:
        """
        Run PrepGen on an AnteChamber output file.

        Args:
            ac_path: Path to .ac file from AnteChamber
            **kwargs: Additional arguments

        Returns:
            Path to generated .prepi file
        """
        ac_path = Path(ac_path)
        name = ac_path.stem
        workdir = self.workdir / name
        workdir.mkdir(parents=True, exist_ok=True)

        # Copy .ac file to workdir
        import shutil

        shutil.copy2(ac_path, workdir / f"{name}.ac")

        prepi_name = f"{name}.prepi"
        prepi_path = workdir / prepi_name

        if not prepi_path.exists():
            submitor = LocalSubmitor("local", {})
            submitor.local_submit(
                job_name="prepgen",
                cmd=f"prepgen -i {name}.ac -o {prepi_name} -f prepi",
                cwd=workdir,
                conda_env=self.conda_env,
                block=True,
            )

        return prepi_path


class Parmchk(AmberToolsComponent):
    """
    ParmChk component for generating frcmod files.

    Usage:
        parmchk = Parmchk(conda_env="AmberTools25")
        frcmod_path = parmchk(prepi_file_path)
    """

    def __call__(self, prepi_path: str | Path, **kwargs) -> Path:
        """
        Run ParmChk on a PrepGen output file.

        Args:
            prepi_path: Path to .prepi file from PrepGen
            **kwargs: Additional arguments

        Returns:
            Path to generated .frcmod file
        """
        prepi_path = Path(prepi_path)
        name = prepi_path.stem
        workdir = self.workdir / name
        workdir.mkdir(parents=True, exist_ok=True)

        # Copy .prepi file to workdir
        import shutil

        shutil.copy2(prepi_path, workdir / f"{name}.prepi")

        frcmod_name = f"{name}.frcmod"
        frcmod_path = workdir / frcmod_name

        if not frcmod_path.exists():
            submitor = LocalSubmitor("local", {})
            submitor.local_submit(
                job_name="parmchk",
                cmd=f"parmchk -i {name}.prepi -f {frcmod_name}",
                cwd=workdir,
                conda_env=self.conda_env,
                block=True,
            )

        return frcmod_path


class TLeap(AmberToolsComponent):
    """
    TLeap component for generating parameter files and coordinates.

    Usage:
        tleap = TLeap(conda_env="AmberTools25")
        prmtop, inpcrd = tleap(name, [prepi_files], [frcmod_files])
    """

    def __call__(
        self,
        name: str,
        prepi_files: list[str | Path] | None = None,
        frcmod_files: list[str | Path] | None = None,
        **kwargs,
    ) -> tuple[Path, Path]:
        """
        Run TLeap to generate parameter files.

        Args:
            name: Name for the output system
            prepi_files: List of .prepi files to load
            frcmod_files: List of .frcmod files to load
            **kwargs: Additional arguments

        Returns:
            Tuple of (prmtop_path, inpcrd_path)
        """
        workdir = self.workdir / name
        workdir.mkdir(parents=True, exist_ok=True)

        # Build leap script
        leap_script = self._build_leap_script(name, prepi_files, frcmod_files)

        # Write leap input file
        leap_in = workdir / f"{name}.in"
        with open(leap_in, "w") as f:
            f.write(leap_script)

        # Run TLeap
        submitor = LocalSubmitor("local", {})
        submitor.local_submit(
            job_name="tleap",
            cmd=f"tleap -f {name}.in",
            cwd=workdir,
            conda_env=self.conda_env,
            block=True,
        )

        return workdir / f"{name}.prmtop", workdir / f"{name}.inpcrd"

    def _build_leap_script(
        self,
        name: str,
        prepi_files: list[str | Path] | None = None,
        frcmod_files: list[str | Path] | None = None,
    ) -> str:
        """Build the TLeap script content."""
        commands = []

        # Source force field libraries
        commands.append("source leaprc.gaff")

        # Load prepi files
        if prepi_files:
            for prepi_file in prepi_files:
                prepi_path = Path(prepi_file)
                commands.append(f"loadprep {prepi_path}")

        # Load frcmod files
        if frcmod_files:
            for frcmod_file in frcmod_files:
                frcmod_path = Path(frcmod_file)
                commands.append(f"loadfrcmod {frcmod_path}")

        # Combine molecules if multiple
        if prepi_files and len(prepi_files) > 1:
            mol_names = [Path(f).stem for f in prepi_files]
            commands.append(f"combine {name} {' '.join(mol_names)}")
        elif prepi_files:
            mol_name = Path(prepi_files[0]).stem
            commands.append(f"combine {name} {mol_name}")

        # Save output files
        commands.append(f"saveamberparm {name} {name}.prmtop {name}.inpcrd")
        commands.append(f"savepdb {name} {name}.pdb")
        commands.append("quit")

        return "\n".join(commands)


class AmberToolsTypifier(AmberToolsComponent):
    """
    Typifier component for automatic atom type assignment using AmberTools.

    Usage:
        typifier = AmberToolsTypifier(conda_env="AmberTools25")
        typed_struct = typifier(struct, forcefield="gaff")
    """

    def __init__(self, conda_env: str = "AmberTools25", workdir: Path | None = None):
        """Initialize AmberTools typifier."""
        super().__init__(conda_env, workdir)

        # Initialize components
        self.antechamber = Antechamber(conda_env, self.workdir)
        self.prepgen = Prepgen(conda_env, self.workdir)
        self.parmchk = Parmchk(conda_env, self.workdir)
        self.tleap = TLeap(conda_env, self.workdir)

    def __call__(
        self,
        struct: Atomistic,
        forcefield: str = "gaff",
        charge_type: str = "bcc",
        net_charge: float = 0.0,
        is_frcmod: bool = True,
        is_prepi: bool = True,
        is_tleap: bool = True,
        **kwargs,
    ) -> Atomistic:
        """
        Typify a molecular structure using AmberTools.

        Args:
            struct: Structure to typify
            forcefield: Force field to use
            charge_type: Charge calculation method
            net_charge: Net charge of the system
            is_frcmod: Whether to generate frcmod file
            is_prepi: Whether to generate prepi file
            is_tleap: Whether to run TLeap
            **kwargs: Additional arguments

        Returns:
            Typified structure
        """
        name = struct.get("name")
        if name is None:
            raise ValueError("Struct must have a name attribute")

        workdir = Path(self.workdir) / name
        workdir.mkdir(parents=True, exist_ok=True)

        # Run AnteChamber
        ac_path, result_frame = self.antechamber(
            struct,
            net_charge=net_charge,
            forcefield=forcefield,
            charge_type=charge_type,
        )
        logger.info(f"AnteChamber step completed for {name}")

        # Run PrepGen if requested
        prepi_path = None
        frcmod_path = None

        if is_prepi:
            prepi_path = self.prepgen(ac_path)
            logger.info(f"PrepGen step completed for {name}")

        # Run ParmChk if requested
        if is_frcmod:
            if prepi_path is None:
                raise ValueError("PrepGen must be run before ParmChk")
            frcmod_path = self.parmchk(prepi_path)
            logger.info(f"ParmChk step completed for {name}")

        # Run TLeap if requested
        if is_tleap:
            if prepi_path is None or frcmod_path is None:
                raise ValueError("PrepGen and ParmChk must be run before TLeap")
            prmtop_path, inpcrd_path = self.tleap(name, [prepi_path], [frcmod_path])

            # Read final structure and update original
            import molpy.io as mp_io

            frame, forcefield = mp_io.read_amber(prmtop_path, inpcrd_path)

            # Update structure with new topology
            self._update_structure_topology(struct, frame)

        else:
            # Just read type information from .ac file
            import molpy.io as mp_io

            frame = mp_io.read_amber_ac(ac_path, frame=None)
            atom_types = frame["atoms"]["type"]
            atom_charges = frame["atoms"]["q"]

            for satom, typ, q in zip(struct.atoms, atom_types, atom_charges):
                satom["type"] = typ.item()
                satom["q"] = q.item()

        return struct

    def _update_structure_topology(self, struct: Atomistic, frame: Frame):
        """Update structure with topology information from TLeap output."""
        # Update atoms
        atoms = []
        atom_ids = frame["atoms"]["id"].tolist()
        atom_names = frame["atoms"]["name"].tolist()
        atom_types = frame["atoms"]["type"].tolist()
        atom_charges = frame["atoms"]["q"].tolist()
        atom_masses = frame["atoms"]["mass"].tolist()
        atom_numbers = frame["atoms"]["atomic_number"].tolist()
        xyz = frame["atoms"]["xyz"].tolist()

        for i in range(len(frame["atoms"]["id"])):
            atoms.append(
                Atom(
                    id=atom_ids[i],
                    name=atom_names[i],
                    type=atom_types[i],
                    element=Element(atom_numbers[i]).symbol,
                    q=atom_charges[i],
                    mass=atom_masses[i],
                    xyz=xyz[i],
                )
            )
        struct["atoms"] = Entities(atoms)

        # Update bonds if present
        if "bonds" in frame:
            bonds = []
            bond_i = frame["bonds"]["i"].tolist()
            bond_j = frame["bonds"]["j"].tolist()
            bond_types = frame["bonds"]["type"].tolist()
            bond_ids = frame["bonds"]["id"].tolist()

            for i in range(len(bond_ids)):
                bonds.append(
                    Bond(
                        atoms[bond_i[i]],
                        atoms[bond_j[i]],
                        type=bond_types[i],
                        id=bond_ids[i],
                    )
                )
            struct["bonds"] = Entities(bonds)

        # Update angles if present
        if "angles" in frame:
            angles = []
            angle_i = frame["angles"]["i"].tolist()
            angle_j = frame["angles"]["j"].tolist()
            angle_k = frame["angles"]["k"].tolist()
            angle_types = frame["angles"]["type"].tolist()
            angle_ids = frame["angles"]["id"].tolist()

            for i in range(len(angle_ids)):
                angles.append(
                    Angle(
                        atoms[angle_i[i]],
                        atoms[angle_j[i]],
                        atoms[angle_k[i]],
                        type=angle_types[i],
                        id=angle_ids[i],
                    )
                )
            struct["angles"] = Entities(angles)

        # Update dihedrals if present
        if "dihedrals" in frame:
            dihedrals = []
            dihedral_i = frame["dihedrals"]["i"].tolist()
            dihedral_j = frame["dihedrals"]["j"].tolist()
            dihedral_k = frame["dihedrals"]["k"].tolist()
            dihedral_l = frame["dihedrals"]["l"].tolist()
            dihedral_types = frame["dihedrals"]["type"].tolist()
            dihedral_ids = frame["dihedrals"]["id"].tolist()

            for i in range(len(dihedral_ids)):
                dihedrals.append(
                    Dihedral(
                        atoms[dihedral_i[i]],
                        atoms[dihedral_j[i]],
                        atoms[dihedral_k[i]],
                        atoms[dihedral_l[i]],
                        type=dihedral_types[i],
                        id=dihedral_ids[i],
                    )
                )
            struct["dihedrals"] = Entities(dihedrals)
