"""
Packmol-based molecular packing for molpy.
Provides clean, composable API for molecular packing workflows using molq for orchestration.
"""

import shutil
import tempfile
from pathlib import Path
from typing import Any, Generator

import molq
import numpy as np

import molpy.pack as mpk
from molpy.core.frame import Frame

from .base import Packer


class PackmolComponent:
    """Base class for all Packmol components."""

    def __init__(self, executable: str = "packmol", workdir: Path | None = None):
        """
        Initialize Packmol component.

        Args:
            executable: Path to packmol executable
            workdir: Working directory for temporary files
        """
        self.executable = executable
        self.workdir = workdir if workdir else Path.cwd() / "packmol_work"
        self.workdir.mkdir(parents=True, exist_ok=True)

        # Check if packmol is available
        self.packmol_path = shutil.which(self.executable)
        if self.packmol_path is None:
            raise FileNotFoundError(
                f"Packmol executable '{self.executable}' not found in PATH"
            )


def map_region_to_packmol_definition(constraint):
    """Convert molpy constraint objects to packmol input format."""

    def box(region, not_flag):
        origin = region.origin
        upper = region.upper.tolist()
        flag = "outside" if not_flag else "inside"
        return f"{flag} box {origin[0]} {origin[1]} {origin[2]} {upper[0]} {upper[1]} {upper[2]}"

    def sphere(region, not_flag):
        origin = region.origin
        flag = "outside" if not_flag else "inside"
        return f"{flag} sphere {origin[0]} {origin[1]} {origin[2]} {region.radius}"

    if isinstance(constraint, mpk.AndConstraint):
        r1_cmd = map_region_to_packmol_definition(constraint.a)
        r2_cmd = map_region_to_packmol_definition(constraint.b)
        return f"  {r1_cmd} \n  {r2_cmd}"

    if isinstance(constraint, mpk.InsideBoxConstraint):
        return box(constraint.region, False)
    elif isinstance(constraint, mpk.InsideSphereConstraint):
        return sphere(constraint.region, False)
    elif isinstance(constraint, mpk.OutsideBoxConstraint):
        return box(constraint.region, True)
    elif isinstance(constraint, mpk.OutsideSphereConstraint):
        return sphere(constraint.region, True)
    elif isinstance(constraint, mpk.MinDistanceConstraint):
        pass
    else:
        raise NotImplementedError(
            f"Packmol does not support constraint type {type(constraint)}"
        )


class Packmol(Packer, PackmolComponent):
    """
    Packmol component for molecular packing.

    Usage:
        packmol = Packmol(executable="packmol")
        result = packmol(targets, max_steps=1000, seed=4628)
    """

    def __init__(self, executable: str = "packmol", workdir: Path | None = None):
        """Initialize Packmol packer."""
        Packer.__init__(self)
        PackmolComponent.__init__(self, executable, workdir)

        self.intermediate_files = [
            ".optimized.pdb",
            ".packmol.inp",
            ".packmol.out",
        ]

    def __call__(
        self,
        targets: list[mpk.Target] | None = None,
        max_steps: int = 1000,
        seed: int | None = None,
        workdir: Path | None = None,
        **kwargs,
    ) -> Frame:
        """
        Pack molecules using Packmol.

        Args:
            targets: List of packing targets
            max_steps: Maximum optimization steps
            seed: Random seed for packing
            workdir: Optional working directory (overrides default)
            **kwargs: Additional arguments

        Returns:
            Packed molecular system as Frame
        """
        # Use provided workdir or default
        if workdir is not None:
            workdir = Path(workdir)
            workdir.mkdir(parents=True, exist_ok=True)
        else:
            workdir = self.workdir

        # Use provided targets or stored targets
        if targets is None:
            targets = self.targets

        if seed is None:
            seed = 4628

        # Generate packmol input
        self._generate_input(targets, max_steps, seed, workdir)

        # Run packmol
        submitor = molq.LocalSubmitor("local", {})
        submitor.local_submit(
            job_name="packmol",
            cmd="packmol < .packmol.inp > .packmol.out",
            cwd=workdir,
            block=True,
        )

        # Read results
        optimized_frame = self._read_packmol_output(workdir)

        # Clean up intermediate files
        self._cleanup_intermediate_files(workdir)

        # Process and return results
        return self._build_final_frame(targets, optimized_frame, workdir)

    def _generate_input(
        self, targets: list[mpk.Target], max_steps: int, seed: int, workdir: Path
    ):
        """Generate Packmol input file."""
        lines = []
        lines.append("tolerance 2.0")
        lines.append("filetype pdb")
        lines.append("output .optimized.pdb")
        lines.append(f"seed {seed}")

        for target in targets:
            frame = target.frame
            number = target.number
            constraint = target.constraint

            # Create temporary PDB file
            tmpfile = tempfile.NamedTemporaryFile(delete=False, suffix=".pdb")
            import molpy.io as mp_io

            mp_io.write_pdb(Path(tmpfile.name), frame)

            lines.append(f"structure {tmpfile.name}")
            lines.append(f"  number {number}")
            lines.append(map_region_to_packmol_definition(constraint))
            lines.append("end structure")

        # Write input file
        input_file = workdir / ".packmol.inp"
        with open(input_file, "w") as f:
            f.write("\n".join(lines))

    def _read_packmol_output(self, workdir: Path) -> Frame:
        """Read Packmol output file."""
        import molpy.io as mp_io

        output_file = workdir / ".optimized.pdb"
        return mp_io.read_pdb(output_file)

    def _cleanup_intermediate_files(self, workdir: Path):
        """Clean up intermediate files."""
        for file in self.intermediate_files:
            file_path = workdir / file
            if file_path.exists():
                file_path.unlink()

    def _build_final_frame(
        self, targets: list[mpk.Target], optimized_frame: Frame, workdir: Path
    ) -> Frame:
        """Build final frame from packing results."""
        # Count atoms per instance
        target_atoms_count = []
        for target in targets:
            atoms_block = target.frame["atoms"]
            n_atoms = len(atoms_block.get("id", atoms_block.get("xyz", [])))
            for _ in range(target.number):
                target_atoms_count.append(n_atoms)

        # Compute cumulative offsets
        atom_offsets = np.concatenate([[0], np.cumsum(target_atoms_count)])

        # Prepare lists for each block
        all_atoms = []
        all_bonds = []
        all_angles = []
        all_dihedrals = []

        current_instance = 0
        bond_id_counter = 0
        angle_id_counter = 0
        dihedral_id_counter = 0

        for target in targets:
            for inst in range(target.number):
                offset = atom_offsets[current_instance]
                frame_dict = {}

                # Process atoms block
                atoms = target.frame["atoms"].copy()
                n = len(atoms.get("id", atoms.get("xyz", [])))

                atoms["mol"] = np.full(n, current_instance + 1, dtype=int)

                # Reassign atom IDs sequentially
                if "id" in atoms:
                    atoms["id"] = np.arange(offset + 1, offset + n + 1, dtype=int)

                frame_dict["atoms"] = atoms

                # Process bonds
                if "bonds" in target.frame:
                    bonds = target.frame["bonds"].copy()
                    m = len(bonds.get("id", bonds.get("i", [])))
                    # IDs
                    if "id" in bonds:
                        bonds["id"] = np.arange(
                            bond_id_counter + 1, bond_id_counter + m + 1, dtype=int
                        )
                        bond_id_counter += m
                    # Endpoints
                    for end in ("i", "j"):
                        if end in bonds:
                            bonds[end] = bonds[end] + offset
                    frame_dict["bonds"] = bonds

                # Process angles
                if "angles" in target.frame:
                    angles = target.frame["angles"].copy()
                    p = len(angles.get("id", angles.get("i", [])))
                    if "id" in angles:
                        angles["id"] = np.arange(
                            angle_id_counter + 1, angle_id_counter + p + 1, dtype=int
                        )
                        angle_id_counter += p
                    for end in ("i", "j", "k"):
                        if end in angles:
                            angles[end] = angles[end] + offset
                    frame_dict["angles"] = angles

                # Process dihedrals
                if "dihedrals" in target.frame:
                    dihedrals = target.frame["dihedrals"].copy()
                    q = len(dihedrals.get("id", dihedrals.get("i", [])))
                    if "id" in dihedrals:
                        dihedrals["id"] = np.arange(
                            dihedral_id_counter + 1,
                            dihedral_id_counter + q + 1,
                            dtype=int,
                        )
                        dihedral_id_counter += q
                    for end in ("i", "j", "k", "l"):
                        if end in dihedrals:
                            dihedrals[end] = dihedrals[end] + offset
                    frame_dict["dihedrals"] = dihedrals

                # Collect instance frame
                all_atoms.append(frame_dict["atoms"])
                if "bonds" in frame_dict:
                    all_bonds.append(frame_dict["bonds"])
                if "angles" in frame_dict:
                    all_angles.append(frame_dict["angles"])
                if "dihedrals" in frame_dict:
                    all_dihedrals.append(frame_dict["dihedrals"])

                current_instance += 1

        # Helper to concatenate dicts of arrays
        def concat_blocks(blocks):
            if not blocks:
                return {}
            keys = blocks[0].keys()
            combined = {}
            for key in keys:
                combined[key] = np.concatenate([blk[key] for blk in blocks], axis=0)
            return combined

        # Build final frame
        final_frame = Frame()
        final_frame["atoms"] = concat_blocks(all_atoms)
        # Overwrite coordinates with optimized positions
        if "xyz" in optimized_frame["atoms"]:
            final_frame["atoms"]["xyz"] = optimized_frame["atoms"]["xyz"]

        final_frame["bonds"] = concat_blocks(all_bonds)
        final_frame["angles"] = concat_blocks(all_angles)
        final_frame["dihedrals"] = concat_blocks(all_dihedrals)

        return final_frame

    def read_packmol_output(self, workdir: Path | str) -> Frame:
        """
        Manually read an existing Packmol output file.

        Args:
            workdir: Working directory containing Packmol output

        Returns:
            Frame object with packed molecular system
        """
        workdir = Path(workdir)
        return self._read_packmol_output(workdir)

    def generate_input_only(
        self,
        targets: list[mpk.Target],
        max_steps: int = 1000,
        seed: int = 4628,
        workdir: Path | None = None,
    ) -> Path:
        """
        Generate Packmol input file without running the packing.

        Args:
            targets: List of packing targets
            max_steps: Maximum optimization steps
            seed: Random seed for packing
            workdir: Optional working directory

        Returns:
            Path to generated input file
        """
        if workdir is None:
            workdir = self.workdir

        workdir = Path(workdir)
        workdir.mkdir(parents=True, exist_ok=True)

        self._generate_input(targets, max_steps, seed, workdir)
        return workdir / ".packmol.inp"

    def pack(
        self,
        targets: list[mpk.Target] | None = None,
        max_steps: int = 1000,
        seed: int | None = None,
    ) -> Frame:
        """Pack molecules using Packmol.
        
        This method implements the abstract pack() method from Packer base class.
        It delegates to __call__() for the actual implementation.
        """
        return self(targets, max_steps=max_steps, seed=seed)
