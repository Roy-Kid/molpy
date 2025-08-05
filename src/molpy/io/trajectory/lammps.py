from io import StringIO
from pathlib import Path
from typing import List, Sequence, Union, Optional
import re
import numpy as np

from molpy.core import Trajectory, Frame, Block, Box

from .base import TrajectoryReader, TrajectoryWriter


class LammpsTrajectoryReader(TrajectoryReader):
    """Reader for LAMMPS trajectory files, supporting multiple files."""

    def __init__(self, fpath: str | Path, trajectory: Trajectory | None = None):
        # Convert fpaths to the expected format
        fpath_converted = Path(fpath)
        if trajectory is None:
            trajectory = Trajectory()
        super().__init__(trajectory, fpath_converted)

    def read_frame(self, index: int) -> "Frame":
        """Read a specific frame from the trajectory."""
        frame_start = self._byte_offsets[index]
        frame_end = self._byte_offsets[index + 1] if index + 1 < len(self._byte_offsets) else None
        frame_bytes = self._mm[frame_start:frame_end]
        frame_lines = frame_bytes.decode().splitlines()
        return self._parse_frame(frame_lines)

    def _parse_trajectory(self):
        """Parse trajectory files and update frame count."""

        # Use regex to find all ITEMs and ATOMS blocks, and compute frame offsets accordingly
        item_re = re.compile(rb"^ITEM:.*", re.MULTILINE)
        atoms_re = re.compile(rb"^ITEM:\s+ATOMS\b")

        # Find all ITEMs and their byte offsets
        item_matches = list(item_re.finditer(self._mm))
        item_offsets = [m.start() for m in item_matches]

        # Find indices of ITEM: ATOMS entries (index into item_offsets)
        atoms_indices = [
            i for i, m in enumerate(item_matches)
            if atoms_re.match(m.group(0))
        ]

        # Compute frame starts: first ITEM, then every ATOMS+1
        frame_offsets = [item_offsets[0]] + [
            item_offsets[i + 1] for i in atoms_indices[:-1]
        ]

        self._byte_offsets = frame_offsets
        self._total_frames = len(self._byte_offsets)

    def _parse_frame(self, frame_lines: Sequence[str]) -> Frame:
        """Parse frame lines into a Frame object."""
        
        header = []
        box_bounds = []

        timestep = None
        data_start = None
        header = []
        for i, line in enumerate(frame_lines):
            if line.startswith("ITEM: TIMESTEP"):
                # The timestep value is on the next line
                timestep = int(frame_lines[i + 1].strip())
            elif line.startswith("ITEM: BOX BOUNDS"):
                periodic = line.split()[-3:]
                for j in range(3):
                    box_bounds.append(
                        list(map(float, frame_lines[i + j + 1].strip().split()))
                    )
            elif line.startswith("ITEM: ATOMS"):
                header = line.split()[2:]
                data_start = i + 1
                break

        # Check if we found atom data
        if not header:
            raise ValueError("No atom data found in trajectory frame")
        
        box_bounds = np.array(box_bounds)

        if box_bounds.shape == (3, 2):
            box_matrix = np.array(
                [
                    [box_bounds[0, 1] - box_bounds[0, 0], 0, 0],
                    [0, box_bounds[1, 1] - box_bounds[1, 0], 0],
                    [0, 0, box_bounds[2, 1] - box_bounds[2, 0]],
                ]
            )
            origin = np.array([box_bounds[0, 0], box_bounds[1, 0], box_bounds[2, 0]])
        elif box_bounds.shape == (3, 3):
            xy, xz, yz = box_bounds[:, 2]
            box_matrix = np.array(
                [
                    [box_bounds[0, 1] - box_bounds[0, 0], xy, xz],
                    [0, box_bounds[1, 1] - box_bounds[1, 0], yz],
                    [0, 0, box_bounds[2, 1] - box_bounds[2, 0]],
                ]
            )
            origin = np.array([box_bounds[0, 0], box_bounds[1, 0], box_bounds[2, 0]])
        else:
            raise ValueError(f"Invalid box bounds shape {box_bounds.shape}")

        box = Box(matrix=box_matrix, origin=origin)
        
        # Create frame with proper structure
        frame = Frame(box=box, timestep=timestep)
        frame["atoms"] = Block.from_csv(StringIO("\n".join(frame_lines[data_start:])), header=header, delimiter=" ")
        
        return frame


class LammpsTrajectoryWriter(TrajectoryWriter):
    """Writer for LAMMPS trajectory files (dump format)."""

    def __init__(self, fpath: Union[str, Path], atom_style: str = "full"):
        super().__init__(fpath)
        self.atom_style = atom_style
        self._frame_count = 0

    def write_frame(self, frame: "Frame", timestep: Optional[int] = None):
        """Write a single frame to the trajectory file.
        
        Args:
            frame: Frame object containing atoms and box information
            timestep: Timestep number (if None, uses auto-increment)
        """
        if timestep is None:
            # Try to get timestep from frame metadata first
            if hasattr(frame, 'metadata') and frame.metadata and 'timestep' in frame.metadata:
                timestep = frame.metadata['timestep']
            else:
                timestep = self._frame_count
        
        self._frame_count += 1
        
        # Get atoms data
        if "atoms" not in frame:
            raise ValueError("Frame must contain atoms data")
        
        atoms = frame["atoms"]
        
        # Get number of atoms
        n_atoms = 0
        if hasattr(atoms, '_vars'):
            # Block object - get count from first variable
            if atoms._vars:
                first_key = next(iter(atoms._vars))
                n_atoms = len(atoms._vars[first_key])
        else:
            # Dict-like - count entries in first field
            if atoms and isinstance(atoms, dict):
                first_key = next(iter(atoms))
                n_atoms = len(atoms[first_key])
        
        # Write timestep header
        self._fp.write(f"ITEM: TIMESTEP\n{timestep}\n".encode())
        
        # Write number of atoms
        self._fp.write(f"ITEM: NUMBER OF ATOMS\n{n_atoms}\n".encode())
        
        # Write box bounds
        box = getattr(frame, 'box', None)
        if box and hasattr(box, 'matrix'):
            matrix = box.matrix
            origin = getattr(box, 'origin', np.zeros(3))
            
            # For orthogonal box
            if np.allclose(matrix, np.diag(np.diag(matrix))):
                xlo, ylo, zlo = origin
                xhi = xlo + matrix[0, 0]
                yhi = ylo + matrix[1, 1]
                zhi = zlo + matrix[2, 2]
                
                self._fp.write(f"ITEM: BOX BOUNDS pp pp pp\n".encode())
                self._fp.write(f"{xlo:.6f} {xhi:.6f}\n".encode())
                self._fp.write(f"{ylo:.6f} {yhi:.6f}\n".encode())
                self._fp.write(f"{zlo:.6f} {zhi:.6f}\n".encode())
            else:
                # For triclinic box
                xlo, ylo, zlo = origin
                xhi = xlo + matrix[0, 0]
                yhi = ylo + matrix[1, 1]
                zhi = zlo + matrix[2, 2]
                xy = matrix[0, 1]
                xz = matrix[0, 2]
                yz = matrix[1, 2]
                
                self._fp.write(f"ITEM: BOX BOUNDS xy xz yz pp pp pp\n".encode())
                self._fp.write(f"{xlo:.6f} {xhi:.6f} {xy:.6f}\n".encode())
                self._fp.write(f"{ylo:.6f} {yhi:.6f} {xz:.6f}\n".encode())
                self._fp.write(f"{zlo:.6f} {zhi:.6f} {yz:.6f}\n".encode())
        else:
            # Default box if none provided
            self._fp.write(f"ITEM: BOX BOUNDS pp pp pp\n".encode())
            self._fp.write(f"0.000000 10.000000\n".encode())
            self._fp.write(f"0.000000 10.000000\n".encode())
            self._fp.write(f"0.000000 10.000000\n".encode())
        
        # Determine atom columns to write
        if self.atom_style == "full":
            cols = ["id", "mol", "type", "q", "x", "y", "z"]
        else:
            # Default atomic style
            cols = ["id", "type", "x", "y", "z"]
        
        # Map frame fields to dump columns
        field_mapping = {
            "mol": "molid",  # molid in frame -> mol in dump
            "q": "q",        # charge field
            "x": "x", "y": "y", "z": "z"  # coordinates (try individual first)
        }
        
        # Write atom header
        self._fp.write(f"ITEM: ATOMS {' '.join(cols)}\n".encode())
        
        # Prepare atom data
        atom_data = {}
        coords = None
        
        # Handle coordinates - try individual x,y,z first, then xyz array
        if all(coord in atoms for coord in ["x", "y", "z"]):
            atom_data["x"] = atoms["x"]
            atom_data["y"] = atoms["y"]
            atom_data["z"] = atoms["z"]
        elif "xyz" in atoms:
            coords = atoms["xyz"]
            coords = np.asarray(coords)
            if coords.ndim == 2 and coords.shape[1] == 3:
                atom_data["x"] = coords[:, 0]
                atom_data["y"] = coords[:, 1]
                atom_data["z"] = coords[:, 2]
            else:
                raise ValueError("xyz coordinates must be Nx3 array")
        else:
            raise ValueError("No coordinate data found in atoms")
        
        # Handle other fields
        for col in cols:
            if col in ["x", "y", "z"]:
                continue  # Already handled
            
            field_name = field_mapping.get(col, col)
            if field_name in atoms:
                atom_data[col] = atoms[field_name]
            else:
                # Provide defaults
                if col == "id":
                    atom_data[col] = np.arange(1, n_atoms + 1)
                elif col == "mol":
                    atom_data[col] = np.ones(n_atoms, dtype=int)
                elif col == "type":
                    atom_data[col] = np.ones(n_atoms, dtype=int)
                elif col == "q":
                    atom_data[col] = np.zeros(n_atoms)
        
        # Write atom data
        for i in range(n_atoms):
            line_parts = []
            for col in cols:
                value = atom_data[col][i]
                if col in ["id", "mol"]:
                    line_parts.append(f"{int(value)}")
                elif col == "type":
                    # Handle both string and numeric types
                    if isinstance(value, str):
                        # Create a simple mapping for common atom types
                        type_map = {'H': 1, 'C': 2, 'N': 3, 'O': 4, 'S': 5, 'P': 6}
                        line_parts.append(f"{type_map.get(value, 1)}")
                    else:
                        line_parts.append(f"{int(value)}")
                else:
                    line_parts.append(f"{float(value):.6f}")
            
            line = " ".join(line_parts) + "\n"
            self._fp.write(line.encode())
