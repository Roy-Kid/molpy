import re
from io import StringIO
from pathlib import Path
from typing import List, Optional, Union

import numpy as np

from molpy.core import Block, Box, Frame

from .base import BaseTrajectoryReader, FrameLocation, TrajectoryWriter


class LammpsTrajectoryReader(BaseTrajectoryReader):
    """Reader for LAMMPS trajectory files, supporting multiple files."""

    def __init__(self, fpath: Union[str, Path, List[str], List[Path]]):
        # Pass the fpath (single or multiple) to the base class
        super().__init__(fpath)

    def _parse_trajectory(self, file_index: int):
        """Parse trajectory file at given index and update frame count."""
        mm = self._get_mmap(file_index)

        # Use regex to find all ITEMs and ATOMS blocks, and compute frame offsets accordingly
        item_re = re.compile(rb"^ITEM:.*", re.MULTILINE)
        atoms_re = re.compile(rb"^ITEM:\s+ATOMS\b")

        # Find all ITEMs and their byte offsets
        item_matches = list(item_re.finditer(mm))
        item_offsets = [m.start() for m in item_matches]

        # Find indices of ITEM: ATOMS entries (index into item_offsets)
        atoms_indices = [
            i for i, m in enumerate(item_matches) if atoms_re.match(m.group(0))
        ]

        # Compute frame starts: first ITEM, then every ATOMS+1
        if item_offsets:
            frame_offsets = [item_offsets[0]] + [
                item_offsets[i + 1] for i in atoms_indices[:-1]
            ]
        else:
            frame_offsets = []

        # Add frame locations to the global list
        for offset in frame_offsets:
            location = FrameLocation(
                file_index=file_index,
                byte_offset=offset,
                file_path=self.fpaths[file_index],
            )
            self._frame_locations.append(location)

        self._total_frames += len(frame_offsets)

    def _parse_frame(self, frame_lines: List[str]) -> Frame:
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
            raise ValueError(f"Invalid box bounds shape: {box_bounds.shape}")

        # Create box
        box = Box(matrix=box_matrix, origin=origin)

        frame = Frame(timestep=timestep)
        frame["atoms"] = Block.from_csv(
            StringIO("\n".join(frame_lines[data_start:])), header=header, delimiter=" "
        )

        # 优化后的列类型映射 - 更加智能的类型推断
        int_columns = ["id", "type", "mol"]
        str_columns = []
        field_type_mappings = {c: int for c in int_columns} | {
            c: str for c in str_columns
        }

        for col_name in header:
            if col_name in field_type_mappings:
                try:
                    frame["atoms"][col_name] = frame["atoms"][col_name].astype(
                        field_type_mappings[col_name]
                    )
                except (ValueError, TypeError):
                    # 如果转换失败，保持原始类型
                    pass
            else:
                # 对于其他列，尝试转换为浮点数，失败则保持原样
                try:
                    frame["atoms"][col_name] = frame["atoms"][col_name].astype(float)
                except (ValueError, TypeError):
                    # 保持原始类型（可能是字符串）
                    pass

        frame.metadata["box"] = box
        return frame


class LammpsTrajectoryWriter(TrajectoryWriter):
    """Writer for LAMMPS trajectory files."""

    def __init__(self, fpath: Union[str, Path], atom_style: str = "full"):
        super().__init__(fpath)
        self.atom_style = atom_style

    def write_frame(self, frame: "Frame", timestep: Optional[int] = None):
        """Write a single frame to the file."""
        if timestep is None:
            timestep = frame.metadata.get("timestep", 0)

        # Write timestep
        self._fp.write(f"ITEM: TIMESTEP\n{timestep}\n".encode())

        # Write number of atoms
        if "atoms" in frame:
            atoms = frame["atoms"]
            first_col = next(iter(atoms.keys()))
            n_atoms = len(atoms[first_col])
            self._fp.write(f"ITEM: NUMBER OF ATOMS\n{n_atoms}\n".encode())

        # Write box bounds
        # Get box from metadata
        box = frame.metadata.get("box")
        if box:
            matrix = box.matrix
            origin = box.origin

            # Check if box is orthogonal
            if np.allclose(matrix, np.diag(np.diag(matrix))):
                # Orthogonal box
                self._fp.write("ITEM: BOX BOUNDS pp pp pp\n".encode())
                for i in range(3):
                    self._fp.write(f"{origin[i]} {origin[i] + matrix[i, i]}\n".encode())
            else:
                # Triclinic box
                self._fp.write("ITEM: BOX BOUNDS pp pp pp xy xz yz\n".encode())
                for i in range(3):
                    if i == 0:
                        self._fp.write(
                            f"{origin[i]} {origin[i] + matrix[i, i]} {matrix[0, 1]}\n".encode()
                        )
                    elif i == 1:
                        self._fp.write(
                            f"{origin[i]} {origin[i] + matrix[i, i]} {matrix[0, 2]}\n".encode()
                        )
                    else:
                        self._fp.write(
                            f"{origin[i]} {origin[i] + matrix[i, i]} {matrix[1, 2]}\n".encode()
                        )

        # Write atoms
        if "atoms" in frame:
            atoms = frame["atoms"]

            # Determine column order based on available data
            columns = []
            if "id" in atoms:
                columns.append("id")
            if "mol_id" in atoms:
                columns.append("mol_id")
            if "type" in atoms:
                columns.append("type")
            if "q" in atoms:
                columns.append("q")
            if "x" in atoms and "y" in atoms and "z" in atoms:
                columns.extend(["x", "y", "z"])
            elif "xu" in atoms and "yu" in atoms and "zu" in atoms:
                columns.extend(["xu", "yu", "zu"])
            elif "xs" in atoms and "ys" in atoms and "zs" in atoms:
                columns.extend(["xs", "ys", "zs"])
            if "vx" in atoms and "vy" in atoms and "vz" in atoms:
                columns.extend(["vx", "vy", "vz"])
            if "fx" in atoms and "fy" in atoms and "fz" in atoms:
                columns.extend(["fx", "fy", "fz"])

            # Write atom header
            self._fp.write(f"ITEM: ATOMS {' '.join(columns)}\n".encode())

            # Write atom data
            n_atoms = len(atoms)

            # 获取第一个可用列来确定实际原子数量
            first_col = next(iter(atoms.keys()))
            actual_n_atoms = len(atoms[first_col])

            for i in range(actual_n_atoms):
                row_data = []
                for col in columns:
                    if col in [
                        "x",
                        "y",
                        "z",
                        "xu",
                        "yu",
                        "zu",
                        "xs",
                        "ys",
                        "zs",
                        "vx",
                        "vy",
                        "vz",
                        "fx",
                        "fy",
                        "fz",
                        "q",
                    ]:
                        row_data.append(f"{atoms[col][i]:.6f}")
                    else:
                        row_data.append(f"{atoms[col][i]}")
                self._fp.write(f"{' '.join(row_data)}\n".encode())

        self._fp.flush()
