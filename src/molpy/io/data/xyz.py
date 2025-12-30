import shlex

import numpy as np

from molpy.core import Block, Box, Frame
from molpy.core.element import Element

from .base import DataReader


class XYZReader(DataReader):
    """
    Parse an XYZ file (single model) into an :class:`Frame`.

    Supports both standard XYZ and extended XYZ (extxyz) formats.

    Standard XYZ Format
    -------------------
        1. integer `N`  - number of atoms
        2. comment line - stored in frame.metadata
        3. N lines: `symbol  x  y  z`

    Extended XYZ Format
    -------------------
        1. integer `N`  - number of atoms
        2. comment line with key=value pairs (e.g., Properties=species:S:1:pos:R:3)
        3. N lines with columns defined by Properties specification
    """

    def read(self, frame: Frame | None = None) -> Frame:
        """
        Parameters
        ----------
        frame
            Optional frame to populate; if *None*, a new one is created.

        Returns
        -------
        Frame
            Frame with:
              * block ``"atoms"``:
                  - ``element``   -> (N,)  <U3   array
                  - ``x``         -> (N,3) float array (Å)
                  - ``z``         -> (N,)  int array (atomic numbers)
                  - additional columns from Properties if extxyz
              * metadata from comment line
        """
        # --- collect lines ------------------------------------------------
        lines: list[str] = self.read_lines()
        if len(lines) < 2:
            raise ValueError("XYZ file too short")

        natoms = int(lines[0])
        if len(lines) < natoms + 2:
            raise ValueError("XYZ record truncated")

        comment = lines[1]
        records = lines[2 : 2 + natoms]

        # --- build / update frame ----------------------------------------
        frame = frame or Frame()

        # Parse comment line for extxyz metadata
        metadata = self._parse_xyz_comment(comment)

        # Check if this is extxyz format with Properties
        if "Properties" in metadata:
            atoms_blk = self._parse_extxyz_atoms(records, metadata["Properties"])
        else:
            atoms_blk = self._parse_standard_xyz_atoms(records)

        # Set box if Lattice is present
        if "Lattice" in metadata:
            lattice_str = metadata.pop("Lattice")
            lattice_values = [float(x) for x in lattice_str.split()]
            frame.box = Box(np.array(lattice_values).reshape(3, 3))

        # Update frame metadata (excluding Properties and Lattice which are structural)
        frame_metadata = {
            k: v for k, v in metadata.items() if k not in ["Properties", "Lattice"]
        }
        frame.metadata.update(frame_metadata)

        frame["atoms"] = atoms_blk
        return frame

    def _parse_standard_xyz_atoms(self, records: list[str]) -> Block:
        """Parse standard XYZ atom records."""
        symbols: list[str] = []
        coords: list[tuple[float, float, float]] = []

        for rec in records:
            parts = rec.split()
            if len(parts) < 4:
                raise ValueError(f"Bad XYZ line: {rec!r}")
            symbols.append(parts[0])
            x, y, z = parts[1:4]
            coords.append((float(x), float(y), float(z)))

        atoms_blk = Block()
        atoms_blk["element"] = np.array(symbols, dtype="U3")

        # Store coordinates as separate x, y, z arrays (standard convention)
        coords_array = np.asarray(coords, dtype=float)
        atoms_blk["x"] = coords_array[:, 0]
        atoms_blk["y"] = coords_array[:, 1]
        atoms_blk["z"] = coords_array[:, 2]

        # Add atomic numbers using "number" field to avoid conflict with "z" coordinate
        number_list = [Element.get_atomic_number(sym) for sym in symbols]
        atoms_blk["number"] = np.array(number_list, dtype=np.int64)

        return atoms_blk

    def _parse_extxyz_atoms(self, records: list[str], properties_spec: list) -> Block:
        """Parse extended XYZ atom records using Properties specification.

        Args:
            records: Atom data lines
            properties_spec: List of (name, type, ncols) tuples from Properties
        """
        # Parse properties specification
        col_specs = []
        col_idx = 0
        for name, dtype, ncols in properties_spec:
            col_specs.append((name, dtype, col_idx, ncols))
            col_idx += ncols

        # Parse all records
        data_arrays = {name: [] for name, _, _, _ in col_specs}

        for rec in records:
            parts = rec.split()
            for name, dtype, start_idx, ncols in col_specs:
                if ncols == 1:
                    data_arrays[name].append(parts[start_idx])
                else:
                    data_arrays[name].append(parts[start_idx : start_idx + ncols])

        # Convert to numpy arrays
        atoms_blk = Block()
        for name, dtype, _, ncols in col_specs:
            arr = data_arrays[name]
            if dtype == "S":  # String
                atoms_blk[name] = np.array(arr, dtype="U3")
            elif dtype == "R":  # Real
                if ncols == 1:
                    atoms_blk[name] = np.array([float(x) for x in arr], dtype=float)
                else:
                    atoms_blk[name] = np.array(
                        [[float(x) for x in row] for row in arr], dtype=float
                    )
            elif dtype == "I":  # Integer
                if ncols == 1:
                    atoms_blk[name] = np.array([int(x) for x in arr], dtype=np.int64)
                else:
                    atoms_blk[name] = np.array(
                        [[int(x) for x in row] for row in arr], dtype=np.int64
                    )

        # Ensure standard fields exist
        # Map common extxyz names to standard names
        if "species" in atoms_blk and "element" not in atoms_blk:
            atoms_blk["element"] = atoms_blk["species"]

        # Split pos (N×3) into separate x, y, z arrays (standard convention)
        if "pos" in atoms_blk and "x" not in atoms_blk:
            pos = atoms_blk["pos"]
            if pos.ndim == 2 and pos.shape[1] == 3:
                atoms_blk["x"] = pos[:, 0]
                atoms_blk["y"] = pos[:, 1]
                atoms_blk["z"] = pos[:, 2]
            else:
                atoms_blk["x"] = pos

        # Add atomic numbers if not present using "number" field to avoid conflict with "z" coordinate
        if "number" not in atoms_blk and "element" in atoms_blk:
            symbols = atoms_blk["element"]
            number_list = [Element.get_atomic_number(str(sym)) for sym in symbols]
            atoms_blk["number"] = np.array(number_list, dtype=np.int64)

        return atoms_blk

    def _parse_xyz_comment(self, comment: str) -> dict:
        """
        Parse an extended XYZ comment line into a dictionary of key-value pairs.

        Args:
            comment (str): The comment line from an XYZ file.

        Returns:
            dict: Parsed key-value pairs.
        """
        result: dict = {}

        for token in shlex.split(comment):
            if "=" in token:
                key, value = token.split("=", 1)
                if key == "Properties":
                    parts = value.split(":")
                    triples = [
                        (parts[i], parts[i + 1], int(parts[i + 2]))
                        for i in range(0, len(parts), 3)
                    ]
                    result[key] = triples
                elif key == "Lattice":
                    # Parse lattice vectors
                    result[key] = value.strip('"')
                else:
                    result[key] = value.strip('"')
            else:
                # Standalone key, treat as bool flag
                result[token] = True

        return result
