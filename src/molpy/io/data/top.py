"""GROMACS topology file reader for Frame objects.

This module provides a reader for GROMACS topology files that extracts
structural information (atoms, bonds, angles, dihedrals, pairs) and creates
Frame objects with Block containers.
"""

from pathlib import Path
from typing import Any

import numpy as np

from molrs import Block, Frame
from molpy._frame_meta import get_frame_meta
from molrs import Element

from .base import DataReader, DataWriter, PathLike


class TopFormatError(Exception):
    """A ``.top`` file is structurally parseable but semantically invalid.

    Deliberately **not** a :class:`ValueError`: the per-line parsers treat
    ``ValueError`` as "this line is not a record, skip it", so a malformed
    *record* raised as ``ValueError`` would be silently dropped instead of
    reported.
    """


def _to_row_index(field: str) -> int:
    """Convert a GROMACS 1-based atom serial to a 0-based ``atoms`` row index.

    GROMACS numbers atoms from 1; canonical endpoint columns (``atomi``,
    ``atomj``, …) are 0-based row indices into the ``atoms`` block. Storing the
    file value unchanged shifts every endpoint by one and makes the last atom
    index out of range — a corruption that is invisible unless a topology
    happens to reference the final atom.

    Raises:
        TopFormatError: If *field* is not a positive integer, i.e. not a valid
            GROMACS serial. Never silently clamps or wraps to a negative index.
    """
    serial = int(field)
    if serial < 1:
        raise TopFormatError(
            f"GROMACS atom serial must be >= 1 (1-based), got {serial}"
        )
    return serial - 1


class TopReader(DataReader):
    """Read GROMACS topology files and create Frame objects.

    This reader parses GROMACS .top files and extracts structural information
    from sections like [atoms], [bonds], [angles], [dihedrals], and [pairs].

    Examples:
        >>> reader = TopReader("molecule.top")
        >>> frame = reader.read()
        >>> frame["atoms"]  # Block with atom data
        >>> frame["bonds"]  # Block with bond data
    """

    def __init__(self, file: PathLike, **open_kwargs):
        """Initialize GROMACS topology reader.

        Args:
            file: Path to GROMACS .top file
            **open_kwargs (Any): Additional arguments passed to file open.
        """
        super().__init__(file, **open_kwargs)
        self._file = Path(file)

    def read(self, frame: Frame | None = None) -> Frame:
        """Read GROMACS topology file and populate Frame.

        Args:
            frame: Optional existing Frame to populate. If None, creates a new one.

        Returns:
            Frame object with atoms, bonds, angles, dihedrals, and pairs blocks.
        """
        if frame is None:
            frame = Frame()

        # Initialize data containers
        atoms_data: list[dict[str, Any]] = []
        bonds_data: list[dict[str, Any]] = []
        pairs_data: list[dict[str, Any]] = []
        angles_data: list[dict[str, Any]] = []
        dihedrals_data: list[dict[str, Any]] = []

        # Read and parse file
        with self._file.open("r", encoding="utf-8") as f:
            section: str | None = None
            for line in f:
                # Remove inline comments first
                if ";" in line:
                    line = line.split(";")[0]
                line = line.strip()

                # Skip empty lines
                if not line:
                    continue

                # Check for section header (must be on its own line)
                if line.startswith("[") and line.endswith("]"):
                    # Normalize: "[ atoms ]", "[atoms]", "[ATOMS]" → "[ atoms ]"
                    inner = line[1:-1].strip().lower()
                    section = f"[ {inner} ]"
                    continue

                # Skip comment-only lines
                if line.startswith("#"):
                    continue

                # Parse content based on current section
                if section == "[ atoms ]":
                    atom_dict = self._parse_atom_line(line)
                    if atom_dict:
                        atoms_data.append(atom_dict)
                elif section == "[ bonds ]":
                    bond_dict = self._parse_bond_line(line)
                    if bond_dict:
                        bonds_data.append(bond_dict)
                elif section == "[ pairs ]":
                    pair_dict = self._parse_pair_line(line)
                    if pair_dict:
                        pairs_data.append(pair_dict)
                elif section == "[ angles ]":
                    angle_dict = self._parse_angle_line(line)
                    if angle_dict:
                        angles_data.append(angle_dict)
                elif section == "[ dihedrals ]":
                    dihedral_dict = self._parse_dihedral_line(line)
                    if dihedral_dict:
                        dihedrals_data.append(dihedral_dict)

        # Convert lists of dicts to Blocks
        if atoms_data:
            frame["atoms"] = self._dicts_to_block(atoms_data)
            # Assign atomic numbers
            self._assign_atomic_numbers(frame["atoms"])

        if bonds_data:
            frame["bonds"] = self._dicts_to_block(bonds_data)

        if pairs_data:
            frame["pairs"] = self._dicts_to_block(pairs_data)

        if angles_data:
            frame["angles"] = self._dicts_to_block(angles_data)

        if dihedrals_data:
            frame["dihedrals"] = self._dicts_to_block(dihedrals_data)

        return frame

    def _parse_atom_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single line from [atoms] section.

        GROMACS atoms format:
        nr type resnr residu atom cgnr charge mass [typeB chargeB massB]

        Args:
            line: Line to parse

        Returns:
            Dictionary with atom data, or None if line is invalid
        """
        # Remove inline comments
        line = line.split(";")[0].strip()
        if not line:
            return None

        parts = line.split()
        if len(parts) < 8:
            return None  # Invalid line

        try:
            return {
                "id": int(parts[0]),
                "type": parts[1],
                "resnr": int(parts[2]),
                "residu": parts[3],
                "name": parts[4],
                "cgnr": int(parts[5]),
                "charge": float(parts[6]),
                "mass": float(parts[7]),
            }
        except (ValueError, IndexError):
            return None

    def _parse_bond_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single line from [bonds] section.

        GROMACS bonds format:
        ai aj funct [parameters...]

        Args:
            line: Line to parse

        Returns:
            Dictionary with bond data, or None if line is invalid
        """
        # Remove inline comments
        line = line.split(";")[0].strip()
        if not line:
            return None

        parts = line.split()
        if len(parts) < 3:
            return None  # Invalid line

        try:
            return {
                "atomi": _to_row_index(parts[0]),
                "atomj": _to_row_index(parts[1]),
                "type": int(parts[2]),
            }
        except (ValueError, IndexError):
            return None

    def _parse_pair_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single line from [pairs] section.

        GROMACS pairs format:
        ai aj funct [parameters...]

        Args:
            line: Line to parse

        Returns:
            Dictionary with pair data, or None if line is invalid
        """
        # Remove inline comments
        line = line.split(";")[0].strip()
        if not line:
            return None

        parts = line.split()
        if len(parts) < 3:
            return None  # Invalid line

        try:
            return {
                "atomi": _to_row_index(parts[0]),
                "atomj": _to_row_index(parts[1]),
                "type": int(parts[2]),
            }
        except (ValueError, IndexError):
            return None

    def _parse_angle_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single line from [angles] section.

        GROMACS angles format:
        ai aj ak funct [parameters...]

        Args:
            line: Line to parse

        Returns:
            Dictionary with angle data, or None if line is invalid
        """
        # Remove inline comments
        line = line.split(";")[0].strip()
        if not line:
            return None

        parts = line.split()
        if len(parts) < 4:
            return None  # Invalid line

        try:
            return {
                "atomi": _to_row_index(parts[0]),
                "atomj": _to_row_index(parts[1]),
                "atomk": _to_row_index(parts[2]),
                "type": int(parts[3]),
            }
        except (ValueError, IndexError):
            return None

    def _parse_dihedral_line(self, line: str) -> dict[str, Any] | None:
        """Parse a single line from [dihedrals] section.

        GROMACS dihedrals format:
        ai aj ak al funct [parameters...]

        Args:
            line: Line to parse

        Returns:
            Dictionary with dihedral data, or None if line is invalid
        """
        # Remove inline comments
        line = line.split(";")[0].strip()
        if not line:
            return None

        parts = line.split()
        if len(parts) < 5:
            return None  # Invalid line

        try:
            return {
                "atomi": _to_row_index(parts[0]),
                "atomj": _to_row_index(parts[1]),
                "atomk": _to_row_index(parts[2]),
                "atoml": _to_row_index(parts[3]),
                "type": int(parts[4]),
            }
        except (ValueError, IndexError):
            return None

    def _dicts_to_block(self, data: list[dict[str, Any]]) -> Block:
        """Convert list of dictionaries to Block.

        Args:
            data: List of dictionaries with consistent keys

        Returns:
            Block object with columns from dictionary keys
        """
        if not data:
            return Block()

        # Get all unique keys
        all_keys = set()
        for d in data:
            all_keys.update(d.keys())

        # Build column arrays
        block_data: dict[str, list[Any]] = {}
        for key in sorted(all_keys):
            block_data[key] = []

        # Fill data
        for d in data:
            for key in all_keys:
                value = d.get(key)
                block_data[key].append(value)

        # Convert to numpy arrays, handling mixed types intelligently
        result: dict[str, np.ndarray] = {}
        for key, values in block_data.items():
            # Check if all values are None
            if all(v is None for v in values):
                result[key] = np.array([None] * len(values), dtype=object)
                continue

            # Try to infer the best dtype
            non_none_values = [v for v in values if v is not None]
            if non_none_values:
                # Only use int if every value is already an int (not a float)
                if all(isinstance(v, int) for v in non_none_values):
                    arr = np.array(
                        [int(v) if v is not None else 0 for v in values],
                        dtype=int,
                    )
                    result[key] = arr
                    continue

                # Try float
                try:
                    arr = np.array(
                        [float(v) if v is not None else np.nan for v in values],
                        dtype=float,
                    )
                    result[key] = arr
                    continue
                except (ValueError, TypeError):
                    pass

            # Fall back to string/object array
            # Convert to strings for consistency
            str_values = [str(v) if v is not None else "" for v in values]
            # Try to use fixed-width unicode string if all are short strings
            max_len = max(len(s) for s in str_values) if str_values else 1
            if max_len <= 32:  # Use fixed-width string for short strings
                result[key] = np.array(str_values, dtype=f"U{max_len}")
            else:
                result[key] = np.array(str_values, dtype=object)

        return Block(result)

    def _assign_atomic_numbers(self, atoms_block: Block) -> None:
        """Assign atomic numbers to atoms based on element names.

        Args:
            atoms_block: Block containing atom data with 'name' and 'type' columns
        """
        if "name" not in atoms_block and "type" not in atoms_block:
            return

        atomic_numbers = []
        names = atoms_block.get("name", atoms_block.get("type", None))
        types = atoms_block.get("type", None)

        if names is None:
            return

        # Convert to list if it's a numpy array
        if isinstance(names, np.ndarray):
            names = names.tolist()
        if isinstance(types, np.ndarray):
            types = types.tolist()

        for i, name in enumerate(names):
            atomic_number = self._guess_atomic_number(name, types[i] if types else None)
            atomic_numbers.append(atomic_number)

        # Add atomic numbers to block
        atoms_block["number"] = np.array(atomic_numbers, dtype=int)

    def _guess_atomic_number(self, name: str, atom_type: str | None = None) -> int:
        """Guess atomic number from element name or type.

        Args:
            name: Atom name
            atom_type: Optional atom type (used as fallback)

        Returns:
            Atomic number (0 if cannot be determined)
        """
        # Try name first
        element_name = "".join(c for c in name if c.isalpha())
        if element_name:
            try:
                element = Element(element_name)
                return element.number
            except (KeyError, AttributeError):
                pass

        # Try first letter + rest lowercase for multi-character names
        if len(element_name) > 1:
            try:
                element = Element(element_name[0].upper() + element_name[1:].lower())
                return element.number
            except (KeyError, AttributeError):
                pass

        # Fallback to atom type
        if atom_type:
            type_name = "".join(c for c in atom_type if c.isalpha())
            if type_name:
                try:
                    element = Element(type_name)
                    return element.number
                except (KeyError, AttributeError):
                    pass

        return 0  # Unknown element


class TopWriter(DataWriter):
    """Write GROMACS topology files from Frame objects.

    Produces a minimal ``.top`` file with ``[ moleculetype ]``,
    ``[ atoms ]``, ``[ bonds ]``, and optional ``[ pairs ]``,
    ``[ angles ]``, ``[ dihedrals ]`` sections, followed by
    ``[ system ]`` and ``[ molecules ]``.

    Examples:
        >>> writer = TopWriter("molecule.top")
        >>> writer.write(frame)
    """

    def __init__(self, file: PathLike, **open_kwargs):
        super().__init__(file, **open_kwargs)
        self._file = Path(file)

    def write(self, frame: Frame) -> None:
        """Write Frame to GROMACS topology file.

        Args:
            frame: Frame object containing ``"atoms"`` and optionally
                ``"bonds"``, ``"pairs"``, ``"angles"``, ``"dihedrals"`` blocks.
        """
        mol_name = get_frame_meta(frame, "name", "MOL")
        lines: list[str] = []

        # Header
        lines.append("; Generated by MolPy TopWriter")
        lines.append("")

        # [ moleculetype ]
        lines.append("[ moleculetype ]")
        lines.append("; name  nrexcl")
        lines.append(f"{mol_name}  3")
        lines.append("")

        # [ atoms ]
        if "atoms" in frame:
            atoms = frame["atoms"]
            # Force-field data must come from the frame. Fabricating a type
            # ("X"), a charge (0.0) or a mass (0.0) produces a topology that
            # loads and runs while being physically meaningless — strictly
            # worse than refusing to write it.
            missing = [c for c in ("type", "charge", "mass") if c not in atoms]
            if missing:
                raise ValueError(
                    f"Cannot write {self._file.name}: the 'atoms' block is missing "
                    f"required column(s) {missing}. GROMACS [ atoms ] records carry "
                    f"force-field data that molpy will not invent; typify or assign "
                    f"these columns before writing."
                )
            lines.append("[ atoms ]")
            lines.append(";  nr  type  resnr  residu  atom  cgnr  charge  mass")
            n_atoms = atoms.nrows
            for i in range(n_atoms):
                row = atoms[i]
                atype = str(row["type"])
                charge = float(row["charge"])
                mass = float(row["mass"])
                # Positional/structural serials, not measurements: GROMACS
                # defines `nr` and `cgnr` as 1-based ordinals, and a single
                # [ moleculetype ] block is one residue by construction. These
                # are format conventions, so deriving them is not fabrication.
                aid = int(row.get("id", i + 1))
                resnr = int(row.get("resnr", 1))
                residu = str(row.get("residu", mol_name))
                name = str(row.get("name", atype))
                cgnr = int(row.get("cgnr", i + 1))
                lines.append(
                    f"  {aid}  {atype}  {resnr}  {residu}  {name}  "
                    f"{cgnr}  {charge:.4f}  {mass:.3f}"
                )
            lines.append("")

        # Helper for index-based sections
        def _write_index_section(
            section_name: str, key: str, columns: list[str]
        ) -> None:
            if key not in frame:
                return
            block = frame[key]
            lines.append(f"[ {section_name} ]")
            header = "  ".join(columns + ["funct"])
            lines.append(f"; {header}")
            for i in range(block.nrows):
                row = block[i]
                # Endpoints are 0-based row indices in the frame; GROMACS
                # serials are 1-based. Mirror of `_to_row_index` on read.
                vals = "  ".join(str(int(row[c]) + 1) for c in columns)
                # `funct` is the GROMACS interaction-function selector, not a
                # molpy field. Frames from other formats carry no equivalent,
                # so 1 (the GROMACS default) is retained here rather than
                # refusing to write. Revisited when `type`/`type_id` split.
                funct = int(row.get("type", 1))
                lines.append(f"  {vals}  {funct}")
            lines.append("")

        _write_index_section("bonds", "bonds", ["atomi", "atomj"])
        _write_index_section("pairs", "pairs", ["atomi", "atomj"])
        _write_index_section("angles", "angles", ["atomi", "atomj", "atomk"])
        _write_index_section(
            "dihedrals", "dihedrals", ["atomi", "atomj", "atomk", "atoml"]
        )

        # [ system ]
        lines.append("[ system ]")
        lines.append(mol_name)
        lines.append("")

        # [ molecules ]
        lines.append("[ molecules ]")
        lines.append(f"{mol_name}  1")
        lines.append("")

        self._file.write_text("\n".join(lines), encoding="utf-8")
