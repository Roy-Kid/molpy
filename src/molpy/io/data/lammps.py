"""
Modern LAMMPS data file I/O using Block.from_csv.

This module provides a clean, imperative approach to reading and writing
LAMMPS data files using the Block.from_csv functionality.
"""

from dataclasses import dataclass
from datetime import datetime
from io import StringIO
from pathlib import Path
from typing import Any, Callable
import re

import numpy as np

from molpy.core.frame import Block
from molpy.core.box import Box
from molpy.core.forcefield import ForceField
from molpy.core.frame import Frame

from .base import DataReader, DataWriter


# Constants
DEFAULT_BOX_SIZE = 10.0
TYPE_KEY_MAPPING = {
    "atoms": "atom_types",
    "bonds": "bond_types",
    "angles": "angle_types",
    "dihedrals": "dihedral_types",
    "impropers": "improper_types",
}

# Connectivity field requirements
CONNECTIVITY_FIELDS = {
    "bonds": ["atomi", "atomj"],
    "angles": ["atomi", "atomj", "atomk"],
    "dihedrals": ["atomi", "atomj", "atomk", "atoml"],
    "impropers": ["atomi", "atomj", "atomk", "atoml"],
}


class TypeSystem:
    """Centralized type mapping for LAMMPS I/O.

    Manages bidirectional mapping between type names (strings) and type IDs (integers).
    Merges metadata type labels with actual types found in data.
    """

    def __init__(
        self,
        metadata_labels: dict[str, list[str]] | None,
        actual_types: dict[str, set[str]],
    ):
        """Initialize type system.

        Args:
            metadata_labels: Optional dict mapping category -> list of type names from metadata
            actual_types: Dict mapping category -> set of type names found in data
        """
        self._mappings = self._build_mappings(metadata_labels, actual_types)
        self._reverse_mappings = self._build_reverse_mappings()

    def _build_mappings(
        self,
        metadata_labels: dict[str, list[str]] | None,
        actual_types: dict[str, set[str]],
    ) -> dict[str, dict[str, int]]:
        """Build type name -> type ID mappings."""
        mappings = {}

        for category in [
            "atom_types",
            "bond_types",
            "angle_types",
            "dihedral_types",
            "improper_types",
        ]:
            metadata_list = metadata_labels.get(category, []) if metadata_labels else []
            actual_set = actual_types.get(category, set())

            # Merge: metadata types + actual types not in metadata
            merged_set = set(metadata_list) | actual_set
            sorted_types = sorted(list(merged_set))

            # Build 1-based mapping
            mappings[category] = {
                type_name: idx + 1 for idx, type_name in enumerate(sorted_types)
            }

        return mappings

    def _build_reverse_mappings(self) -> dict[str, dict[int, str]]:
        """Build type ID -> type name reverse mappings."""
        reverse = {}
        for category, mapping in self._mappings.items():
            reverse[category] = {
                type_id: type_name for type_name, type_id in mapping.items()
            }
        return reverse

    def get_id(self, type_name: str, category: str) -> int:
        """Get type ID for a type name."""
        return self._mappings.get(category, {}).get(type_name, 1)

    def get_label(self, type_id: int, category: str) -> str:
        """Get type name for a type ID."""
        return self._reverse_mappings.get(category, {}).get(type_id, str(type_id))

    def get_type_list(self, category: str) -> list[str]:
        """Get sorted list of type names for a category."""
        mapping = self._mappings.get(category, {})
        return sorted(mapping.keys(), key=lambda x: mapping[x])

    def has_types(self, category: str) -> bool:
        """Check if category has any types."""
        return bool(self._mappings.get(category))


@dataclass
class HeaderField:
    """Configuration for a header field."""

    pattern: str
    key: str
    parser: Callable[[tuple], Any] = lambda x: int(x[0])


@dataclass
class CoeffConfig:
    """Configuration for force field coefficient parsing."""

    style_name: str
    style_def: Callable[[ForceField], None]
    min_fields: int


def parse_bounds(parts: tuple[str, str]) -> tuple[float, float]:
    """Parse box bounds from string parts."""
    return (float(parts[0]), float(parts[1]))


def build_id_to_index_mapping(atoms_block: Block) -> dict[int, int]:
    """Build mapping from atom IDs to 0-based indices."""
    mapping = {}
    for idx, atom_id in enumerate(atoms_block["id"]):
        mapping[int(atom_id)] = idx
    return mapping


# Header field configurations
HEADER_FIELDS = [
    # Count fields
    HeaderField(r"^(\d+)\s+atoms$", "atoms"),
    HeaderField(r"^(\d+)\s+bonds$", "bonds"),
    HeaderField(r"^(\d+)\s+angles$", "angles"),
    HeaderField(r"^(\d+)\s+dihedrals$", "dihedrals"),
    HeaderField(r"^(\d+)\s+impropers$", "impropers"),
    HeaderField(r"^(\d+)\s+atom\s+types$", "atom_types"),
    HeaderField(r"^(\d+)\s+bond\s+types$", "bond_types"),
    HeaderField(r"^(\d+)\s+angle\s+types$", "angle_types"),
    HeaderField(r"^(\d+)\s+dihedral\s+types$", "dihedral_types"),
    HeaderField(r"^(\d+)\s+improper\s+types$", "improper_types"),
    # Box bounds
    HeaderField(
        r"^([\d.\-+eE]+)\s+([\d.\-+eE]+)\s+xlo\s+xhi$", "x_bounds", parse_bounds
    ),
    HeaderField(
        r"^([\d.\-+eE]+)\s+([\d.\-+eE]+)\s+ylo\s+yhi$", "y_bounds", parse_bounds
    ),
    HeaderField(
        r"^([\d.\-+eE]+)\s+([\d.\-+eE]+)\s+zlo\s+zhi$", "z_bounds", parse_bounds
    ),
]


# Force field coefficient configurations
COEFF_CONFIGS = {
    "PairCoeffs": CoeffConfig("lj/cut", lambda ff: ff.def_pairstyle("lj/cut"), 3),
    "BondCoeffs": CoeffConfig("harmonic", lambda ff: ff.def_bondstyle("harmonic"), 3),
    "AngleCoeffs": CoeffConfig("harmonic", lambda ff: ff.def_anglestyle("harmonic"), 3),
    "DihedralCoeffs": CoeffConfig(
        "harmonic", lambda ff: ff.def_dihedralstyle("harmonic"), 4
    ),
    "ImproperCoeffs": CoeffConfig(
        "harmonic", lambda ff: ff.def_improperstyle("harmonic"), 4
    ),
}


class LammpsDataReader(DataReader):
    """Modern LAMMPS data file reader using Block.from_csv."""

    def __init__(self, path: str | Path, atom_style: str = "full") -> None:
        super().__init__(Path(path))
        self.atom_style = atom_style

    def read(self, frame: Frame | None = None) -> Frame:
        """Read LAMMPS data file into a Frame."""
        frame = frame or Frame()

        # Read and parse the file
        lines = self._read_lines()
        sections = self._extract_sections(lines)

        # Parse header and set up box
        header_info = self._parse_header(sections.get("header", []))
        frame.metadata["box"] = self._create_box(header_info["box_bounds"])

        # Parse masses if present
        masses = self._parse_masses(sections.get("Masses", []))

        # Parse type labels
        type_labels = self._parse_type_labels(sections)

        # Parse force field parameters
        forcefield = self._parse_force_field(sections)

        # Parse atoms section
        if "Atoms" in sections:
            frame["atoms"] = self._parse_atoms_section(
                sections["Atoms"], masses, type_labels.get("atom", {})
            )

        # Build id to index mapping for connectivity
        id_to_idx = {}
        if "atoms" in frame:
            for i, atom_id in enumerate(frame["atoms"]["id"]):
                id_to_idx[int(atom_id)] = i

        # Parse connectivity sections
        if "Bonds" in sections and header_info["counts"].get("bonds", 0) > 0:
            frame["bonds"] = self._parse_connectivity_section(
                sections["Bonds"], "bond", type_labels.get("bond", {}), id_to_idx
            )

        if "Angles" in sections and header_info["counts"].get("angles", 0) > 0:
            frame["angles"] = self._parse_connectivity_section(
                sections["Angles"], "angle", type_labels.get("angle", {}), id_to_idx
            )

        if "Dihedrals" in sections and header_info["counts"].get("dihedrals", 0) > 0:
            frame["dihedrals"] = self._parse_connectivity_section(
                sections["Dihedrals"],
                "dihedral",
                type_labels.get("dihedral", {}),
                id_to_idx,
            )

        if "Impropers" in sections and header_info["counts"].get("impropers", 0) > 0:
            frame["impropers"] = self._parse_connectivity_section(
                sections["Impropers"],
                "improper",
                type_labels.get("improper", {}),
                id_to_idx,
            )

        # Store metadata
        frame.metadata.update(
            {
                "format": "lammps_data",
                "atom_style": self.atom_style,
                "counts": header_info["counts"],
                "source_file": str(self._path),
                "forcefield": forcefield,
            }
        )

        return frame

    def _read_lines(self) -> list[str]:
        """Read file and return non-empty, non-comment lines."""
        with open(self._path) as f:
            return [
                line.strip()
                for line in f
                if line.strip() and not line.strip().startswith("#")
            ]

    def _extract_sections(self, lines: list[str]) -> dict[str, list[str]]:
        """Extract sections from LAMMPS data file."""
        sections = {"header": []}
        current_section = "header"

        section_keywords = [
            "atoms",
            "masses",
            "bonds",
            "angles",
            "dihedrals",
            "impropers",
            "atom type labels",
            "bond type labels",
            "angle type labels",
            "dihedral type labels",
            "improper type labels",
            "pair coeffs",
            "bond coeffs",
            "angle coeffs",
            "dihedral coeffs",
            "improper coeffs",
        ]

        for line in lines:
            line_lower = line.lower()

            # Check if this line starts a new section
            if any(line_lower.startswith(keyword) for keyword in section_keywords):
                if "type labels" in line_lower:
                    # Handle type labels sections
                    section_name = line.replace("Type Labels", "TypeLabels").replace(
                        " ", ""
                    )
                elif "coeffs" in line_lower:
                    # Handle force field coefficients sections
                    section_name = line.replace(" ", "")
                else:
                    # Handle regular sections
                    section_name = line.split()[0].capitalize()
                current_section = section_name
                sections[current_section] = []
            else:
                # Add line to current section
                if current_section not in sections:
                    sections[current_section] = []
                sections[current_section].append(line)

        return sections

    def _parse_header(self, header_lines: list[str]) -> dict[str, Any]:
        """Parse header information using data-driven configuration."""
        counts = {}
        box_bounds = {}

        for line in header_lines:
            # Try to match against configured patterns
            matched = False
            for field in HEADER_FIELDS:
                match = re.match(field.pattern, line)
                if match:
                    value = field.parser(match.groups())
                    if field.key.endswith("_bounds"):
                        box_bounds[field.key[0]] = value  # x_bounds -> x
                    else:
                        counts[field.key] = value
                    matched = True
                    break

            # If line looks like it should be a header field but didn't match, validate it
            if not matched and line.strip():
                # Check for malformed count lines
                if any(
                    keyword in line.lower()
                    for keyword in [
                        "atoms",
                        "bonds",
                        "angles",
                        "dihedrals",
                        "impropers",
                        "types",
                    ]
                ):
                    parts = line.split()
                    if len(parts) >= 2:
                        # Try to parse the first part as int to trigger error if malformed
                        int(parts[0])

                # Check for malformed box bounds
                if any(
                    keyword in line.lower()
                    for keyword in ["xlo", "ylo", "zlo", "xhi", "yhi", "zhi"]
                ):
                    parts = line.split()
                    if len(parts) >= 2:
                        # Try to parse first two parts as floats to trigger error if malformed
                        float(parts[0])
                        float(parts[1])

        return {"counts": counts, "box_bounds": box_bounds if box_bounds else None}

    def _create_box(self, box_bounds: dict[str, tuple[float, float]] | None) -> Box:
        """Create Box from bounds."""
        if not box_bounds:
            return Box(np.array([DEFAULT_BOX_SIZE] * 3))

        # Use provided bounds or default for each dimension
        default_bounds = (0.0, DEFAULT_BOX_SIZE)
        x_bounds = box_bounds.get("x", default_bounds)
        y_bounds = box_bounds.get("y", default_bounds)
        z_bounds = box_bounds.get("z", default_bounds)

        lengths = np.array(
            [
                x_bounds[1] - x_bounds[0],
                y_bounds[1] - y_bounds[0],
                z_bounds[1] - z_bounds[0],
            ]
        )
        origin = np.array(
            [
                x_bounds[0],
                y_bounds[0],
                z_bounds[0],
            ]
        )
        return Box(lengths, origin=origin)

    def _parse_masses(self, mass_lines: list[str]) -> dict[str, float]:
        """Parse mass section."""
        masses = {}
        for line in mass_lines:
            parts = line.split()
            if len(parts) >= 2:
                type_str = parts[0]
                mass = float(parts[1])
                masses[type_str] = mass
        return masses

    def _parse_type_labels(
        self, sections: dict[str, list[str]]
    ) -> dict[str, dict[int, str]]:
        """Parse all type labels sections."""
        type_labels = {}

        label_sections = {
            "atom": "AtomTypeLabels",
            "bond": "BondTypeLabels",
            "angle": "AngleTypeLabels",
            "dihedral": "DihedralTypeLabels",
            "improper": "ImproperTypeLabels",
        }

        for label_type, section_name in label_sections.items():
            if section_name in sections:
                id_to_label = {}
                for line in sections[section_name]:
                    parts = line.split()
                    if len(parts) >= 2:
                        type_id = int(parts[0])
                        label = parts[1]
                        id_to_label[type_id] = label
                type_labels[label_type] = id_to_label

        return type_labels

    def _parse_force_field(self, sections: dict[str, list[str]]) -> ForceField:
        """Parse force field parameters using data-driven configuration."""
        forcefield = ForceField()

        for section_name, config in COEFF_CONFIGS.items():
            if section_name in sections:
                config.style_def(forcefield)
                for line in sections[section_name]:
                    parts = line.split()
                    if len(parts) >= config.min_fields:
                        # Validate parsing (will raise ValueError if malformed)
                        int(parts[0])
                        for i in range(1, config.min_fields):
                            float(parts[i])
                        # TODO: Store in forcefield when implementation ready

        return forcefield

    def _parse_atoms_section(
        self,
        atom_lines: list[str],
        masses: dict[str, float],
        type_labels: dict[int, str],
    ) -> Block:
        """Parse atoms section using Block.from_csv with space delimiter."""
        if not atom_lines:
            return Block()

        # Create space-separated string for Block.from_csv
        csv_lines = []

        # Add header based on atom style
        if self.atom_style == "full":
            header = ["id", "mol", "type", "q", "x", "y", "z"]
        elif self.atom_style == "charge":
            header = ["id", "type", "q", "x", "y", "z"]
        else:  # atomic
            header = ["id", "type", "x", "y", "z"]

        csv_lines.append(" ".join(header))

        # Add data lines directly
        for line in atom_lines:
            parts = line.split()
            if len(parts) >= len(header):
                csv_lines.append(line)

        # Parse using Block.from_csv with space delimiter
        csv_string = "\n".join(csv_lines)
        block = Block.from_csv(
            StringIO(csv_string), delimiter=" ", skipinitialspace=True
        )

        # Add mass information
        if block.nrows > 0:
            mass_values = []
            for type_str in block["type"]:
                mass_values.append(masses.get(str(type_str), 1.0))
            block["mass"] = np.array(mass_values)

            # Convert numeric types back to string types using type labels
            if type_labels:
                converted_types = []
                for type_id in block["type"]:
                    # Check if type_id is numeric (needs conversion) or already a string label
                    type_id_str = str(type_id)
                    # Try to convert to int - if it works, it's a numeric ID that needs label lookup
                    # If it fails, it's already a string label
                    try:
                        type_id_int = int(type_id_str)
                        # It's numeric, look up the label
                        converted_type = type_labels.get(type_id_int, type_id_str)
                        converted_types.append(converted_type)
                    except ValueError:
                        # It's already a string label, keep it as is
                        converted_types.append(type_id_str)
                block["type"] = np.array(converted_types)

        return block

    def _parse_connectivity_section(
        self,
        lines: list[str],
        section_type: str,
        type_labels: dict[int, str],
        id_to_idx: dict[int, int],
    ) -> Block:
        """Parse connectivity sections (bonds, angles, dihedrals, impropers)."""
        if not lines:
            return Block()

        # Define temporary headers for parsing (using atom1, atom2 etc. as placeholders)
        headers = {
            "bond": ["id", "type", "atom1", "atom2"],
            "angle": ["id", "type", "atom1", "atom2", "atom3"],
            "dihedral": ["id", "type", "atom1", "atom2", "atom3", "atom4"],
            "improper": ["id", "type", "atom1", "atom2", "atom3", "atom4"],
        }

        header = headers[section_type]
        csv_lines = [" ".join(header)]

        # Add data lines directly
        for line in lines:
            parts = line.split()
            if len(parts) >= len(header):
                csv_lines.append(line)

        # Parse using Block.from_csv
        csv_string = "\n".join(csv_lines)
        block = Block.from_csv(
            StringIO(csv_string), delimiter=" ", skipinitialspace=True
        )

        # Map atom IDs to 0-based indices and rename columns
        atom_key_map = {
            "atom1": "atomi",
            "atom2": "atomj",
            "atom3": "atomk",
            "atom4": "atoml",
        }

        for old_key, new_key in atom_key_map.items():
            if old_key in block:
                # Convert IDs to indices using the provided mapping
                ids = block[old_key].astype(int)
                indices = np.array([id_to_idx.get(id_val, -1) for id_val in ids])

                # Check for unmapped IDs
                if np.any(indices == -1):
                    unmapped = ids[indices == -1]
                    import warnings

                    warnings.warn(
                        f"Found {len(unmapped)} atom IDs in {section_type} section "
                        f"that could not be mapped to atom indices: {unmapped[:5]}..."
                    )

                # Add new column and remove old one
                block[new_key] = indices
                del block[old_key]

        return block


class LammpsDataWriter(DataWriter):
    """Modern LAMMPS data file writer using Block.to_csv approach.

    **Important Requirements:**
    - Atoms in the frame must have an 'id' field. This field is required
      to map atom indices to atom IDs for LAMMPS output.
    - Connectivity data (bonds, angles, dihedrals) in the frame uses atom
      indices (0-based from to_frame()). The writer automatically converts
      these indices to atom IDs using the index->ID mapping from the atoms
      'id' field.

    **Frame Structure:**
    - Atoms: Must include 'id' field. Other required fields depend on atom_style.
    - Bonds/Angles/Dihedrals: Use atom indices in 'atomi', 'atomj', 'atomk', 'atoml'
      (from to_frame()). These are 0-based indices that will be converted to 1-based atom IDs.
    """

    def __init__(self, path: str | Path, atom_style: str = "full") -> None:
        super().__init__(Path(path))
        self.atom_style = atom_style

    def write(self, frame: Frame) -> None:
        """Write Frame to LAMMPS data file.

        Args:
            frame: Frame containing atoms and optionally bonds/angles/dihedrals.
                  Atoms must have 'id' field.

        Raises:
            ValueError: If atoms are missing 'id' field.
        """
        lines = []

        # Build TypeSystem once for entire write operation
        metadata_type_labels = frame.metadata.get("type_labels")
        actual_types = self._collect_actual_types(frame)
        type_system = TypeSystem(metadata_type_labels, actual_types)

        # Header
        lines.append(
            f"# LAMMPS data file written by molpy on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}"
        )
        lines.append("")

        # Count sections
        counts = self._get_counts(frame)
        self._write_counts(lines, counts)
        lines.append("")

        # Type counts
        self._write_type_counts(lines, type_system)
        lines.append("")

        # Box bounds
        self._write_box_bounds(lines, frame)
        lines.append("")

        # Type labels sections (must come before Masses for LAMMPS)
        self._write_type_labels_sections(lines, type_system)

        # Masses section
        if "atoms" in frame:
            self._write_masses_section(lines, frame, type_system)

        # Force field coefficients sections
        self._write_force_field_coeffs_sections(lines, frame)

        # Data sections
        if "atoms" in frame:
            self._write_atoms_section(lines, frame, type_system)

        if "bonds" in frame and counts.get("bonds", 0) > 0:
            self._write_connectivity_section(lines, frame, "bonds", type_system)

        if "angles" in frame and counts.get("angles", 0) > 0:
            self._write_connectivity_section(lines, frame, "angles", type_system)

        if "dihedrals" in frame and counts.get("dihedrals", 0) > 0:
            self._write_connectivity_section(lines, frame, "dihedrals", type_system)

        if "impropers" in frame and counts.get("impropers", 0) > 0:
            self._write_connectivity_section(lines, frame, "impropers", type_system)

        # Write to file
        with open(self._path, "w") as f:
            f.write("\n".join(lines))

    def _get_counts(self, frame: Frame) -> dict[str, int]:
        """Get counts from frame."""
        block_names = ["atoms", "bonds", "angles", "dihedrals", "impropers"]
        return {name: frame[name].nrows for name in block_names if name in frame}

    def _collect_actual_types(self, frame: Frame) -> dict[str, set[str]]:
        """Collect actual types used in frame blocks."""
        actual_types = {}

        for block_name, type_key in TYPE_KEY_MAPPING.items():
            if block_name in frame and frame[block_name].nrows > 0:
                types = frame[block_name]["type"]
                # Validate types
                for t in types:
                    if t is None:
                        raise ValueError(
                            f"Found None type in {block_name} block. "
                            f"All {block_name} must have valid type values."
                        )
                    if isinstance(t, str) and not t.strip():
                        raise ValueError(
                            f"Found empty string type in {block_name} block. "
                            f"All {block_name} must have non-empty type values."
                        )
                # Convert to string set
                actual_types[type_key] = set(str(t) for t in np.unique(types))

        return actual_types

    def _write_counts(self, lines: list[str], counts: dict[str, int]) -> None:
        """Write count lines."""
        if "atoms" in counts:
            lines.append(f"{counts['atoms']} atoms")
        if "bonds" in counts and counts["bonds"] > 0:
            lines.append(f"{counts['bonds']} bonds")
        if "angles" in counts and counts["angles"] > 0:
            lines.append(f"{counts['angles']} angles")
        if "dihedrals" in counts and counts["dihedrals"] > 0:
            lines.append(f"{counts['dihedrals']} dihedrals")
        if "impropers" in counts and counts["impropers"] > 0:
            lines.append(f"{counts['impropers']} impropers")

    def _write_type_counts(self, lines: list[str], type_system: TypeSystem) -> None:
        """Write type count lines using TypeSystem."""
        for category, label in [
            ("atom_types", "atom types"),
            ("bond_types", "bond types"),
            ("angle_types", "angle types"),
            ("dihedral_types", "dihedral types"),
            ("improper_types", "improper types"),
        ]:
            if type_system.has_types(category):
                type_list = type_system.get_type_list(category)
                lines.append(f"{len(type_list)} {label}")

    def _write_box_bounds(self, lines: list[str], frame: Frame) -> None:
        """Write box bounds."""
        if frame.metadata.get("box") is not None:
            box = frame.metadata["box"]
            lines.append(
                f"{box.origin[0]:.6f} {box.origin[0] + box.lengths[0]:.6f} xlo xhi"
            )
            lines.append(
                f"{box.origin[1]:.6f} {box.origin[1] + box.lengths[1]:.6f} ylo yhi"
            )
            lines.append(
                f"{box.origin[2]:.6f} {box.origin[2] + box.lengths[2]:.6f} zlo zhi"
            )
            # Add tilt factors for triclinic boxes
            if box.style == box.Style.TRICLINIC:
                lines.append(f"{box.xy:.6f} {box.xz:.6f} {box.yz:.6f} xy xz yz")
        else:
            lines.append("0.0 10.0 xlo xhi")
            lines.append("0.0 10.0 ylo yhi")
            lines.append("0.0 10.0 zlo zhi")

    def _write_masses_section(
        self, lines: list[str], frame: Frame, type_system: TypeSystem
    ) -> None:
        """Write masses section using TypeSystem."""
        lines.append("Masses")
        lines.append("")

        atoms_data = frame["atoms"]

        # Require mass field
        if "mass" not in atoms_data:
            raise ValueError("Atoms must have 'mass' field for LAMMPS output")

        atom_type_list = type_system.get_type_list("atom_types")

        # Build type -> mass mapping from actual atoms
        type_to_mass = {}
        for atom_type in np.unique(atoms_data["type"]):
            atom_type_str = str(atom_type)
            mask = atoms_data["type"] == atom_type
            mass = atoms_data["mass"][mask][0]
            type_to_mass[atom_type_str] = mass

        # Write masses for all types
        for atom_type in atom_type_list:
            type_id = type_system.get_id(atom_type, "atom_types")
            mass = type_to_mass.get(atom_type, 1.0)
            lines.append(f"{type_id} {mass:.6f}")

        lines.append("")

    def _write_type_labels_sections(
        self, lines: list[str], type_system: TypeSystem
    ) -> None:
        """Write type labels sections using TypeSystem."""

        # Write sections
        section_configs = [
            ("atom_types", "Atom Type Labels"),
            ("bond_types", "Bond Type Labels"),
            ("angle_types", "Angle Type Labels"),
            ("dihedral_types", "Dihedral Type Labels"),
            ("improper_types", "Improper Type Labels"),
        ]

        for type_key, section_name in section_configs:
            if type_system.has_types(type_key):
                type_list = type_system.get_type_list(type_key)
                if self._needs_type_labels(type_list):
                    lines.append(section_name)
                    lines.append("")
                    for type_idx, type_name in enumerate(type_list):
                        type_id = type_idx + 1
                        lines.append(f"{type_id} {type_name}")
                    lines.append("")

    def _needs_type_labels(self, types: np.ndarray) -> bool:
        """Check if type labels section is needed."""
        # Only write type labels if types are non-numeric (strings)
        # Numeric types (integers) don't need labels
        if len(types) == 0:
            return False
        # Check if any type is a string (not numeric)
        if isinstance(types, list):
            # For lists from TypeSystem
            for t in types:
                try:
                    int(t)
                except (ValueError, TypeError):
                    return True
            return False
        # For numpy arrays
        return types.dtype.kind in ("U", "S", "O")  # Unicode, byte string, or object

    def _write_force_field_coeffs_sections(
        self, lines: list[str], frame: Frame
    ) -> None:
        """Write force field coefficients sections."""
        forcefield = frame.metadata.get("forcefield")
        if not forcefield:
            return

        # Import style classes
        from molpy import (
            AngleStyle,
            BondStyle,
            DihedralStyle,
            ImproperStyle,
            PairStyle,
            Type,
        )

        # Write pair coefficients
        pair_styles = forcefield.get_styles(PairStyle)
        if pair_styles:
            lines.append("Pair Coeffs")
            lines.append("")
            for style in pair_styles:
                for type_obj in style.types.bucket(Type):
                    type_id = int(type_obj.name.split("_")[1])
                    epsilon = type_obj.get("epsilon", 0.0)
                    sigma = type_obj.get("sigma", 1.0)
                    lines.append(f"{type_id} {epsilon:.6f} {sigma:.6f}")
            lines.append("")

        # Write bond coefficients
        bond_styles = forcefield.get_styles(BondStyle)
        if bond_styles:
            lines.append("Bond Coeffs")
            lines.append("")
            for style in bond_styles:
                for type_obj in style.types.bucket(Type):
                    type_id = int(type_obj.name.split("_")[1])
                    k = type_obj.get("k", 0.0)
                    r0 = type_obj.get("r0", 1.0)
                    lines.append(f"{type_id} {k:.6f} {r0:.6f}")
            lines.append("")

        # Write angle coefficients
        angle_styles = forcefield.get_styles(AngleStyle)
        if angle_styles:
            lines.append("Angle Coeffs")
            lines.append("")
            for style in angle_styles:
                for type_obj in style.types.bucket(Type):
                    type_id = int(type_obj.name.split("_")[1])
                    k = type_obj.get("k", 0.0)
                    theta0 = type_obj.get("theta0", 0.0)
                    lines.append(f"{type_id} {k:.6f} {theta0:.6f}")
            lines.append("")

        # Write dihedral coefficients
        dihedral_styles = forcefield.get_styles(DihedralStyle)
        if dihedral_styles:
            lines.append("Dihedral Coeffs")
            lines.append("")
            for style in dihedral_styles:
                for type_obj in style.types.bucket(Type):
                    type_id = int(type_obj.name.split("_")[1])
                    k = type_obj.get("k", 0.0)
                    d = type_obj.get("d", 1)
                    n = type_obj.get("n", 1)
                    lines.append(f"{type_id} {k:.6f} {d} {n}")
            lines.append("")

        # Write improper coefficients
        improper_styles = forcefield.get_styles(ImproperStyle)
        if improper_styles:
            lines.append("Improper Coeffs")
            lines.append("")
            for style in improper_styles:
                for type_obj in style.types.bucket(Type):
                    type_id = int(type_obj.name.split("_")[1])
                    k = type_obj.get("k", 0.0)
                    d = type_obj.get("d", 1)
                    n = type_obj.get("n", 1)
                    lines.append(f"{type_id} {k:.6f} {d} {n}")
            lines.append("")

    def _write_atoms_section(
        self, lines: list[str], frame: Frame, type_system: TypeSystem
    ) -> None:
        """Write atoms section.

        Uses merged type labels to ensure type_id consistency.
        Requires that atoms have an 'id' field.

        Args:
            lines: List of lines to append to
            frame: Frame containing atoms data

        Raises:
            ValueError: If atoms are missing 'id' field
        """
        lines.append("Atoms")
        lines.append("")

        atoms_data = frame["atoms"]

        # Require that all atoms have an 'id' field
        if "id" not in atoms_data:
            raise ValueError(
                "Atoms in frame must have 'id' field. "
                "This field is required for LAMMPS output to map indices to atom IDs."
            )

        n_atoms = len(atoms_data["type"])

        for idx in range(n_atoms):
            # Use atom ID from the 'id' field
            atom_id = int(atoms_data["id"][idx])
            atom_type_str = str(atoms_data["type"][idx])
            atom_type = type_system.get_id(atom_type_str, "atom_types")

            # Get coordinates - must use separate x, y, z fields
            x = float(atoms_data["x"][idx])
            y = float(atoms_data["y"][idx])
            z = float(atoms_data["z"][idx])

            if self.atom_style == "full":
                mol_id = int(atoms_data["mol"][idx])
                charge = float(atoms_data["q"][idx])
                lines.append(
                    f"{atom_id} {mol_id} {atom_type} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}"
                )
            elif self.atom_style == "charge":
                charge = float(atoms_data["q"][idx])
                lines.append(
                    f"{atom_id} {atom_type} {charge:.6f} {x:.6f} {y:.6f} {z:.6f}"
                )
            else:  # atomic
                lines.append(f"{atom_id} {atom_type} {x:.6f} {y:.6f} {z:.6f}")

        lines.append("")

    def _validate_and_convert_indices(
        self,
        data: Block,
        section_name: str,
        idx: int,
        index_to_id: dict[int, int],
        id_to_index: dict[int, int],
    ) -> list[int]:
        """Validate and convert atom indices to IDs for connectivity.

        Args:
            data: Connectivity data block
            section_name: Name of section (bonds, angles, etc.)
            idx: Current item index
            index_to_id: Mapping from 0-based index to atom ID
            id_to_index: Mapping from atom ID to 0-based index

        Returns:
            List of atom IDs

        Raises:
            ValueError: If indices are out of range
        """
        fields = CONNECTIVITY_FIELDS[section_name]
        atom_ids = []

        for field in fields:
            atom_val = int(data[field][idx])

            # Special handling for bonds: support both 0-based indices and 1-based IDs
            if section_name == "bonds":
                if atom_val in index_to_id:
                    # 0-based index: convert to ID
                    atom_ids.append(index_to_id[atom_val])
                elif atom_val in id_to_index:
                    # 1-based ID: use directly (non-standard but supported)
                    atom_ids.append(atom_val)
                else:
                    raise ValueError(
                        f"Bond {idx + 1}: {field} value {atom_val} is neither a valid "
                        f"0-based index (0-{len(index_to_id) - 1}) nor a valid 1-based atom ID."
                    )
            else:
                # Other connectivity types: only support 0-based indices
                if atom_val not in index_to_id:
                    raise ValueError(
                        f"{section_name.capitalize()} {idx + 1}: {field} index {atom_val} is out of range. "
                        f"Valid indices: 0-{len(index_to_id) - 1}"
                    )
                atom_ids.append(index_to_id[atom_val])

        return atom_ids

    def _write_connectivity_section(
        self, lines: list[str], frame: Frame, section_name: str, type_system: TypeSystem
    ) -> None:
        """Write connectivity section (bonds, angles, dihedrals, impropers).

        Uses merged type labels to ensure type_id consistency.
        Converts atom indices to atom IDs using the index->ID mapping from atoms.

        Args:
            lines: List of lines to append to
            frame: Frame containing connectivity data
            section_name: Name of the connectivity section (bonds, angles, etc.)

        Raises:
            ValueError: If atoms are missing 'id' field
        """
        lines.append(section_name.capitalize())
        lines.append("")

        atoms_data = frame["atoms"]

        # Require that all atoms have an 'id' field
        if "id" not in atoms_data:
            raise ValueError(
                "Atoms in frame must have 'id' field. "
                "This field is required for LAMMPS output to map indices to atom IDs."
            )

        # Build index to ID mapping
        # Frame uses index (0-based), we need to map to atom ID
        index_to_id = {}
        # Also build ID to index mapping for reverse lookup
        # Bonds may contain either 0-based indices or 1-based IDs
        id_to_index = {}
        # Also build ID set for validation
        atom_ids_set = set()
        for idx in range(len(atoms_data["type"])):
            atom_id = int(atoms_data["id"][idx])
            index_to_id[idx] = atom_id
            id_to_index[atom_id] = idx
            atom_ids_set.add(atom_id)

        data = frame[section_name]

        # Validate that 'type' field exists
        if "type" not in data:
            raise ValueError(
                f"{section_name.capitalize()} data must have 'type' field. "
                f"Available fields: {list(data.keys())}"
            )

        # Get number of items from the data block
        n_items = data.nrows

        # Validate that all required atom index fields exist
        required_fields = CONNECTIVITY_FIELDS[section_name]
        missing_fields = [f for f in required_fields if f not in data]
        if missing_fields:
            raise ValueError(
                f"{section_name.capitalize()} must have {', '.join(required_fields)} fields (0-based atom indices). "
                f"Missing: {', '.join(missing_fields)}"
            )

        # Validate that type field has the same length as atom index fields
        if len(data["type"]) != n_items:
            raise ValueError(
                f"{section_name.capitalize()} 'type' field has {len(data['type'])} values, "
                f"but expected {n_items} (based on atom index fields)"
            )

        for idx in range(n_items):
            item_id = idx + 1
            item_type_str = str(data["type"][idx])
            category = TYPE_KEY_MAPPING[section_name]
            item_type = type_system.get_id(item_type_str, category)

            # Validate and convert indices to IDs
            atom_ids = self._validate_and_convert_indices(
                data, section_name, idx, index_to_id, id_to_index
            )

            # Format output line
            atom_ids_str = " ".join(str(aid) for aid in atom_ids)
            lines.append(f"{item_id} {item_type} {atom_ids_str}")

        lines.append("")
