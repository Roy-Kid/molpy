import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any, Callable, Dict, List, Optional

import molpy as mp

logger = logging.getLogger("molpy")


class XMLForceFieldReader:
    """
    Universal XML force field parser

    Features:
    1. Manual handling of special tags (e.g., <Author>)
    2. Universal potential definition processing (e.g., <vdW type="pair">)
    3. Metadata field exclusion (e.g., smirks attribute)
    4. Interchangeable with LAMMPS force field definitions
    """

    def __init__(self, file: Path):
        self._file = file
        self._metadata_fields = {
            "smirks",
            "name",
            "id",
            "type",
            "type1",
            "type2",
            "type3",
            "type4",
        }  # Metadata fields, not force field parameters
        self._special_handlers = {}  # Special tag handlers
        self._potential_handlers = {}  # Potential handlers
        self._atomtypes = {}  # Cache for atom types

        # Register default handlers
        self._register_default_handlers()

    def _register_default_handlers(self):
        """Register default special tag handlers"""
        self._special_handlers = {
            "Author": self._parse_author_section,
            "Version": self._parse_version_section,
            "Date": self._parse_date_section,
            "Description": self._parse_description_section,
            "Units": self._parse_units_section,
        }

        # Register default potential handlers
        self._potential_handlers = {
            "pair": self._parse_pair_potential,
            "bond": self._parse_bond_potential,
            "angle": self._parse_angle_potential,
            "dihedral": self._parse_dihedral_potential,
            "improper": self._parse_improper_potential,
        }

    def register_special_handler(self, tag: str, handler: Callable):
        """Register special tag handler"""
        self._special_handlers[tag] = handler

    def register_potential_handler(self, potential_type: str, handler: Callable):
        """Register potential handler"""
        self._potential_handlers[potential_type] = handler

    def set_metadata_fields(self, fields: set):
        """Set metadata fields"""
        self._metadata_fields.update(fields)

    def read(self, system):
        """Read XML file and parse force field definitions"""
        ff = system.forcefield

        tree = ET.parse(self._file)
        root = tree.getroot()

        logger.info(f"Starting to parse XML force field file: {self._file}")

        # Parse all child elements
        for child in root:
            self._parse_element(child, ff)

        logger.info(f"Completed parsing XML force field file: {self._file}")
        return system

    def _parse_element(self, element: ET.Element, ff: mp.ForceField):
        """Parse single XML element"""
        tag = element.tag

        # Check if it's a special tag
        if tag in self._special_handlers:
            self._special_handlers[tag](element, ff)
            return

        # Check if it's atom types section
        if tag == "AtomTypes":
            self._parse_atomtypes_section(element, ff)
            return

        # Check if it's residues section
        if tag == "Residues":
            self._parse_residues_section(element, ff)
            return

        # Check if it's a potential definition
        potential_type = self._detect_potential_type(element)
        if potential_type and potential_type in self._potential_handlers:
            self._potential_handlers[potential_type](element, ff)
            return

        # Recursively process child elements
        for child in element:
            self._parse_element(child, ff)

    def _detect_potential_type(self, element: ET.Element) -> Optional[str]:
        """Detect potential type based on content analysis"""
        tag = element.tag.lower()
        attrib = element.attrib

        # Check if type attribute is explicitly set
        potential_type = attrib.get("type", "").lower()
        if potential_type in self._potential_handlers:
            return potential_type

        # Analyze child elements to guess potential type
        child_tags = [child.tag.lower() for child in element]

        # Check for pair potential indicators
        if any(keyword in tag for keyword in ["vdw", "nonbond", "pair", "lj"]):
            return "pair"
        if any(child in ["atom", "pair"] for child in child_tags):
            return "pair"

        # Check for bond potential indicators
        if "bond" in tag:
            return "bond"
        if any(child in ["bond"] for child in child_tags):
            return "bond"

        # Check for angle potential indicators
        if "angle" in tag:
            return "angle"
        if any(child in ["angle"] for child in child_tags):
            return "angle"

        # Check for dihedral potential indicators
        if any(keyword in tag for keyword in ["dihedral", "torsion"]):
            return "dihedral"
        if any(child in ["dihedral", "torsion", "proper"] for child in child_tags):
            return "dihedral"

        # Check for improper potential indicators
        if "improper" in tag:
            return "improper"
        if any(child in ["improper"] for child in child_tags):
            return "improper"

        return None

    def _extract_parameters(self, element: ET.Element) -> Dict[str, Any]:
        """Extract force field parameters, excluding metadata fields"""
        params = {}
        for key, value in element.attrib.items():
            if key not in self._metadata_fields:
                # Try to convert to numeric value
                try:
                    params[key] = float(value)
                except ValueError:
                    params[key] = value
        return params

    def _extract_metadata(self, element: ET.Element) -> Dict[str, Any]:
        """Extract metadata"""
        metadata = {}
        for key, value in element.attrib.items():
            if key in self._metadata_fields:
                metadata[key] = value
        return metadata

    # Special tag handlers
    def _parse_author_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse author information"""
        author = element.text.strip() if element.text else ""
        logger.info(f"Force field author: {author}")
        # Store author information in force field metadata
        if not hasattr(ff, "_metadata"):
            ff._metadata = {}
        ff._metadata["author"] = author

    def _parse_version_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse version information"""
        version = element.text.strip() if element.text else ""
        logger.info(f"Force field version: {version}")
        if not hasattr(ff, "_metadata"):
            ff._metadata = {}
        ff._metadata["version"] = version

    def _parse_date_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse date information"""
        date = element.text.strip() if element.text else ""
        logger.info(f"Force field date: {date}")
        if not hasattr(ff, "_metadata"):
            ff._metadata = {}
        ff._metadata["date"] = date

    def _parse_description_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse description information"""
        description = element.text.strip() if element.text else ""
        logger.info(f"Force field description: {description}")
        if not hasattr(ff, "_metadata"):
            ff._metadata = {}
        ff._metadata["description"] = description

    def _parse_units_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse units information"""
        units = element.text.strip() if element.text else ""
        logger.info(f"Force field units: {units}")
        ff.units = units

    # Potential handlers
    def _parse_pair_potential(self, element: ET.Element, ff: mp.ForceField):
        """Parse pair potential parameters"""
        # Get or create pair style
        style_name = element.tag
        pairstyle = ff.def_pairstyle(style_name, element.attrib)

        # Process child elements
        for child in element:
            if child.tag == "Atom":
                self._parse_atom_pair(child, pairstyle)
            elif child.tag == "Pair":
                self._parse_pair_definition(child, pairstyle)
            elif child.tag == "UseAttributeFromResidue":
                # Handle residue charge usage
                logger.info(f"Using residue charges for {element.tag}")

        logger.info(f"Parsed {len(element)} pair potential parameters")

    def _parse_bond_potential(self, element: ET.Element, ff: mp.ForceField):
        """Parse bond potential parameters"""
        style_name = element.tag
        bondstyle = ff.def_bondstyle(style_name, element.attrib)

        for child in element:
            if child.tag == "Bond":
                self._parse_bond_definition(child, bondstyle)
            else:
                # Handle direct bond definitions (like in tip3p.xml)
                self._parse_bond_definition(child, bondstyle)

        logger.info(f"Parsed {len(element)} bond potential parameters")

    def _parse_angle_potential(self, element: ET.Element, ff: mp.ForceField):
        """Parse angle potential parameters"""
        style_name = element.tag
        anglestyle = ff.def_anglestyle(style_name, element.attrib)

        for child in element:
            if child.tag == "Angle":
                self._parse_angle_definition(child, anglestyle)
            else:
                # Handle direct angle definitions (like in tip3p.xml)
                self._parse_angle_definition(child, anglestyle)

        logger.info(f"Parsed {len(element)} angle potential parameters")

    def _parse_dihedral_potential(self, element: ET.Element, ff: mp.ForceField):
        """Parse dihedral potential parameters"""
        style_name = element.tag
        dihedralstyle = ff.def_dihedralstyle(style_name, element.attrib)

        for child in element:
            if child.tag in ["Dihedral", "Torsion", "Proper"]:
                self._parse_dihedral_definition(child, dihedralstyle)

        logger.info(f"Parsed {len(element)} dihedral potential parameters")

    def _parse_improper_potential(self, element: ET.Element, ff: mp.ForceField):
        """Parse improper potential parameters"""
        style_name = element.tag
        improperstyle = ff.def_improperstyle(style_name, element.attrib)

        for child in element:
            if child.tag == "Improper":
                self._parse_improper_definition(child, improperstyle)

        logger.info(f"Parsed {len(element)} improper potential parameters")

    # Specific definition parsers
    def _parse_atom_pair(self, element: ET.Element, pairstyle):
        """Parse atom pair definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Get atom type identifier
        atom_id = self._get_atom_identifier(metadata)
        if atom_id:
            # Get or create atom type
            atomtype = self._get_or_create_atomtype(atom_id)
            pairstyle.def_type(atomtype, atomtype, **params)

    def _parse_pair_definition(self, element: ET.Element, pairstyle):
        """Parse pair definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Get two atom types
        type1 = metadata.get("type1", metadata.get("class1"))
        type2 = metadata.get("type2", metadata.get("class2"))

        if type1 and type2:
            atomtype1 = self._get_or_create_atomtype(type1)
            atomtype2 = self._get_or_create_atomtype(type2)
            pairstyle.def_type(atomtype1, atomtype2, **params)

    def _parse_bond_definition(self, element: ET.Element, bondstyle):
        """Parse bond definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Extract type information from metadata
        type1 = metadata.pop("type1", metadata.pop("class1", None))
        type2 = metadata.pop("type2", metadata.pop("class2", None))

        # Map parameter names if needed
        if "length" in params:
            params["r0"] = params.pop("length")

        if type1 and type2:
            atomtype1 = self._get_or_create_atomtype(type1)
            atomtype2 = self._get_or_create_atomtype(type2)
            bondstyle.def_type(atomtype1, atomtype2, **params)

    def _parse_angle_definition(self, element: ET.Element, anglestyle):
        """Parse angle definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Extract type information from metadata
        type1 = metadata.pop("type1", metadata.pop("class1", None))
        type2 = metadata.pop("type2", metadata.pop("class2", None))
        type3 = metadata.pop("type3", metadata.pop("class3", None))

        # Map parameter names if needed
        if "angle" in params:
            params["theta0"] = params.pop("angle")

        if type1 and type2 and type3:
            atomtype1 = self._get_or_create_atomtype(type1)
            atomtype2 = self._get_or_create_atomtype(type2)
            atomtype3 = self._get_or_create_atomtype(type3)
            anglestyle.def_type(atomtype1, atomtype2, atomtype3, **params)

    def _parse_dihedral_definition(self, element: ET.Element, dihedralstyle):
        """Parse dihedral definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Extract type information from metadata
        type1 = metadata.pop("type1", metadata.pop("class1", None))
        type2 = metadata.pop("type2", metadata.pop("class2", None))
        type3 = metadata.pop("type3", metadata.pop("class3", None))
        type4 = metadata.pop("type4", metadata.pop("class4", None))

        if type1 and type2 and type3 and type4:
            atomtype1 = self._get_or_create_atomtype(type1)
            atomtype2 = self._get_or_create_atomtype(type2)
            atomtype3 = self._get_or_create_atomtype(type3)
            atomtype4 = self._get_or_create_atomtype(type4)
            dihedralstyle.def_type(atomtype1, atomtype2, atomtype3, atomtype4, **params)

    def _parse_improper_definition(self, element: ET.Element, improperstyle):
        """Parse improper definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Extract type information from metadata
        type1 = metadata.pop("type1", metadata.pop("class1", None))
        type2 = metadata.pop("type2", metadata.pop("class2", None))
        type3 = metadata.pop("type3", metadata.pop("class3", None))
        type4 = metadata.pop("type4", metadata.pop("class4", None))

        if type1 and type2 and type3 and type4:
            atomtype1 = self._get_or_create_atomtype(type1)
            atomtype2 = self._get_or_create_atomtype(type2)
            atomtype3 = self._get_or_create_atomtype(type3)
            atomtype4 = self._get_or_create_atomtype(type4)
            improperstyle.def_type(atomtype1, atomtype2, atomtype3, atomtype4, **params)

    def _get_atom_identifier(self, metadata: Dict[str, Any]) -> Optional[str]:
        """Get atom type identifier"""
        # Try different identifier fields in order of priority
        for key in ["name", "type", "id", "class", "smirks"]:
            if key in metadata:
                return metadata[key]
        return None

    def _get_or_create_atomtype(self, atom_id: str) -> mp.AtomType:
        """Get or create atom type"""
        if atom_id in self._atomtypes:
            return self._atomtypes[atom_id]

        # Create new atom type
        atomtype = mp.AtomType(atom_id)
        self._atomtypes[atom_id] = atomtype
        return atomtype

    def _parse_atomtypes_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse atom types section"""
        atomstyle = ff.def_atomstyle("default")

        for child in element:
            if child.tag in ["Type", "Atom"]:
                self._parse_atomtype_definition(child, atomstyle)

        logger.info(f"Parsed {len(element)} atom type definitions")

    def _parse_atomtype_definition(self, element: ET.Element, atomstyle):
        """Parse atom type definition"""
        metadata = self._extract_metadata(element)
        params = self._extract_parameters(element)

        # Get atom type identifier
        atom_id = self._get_atom_identifier(metadata)
        if atom_id:
            # Create atom type
            atomtype = atomstyle.def_type(atom_id, **params)
            self._atomtypes[atom_id] = atomtype

    def _parse_residues_section(self, element: ET.Element, ff: mp.ForceField):
        """Parse residues section"""
        for child in element:
            if child.tag == "Residue":
                self._parse_residue_definition(child, ff)

        logger.info(f"Parsed {len(element)} residue definitions")

    def _parse_residue_definition(self, element: ET.Element, ff: mp.ForceField):
        """Parse residue definition"""
        residue_name = element.attrib.get("name", "UNK")
        logger.info(f"Parsing residue: {residue_name}")

        # Parse atoms in the residue
        for atom in element.findall("Atom"):
            atom_name = atom.attrib.get("name")
            atom_type = atom.attrib.get("type")
            charge = atom.attrib.get("charge")

            # Store residue-specific information
            if atom_type and atom_type in self._atomtypes:
                atom_obj = self._atomtypes[atom_type]
                if charge is not None:
                    # Store charge in residue context
                    if not hasattr(atom_obj, "_residue_charges"):
                        atom_obj._residue_charges = {}
                    atom_obj._residue_charges[residue_name] = float(charge)
