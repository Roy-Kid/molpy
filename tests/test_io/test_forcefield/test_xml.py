"""
Unit tests for XML force field parser.

Tests the XMLForceFieldReader with the OPLS-AA force field.
"""

import sys
from pathlib import Path

import pytest

# Add src to path
sys.path.insert(0, str(Path(__file__).parents[3] / "src"))

# Now we can import directly
from molpy.io.forcefield.xml import XMLForceFieldReader, read_xml_forcefield
from molpy.core.forcefield import (
    AtomisticForcefield,
    AtomType,
    BondType,
    AngleType,
    DihedralType,
    PairType,
)


# Path to the OPLS-AA force field file
OPLSAA_PATH = Path(__file__).parents[3] / "src" / "molpy" / "data" / "forcefield" / "oplsaa.xml"


class TestXMLForceFieldReader:
    """Test XML force field reader basic functionality."""

    def test_reader_initialization(self):
        """Test creating a reader instance."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        assert reader._file == OPLSAA_PATH
        assert len(reader._atomtype_cache) == 0
        assert reader._ff is None

    def test_file_not_found(self):
        """Test error handling for missing file."""
        reader = XMLForceFieldReader(Path("nonexistent.xml"))
        with pytest.raises(FileNotFoundError):
            reader.read()

    def test_read_creates_forcefield(self):
        """Test that reading creates a force field."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read()
        assert isinstance(ff, AtomisticForcefield)
        assert ff.name == "OPLS-AA"

    def test_read_with_existing_forcefield(self):
        """Test reading into an existing force field."""
        existing_ff = AtomisticForcefield(name="Custom", units="real")
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read(forcefield=existing_ff)
        assert ff is existing_ff
        # Should have populated the existing FF
        assert len(ff.get_atomtypes()) > 0


class TestOPLSAAAtomTypes:
    """Test parsing of OPLS-AA atom types."""

    @pytest.fixture
    def oplsaa_ff(self):
        """Load OPLS-AA force field."""
        return read_xml_forcefield(OPLSAA_PATH)

    def test_atomtypes_loaded(self, oplsaa_ff):
        """Test that atom types are loaded."""
        atom_types = oplsaa_ff.get_atomtypes()
        assert len(atom_types) > 0, "No atom types loaded"

    def test_specific_atomtype(self, oplsaa_ff):
        """Test that specific atom types exist."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read()
        
        # Check cache for specific types
        assert "opls_001" in reader._atomtype_cache
        assert "opls_135" in reader._atomtype_cache
        
        # Test atom type properties
        opls_001 = reader._atomtype_cache["opls_001"]
        assert opls_001.name == "opls_001"
        assert opls_001.params.kwargs.get("element") == "C"
        assert opls_001.params.kwargs.get("mass") == 12.011

    def test_atomtype_with_class(self, oplsaa_ff):
        """Test atom types with class names."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read()
        
        # opls_135 has class="CT"
        if "opls_135" in reader._atomtype_cache:
            opls_135 = reader._atomtype_cache["opls_135"]
            assert opls_135.params.kwargs.get("class_name") == "CT"
            # Should also be accessible by class name
            assert "CT" in reader._atomtype_cache

    def test_multiple_elements(self, oplsaa_ff):
        """Test that multiple element types are represented."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read()
        
        elements = set()
        for atomtype in reader._atomtype_cache.values():
            elem = atomtype.params.kwargs.get("element")
            if elem:
                elements.add(elem)
        
        # OPLS-AA should have C, H, O, N, S, etc.
        assert "C" in elements
        assert "H" in elements
        assert "O" in elements
        assert "N" in elements


class TestOPLSAABonds:
    """Test parsing of OPLS-AA bond parameters."""

    @pytest.fixture
    def oplsaa_ff(self):
        """Load OPLS-AA force field."""
        return read_xml_forcefield(OPLSAA_PATH)

    def test_bonds_loaded(self, oplsaa_ff):
        """Test that bond types are loaded."""
        bond_types = oplsaa_ff.get_bondtypes()
        assert len(bond_types) > 0, "No bond types loaded"

    def test_bond_parameters(self, oplsaa_ff):
        """Test that bonds have correct parameters."""
        bond_types = oplsaa_ff.get_bondtypes()
        
        # Check that bonds have r0 and k parameters
        for bond in bond_types[:10]:  # Check first 10
            assert "r0" in bond.params.kwargs or "k" in bond.params.kwargs
            
    def test_bond_atomtypes(self, oplsaa_ff):
        """Test that bonds reference atom types."""
        bond_types = oplsaa_ff.get_bondtypes()
        
        # Bonds should have atom types in args
        for bond in bond_types[:5]:
            assert len(bond.params.args) >= 2
            assert all(isinstance(arg, AtomType) for arg in bond.params.args[:2])


class TestOPLSAAAngles:
    """Test parsing of OPLS-AA angle parameters."""

    @pytest.fixture
    def oplsaa_ff(self):
        """Load OPLS-AA force field."""
        return read_xml_forcefield(OPLSAA_PATH)

    def test_angles_loaded(self, oplsaa_ff):
        """Test that angle types are loaded."""
        angle_types = oplsaa_ff.get_angletypes()
        assert len(angle_types) > 0, "No angle types loaded"

    def test_angle_parameters(self, oplsaa_ff):
        """Test that angles have correct parameters."""
        angle_types = oplsaa_ff.get_angletypes()
        
        # Check that angles have theta0 and k parameters
        for angle in angle_types[:10]:
            assert "theta0" in angle.params.kwargs or "k" in angle.params.kwargs

    def test_angle_atomtypes(self, oplsaa_ff):
        """Test that angles reference three atom types."""
        angle_types = oplsaa_ff.get_angletypes()
        
        for angle in angle_types[:5]:
            assert len(angle.params.args) >= 3
            assert all(isinstance(arg, AtomType) for arg in angle.params.args[:3])


class TestOPLSAADihedrals:
    """Test parsing of OPLS-AA dihedral parameters."""

    @pytest.fixture
    def oplsaa_ff(self):
        """Load OPLS-AA force field."""
        return read_xml_forcefield(OPLSAA_PATH)

    def test_dihedrals_loaded(self, oplsaa_ff):
        """Test that dihedral types are loaded."""
        from molpy.core.forcefield import DihedralType, DihedralStyle
        
        dihedral_styles = oplsaa_ff.get_styles(DihedralStyle)
        assert len(dihedral_styles) > 0, "No dihedral styles found"
        
        # Get dihedrals from the first style
        dihedralstyle = dihedral_styles[0]
        dihedrals = dihedralstyle.types.bucket(DihedralType)
        assert len(dihedrals) > 0, "No dihedral types loaded"

    def test_dihedral_rb_coefficients(self, oplsaa_ff):
        """Test that dihedrals have RB coefficients (c0-c5)."""
        from molpy.core.forcefield import DihedralType, DihedralStyle
        
        dihedral_styles = oplsaa_ff.get_styles(DihedralStyle)
        dihedralstyle = dihedral_styles[0]
        dihedrals = dihedralstyle.types.bucket(DihedralType)
        
        # Check first few dihedrals for RB coefficients
        for dihedral in list(dihedrals)[:5]:
            params = dihedral.params.kwargs
            # Should have at least some c coefficients
            has_c_params = any(key.startswith("c") for key in params.keys())
            assert has_c_params, f"Dihedral {dihedral.name} missing RB coefficients"

    def test_dihedral_atomtypes(self, oplsaa_ff):
        """Test that dihedrals reference four atom types."""
        from molpy.core.forcefield import DihedralType, DihedralStyle
        
        dihedral_styles = oplsaa_ff.get_styles(DihedralStyle)
        dihedralstyle = dihedral_styles[0]
        dihedrals = dihedralstyle.types.bucket(DihedralType)
        
        for dihedral in list(dihedrals)[:5]:
            assert len(dihedral.params.args) >= 4
            assert all(isinstance(arg, AtomType) for arg in dihedral.params.args[:4])


class TestOPLSAANonbonded:
    """Test parsing of OPLS-AA nonbonded parameters."""

    @pytest.fixture
    def oplsaa_ff(self):
        """Load OPLS-AA force field."""
        return read_xml_forcefield(OPLSAA_PATH)

    def test_nonbonded_loaded(self, oplsaa_ff):
        """Test that nonbonded pair types are loaded."""
        from molpy.core.forcefield import PairStyle, PairType
        
        pair_styles = oplsaa_ff.get_styles(PairStyle)
        assert len(pair_styles) > 0, "No pair styles found"
        
        pairstyle = pair_styles[0]
        pairs = pairstyle.types.bucket(PairType)
        assert len(pairs) > 0, "No pair types loaded"

    def test_nonbonded_parameters(self, oplsaa_ff):
        """Test that nonbonded types have LJ parameters."""
        from molpy.core.forcefield import PairStyle, PairType
        
        pair_styles = oplsaa_ff.get_styles(PairStyle)
        pairstyle = pair_styles[0]
        pairs = pairstyle.types.bucket(PairType)
        
        # Check first few pairs for parameters
        for pair in list(pairs)[:10]:
            params = pair.params.kwargs
            # Should have sigma, epsilon, and charge
            assert "sigma" in params or "epsilon" in params or "charge" in params

    def test_nonbonded_scaling_factors(self, oplsaa_ff):
        """Test that nonbonded force has scaling factors."""
        from molpy.core.forcefield import PairStyle
        
        pair_styles = oplsaa_ff.get_styles(PairStyle)
        pairstyle = pair_styles[0]
        
        # OPLS-AA uses 0.5 scaling for 1-4 interactions
        params = pairstyle.params.kwargs
        assert "coulomb14scale" in params
        assert "lj14scale" in params
        assert params["coulomb14scale"] == 0.5
        assert params["lj14scale"] == 0.5

    def test_specific_nonbonded_parameters(self, oplsaa_ff):
        """Test specific nonbonded parameters from OPLS-AA."""
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff = reader.read()
        
        from molpy.core.forcefield import PairStyle, PairType
        
        pair_styles = ff.get_styles(PairStyle)
        pairstyle = pair_styles[0]
        pairs = pairstyle.types.bucket(PairType)
        
        # Find opls_001 pair parameters (if exists)
        opls_001_pair = None
        for pair in pairs:
            # Check if this is the opls_001 self-interaction
            if len(pair.params.args) >= 2:
                at1, at2 = pair.params.args[:2]
                if at1.name == "opls_001" and at2.name == "opls_001":
                    opls_001_pair = pair
                    break
        
        if opls_001_pair:
            # opls_001: charge=0.5, sigma=0.375, epsilon=0.43932
            assert opls_001_pair.params.kwargs.get("charge") == 0.5
            assert opls_001_pair.params.kwargs.get("sigma") == 0.375
            assert opls_001_pair.params.kwargs.get("epsilon") == 0.43932


class TestOPLSAACompleteness:
    """Test overall completeness of OPLS-AA parsing."""

    def test_full_parse(self):
        """Test that a complete parse works without errors."""
        ff = read_xml_forcefield(OPLSAA_PATH)
        
        # Should have all major components
        assert len(ff.get_atomtypes()) > 100, "Too few atom types"
        assert len(ff.get_bondtypes()) > 10, "Too few bond types"
        assert len(ff.get_angletypes()) > 10, "Too few angle types"

    def test_force_field_name_and_version(self):
        """Test that force field metadata is correct."""
        ff = read_xml_forcefield(OPLSAA_PATH)
        assert ff.name == "OPLS-AA"
        assert ff.units == "real"

    def test_no_duplicate_atomtypes(self):
        """Test that atom types are unique."""
        ff = read_xml_forcefield(OPLSAA_PATH)
        atom_types = ff.get_atomtypes()
        
        names = [at.name for at in atom_types]
        assert len(names) == len(set(names)), "Duplicate atom type names found"

    def test_reader_convenience_function(self):
        """Test the convenience function works."""
        ff1 = read_xml_forcefield(OPLSAA_PATH)
        
        reader = XMLForceFieldReader(OPLSAA_PATH)
        ff2 = reader.read()
        
        # Both should have same number of atom types
        assert len(ff1.get_atomtypes()) == len(ff2.get_atomtypes())


class TestXMLParserEdgeCases:
    """Test edge cases and error handling."""

    def test_missing_type_name(self, tmp_path):
        """Test handling of atom type without name."""
        xml_content = """<?xml version="1.0"?>
        <ForceField name="Test">
            <AtomTypes>
                <Type element="C" mass="12.011"/>
            </AtomTypes>
        </ForceField>
        """
        xml_file = tmp_path / "test.xml"
        xml_file.write_text(xml_content)
        
        ff = read_xml_forcefield(xml_file)
        # Should skip atom type without name
        assert len(ff.get_atomtypes()) == 0

    def test_missing_bond_types(self, tmp_path):
        """Test handling of bonds without type information."""
        xml_content = """<?xml version="1.0"?>
        <ForceField name="Test">
            <HarmonicBondForce>
                <Bond length="1.5" k="1000.0"/>
            </HarmonicBondForce>
        </ForceField>
        """
        xml_file = tmp_path / "test.xml"
        xml_file.write_text(xml_content)
        
        ff = read_xml_forcefield(xml_file)
        # Should skip bond without type information
        assert len(ff.get_bondtypes()) == 0

    def test_on_the_fly_atomtype_creation(self, tmp_path):
        """Test that missing atom types are created on-the-fly."""
        xml_content = """<?xml version="1.0"?>
        <ForceField name="Test">
            <HarmonicBondForce>
                <Bond class1="CA" class2="CB" length="1.5" k="1000.0"/>
            </HarmonicBondForce>
        </ForceField>
        """
        xml_file = tmp_path / "test.xml"
        xml_file.write_text(xml_content)
        
        reader = XMLForceFieldReader(xml_file)
        ff = reader.read()
        
        # Should have created CA and CB atom types
        assert "CA" in reader._atomtype_cache
        assert "CB" in reader._atomtype_cache
        assert len(ff.get_bondtypes()) == 1
