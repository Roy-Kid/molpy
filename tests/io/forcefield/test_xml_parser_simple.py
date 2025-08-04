#!/usr/bin/env python3
"""
Simple tests for XML force field parser
"""

import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../'))

import sys
sys.path.insert(0, 'molpy/src')
import molpy as mp
from molpy.io.forcefield.xml import XMLForceFieldReader
from pathlib import Path


def test_basic_initialization():
    """Test basic initialization of XMLForceFieldReader."""
    print("Testing basic initialization...")
    xml_file = Path("molpy/src/molpy/data/forcefield/tip3p.xml")
    reader = XMLForceFieldReader(xml_file)
    
    assert reader._file == xml_file
    assert 'smirks' in reader._metadata_fields
    assert 'type1' in reader._metadata_fields
    assert 'Author' in reader._special_handlers
    assert 'pair' in reader._potential_handlers
    print("✓ Basic initialization test passed")


def test_metadata_field_exclusion():
    """Test metadata field exclusion functionality."""
    print("Testing metadata field exclusion...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    # Test default metadata fields
    assert 'smirks' in reader._metadata_fields
    assert 'type1' in reader._metadata_fields
    
    # Test adding custom metadata fields
    reader.set_metadata_fields({'charge', 'mass'})
    assert 'charge' in reader._metadata_fields
    assert 'mass' in reader._metadata_fields
    print("✓ Metadata field exclusion test passed")


def test_special_handler_registration():
    """Test special handler registration."""
    print("Testing special handler registration...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    def custom_handler(element, ff):
        pass
    
    reader.register_special_handler('CustomTag', custom_handler)
    assert 'CustomTag' in reader._special_handlers
    assert reader._special_handlers['CustomTag'] == custom_handler
    print("✓ Special handler registration test passed")


def test_potential_handler_registration():
    """Test potential handler registration."""
    print("Testing potential handler registration...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    def custom_potential_handler(element, ff):
        pass
    
    reader.register_potential_handler('custom_type', custom_potential_handler)
    assert 'custom_type' in reader._potential_handlers
    assert reader._potential_handlers['custom_type'] == custom_potential_handler
    print("✓ Potential handler registration test passed")


def test_potential_type_detection():
    """Test potential type detection logic."""
    print("Testing potential type detection...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    # Test type attribute detection
    import xml.etree.ElementTree as ET
    element = ET.Element('Test', type='pair')
    assert reader._detect_potential_type(element) == 'pair'
    
    # Test tag name detection
    element = ET.Element('vdW')
    assert reader._detect_potential_type(element) == 'pair'
    
    element = ET.Element('Bonds')
    assert reader._detect_potential_type(element) == 'bond'
    
    element = ET.Element('Angles')
    assert reader._detect_potential_type(element) == 'angle'
    
    element = ET.Element('Dihedrals')
    assert reader._detect_potential_type(element) == 'dihedral'
    
    element = ET.Element('Impropers')
    assert reader._detect_potential_type(element) == 'improper'
    print("✓ Potential type detection test passed")


def test_tip3p_parsing():
    """Test parsing of tip3p.xml file."""
    print("Testing tip3p.xml parsing...")
    xml_file = Path("molpy/src/molpy/data/forcefield/tip3p.xml")
    reader = XMLForceFieldReader(xml_file)
    
    # Create system
    system = mp.FrameSystem()
    
    # Parse XML file
    system = reader.read(system)
    
    # Check force field
    ff = system.forcefield
    assert ff.units == "real"
    
    # Check atom types
    assert len(ff.atomstyles) == 1
    atomstyle = ff.atomstyles[0]
    assert len(atomstyle.types) == 2
    
    # Check for specific atom types
    atom_types = [at.name for at in atomstyle.types]
    assert 'tip3p-O' in atom_types
    assert 'tip3p-H' in atom_types
    
    # Check pair styles
    assert len(ff.pairstyles) == 1
    pairstyle = ff.pairstyles[0]
    assert pairstyle.name == "NonbondedForce"
    assert len(pairstyle.types) == 2
    
    # Check bond styles
    assert len(ff.bondstyles) == 1
    bondstyle = ff.bondstyles[0]
    assert bondstyle.name == "HarmonicBondForce"
    assert len(bondstyle.types) == 1
    
    # Check angle styles
    assert len(ff.anglestyles) == 1
    anglestyle = ff.anglestyles[0]
    assert anglestyle.name == "HarmonicAngleForce"
    assert len(anglestyle.types) == 1
    
    # Check dihedral and improper styles (should be empty for tip3p)
    assert len(ff.dihedralstyles) == 0
    assert len(ff.improperstyles) == 0
    print("✓ tip3p.xml parsing test passed")


def test_parameter_mapping():
    """Test parameter name mapping functionality."""
    print("Testing parameter mapping...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    # Test bond parameter mapping
    import xml.etree.ElementTree as ET
    bond_element = ET.Element('Bond', type1='C', type2='H', length='1.0', k='100.0')
    
    metadata = reader._extract_metadata(bond_element)
    params = reader._extract_parameters(bond_element)
    
    # Check that type information is in metadata
    assert 'type1' in metadata
    assert 'type2' in metadata
    
    # Check that parameters are extracted correctly
    assert 'length' in params
    assert 'k' in params
    
    # Test angle parameter mapping
    angle_element = ET.Element('Angle', type1='C', type2='O', type3='H', angle='120.0', k='50.0')
    
    metadata = reader._extract_metadata(angle_element)
    params = reader._extract_parameters(angle_element)
    
    # Check that type information is in metadata
    assert 'type1' in metadata
    assert 'type2' in metadata
    assert 'type3' in metadata
    
    # Check that parameters are extracted correctly
    assert 'angle' in params
    assert 'k' in params
    print("✓ Parameter mapping test passed")


def test_atom_type_creation():
    """Test atom type creation and caching."""
    print("Testing atom type creation...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    # Test creating new atom type
    atomtype1 = reader._get_or_create_atomtype("C")
    assert atomtype1.name == "C"
    assert "C" in reader._atomtypes
    
    # Test retrieving existing atom type
    atomtype2 = reader._get_or_create_atomtype("C")
    assert atomtype1 is atomtype2  # Should be the same object
    
    # Test creating different atom type
    atomtype3 = reader._get_or_create_atomtype("H")
    assert atomtype3.name == "H"
    assert atomtype3 is not atomtype1
    print("✓ Atom type creation test passed")


def test_metadata_extraction():
    """Test metadata and parameter extraction."""
    print("Testing metadata extraction...")
    reader = XMLForceFieldReader(Path("dummy.xml"))
    
    import xml.etree.ElementTree as ET
    element = ET.Element('Test', 
                       type1='C', type2='H', 
                       k='100.0', r0='1.0',
                       smirks='[#6]', name='test')
    
    metadata = reader._extract_metadata(element)
    params = reader._extract_parameters(element)
    
    # Check metadata extraction
    assert 'type1' in metadata
    assert 'type2' in metadata
    assert 'smirks' in metadata
    assert 'name' in metadata
    
    # Check parameter extraction
    assert 'k' in params
    assert 'r0' in params
    assert params['k'] == 100.0
    assert params['r0'] == 1.0
    
    # Check that metadata fields are not in parameters
    assert 'type1' not in params
    assert 'type2' not in params
    assert 'smirks' not in params
    assert 'name' not in params
    print("✓ Metadata extraction test passed")


def run_all_tests():
    """Run all tests."""
    print("=== Running XML Force Field Parser Tests ===\n")
    
    try:
        test_basic_initialization()
        test_metadata_field_exclusion()
        test_special_handler_registration()
        test_potential_handler_registration()
        test_potential_type_detection()
        test_tip3p_parsing()
        test_parameter_mapping()
        test_atom_type_creation()
        test_metadata_extraction()
        
        print("\n=== All tests passed! ===")
        return True
    except Exception as e:
        print(f"\n❌ Test failed: {e}")
        import traceback
        traceback.print_exc()
        return False


if __name__ == "__main__":
    success = run_all_tests()
    sys.exit(0 if success else 1) 