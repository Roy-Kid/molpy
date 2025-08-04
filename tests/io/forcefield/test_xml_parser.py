#!/usr/bin/env python3
"""
Tests for XML force field parser
"""

import pytest
import sys
import os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../../'))

import molpy as mp
from molpy.io.forcefield.xml import XMLForceFieldReader
from pathlib import Path


class TestXMLForceFieldReader:
    """Test cases for XMLForceFieldReader class."""
    
    def test_basic_initialization(self):
        """Test basic initialization of XMLForceFieldReader."""
        xml_file = Path("molpy/src/molpy/data/forcefield/tip3p.xml")
        reader = XMLForceFieldReader(xml_file)
        
        assert reader._file == xml_file
        assert 'smirks' in reader._metadata_fields
        assert 'type1' in reader._metadata_fields
        assert 'Author' in reader._special_handlers
        assert 'pair' in reader._potential_handlers
    
    def test_metadata_field_exclusion(self):
        """Test metadata field exclusion functionality."""
        reader = XMLForceFieldReader(Path("dummy.xml"))
        
        # Test default metadata fields
        assert 'smirks' in reader._metadata_fields
        assert 'type1' in reader._metadata_fields
        
        # Test adding custom metadata fields
        reader.set_metadata_fields({'charge', 'mass'})
        assert 'charge' in reader._metadata_fields
        assert 'mass' in reader._metadata_fields
    
    def test_special_handler_registration(self):
        """Test special handler registration."""
        reader = XMLForceFieldReader(Path("dummy.xml"))
        
        def custom_handler(element, ff):
            pass
        
        reader.register_special_handler('CustomTag', custom_handler)
        assert 'CustomTag' in reader._special_handlers
        assert reader._special_handlers['CustomTag'] == custom_handler
    
    def test_potential_handler_registration(self):
        """Test potential handler registration."""
        reader = XMLForceFieldReader(Path("dummy.xml"))
        
        def custom_potential_handler(element, ff):
            pass
        
        reader.register_potential_handler('custom_type', custom_potential_handler)
        assert 'custom_type' in reader._potential_handlers
        assert reader._potential_handlers['custom_type'] == custom_potential_handler
    
    def test_potential_type_detection(self):
        """Test potential type detection logic."""
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
    
    def test_tip3p_parsing(self):
        """Test parsing of tip3p.xml file."""
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
    
    def test_parameter_mapping(self):
        """Test parameter name mapping functionality."""
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
    
    def test_atom_type_creation(self):
        """Test atom type creation and caching."""
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
    
    def test_metadata_extraction(self):
        """Test metadata and parameter extraction."""
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


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v"]) 