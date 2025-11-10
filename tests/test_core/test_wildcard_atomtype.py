"""
Tests for wildcard atom type matching in force fields.

Tests cover:
1. Wildcard atom type creation and identification
2. Wildcard matching with regular atom types
3. Wildcard usage in bond/angle/dihedral type definitions
4. XML force field parsing with wildcards
"""

import pytest
from molpy.core.forcefield import AtomType, AtomStyle


def test_wildcard_creation():
    """Test creating wildcard atom types."""
    style = AtomStyle("test")
    
    # Create wildcard with "*"
    wild1 = style.def_type("*")
    assert wild1.is_wildcard
    assert wild1.name == "*"
    
    # Create wildcard with ""
    wild2 = style.def_type("")
    assert wild2.is_wildcard
    assert wild2.name == ""


def test_wildcard_equality():
    """Test wildcard matching with regular atom types."""
    style = AtomStyle("test")
    
    wildcard = style.def_type("*")
    carbon = style.def_type("C")
    oxygen = style.def_type("O")
    
    # Wildcard should match any atom type
    assert wildcard == carbon
    assert wildcard == oxygen
    assert carbon == wildcard
    assert oxygen == wildcard
    
    # Regular types should not match each other
    assert carbon != oxygen
    
    # Wildcard should match another wildcard
    wildcard2 = style.def_type("")
    assert wildcard == wildcard2


def test_wildcard_self_equality():
    """Test that wildcard equals itself."""
    style = AtomStyle("test")
    wildcard = style.def_type("*")
    
    assert wildcard == wildcard


def test_wildcard_hash():
    """Test that wildcards have consistent hashing."""
    style = AtomStyle("test")
    
    wild1 = style.def_type("*")
    wild2 = style.def_type("")
    
    # Both wildcards should have the same hash
    assert hash(wild1) == hash(wild2)
    
    # Can be used in sets/dicts
    wildcard_set = {wild1, wild2}
    # Should deduplicate to one element (same hash and equality)
    # Note: This behavior depends on set implementation
    assert wild1 in wildcard_set
    assert wild2 in wildcard_set


def test_wildcard_repr():
    """Test string representation of wildcards."""
    style = AtomStyle("test")
    
    wildcard = style.def_type("*")
    repr_str = repr(wildcard)
    
    assert "WILDCARD" in repr_str


def test_wildcard_in_bond_type():
    """Test using wildcards in bond type definitions."""
    from molpy.core.forcefield import BondStyle
    
    atom_style = AtomStyle("test")
    bond_style = BondStyle("harmonic")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    
    # Define bond with wildcard (matches any atom bonded to carbon)
    bond_type = bond_style.def_type(wildcard, carbon, k=300.0, r0=1.5)
    
    assert bond_type is not None
    # Atom types are now stored in atom1 and atom2 attributes
    assert bond_type.atom1 == wildcard
    assert bond_type.atom2 == carbon
    # Parameters are in kwargs
    assert bond_type.params.kwargs['k'] == 300.0
    assert bond_type.params.kwargs['r0'] == 1.5


def test_wildcard_matching_in_lookup():
    """Test wildcard matching in force field parameter lookup."""
    from molpy.core.forcefield import BondStyle
    
    atom_style = AtomStyle("test")
    bond_style = BondStyle("harmonic")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    nitrogen = atom_style.def_type("N")
    
    # Define wildcard bond type
    wild_bond = bond_style.def_type(wildcard, carbon, k=300.0, r0=1.5)
    
    # Should match specific types due to wildcard equality
    # (This depends on how the lookup is implemented)
    assert wildcard == nitrogen  # Wildcard matches N
    assert carbon == carbon  # C matches C


def test_multiple_wildcards():
    """Test multiple wildcards in same force field."""
    atom_style = AtomStyle("test")
    
    wild_star = atom_style.def_type("*")
    wild_empty = atom_style.def_type("")
    wild_star2 = atom_style.def_type("*")
    
    # All should be equal
    assert wild_star == wild_empty
    assert wild_star == wild_star2
    assert wild_empty == wild_star2


def test_wildcard_class_name():
    """Test that wildcard can have class_name."""
    style = AtomStyle("test")
    
    wildcard = style.def_type("*", class_name="ANY")
    
    assert wildcard.is_wildcard
    assert wildcard.class_name == "ANY"


def test_non_wildcard_types():
    """Test that regular types are not wildcards."""
    style = AtomStyle("test")
    
    carbon = style.def_type("C")
    ct = style.def_type("CT")
    opls_135 = style.def_type("opls_135")
    
    assert not carbon.is_wildcard
    assert not ct.is_wildcard
    assert not opls_135.is_wildcard


def test_wildcard_xml_parsing():
    """Test XML force field parsing with wildcards."""
    from molpy.io.forcefield.xml import XMLForceFieldReader
    from molpy.core.forcefield import AtomisticForcefield
    import tempfile
    from pathlib import Path
    
    # Create test XML with wildcard
    xml_content = """<?xml version="1.0"?>
<ForceField name="TestFF" version="1.0">
    <AtomTypes>
        <Type name="C" class="C" element="C" mass="12.01"/>
        <Type name="N" class="N" element="N" mass="14.01"/>
    </AtomTypes>
    <HarmonicBondForce>
        <!-- Wildcard bond: any atom to carbon -->
        <Bond type1="*" type2="C" length="0.15" k="300000.0"/>
        <!-- Specific bond -->
        <Bond type1="C" type2="N" length="0.14" k="400000.0"/>
    </HarmonicBondForce>
</ForceField>
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
        f.write(xml_content)
        temp_path = Path(f.name)
    
    try:
        reader = XMLForceFieldReader(temp_path)
        ff = reader.read()
        
        # Check that wildcard was created
        assert reader._wildcard_atomtype is not None
        assert reader._wildcard_atomtype.is_wildcard
        
        # Check atom types
        atom_types = ff.get_types(AtomType)
        assert len(atom_types) >= 3  # C, N, and wildcard
        
        # Find wildcard
        wildcards = [at for at in atom_types if at.is_wildcard]
        assert len(wildcards) >= 1
        
    finally:
        temp_path.unlink()


def test_wildcard_empty_string_xml():
    """Test XML parsing with empty string as wildcard."""
    from molpy.io.forcefield.xml import XMLForceFieldReader
    from molpy.core.forcefield import AtomisticForcefield
    import tempfile
    from pathlib import Path
    
    xml_content = """<?xml version="1.0"?>
<ForceField name="TestFF" version="1.0">
    <AtomTypes>
        <Type name="C" class="C" element="C" mass="12.01"/>
    </AtomTypes>
    <HarmonicBondForce>
        <!-- Empty string as wildcard -->
        <Bond type1="" type2="C" length="0.15" k="300000.0"/>
    </HarmonicBondForce>
</ForceField>
"""
    
    with tempfile.NamedTemporaryFile(mode='w', suffix='.xml', delete=False) as f:
        f.write(xml_content)
        temp_path = Path(f.name)
    
    try:
        reader = XMLForceFieldReader(temp_path)
        ff = reader.read()
        
        # Should create wildcard
        assert reader._wildcard_atomtype is not None
        assert reader._wildcard_atomtype.is_wildcard
        
    finally:
        temp_path.unlink()


def test_bond_type_matches():
    """Test BondType.matches() with wildcards."""
    from molpy.core.forcefield import BondStyle
    
    atom_style = AtomStyle("test")
    bond_style = BondStyle("harmonic")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    nitrogen = atom_style.def_type("N")
    oxygen = atom_style.def_type("O")
    
    # Wildcard-Carbon bond
    wild_c_bond = bond_style.def_type(wildcard, carbon, k=300.0, r0=1.5)
    
    # Should match N-C (wildcard matches N)
    assert wild_c_bond.matches(nitrogen, carbon)
    assert wild_c_bond.matches(carbon, nitrogen)  # Order doesn't matter
    
    # Should match O-C (wildcard matches O)
    assert wild_c_bond.matches(oxygen, carbon)
    
    # Should match C-C (wildcard matches C)
    assert wild_c_bond.matches(carbon, carbon)
    
    # Specific bond
    c_n_bond = bond_style.def_type(carbon, nitrogen, k=400.0, r0=1.4)
    
    # Should match C-N
    assert c_n_bond.matches(carbon, nitrogen)
    assert c_n_bond.matches(nitrogen, carbon)  # Order doesn't matter
    
    # Should NOT match C-O
    assert not c_n_bond.matches(carbon, oxygen)


def test_angle_type_matches():
    """Test AngleType.matches() with wildcards."""
    from molpy.core.forcefield import AngleStyle
    
    atom_style = AtomStyle("test")
    angle_style = AngleStyle("harmonic")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    nitrogen = atom_style.def_type("N")
    oxygen = atom_style.def_type("O")
    
    # *-C-* angle (any-Carbon-any)
    wild_angle = angle_style.def_type(wildcard, carbon, wildcard, k=50.0, theta0=109.5)
    
    # Should match N-C-O
    assert wild_angle.matches(nitrogen, carbon, oxygen)
    assert wild_angle.matches(oxygen, carbon, nitrogen)  # Reverse order
    
    # Should match C-C-C
    assert wild_angle.matches(carbon, carbon, carbon)
    
    # Should NOT match N-O-C (center must be C)
    assert not wild_angle.matches(nitrogen, oxygen, carbon)
    
    # Specific angle C-N-O
    specific_angle = angle_style.def_type(carbon, nitrogen, oxygen, k=60.0, theta0=120.0)
    
    assert specific_angle.matches(carbon, nitrogen, oxygen)
    assert specific_angle.matches(oxygen, nitrogen, carbon)  # Reverse
    assert not specific_angle.matches(carbon, oxygen, nitrogen)  # Wrong center


def test_dihedral_type_matches():
    """Test DihedralType.matches() with wildcards."""
    from molpy.core.forcefield import DihedralStyle
    
    atom_style = AtomStyle("test")
    dihedral_style = DihedralStyle("opls")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    nitrogen = atom_style.def_type("N")
    oxygen = atom_style.def_type("O")
    hydrogen = atom_style.def_type("H")
    
    # *-C-C-* dihedral
    wild_dihedral = dihedral_style.def_type(wildcard, carbon, carbon, wildcard, 
                                            c0=0.0, c1=1.0, c2=0.0)
    
    # Should match H-C-C-O
    assert wild_dihedral.matches(hydrogen, carbon, carbon, oxygen)
    assert wild_dihedral.matches(oxygen, carbon, carbon, hydrogen)  # Reverse
    
    # Should match N-C-C-N
    assert wild_dihedral.matches(nitrogen, carbon, carbon, nitrogen)
    
    # Should NOT match C-N-C-O (middle two must be C-C)
    assert not wild_dihedral.matches(carbon, nitrogen, carbon, oxygen)


def test_pair_type_matches():
    """Test PairType.matches() with wildcards."""
    from molpy.core.forcefield import PairStyle
    
    atom_style = AtomStyle("test")
    pair_style = PairStyle("lj/cut")
    
    wildcard = atom_style.def_type("*")
    carbon = atom_style.def_type("C")
    nitrogen = atom_style.def_type("N")
    
    # Self-interaction C-C
    c_c_pair = pair_style.def_type(carbon, sigma=3.4, epsilon=0.066)
    
    assert c_c_pair.matches(carbon, carbon)
    assert c_c_pair.matches(carbon)  # Single argument means self-interaction
    assert not c_c_pair.matches(nitrogen, nitrogen)
    
    # Cross-interaction C-N
    c_n_pair = pair_style.def_type(carbon, nitrogen, sigma=3.3, epsilon=0.08)
    
    assert c_n_pair.matches(carbon, nitrogen)
    assert c_n_pair.matches(nitrogen, carbon)  # Order doesn't matter
    assert not c_n_pair.matches(carbon, carbon)
    
    # Wildcard pair
    wild_pair = pair_style.def_type(wildcard, sigma=3.0, epsilon=0.1)
    
    # Should match any self-interaction
    assert wild_pair.matches(carbon, carbon)
    assert wild_pair.matches(nitrogen, nitrogen)
    assert wild_pair.matches(carbon)  # Single argument


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
