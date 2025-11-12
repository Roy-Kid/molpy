"""
Test OPLS type dependencies (e.g., opls_155 depends on opls_154).

This tests the layered typing engine's ability to handle type references
in SMARTS patterns.
"""

import pytest
from molpy.core.atomistic import Atomistic
from molpy.io.forcefield.xml import read_xml_forcefield
from molpy.typifier.atomistic import OplsAtomTypifier


@pytest.fixture(scope="module")
def oplsaa_forcefield():
    """Load OPLS-AA force field once for all tests."""
    return read_xml_forcefield("oplsaa.xml")


@pytest.fixture(scope="module")
def oplsaa_typifier(oplsaa_forcefield):
    """Create OPLS-AA typifier once for all tests."""
    return OplsAtomTypifier(oplsaa_forcefield)


def test_opls_154_155_dependency(oplsaa_typifier):
    """Test that opls_155 (H on alcohol O) depends on opls_154 (alcohol O).
    
    opls_154: [O;X2](H)([!H]) - alcohol oxygen
    opls_155: H[O;%opls_154] - hydrogen on opls_154 oxygen
    
    This tests the dependency resolution: opls_154 must be assigned
    before opls_155 can be matched.
    """
    # Check that patterns exist
    pattern_dict = oplsaa_typifier.pattern_dict
    assert "opls_154" in pattern_dict, "opls_154 pattern should exist"
    assert "opls_155" in pattern_dict, "opls_155 pattern should exist"
    
    # Check dependency
    opls_154_pattern = pattern_dict["opls_154"]
    opls_155_pattern = pattern_dict["opls_155"]
    
    assert "opls_154" in opls_155_pattern.dependencies, \
        f"opls_155 should depend on opls_154, got dependencies: {opls_155_pattern.dependencies}"
    
    # Check levels
    assert opls_154_pattern.level is not None, "opls_154 should have a level"
    assert opls_155_pattern.level is not None, "opls_155 should have a level"
    assert opls_155_pattern.level > opls_154_pattern.level, \
        f"opls_155 (level {opls_155_pattern.level}) should be after opls_154 (level {opls_154_pattern.level})"


def test_ethanol_opls_154_155(oplsaa_typifier):
    """Test typing ethanol molecule to get opls_154 and opls_155.
    
    Ethanol structure: CH3-CH2-OH
    - O should be typed as opls_154 (alcohol oxygen)
    - H (on O) should be typed as opls_155 (alcohol hydrogen)
    """
    # Build ethanol: CH3-CH2-OH
    mol = Atomistic()
    
    c1 = mol.add_atom(symbol="C", number=6)  # CH3
    c2 = mol.add_atom(symbol="C", number=6)  # CH2  
    o = mol.add_atom(symbol="O", number=8)   # OH oxygen
    h = mol.add_atom(symbol="H", number=1)   # OH hydrogen
    
    mol.add_bond(c1, c2, order=1)
    mol.add_bond(c2, o, order=1)
    mol.add_bond(o, h, order=1)
    
    # Apply typifier
    oplsaa_typifier.typify(mol)
    
    # Check results
    o_type = o.get("type")
    h_type = h.get("type")
    
    assert o_type == "opls_154", \
        f"Oxygen should be opls_154, got {o_type}"
    assert h_type == "opls_155", \
        f"Hydrogen on alcohol O should be opls_155, got {h_type}"


def test_methanol_opls_154_155(oplsaa_typifier):
    """Test typing methanol molecule: CH3-OH
    
    Simpler case than ethanol, should also give opls_154/155.
    """
    mol = Atomistic()
    
    c = mol.add_atom(symbol="C", number=6)   # CH3
    o = mol.add_atom(symbol="O", number=8)   # OH oxygen
    h = mol.add_atom(symbol="H", number=1)   # OH hydrogen
    
    mol.add_bond(c, o, order=1)
    mol.add_bond(o, h, order=1)
    
    # Apply typifier
    oplsaa_typifier.typify(mol)
    
    # Check results
    o_type = o.get("type")
    h_type = h.get("type")
    
    assert o_type == "opls_154", \
        f"Oxygen should be opls_154, got {o_type}"
    assert h_type == "opls_155", \
        f"Hydrogen on alcohol O should be opls_155, got {h_type}"


def test_opls_157_depends_on_154(oplsaa_typifier):
    """Test that opls_157 (alcohol carbon) also depends on opls_154.
    
    opls_157: [C;X4]([H])([H])([*])[O;%opls_154] - CH2 or CH3 connected to alcohol O
    
    This pattern references opls_154 and should be at a higher level.
    """
    pattern_dict = oplsaa_typifier.pattern_dict
    
    if "opls_157" in pattern_dict:
        opls_154_pattern = pattern_dict["opls_154"]
        opls_157_pattern = pattern_dict["opls_157"]
        
        # Check dependency
        assert "opls_154" in opls_157_pattern.dependencies, \
            f"opls_157 should depend on opls_154, got: {opls_157_pattern.dependencies}"
        
        # Check levels
        assert opls_157_pattern.level > opls_154_pattern.level, \
            f"opls_157 (level {opls_157_pattern.level}) should be after opls_154 (level {opls_154_pattern.level})"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
