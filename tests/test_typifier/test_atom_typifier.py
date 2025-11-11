"""Tests for OplsAtomTypifier with SMARTS matcher."""

import pytest
from molpy.core.atomistic import Atomistic
from molpy.core.forcefield import ForceField, AtomType, Style
from molpy.typifier.atomistic import OplsAtomTypifier
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import is_element, bond_order


class CustomOplsAtomTypifier(OplsAtomTypifier):
    """Custom atom typifier that uses pre-defined patterns."""
    
    def _extract_patterns(self):
        """Override to use custom patterns if available."""
        if hasattr(self.ff, '_atom_patterns'):
            return self.ff._atom_patterns
        return super()._extract_patterns()


@pytest.fixture
def ethanol_forcefield():
    """Create a forcefield for ethanol typing."""
    ff = ForceField(name="ethanol_ff")
    
    # Define atom types
    at_ct = AtomType("CT", element="C", priority=0, mass=12.01, charge=0.0)
    at_oh = AtomType("OH", element="O", priority=5, mass=16.00, charge=-0.5)
    at_c_oh = AtomType("C_OH", element="C", priority=2, mass=12.01, charge=0.2)
    
    # Add to forcefield
    atom_style = Style("atom")
    atom_style.types.add(at_ct)
    atom_style.types.add(at_oh)
    atom_style.types.add(at_c_oh)
    ff.styles.add(atom_style)
    
    # Create patterns
    patterns = []
    
    # Generic carbon (lowest priority)
    p1 = quick_pattern("CT", "C", priority=0)
    patterns.append(p1)
    
    # Hydroxyl oxygen (high priority)
    pb2 = PatternBuilder("OH", priority=5, source="hydroxyl_pattern")
    v_c = pb2.add_vertex(preds=[is_element("C")])
    v_o = pb2.add_vertex(preds=[is_element("O")])
    pb2.add_edge(v_c, v_o, preds=[bond_order(1)])
    pb2.set_target_vertices([v_o])
    p2 = pb2.build()
    patterns.append(p2)
    
    # Alcohol carbon (medium priority)
    pb3 = PatternBuilder("C_OH", priority=2, source="alcohol_pattern")
    v_c3 = pb3.add_vertex(preds=[is_element("C")])
    v_o3 = pb3.add_vertex(preds=[is_element("O")])
    pb3.add_edge(v_c3, v_o3, preds=[bond_order(1)])
    pb3.set_target_vertices([v_c3])
    p3 = pb3.build()
    patterns.append(p3)
    
    ff._atom_patterns = patterns  # type: ignore
    
    return ff


@pytest.fixture
def ethanol_molecule():
    """Create ethanol molecule: CH3-CH2-OH"""
    mol = Atomistic()
    
    # Atoms
    c1 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    c2 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    o = mol.add_atom(element="O", is_aromatic=False, charge=0)
    
    # Bonds
    mol.add_bond(c1, c2, order=1)
    mol.add_bond(c2, o, order=1)
    
    return mol, {"C1": c1, "C2": c2, "O": o}


def test_typifier_initialization(ethanol_forcefield):
    """Test that typifier initializes correctly."""
    typifier = CustomOplsAtomTypifier(ethanol_forcefield)
    
    assert typifier.ff == ethanol_forcefield
    assert len(typifier.patterns) == 3
    assert typifier.matcher is not None


def test_ethanol_typing(ethanol_forcefield, ethanol_molecule):
    """Test that ethanol is typed correctly."""
    typifier = CustomOplsAtomTypifier(ethanol_forcefield)
    mol, atoms = ethanol_molecule
    
    # Before typing
    for atom in mol.atoms:
        assert atom.get("type") is None
    
    # Apply typifier
    typifier.typify(mol)
    
    # After typing
    assert atoms["C1"].get("type") == "CT"
    assert atoms["C2"].get("type") == "C_OH"
    assert atoms["O"].get("type") == "OH"


def test_atom_parameters_assigned(ethanol_forcefield, ethanol_molecule):
    """Test that atom parameters are assigned from forcefield."""
    typifier = CustomOplsAtomTypifier(ethanol_forcefield)
    mol, atoms = ethanol_molecule
    
    typifier.typify(mol)
    
    # Check that parameters are assigned
    assert atoms["C1"].get("mass") == 12.01
    assert atoms["C1"].get("charge") == 0.0
    
    assert atoms["C2"].get("mass") == 12.01
    assert atoms["C2"].get("charge") == 0.2
    
    assert atoms["O"].get("mass") == 16.00
    assert atoms["O"].get("charge") == -0.5


def test_priority_resolution(ethanol_forcefield):
    """Test that higher priority patterns win."""
    typifier = CustomOplsAtomTypifier(ethanol_forcefield)
    
    # Create a molecule with C-O bond
    mol = Atomistic()
    c = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    o = mol.add_atom(element="O", is_aromatic=False, charge=0)
    mol.add_bond(c, o, order=1)
    
    typifier.typify(mol)
    
    # C should be C_OH (priority 2) not CT (priority 0)
    assert c.get("type") == "C_OH"
    
    # O should be OH (priority 5)
    assert o.get("type") == "OH"


def test_simple_molecule():
    """Test typing a simple methane molecule."""
    ff = ForceField(name="simple_ff")
    
    # Just one atom type: generic carbon
    at_c = AtomType("C", element="C", priority=0, mass=12.01)
    atom_style = Style("atom")
    atom_style.types.add(at_c)
    ff.styles.add(atom_style)
    
    # Simple pattern
    pattern = quick_pattern("C", "C", priority=0)
    ff._atom_patterns = [pattern]  # type: ignore
    
    # Create methane (just C atom)
    mol = Atomistic()
    c = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    
    # Apply typifier
    typifier = CustomOplsAtomTypifier(ff)
    typifier.typify(mol)
    
    assert c.get("type") == "C"
    assert c.get("mass") == 12.01


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
