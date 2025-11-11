"""
Unit tests for generate_3d_coords function in rdkit adapter.
"""

import pytest
from molpy.parser.smiles import SmilesParser
from molpy.core.wrappers.monomer import Monomer
from molpy.adapter.rdkit_adapter import generate_3d_coords
from molpy.adapter.converter import convert
from rdkit import Chem
from molpy.core.atomistic import Atomistic


@pytest.fixture
def ethanol_monomer():
    """Create ethanol monomer without 3D coordinates."""
    parser = SmilesParser()
    smiles_ir = parser.parse_smiles("CCO")
    mol = convert(smiles_ir, Chem.Mol)
    atomistic = convert(mol, Atomistic)
    return Monomer(atomistic)


@pytest.fixture
def benzene_monomer():
    """Create benzene monomer without 3D coordinates."""
    parser = SmilesParser()
    smiles_ir = parser.parse_smiles("c1ccccc1")
    mol = convert(smiles_ir, Chem.Mol)
    atomistic = convert(mol, Atomistic)
    return Monomer(atomistic)


def test_generate_3d_coords_basic(ethanol_monomer):
    """Test basic 3D coordinate generation."""
    # Check initial state - no coordinates
    atoms = list(ethanol_monomer.unwrap().atoms)
    assert len(atoms) == 3
    assert all(atom.data.get("pos") is None for atom in atoms)
    
    # Generate 3D coordinates
    result = generate_3d_coords(ethanol_monomer)
    
    # Should return same object
    assert result is ethanol_monomer
    
    # Should have coordinates now
    atoms = list(ethanol_monomer.unwrap().atoms)
    assert all(atom.data.get("pos") is not None for atom in atoms)
    assert all(atom.data.get("xyz") is not None for atom in atoms)
    
    # Coordinates should be 3D lists
    for atom in atoms:
        pos = atom.data["pos"]
        assert isinstance(pos, list)
        assert len(pos) == 3
        assert all(isinstance(x, float) for x in pos)


def test_generate_3d_coords_without_optimization(ethanol_monomer):
    """Test 3D generation without MMFF optimization."""
    result = generate_3d_coords(ethanol_monomer, optimize=False)
    
    # Should still have coordinates
    atoms = list(result.unwrap().atoms)
    assert all(atom.data.get("pos") is not None for atom in atoms)


def test_generate_3d_coords_benzene(benzene_monomer):
    """Test 3D generation for aromatic molecule."""
    result = generate_3d_coords(benzene_monomer)
    
    # Benzene should have 6 carbons
    atoms = list(result.unwrap().atoms)
    assert len(atoms) == 6
    assert all(atom.data.get("symbol") == "C" for atom in atoms)
    assert all(atom.data.get("pos") is not None for atom in atoms)


def test_generate_3d_coords_with_ports(ethanol_monomer):
    """Test that ports are preserved after 3D generation."""
    # Add a port before 3D generation
    atoms = list(ethanol_monomer.unwrap().atoms)
    ethanol_monomer.set_port("test_port", atoms[0], role="left")
    
    # Generate 3D coordinates
    generate_3d_coords(ethanol_monomer)
    
    # Port should still be there
    assert "test_port" in ethanol_monomer.ports
    port = ethanol_monomer.get_port("test_port")
    assert port.role == "left"
    assert port.target.data.get("pos") is not None


def test_generate_3d_coords_random_seed(ethanol_monomer):
    """Test that random seed makes results reproducible."""
    # Generate with seed 42
    generate_3d_coords(ethanol_monomer, random_seed=42)
    atoms1 = list(ethanol_monomer.unwrap().atoms)
    coords1 = [atom.data["pos"] for atom in atoms1]
    
    # Create fresh monomer
    parser = SmilesParser()
    smiles_ir = parser.parse_smiles("CCO")
    mol = convert(smiles_ir, Chem.Mol)
    atomistic = convert(mol, Atomistic)
    ethanol_monomer2 = Monomer(atomistic)
    
    # Generate with same seed
    generate_3d_coords(ethanol_monomer2, random_seed=42)
    atoms2 = list(ethanol_monomer2.unwrap().atoms)
    coords2 = [atom.data["pos"] for atom in atoms2]
    
    # Should be identical
    for c1, c2 in zip(coords1, coords2):
        assert c1 == pytest.approx(c2, abs=1e-6)


def test_generate_3d_coords_type_error():
    """Test that TypeError is raised for non-Monomer input."""
    with pytest.raises(TypeError, match="Expected Monomer"):
        generate_3d_coords("not a monomer")
