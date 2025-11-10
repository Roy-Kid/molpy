"""Unit tests for molecule adapter.

Tests conversion from Atomistic to igraph.Graph with proper attributes.
"""

import pytest
from molpy.core.atomistic import Atomistic
from molpy.typifier.adapter import (
    build_mol_graph,
    _extract_atom_attributes,
    _extract_bond_attributes,
    _detect_rings,
    _find_cycles_containing_vertex,
)


class TestBuildMolGraph:
    """Test build_mol_graph function."""
    
    def test_single_atom(self):
        """Test converting single atom."""
        mol = Atomistic()
        c = mol.add_atom(element="C", charge=0)
        
        graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)
        
        assert graph.vcount() == 1
        assert graph.ecount() == 0
        assert len(vs_to_atomid) == 1
        assert len(atomid_to_vs) == 1
    
    def test_two_atoms_with_bond(self):
        """Test converting two atoms with bond."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        o = mol.add_atom(element="O")
        mol.add_bond(c, o, order=1)
        
        graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)
        
        assert graph.vcount() == 2
        assert graph.ecount() == 1
    
    def test_bidirectional_mapping(self):
        """Test bidirectional atom ID mapping."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        
        graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)
        
        atom_id = id(c)
        vs_idx = atomid_to_vs[atom_id]
        assert vs_to_atomid[vs_idx] == atom_id
    
    def test_deterministic_ordering(self):
        """Test that ordering is deterministic."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(5)]
        
        graph1, vs_to_id1, _ = build_mol_graph(mol)
        graph2, vs_to_id2, _ = build_mol_graph(mol)
        
        # Should produce same ordering
        assert vs_to_id1 == vs_to_id2


class TestExtractAtomAttributes:
    """Test _extract_atom_attributes function."""
    
    def test_basic_attributes(self):
        """Test extracting basic atom attributes."""
        mol = Atomistic()
        atom = mol.add_atom(element="C", charge=0, is_aromatic=False)
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["element"] == "C"
        assert attrs["charge"] == 0
        assert attrs["is_aromatic"] is False
    
    def test_element_from_symbol(self):
        """Test extracting element from symbol."""
        mol = Atomistic()
        atom = mol.add_atom(element="N")
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["element"] == "N"
        assert attrs["number"] == 7
    
    def test_element_from_number(self):
        """Test extracting element from atomic number."""
        mol = Atomistic()
        atom = mol.add_atom(number=8)  # Oxygen
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["element"] == "O"
        assert attrs["number"] == 8
    
    def test_hybridization(self):
        """Test extracting hybridization."""
        mol = Atomistic()
        atom = mol.add_atom(element="C", hyb=2)
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["hyb"] == 2
    
    def test_hybridization_string(self):
        """Test extracting hybridization from string."""
        mol = Atomistic()
        atom = mol.add_atom(element="C", hyb="sp2")
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["hyb"] == 2
    
    def test_default_values(self):
        """Test default values for missing attributes."""
        mol = Atomistic()
        atom = mol.add_atom(element="C")
        
        attrs = _extract_atom_attributes(atom, mol)
        
        assert attrs["charge"] == 0
        assert attrs["is_aromatic"] is False
        assert attrs["degree"] == 0  # Will be computed later
        assert attrs["in_ring"] is False  # Will be computed later


class TestExtractBondAttributes:
    """Test _extract_bond_attributes function."""
    
    def test_single_bond(self):
        """Test extracting single bond attributes."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        o = mol.add_atom(element="O")
        bond = mol.add_bond(c, o, order=1)
        
        attrs = _extract_bond_attributes(bond, mol)
        
        assert attrs["order"] == 1
        assert attrs["is_aromatic"] is False
    
    def test_double_bond(self):
        """Test extracting double bond attributes."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        o = mol.add_atom(element="O")
        bond = mol.add_bond(c, o, order=2)
        
        attrs = _extract_bond_attributes(bond, mol)
        
        assert attrs["order"] == 2
    
    def test_aromatic_bond(self):
        """Test extracting aromatic bond attributes."""
        mol = Atomistic()
        c1 = mol.add_atom(element="C", is_aromatic=True)
        c2 = mol.add_atom(element="C", is_aromatic=True)
        bond = mol.add_bond(c1, c2, order=":", is_aromatic=True)
        
        attrs = _extract_bond_attributes(bond, mol)
        
        assert attrs["order"] == ":"
        assert attrs["is_aromatic"] is True
    
    def test_bond_order_string(self):
        """Test bond order as string."""
        mol = Atomistic()
        c = mol.add_atom(element="C")
        o = mol.add_atom(element="O")
        bond = mol.add_bond(c, o, order="1")
        
        attrs = _extract_bond_attributes(bond, mol)
        
        assert attrs["order"] == 1


class TestRingDetection:
    """Test ring detection algorithms."""
    
    def test_no_rings(self):
        """Test linear chain has no rings."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(4)]
        for i in range(3):
            mol.add_bond(atoms[i], atoms[i+1], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # No atoms should be in rings
        for v in graph.vs:
            assert v["in_ring"] is False
            assert len(v["cycles"]) == 0
    
    def test_three_membered_ring(self):
        """Test 3-membered ring detection."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(3)]
        mol.add_bond(atoms[0], atoms[1], order=1)
        mol.add_bond(atoms[1], atoms[2], order=1)
        mol.add_bond(atoms[2], atoms[0], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # All atoms should be in ring
        for v in graph.vs:
            assert v["in_ring"] is True
            assert len(v["cycles"]) == 1
    
    def test_six_membered_ring(self):
        """Test 6-membered ring detection (cyclohexane)."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(6)]
        for i in range(6):
            j = (i + 1) % 6
            mol.add_bond(atoms[i], atoms[j], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # All atoms should be in ring
        for v in graph.vs:
            assert v["in_ring"] is True
            assert len(v["cycles"]) == 1
    
    def test_fused_rings(self):
        """Test fused ring system detection."""
        # Create naphthalene-like structure (two fused 6-rings)
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(10)]
        
        # First ring: 0-1-2-3-4-5-0
        bonds_ring1 = [(0,1), (1,2), (2,3), (3,4), (4,5), (5,0)]
        
        # Second ring shares edge 4-5: 4-5-6-7-8-9-4
        bonds_ring2 = [(5,6), (6,7), (7,8), (8,9), (9,4)]
        
        for i, j in bonds_ring1 + bonds_ring2:
            mol.add_bond(atoms[i], atoms[j], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # Shared atoms (4, 5) should be in 2 rings
        assert graph.vs[4]["in_ring"] is True
        assert graph.vs[5]["in_ring"] is True
        # They should have 2 cycles
        assert len(graph.vs[4]["cycles"]) == 2
        assert len(graph.vs[5]["cycles"]) == 2


class TestDegreeComputation:
    """Test degree computation."""
    
    def test_single_bond_degree(self):
        """Test degree computation for single bond."""
        mol = Atomistic()
        c1 = mol.add_atom(element="C")
        c2 = mol.add_atom(element="C")
        mol.add_bond(c1, c2, order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        assert graph.vs[0]["degree"] == 1
        assert graph.vs[1]["degree"] == 1
    
    def test_branched_degree(self):
        """Test degree for branched structure."""
        mol = Atomistic()
        c_center = mol.add_atom(element="C")
        c1 = mol.add_atom(element="C")
        c2 = mol.add_atom(element="C")
        c3 = mol.add_atom(element="C")
        
        mol.add_bond(c_center, c1, order=1)
        mol.add_bond(c_center, c2, order=1)
        mol.add_bond(c_center, c3, order=1)
        
        graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)
        
        # Central carbon should have degree 3
        center_idx = atomid_to_vs[id(c_center)]
        assert graph.vs[center_idx]["degree"] == 3


class TestEdgeInRing:
    """Test edge in ring detection."""
    
    def test_edge_in_ring(self):
        """Test edges in ring are marked."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(6)]
        for i in range(6):
            j = (i + 1) % 6
            mol.add_bond(atoms[i], atoms[j], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # All edges should be in ring
        for e in graph.es:
            assert e["is_in_ring"] is True
    
    def test_edge_not_in_ring(self):
        """Test edges not in ring."""
        mol = Atomistic()
        atoms = [mol.add_atom(element="C") for _ in range(3)]
        mol.add_bond(atoms[0], atoms[1], order=1)
        mol.add_bond(atoms[1], atoms[2], order=1)
        
        graph, _, _ = build_mol_graph(mol)
        
        # No edges should be in ring
        for e in graph.es:
            assert e["is_in_ring"] is False


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
