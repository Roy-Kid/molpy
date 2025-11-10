"""Unit tests for PatternBuilder class.

Tests the fluent API for building SMARTS patterns programmatically.
"""

import pytest
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import (
    is_element, is_aromatic, charge, degree, bond_order
)


class TestPatternBuilder:
    """Test PatternBuilder class."""
    
    def test_init(self):
        """Test PatternBuilder initialization."""
        pb = PatternBuilder("test_type", priority=5, source="test_source")
        assert pb.atomtype_name == "test_type"
        assert pb.priority == 5
        assert pb.source == "test_source"
    
    def test_add_single_vertex(self):
        """Test adding a single vertex."""
        pb = PatternBuilder("t_C")
        v = pb.add_vertex(preds=[is_element("C")])
        
        assert v == 0  # First vertex has index 0
        assert len(pb._vertices) == 1
    
    def test_add_multiple_vertices(self):
        """Test adding multiple vertices."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        
        assert v1 == 0
        assert v2 == 1
        assert len(pb._vertices) == 2
    
    def test_add_edge(self):
        """Test adding an edge."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(v1, v2, preds=[bond_order(1)])
        
        assert len(pb._edges) == 1
    
    def test_add_edge_invalid_vertex(self):
        """Test adding edge with invalid vertex raises error."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        
        with pytest.raises(ValueError, match="Vertex .* not in pattern"):
            pb.add_edge(v1, 999, preds=[])
    
    def test_set_target_vertices(self):
        """Test setting target vertices."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(v1, v2)
        
        pb.set_target_vertices([v2])
        assert pb._target_vertices == [v2]
    
    def test_set_target_vertices_invalid(self):
        """Test setting invalid target vertices raises error."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        
        with pytest.raises(ValueError, match="Vertex .* not in pattern"):
            pb.set_target_vertices([999])
    
    def test_build_single_vertex(self):
        """Test building single-vertex pattern."""
        pb = PatternBuilder("t_C", priority=0)
        v = pb.add_vertex(preds=[is_element("C")])
        pattern = pb.build()
        
        assert pattern.atomtype_name == "t_C"
        assert pattern.vcount() == 1
        assert pattern.ecount() == 0
    
    def test_build_two_vertices_with_edge(self):
        """Test building two-vertex pattern with edge."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(v1, v2, preds=[bond_order(1)])
        pattern = pb.build()
        
        assert pattern.vcount() == 2
        assert pattern.ecount() == 1
    
    def test_build_preserves_predicates(self):
        """Test that build preserves predicates."""
        pb = PatternBuilder("t_C")
        v = pb.add_vertex(preds=[is_element("C"), charge(-1)])
        pattern = pb.build()
        
        # Check predicates are attached
        preds = pattern.vs[0]["preds"]
        assert len(preds) == 2
    
    def test_build_with_target_vertices(self):
        """Test building with target vertices."""
        pb = PatternBuilder("t_CO")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(v1, v2)
        pb.set_target_vertices([v2])
        
        pattern = pb.build()
        assert pattern.target_vertices == [1]  # v2 maps to index 1


class TestQuickPattern:
    """Test quick_pattern helper function."""
    
    def test_simple_element(self):
        """Test simple element pattern."""
        pattern = quick_pattern("t_C", "C", priority=0)
        
        assert pattern.atomtype_name == "t_C"
        assert pattern.vcount() == 1
        assert pattern.ecount() == 0
    
    def test_with_aromatic(self):
        """Test pattern with aromatic constraint."""
        pattern = quick_pattern("t_ar_C", "C", priority=0, is_aromatic=True)
        
        preds = pattern.vs[0]["preds"]
        assert len(preds) == 2  # element + aromatic
    
    def test_with_charge(self):
        """Test pattern with charge constraint."""
        pattern = quick_pattern("t_O_neg", "O", priority=0, charge=-1)
        
        preds = pattern.vs[0]["preds"]
        assert len(preds) == 2  # element + charge
    
    def test_with_multiple_constraints(self):
        """Test pattern with multiple constraints."""
        pattern = quick_pattern(
            "t_C_complex", "C", priority=5,
            is_aromatic=False, charge=0, degree=4, hyb=3
        )
        
        preds = pattern.vs[0]["preds"]
        assert len(preds) == 5  # element + 4 constraints


class TestPatternBuilderComplexPatterns:
    """Test building complex patterns."""
    
    def test_linear_chain(self):
        """Test building linear chain C-C-C."""
        pb = PatternBuilder("t_chain")
        c1 = pb.add_vertex(preds=[is_element("C")])
        c2 = pb.add_vertex(preds=[is_element("C")])
        c3 = pb.add_vertex(preds=[is_element("C")])
        pb.add_edge(c1, c2, preds=[bond_order(1)])
        pb.add_edge(c2, c3, preds=[bond_order(1)])
        
        pattern = pb.build()
        assert pattern.vcount() == 3
        assert pattern.ecount() == 2
    
    def test_branched_pattern(self):
        """Test building branched pattern."""
        pb = PatternBuilder("t_branch")
        c_center = pb.add_vertex(preds=[is_element("C")])
        c1 = pb.add_vertex(preds=[is_element("C")])
        c2 = pb.add_vertex(preds=[is_element("C")])
        c3 = pb.add_vertex(preds=[is_element("C")])
        
        pb.add_edge(c_center, c1, preds=[bond_order(1)])
        pb.add_edge(c_center, c2, preds=[bond_order(1)])
        pb.add_edge(c_center, c3, preds=[bond_order(1)])
        
        pattern = pb.build()
        assert pattern.vcount() == 4
        assert pattern.ecount() == 3
    
    def test_carbonyl_pattern(self):
        """Test building carbonyl C=O pattern."""
        pb = PatternBuilder("t_carbonyl")
        c = pb.add_vertex(preds=[is_element("C")])
        o = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(c, o, preds=[bond_order(2)])
        
        pattern = pb.build()
        assert pattern.vcount() == 2
        assert pattern.ecount() == 1
        
        # Check edge predicates
        edge_preds = pattern.es[0]["preds"]
        assert len(edge_preds) == 1


class TestPatternBuilderEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_pattern(self):
        """Test building pattern with no vertices."""
        pb = PatternBuilder("t_empty")
        pattern = pb.build()
        
        assert pattern.vcount() == 0
        assert pattern.ecount() == 0
    
    def test_single_vertex_no_predicates(self):
        """Test vertex with no predicates (defaults to wildcard)."""
        pb = PatternBuilder("t_any")
        v = pb.add_vertex()  # No predicates
        pattern = pb.build()
        
        preds = pattern.vs[0]["preds"]
        assert len(preds) == 1  # Default wildcard
    
    def test_edge_no_predicates(self):
        """Test edge with no predicates."""
        pb = PatternBuilder("t_any_bond")
        v1 = pb.add_vertex(preds=[is_element("C")])
        v2 = pb.add_vertex(preds=[is_element("O")])
        pb.add_edge(v1, v2)  # No predicates
        
        pattern = pb.build()
        edge_preds = pattern.es[0]["preds"]
        assert len(edge_preds) == 0


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
