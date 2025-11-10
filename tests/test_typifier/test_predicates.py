"""Unit tests for predicates module.

Tests all predicate factories and their matching logic.
"""

import pytest
from molpy.typifier.predicates import (
    # Vertex predicates
    is_element, is_aromatic, charge, degree, hyb, in_ring,
    ring_size, ring_count, custom_vertex, wildcard,
    ElementPredicate, AromaticPredicate, ChargePredicate,
    DegreePredicate, HybridizationPredicate, InRingPredicate,
    RingSizePredicate, RingCountPredicate, CustomPredicate,
    # Edge predicates
    bond_order, bond_aromatic, bond_in_ring,
    BondOrderPredicate, BondAromaticPredicate, BondInRingPredicate,
)


# ===================================================================
#                     Test Vertex Predicates
# ===================================================================

class TestElementPredicate:
    """Test ElementPredicate class."""
    
    def test_element_match(self):
        """Test element matching."""
        pred = is_element("C")
        assert pred({"element": "C"}) is True
        assert pred({"element": "N"}) is False
    
    def test_case_insensitive(self):
        """Test case insensitivity."""
        pred = is_element("c")
        assert pred({"element": "C"}) is True
        assert pred({"element": "c"}) is True
    
    def test_wildcard(self):
        """Test wildcard matching."""
        pred = is_element("*")
        assert pred({"element": "C"}) is True
        assert pred({"element": "N"}) is True
        assert pred({"element": "O"}) is True
    
    def test_weight(self):
        """Test predicate weight."""
        pred = is_element("C")
        assert pred.meta.weight == 0  # Baseline


class TestAromaticPredicate:
    """Test AromaticPredicate class."""
    
    def test_aromatic_true(self):
        """Test aromatic matching."""
        pred = is_aromatic(True)
        assert pred({"is_aromatic": True}) is True
        assert pred({"is_aromatic": False}) is False
    
    def test_aromatic_false(self):
        """Test non-aromatic matching."""
        pred = is_aromatic(False)
        assert pred({"is_aromatic": False}) is True
        assert pred({"is_aromatic": True}) is False
    
    def test_default_value(self):
        """Test default aromatic=True."""
        pred = is_aromatic()
        assert pred({"is_aromatic": True}) is True
    
    def test_weight(self):
        """Test predicate weight."""
        pred = is_aromatic()
        assert pred.meta.weight == 2


class TestChargePredicate:
    """Test ChargePredicate class."""
    
    def test_positive_charge(self):
        """Test positive charge matching."""
        pred = charge(+1)
        assert pred({"charge": +1}) is True
        assert pred({"charge": 0}) is False
    
    def test_negative_charge(self):
        """Test negative charge matching."""
        pred = charge(-1)
        assert pred({"charge": -1}) is True
        assert pred({"charge": 0}) is False
    
    def test_neutral(self):
        """Test neutral charge."""
        pred = charge(0)
        assert pred({"charge": 0}) is True
        assert pred({}) is True  # Missing charge defaults to 0


class TestDegreePredicate:
    """Test DegreePredicate class."""
    
    def test_degree_match(self):
        """Test degree matching."""
        pred = degree(3)
        assert pred({"degree": 3}) is True
        assert pred({"degree": 4}) is False
    
    def test_terminal(self):
        """Test terminal atom (degree=1)."""
        pred = degree(1)
        assert pred({"degree": 1}) is True


class TestHybridizationPredicate:
    """Test HybridizationPredicate class."""
    
    def test_sp(self):
        """Test sp hybridization."""
        pred = hyb(1)
        assert pred({"hyb": 1}) is True
        assert pred({"hyb": 2}) is False
    
    def test_sp2(self):
        """Test sp2 hybridization."""
        pred = hyb(2)
        assert pred({"hyb": 2}) is True
    
    def test_sp3(self):
        """Test sp3 hybridization."""
        pred = hyb(3)
        assert pred({"hyb": 3}) is True


class TestRingPredicates:
    """Test ring-related predicates."""
    
    def test_in_ring(self):
        """Test in_ring predicate."""
        pred = in_ring(True)
        assert pred({"in_ring": True}) is True
        assert pred({"in_ring": False}) is False
    
    def test_not_in_ring(self):
        """Test not in ring."""
        pred = in_ring(False)
        assert pred({"in_ring": False}) is True
    
    def test_ring_size(self):
        """Test ring size matching."""
        pred = ring_size(6)
        
        # Has 6-membered ring
        assert pred({"cycles": {(0, 1, 2, 3, 4, 5)}}) is True
        
        # Has 5-membered ring only
        assert pred({"cycles": {(0, 1, 2, 3, 4)}}) is False
    
    def test_ring_count(self):
        """Test ring count matching."""
        pred = ring_count(2)
        
        # In 2 rings
        assert pred({"cycles": {(0, 1, 2), (0, 3, 4)}}) is True
        
        # In 1 ring
        assert pred({"cycles": {(0, 1, 2)}}) is False


class TestCustomPredicate:
    """Test CustomPredicate class."""
    
    def test_custom_function(self):
        """Test custom predicate function."""
        # Heavy atom (not hydrogen)
        def is_heavy(attrs):
            return attrs.get("element", "H") != "H"
        
        pred = custom_vertex(is_heavy, name="is_heavy", weight=2)
        
        assert pred({"element": "C"}) is True
        assert pred({"element": "H"}) is False
        assert pred.meta.name == "is_heavy"
        assert pred.meta.weight == 2


# ===================================================================
#                     Test Edge Predicates
# ===================================================================

class TestBondOrderPredicate:
    """Test BondOrderPredicate class."""
    
    def test_single_bond(self):
        """Test single bond."""
        pred = bond_order(1)
        assert pred({"order": 1}) is True
        assert pred({"order": 2}) is False
    
    def test_double_bond(self):
        """Test double bond."""
        pred = bond_order(2)
        assert pred({"order": 2}) is True
    
    def test_triple_bond(self):
        """Test triple bond."""
        pred = bond_order(3)
        assert pred({"order": 3}) is True
    
    def test_aromatic_bond(self):
        """Test aromatic bond."""
        pred = bond_order(":")
        assert pred({"order": ":"}) is True
        assert pred({"order": 1}) is False
    
    def test_mixed_types(self):
        """Test int vs string order."""
        pred = bond_order(1)
        assert pred({"order": "1"}) is True  # String "1" matches int 1


class TestBondAromaticPredicate:
    """Test BondAromaticPredicate class."""
    
    def test_aromatic_bond(self):
        """Test aromatic bond matching."""
        pred = bond_aromatic(True)
        assert pred({"is_aromatic": True}) is True
        assert pred({"is_aromatic": False}) is False


class TestBondInRingPredicate:
    """Test BondInRingPredicate class."""
    
    def test_bond_in_ring(self):
        """Test bond in ring matching."""
        pred = bond_in_ring(True)
        assert pred({"is_in_ring": True}) is True
        assert pred({"is_in_ring": False}) is False


# ===================================================================
#                     Integration Tests
# ===================================================================

class TestPredicateCombinations:
    """Test combining multiple predicates."""
    
    def test_aromatic_carbon(self):
        """Test aromatic carbon (multiple predicates)."""
        elem_pred = is_element("C")
        arom_pred = is_aromatic(True)
        
        attrs = {"element": "C", "is_aromatic": True}
        assert elem_pred(attrs) is True
        assert arom_pred(attrs) is True
    
    def test_charged_oxygen(self):
        """Test charged oxygen."""
        elem_pred = is_element("O")
        charge_pred = charge(-1)
        
        attrs = {"element": "O", "charge": -1}
        assert elem_pred(attrs) is True
        assert charge_pred(attrs) is True
    
    def test_sp3_carbon_in_ring(self):
        """Test sp3 carbon in ring."""
        elem_pred = is_element("C")
        hyb_pred = hyb(3)
        ring_pred = in_ring(True)
        
        attrs = {"element": "C", "hyb": 3, "in_ring": True}
        assert elem_pred(attrs) is True
        assert hyb_pred(attrs) is True
        assert ring_pred(attrs) is True


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
