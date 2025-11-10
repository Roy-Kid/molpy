"""
Tests for __repr__ and tree_repr() methods.

This module tests the string representations of:
1. Atomistic structures
2. Wrapper objects
3. Multi-layer wrapper stacks
4. Monomer with ports
"""

import pytest
from molpy.core.atomistic import Atom, Bond, Atomistic
from molpy.core.wrappers.base import Wrapper
from molpy.core.wrappers.monomer import Monomer, Port
from molpy.core.entity import Assembly


# ============================================================================
# Fixtures
# ============================================================================


@pytest.fixture
def simple_atomistic():
    """Create simple atomistic with 3 atoms."""
    asm = Atomistic()
    a1 = Atom(symbol="C", pos=[0, 0, 0])
    a2 = Atom(symbol="C", pos=[1, 0, 0])
    a3 = Atom(symbol="O", pos=[2, 0, 0])
    b1 = Bond(a1, a2, order=1)
    b2 = Bond(a2, a3, order=1)
    asm.add_entity(a1, a2, a3)
    asm.add_link(b1, b2)
    return asm


@pytest.fixture
def large_atomistic():
    """Create larger atomistic with multiple element types."""
    asm = Atomistic()
    # Add 10 carbons, 5 oxygens, 3 nitrogens
    for i in range(10):
        asm.add_atom(symbol="C")
    for i in range(5):
        asm.add_atom(symbol="O")
    for i in range(3):
        asm.add_atom(symbol="N")
    return asm


@pytest.fixture
def monomer_with_ports(simple_atomistic):
    """Create monomer with ports."""
    mon = Monomer(simple_atomistic)
    atoms = simple_atomistic.atoms
    mon.set_port("in", atoms[0])
    mon.set_port("out", atoms[2])
    return mon


# Custom wrappers for testing
class TaggedWrapper(Wrapper):
    """Wrapper that adds tags."""
    
    def __post_init__(self, **props):
        self._tags = props.pop('tags', [])
        return props
    
    @property
    def tags(self):
        return self._tags
    
    def _repr_info(self) -> str:
        if self._tags:
            return f"tags={self._tags}"
        return "no tags"


class PropertyWrapper(Wrapper):
    """Wrapper that adds properties."""
    
    def __post_init__(self, **props):
        self._properties = props.pop('properties', {})
        return props
    
    def _repr_info(self) -> str:
        if self._properties:
            items = ", ".join(f"{k}={v}" for k, v in self._properties.items())
            return f"props: {items}"
        return "no props"


# ============================================================================
# Test Atom and Bond repr
# ============================================================================


class TestAtomBondRepr:
    """Test Atom and Bond __repr__ methods."""
    
    def test_atom_repr_with_symbol(self):
        """Test atom repr shows symbol."""
        atom = Atom(symbol="C")
        repr_str = repr(atom)
        assert "Atom" in repr_str
        assert "C" in repr_str
    
    def test_atom_repr_with_type(self):
        """Test atom repr falls back to type."""
        atom = Atom(type="CT")
        repr_str = repr(atom)
        assert "Atom" in repr_str
        assert "CT" in repr_str
    
    def test_atom_repr_without_symbol(self):
        """Test atom repr without symbol uses id."""
        atom = Atom()
        repr_str = repr(atom)
        assert "Atom" in repr_str
    
    def test_bond_repr(self):
        """Test bond repr shows connected atoms."""
        a1 = Atom(symbol="C")
        a2 = Atom(symbol="O")
        bond = Bond(a1, a2, order=2)
        repr_str = repr(bond)
        assert "Bond" in repr_str
        assert "C" in repr_str
        assert "O" in repr_str


# ============================================================================
# Test Atomistic repr
# ============================================================================


class TestAtomisticRepr:
    """Test Atomistic __repr__ method."""
    
    def test_simple_atomistic_repr(self, simple_atomistic):
        """Test repr of simple atomistic."""
        repr_str = repr(simple_atomistic)
        
        print(f"\nSimple atomistic repr:\n{repr_str}")
        
        assert "Atomistic" in repr_str
        assert "3 atoms" in repr_str
        assert "2 bonds" in repr_str
        # Should show composition
        assert "C:2" in repr_str or "C" in repr_str
        assert "O:1" in repr_str or "O" in repr_str
        assert "with coords" in repr_str
    
    def test_large_atomistic_repr(self, large_atomistic):
        """Test repr of larger atomistic."""
        repr_str = repr(large_atomistic)
        
        print(f"\nLarge atomistic repr:\n{repr_str}")
        
        assert "Atomistic" in repr_str
        assert "18 atoms" in repr_str
        # Should show element counts
        assert "C:10" in repr_str
        assert "O:5" in repr_str
        assert "N:3" in repr_str
    
    def test_atomistic_without_coords(self):
        """Test repr without coordinates."""
        asm = Atomistic()
        asm.add_atom(symbol="C")
        asm.add_atom(symbol="C")
        
        repr_str = repr(asm)
        
        print(f"\nAtomistic without coords:\n{repr_str}")
        
        assert "Atomistic" in repr_str
        assert "2 atoms" in repr_str
        assert "with coords" not in repr_str
    
    def test_empty_atomistic_repr(self):
        """Test repr of empty atomistic."""
        asm = Atomistic()
        repr_str = repr(asm)
        
        print(f"\nEmpty atomistic repr:\n{repr_str}")
        
        assert "Atomistic" in repr_str
        assert "0 atoms" in repr_str


# ============================================================================
# Test Port repr
# ============================================================================


class TestPortRepr:
    """Test Port __repr__ method."""
    
    def test_port_repr(self):
        """Test port repr shows name and target."""
        atom = Atom(symbol="C")
        port = Port("in", atom)
        
        repr_str = repr(port)
        
        print(f"\nPort repr:\n{repr_str}")
        
        assert "Port" in repr_str
        assert "in" in repr_str
        assert "Atom" in repr_str


# ============================================================================
# Test Wrapper repr
# ============================================================================


class TestWrapperRepr:
    """Test Wrapper __repr__ method."""
    
    def test_simple_wrapper_repr(self, simple_atomistic):
        """Test basic wrapper repr."""
        wrapper = Wrapper(simple_atomistic)
        repr_str = repr(wrapper)
        
        print(f"\nSimple wrapper repr:\n{repr_str}")
        
        assert "Wrapper" in repr_str
        assert "Atomistic" in repr_str
    
    def test_monomer_repr(self, monomer_with_ports):
        """Test monomer repr."""
        repr_str = repr(monomer_with_ports)
        
        print(f"\nMonomer repr:\n{repr_str}")
        
        assert "Monomer" in repr_str
        assert "Atomistic" in repr_str
    
    def test_nested_wrapper_repr(self, simple_atomistic):
        """Test nested wrapper repr."""
        inner = TaggedWrapper(simple_atomistic, tags=["polymer"])
        outer = PropertyWrapper(inner, properties={"mass": 100})
        
        repr_str = repr(outer)
        
        print(f"\nNested wrapper repr:\n{repr_str}")
        
        # Should show outermost type and innermost wrapped type
        assert "PropertyWrapper" in repr_str
        assert "TaggedWrapper" in repr_str


# ============================================================================
# Test tree_repr
# ============================================================================


class TestTreeRepr:
    """Test Wrapper.tree_repr() method."""
    
    def test_simple_wrapper_tree(self, simple_atomistic):
        """Test tree repr of simple wrapper."""
        wrapper = Wrapper(simple_atomistic)
        tree = wrapper.tree_repr()
        
        print(f"\nSimple wrapper tree:\n{tree}")
        
        assert "Wrapper" in tree
        assert "Atomistic" in tree
        assert "└──" in tree  # Tree connector
    
    def test_monomer_tree(self, monomer_with_ports):
        """Test tree repr of monomer."""
        tree = monomer_with_ports.tree_repr()
        
        print(f"\nMonomer tree:\n{tree}")
        
        assert "Monomer" in tree
        assert "port" in tree
        assert "in" in tree
        assert "out" in tree
        assert "Atomistic" in tree
    
    def test_two_layer_tree(self, simple_atomistic):
        """Test tree repr with two wrapper layers."""
        inner = TaggedWrapper(simple_atomistic, tags=["polymer", "reactive"])
        outer = Monomer(inner)  # type: ignore
        outer.set_port("test", simple_atomistic.atoms[0])
        
        tree = outer.tree_repr()
        
        print(f"\nTwo-layer tree:\n{tree}")
        
        # Check structure
        lines = tree.split('\n')
        assert len(lines) >= 3  # At least 3 layers
        
        # Check content
        assert "Monomer" in tree
        assert "TaggedWrapper" in tree
        assert "Atomistic" in tree
        
        # Check tree structure
        assert "└──" in tree
    
    def test_three_layer_tree(self, simple_atomistic):
        """Test tree repr with three wrapper layers."""
        layer1 = TaggedWrapper(simple_atomistic, tags=["A"])
        layer2 = PropertyWrapper(layer1, properties={"x": 1, "y": 2})
        layer3 = Monomer(layer2)  # type: ignore
        layer3.set_port("in", simple_atomistic.atoms[0])
        layer3.set_port("out", simple_atomistic.atoms[1])
        
        tree = layer3.tree_repr()
        
        print(f"\nThree-layer tree:\n{tree}")
        
        # Check all layers present
        assert "Monomer" in tree
        assert "PropertyWrapper" in tree
        assert "TaggedWrapper" in tree
        assert "Atomistic" in tree
        
        # Check specific info
        assert "port" in tree
        assert "in" in tree
        assert "out" in tree
        assert "props" in tree or "x=1" in tree
        assert "tags" in tree or "A" in tree
    
    def test_deep_stack_tree(self, simple_atomistic):
        """Test tree repr with deep wrapper stack."""
        # Create 5-layer stack
        w1 = TaggedWrapper(simple_atomistic, tags=["layer1"])
        w2 = PropertyWrapper(w1, properties={"level": 2})
        w3 = TaggedWrapper(w2, tags=["layer3"])
        w4 = PropertyWrapper(w3, properties={"level": 4})
        w5 = Monomer(w4)  # type: ignore
        w5.set_port("top", simple_atomistic.atoms[0])
        
        tree = w5.tree_repr()
        
        print(f"\nDeep stack tree:\n{tree}")
        
        # Count lines (should have 6: 5 wrappers + 1 atomistic)
        lines = tree.split('\n')
        assert len(lines) == 6
        
        # Verify all layers
        assert "Monomer" in tree
        assert tree.count("PropertyWrapper") == 2
        assert tree.count("TaggedWrapper") == 2
        assert "Atomistic" in tree


# ============================================================================
# Test repr with real molecular structures
# ============================================================================


class TestReprWithMolecules:
    """Test repr methods with realistic molecular structures."""
    
    def test_ethane_repr(self):
        """Test repr of ethane molecule."""
        # Build C2H6
        asm = Atomistic()
        c1 = asm.add_atom(symbol="C", pos=[0, 0, 0])
        c2 = asm.add_atom(symbol="C", pos=[1.5, 0, 0])
        
        # Add hydrogens
        for i in range(3):
            h = asm.add_atom(symbol="H")
        for i in range(3):
            h = asm.add_atom(symbol="H")
        
        # Add C-C bond
        asm.add_bond(c1, c2)
        
        repr_str = repr(asm)
        
        print(f"\nEthane repr:\n{repr_str}")
        
        assert "8 atoms" in repr_str
        assert "C:2" in repr_str
        assert "H:6" in repr_str
        assert "1 bond" in repr_str
    
    def test_monomer_molecule_tree(self):
        """Test tree repr of monomer wrapping molecule."""
        # Build simple monomer unit
        asm = Atomistic()
        c1 = asm.add_atom(symbol="C", pos=[0, 0, 0])
        c2 = asm.add_atom(symbol="C", pos=[1, 0, 0])
        o1 = asm.add_atom(symbol="O", pos=[2, 0, 0])
        asm.add_bond(c1, c2)
        asm.add_bond(c2, o1)
        
        mon = Monomer(asm)
        mon.set_port("left", c1)
        mon.set_port("right", o1)
        
        tree = mon.tree_repr()
        
        print(f"\nMonomer molecule tree:\n{tree}")
        
        assert "Monomer" in tree
        assert "2 port" in tree
        assert "left" in tree
        assert "right" in tree
        assert "3 atoms" in tree
        assert "C:2" in tree
        assert "O:1" in tree


# ============================================================================
# Test repr customization
# ============================================================================


class TestReprCustomization:
    """Test custom _repr_info() methods."""
    
    def test_custom_wrapper_info(self, simple_atomistic):
        """Test custom wrapper provides custom info."""
        wrapper = TaggedWrapper(simple_atomistic, tags=["test", "demo"])
        tree = wrapper.tree_repr()
        
        print(f"\nCustom wrapper tree:\n{tree}")
        
        assert "tags" in tree.lower()
        assert "test" in tree or "demo" in tree
    
    def test_property_wrapper_info(self, simple_atomistic):
        """Test property wrapper shows properties."""
        wrapper = PropertyWrapper(
            simple_atomistic, 
            properties={"temperature": 300, "pressure": 1}
        )
        tree = wrapper.tree_repr()
        
        print(f"\nProperty wrapper tree:\n{tree}")
        
        assert "prop" in tree.lower()
        assert "temperature" in tree.lower() or "300" in tree
