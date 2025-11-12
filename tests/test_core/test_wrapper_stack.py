"""
Tests for Wrapper design and multi-layer wrapper stacks.

This module tests:
1. Basic Wrapper behavior (delegation, attribute access)
2. Multi-layer wrapper stacks
3. Monomer as a Wrapper with Port management
4. Copy semantics through wrapper layers
"""

import pytest
from molpy.core.wrappers.base import Wrapper
from molpy.core.wrappers.monomer import Monomer, Port
from molpy.core.entity import Entity, Assembly
from molpy.core.atomistic import Atom, Bond, Atomistic


# ============================================================================
# Test Fixtures
# ============================================================================


@pytest.fixture
def simple_assembly():
    """Create a simple assembly with 2 atoms and 1 bond."""
    asm = Atomistic()
    a1 = Atom(symbol="C", pos=[0, 0, 0])
    a2 = Atom(symbol="C", pos=[1, 0, 0])
    b = Bond(a1, a2, order=1)
    # add_link with include_endpoints=True (default) will add the atoms
    asm.add_link(b)
    return asm


@pytest.fixture
def simple_monomer(simple_assembly):
    """Create a monomer wrapping simple_assembly."""
    mon = Monomer(simple_assembly)
    atoms = simple_assembly.atoms
    mon.set_port("in", atoms[0])
    mon.set_port("out", atoms[1])
    return mon


# ============================================================================
# Custom Wrapper for Testing
# ============================================================================


class TaggedWrapper(Wrapper[Assembly]):
    """Wrapper that adds tags to entities."""
    
    def __post_init__(self, **props):
        """Extract tags from props."""
        self._tags = props.pop('tags', [])
        return props
    
    @property
    def tags(self) -> list[str]:
        return self._tags
    
    def add_tag(self, tag: str) -> None:
        if tag not in self._tags:
            self._tags.append(tag)


class PropertyWrapper(Wrapper[Assembly]):
    """Wrapper that adds arbitrary properties."""
    
    def __post_init__(self, **props):
        """Extract properties."""
        self._properties = props.pop('properties', {})
        return props
    
    def set_property(self, key: str, value) -> None:
        self._properties[key] = value
    
    def get_property(self, key: str, default=None):
        return self._properties.get(key, default)


# ============================================================================
# Test Basic Wrapper Behavior
# ============================================================================


class TestBasicWrapper:
    """Test fundamental wrapper behavior."""
    
    def test_wrapper_delegates_to_wrapped(self, simple_assembly):
        """Test that wrapper delegates attribute access to wrapped object."""
        wrapper = Wrapper(simple_assembly)
        
        # Should access wrapped assembly's attributes
        assert wrapper.entities is simple_assembly.entities
        assert wrapper.links is simple_assembly.links
    
    def test_wrapper_unwrap(self, simple_assembly):
        """Test unwrap() returns wrapped object."""
        wrapper = Wrapper(simple_assembly)
        assert wrapper.unwrap() is simple_assembly
    
    def test_wrapper_get_stack(self, simple_assembly):
        """Test get_stack() returns wrapper chain."""
        wrapper = Wrapper(simple_assembly)
        stack = wrapper.get_stack()
        
        assert len(stack) == 2
        assert stack[0] is wrapper
        assert stack[1] is simple_assembly
    
    def test_wrapper_types(self, simple_assembly):
        """Test wrapper_types() returns class names."""
        wrapper = Wrapper(simple_assembly)
        types = wrapper.wrapper_types()
        
        assert types == ["Wrapper", "Atomistic"]
    
    def test_wrapper_depth(self, simple_assembly):
        """Test wrapper_depth() counts wrapper layers."""
        wrapper = Wrapper(simple_assembly)
        assert wrapper.wrapper_depth() == 1


# ============================================================================
# Test Multi-Layer Wrapper Stack
# ============================================================================


class TestMultiLayerStack:
    """Test multi-layer wrapper stacks."""
    
    def test_two_layer_stack(self, simple_assembly):
        """Test wrapping a wrapper."""
        inner = TaggedWrapper(simple_assembly, tags=["polymer"])
        outer = PropertyWrapper(inner, properties={"mass": 100})
        
        # Check stack
        stack = outer.get_stack()
        assert len(stack) == 3
        assert stack[0] is outer
        assert stack[1] is inner
        assert stack[2] is simple_assembly
        
        # Check types
        assert outer.wrapper_types() == ["PropertyWrapper", "TaggedWrapper", "Atomistic"]
        assert outer.wrapper_depth() == 2
    
    def test_three_layer_stack(self, simple_assembly):
        """Test three-layer wrapper stack."""
        layer1 = TaggedWrapper(simple_assembly, tags=["A"])
        layer2 = PropertyWrapper(layer1, properties={"x": 1})
        layer3 = Wrapper(layer2)
        
        stack = layer3.get_stack()
        assert len(stack) == 4
        assert layer3.wrapper_depth() == 3
    
    def test_delegation_through_stack(self, simple_assembly):
        """Test attribute access delegates through entire stack."""
        inner = TaggedWrapper(simple_assembly, tags=["test"])
        outer = PropertyWrapper(inner, properties={"mass": 100})
        
        # Should access inner wrapper's method
        assert outer.tags == ["test"]
        
        # Should access outer wrapper's method
        assert outer.get_property("mass") == 100
        
        # Should access innermost assembly
        assert len(outer.atoms) == 2
    
    def test_method_calls_through_stack(self, simple_assembly):
        """Test method calls work through wrapper stack."""
        inner = TaggedWrapper(simple_assembly, tags=[])
        outer = PropertyWrapper(inner, properties={})
        
        # Modify through outer wrapper
        outer.add_tag("polymer")
        outer.set_property("temp", 300)
        
        # Check state
        assert inner.tags == ["polymer"]
        assert outer.get_property("temp") == 300


# ============================================================================
# Test Port Class
# ============================================================================


class TestPort:
    """Test Port entity behavior."""
    
    def test_port_creation(self):
        """Test creating a port."""
        atom = Atom(symbol="C")
        port = Port("in", atom)
        
        assert port.name == "in"
        assert port.target is atom
    
    def test_port_is_entity(self):
        """Test that Port is an Entity."""
        atom = Atom(symbol="O")
        port = Port("out", atom)
        
        assert isinstance(port, Entity)
        assert "name" in port
        assert "target" in port
    
    def test_port_dict_access(self):
        """Test port behaves like a dict."""
        atom = Atom(symbol="N")
        port = Port("branch", atom)
        
        assert port["name"] == "branch"
        assert port["target"] is atom
        
        # Can add custom properties
        port["label"] = "reactive"
        assert port["label"] == "reactive"


# ============================================================================
# Test Monomer as Wrapper
# ============================================================================


class TestMonomerWrapper:
    """Test Monomer wrapper behavior."""
    
    def test_monomer_wraps_assembly(self, simple_assembly):
        """Test monomer wraps assembly correctly."""
        mon = Monomer(simple_assembly)
        
        assert mon.unwrap() is simple_assembly
        assert isinstance(mon, Wrapper)
    
    def test_monomer_delegates_to_assembly(self, simple_assembly):
        """Test monomer delegates to wrapped assembly."""
        mon = Monomer(simple_assembly)
        
        # Should access assembly's entities
        assert mon.entities is simple_assembly.entities
        assert mon.atoms == simple_assembly.atoms
    
    def test_monomer_port_management(self, simple_monomer):
        """Test monomer manages ports correctly."""
        assert set(simple_monomer.port_names()) == {"in", "out"}
        
        # Check ports
        in_port = simple_monomer.get_port("in")
        out_port = simple_monomer.get_port("out")
        
        assert isinstance(in_port, Port)
        assert isinstance(out_port, Port)
        assert in_port.name == "in"
        assert out_port.name == "out"
    
    def test_monomer_port_targets(self, simple_monomer):
        """Test port targets point to correct entities."""
        atoms = simple_monomer.unwrap().atoms
        
        in_port = simple_monomer.get_port("in")
        out_port = simple_monomer.get_port("out")
        
        assert in_port.target is atoms[0]
        assert out_port.target is atoms[1]
    
    def test_monomer_add_port(self, simple_assembly):
        """Test adding ports dynamically."""
        mon = Monomer(simple_assembly)
        atoms = simple_assembly.atoms
        
        mon.set_port("port_1", atoms[0])
        mon.set_port("port_2", atoms[1])
        
        assert set(mon.port_names()) == {"port_1", "port_2"}
    
    def test_monomer_update_port(self, simple_assembly):
        """Test updating existing port."""
        mon = Monomer(simple_assembly)
        atoms = simple_assembly.atoms
        
        mon.set_port("test", atoms[0])
        assert mon.get_port("test").target is atoms[0]
        
        # Update to point to different atom
        mon.set_port("test", atoms[1])
        assert mon.get_port("test").target is atoms[1]


# ============================================================================
# Test Monomer Copy Behavior
# ============================================================================


class TestMonomerCopy:
    """Test monomer copy semantics."""
    
    def test_monomer_copy_creates_new_instance(self, simple_monomer):
        """Test copy creates a new monomer instance."""
        copy = simple_monomer.copy()
        
        assert copy is not simple_monomer
        assert isinstance(copy, Monomer)
    
    def test_monomer_copy_deep_copies_assembly(self, simple_monomer):
        """Test copy deep-copies wrapped assembly."""
        copy = simple_monomer.copy()
        
        original_asm = simple_monomer.unwrap()
        copied_asm = copy.unwrap()
        
        assert copied_asm is not original_asm
        assert copied_asm.atoms[0] is not original_asm.atoms[0]
    
    def test_monomer_copy_remaps_ports(self, simple_monomer):
        """Test copy remaps ports to new entities."""
        copy = simple_monomer.copy()
        
        # Original and copied ports should have same names
        assert set(copy.port_names()) == set(simple_monomer.port_names())
        
        # But port targets should be different objects
        original_in = simple_monomer.get_port("in").target
        copied_in = copy.get_port("in").target
        
        assert copied_in is not original_in
        assert copied_in["symbol"] == original_in["symbol"]
    
    def test_monomer_copy_preserves_structure(self, simple_monomer):
        """Test copy preserves molecular structure."""
        copy = simple_monomer.copy()
        
        # Check atom count
        assert len(copy.atoms) == len(simple_monomer.atoms)
        
        # Check bond count
        assert len(copy.bonds) == len(simple_monomer.bonds)


# ============================================================================
# Test Monomer in Wrapper Stack
# ============================================================================


class TestMonomerInStack:
    """Test monomer as part of multi-layer wrapper stack."""
    
    def test_monomer_wrapped_by_tagged(self, simple_monomer):
        """Test wrapping monomer with another wrapper."""
        tagged = TaggedWrapper(simple_monomer, tags=["reactive"])
        
        # Check stack
        stack = tagged.get_stack()
        assert len(stack) == 3
        assert isinstance(stack[0], TaggedWrapper)
        assert isinstance(stack[1], Monomer)
        assert isinstance(stack[2], Atomistic)
    
    def test_access_ports_through_wrapper(self, simple_monomer):
        """Test accessing monomer ports through outer wrapper."""
        tagged = TaggedWrapper(simple_monomer, tags=["test"])
        
        # Should delegate to monomer's port_names
        assert set(tagged.port_names()) == {"in", "out"}
        
        # Should delegate to monomer's get_port
        port = tagged.get_port("in")
        assert isinstance(port, Port)
    
    def test_access_assembly_through_stack(self, simple_monomer):
        """Test accessing assembly through multi-layer stack."""
        tagged = TaggedWrapper(simple_monomer, tags=["test"])
        prop = PropertyWrapper(tagged, properties={"temp": 300})
        
        # Should access through entire stack
        assert len(prop.atoms) == 2
        assert prop.tags == ["test"]
        assert prop.get_property("temp") == 300
        assert set(prop.port_names()) == {"in", "out"}
    
    def test_stack_depth_with_monomer(self, simple_monomer):
        """Test wrapper depth with monomer in stack."""
        tagged = TaggedWrapper(simple_monomer, tags=[])
        
        assert tagged.wrapper_depth() == 2  # Tagged + Monomer
        assert simple_monomer.wrapper_depth() == 1  # Just Monomer


# ============================================================================
# Test Wrapper Consistency
# ============================================================================


class TestWrapperConsistency:
    """Test that wrapped objects behave consistently with unwrapped."""
    
    def test_wrapped_behaves_like_assembly(self, simple_assembly):
        """Test wrapped assembly behaves like original assembly."""
        wrapper = Wrapper(simple_assembly)
        
        # Should have same entities
        assert len(wrapper.atoms) == len(simple_assembly.atoms)
        assert len(wrapper.bonds) == len(simple_assembly.bonds)
    
    def test_monomer_behaves_like_assembly(self, simple_assembly):
        """Test monomer behaves like assembly with added port features."""
        mon = Monomer(simple_assembly)
        
        # Assembly behavior
        assert len(mon.atoms) == 2
        assert len(mon.bonds) == 1
        
        # Extended monomer behavior
        assert hasattr(mon, 'port_names')
        assert hasattr(mon, 'set_port')
    
    def test_modifications_through_wrapper(self, simple_assembly):
        """Test modifying assembly through wrapper."""
        mon = Monomer(simple_assembly)
        
        # Add new atom through monomer
        new_atom = Atom(symbol="O")
        mon.add_entity(new_atom)
        
        # Should reflect in wrapped assembly
        assert new_atom in simple_assembly.atoms
        assert len(mon.atoms) == 3


# ============================================================================
# Test Resolve Ports
# ============================================================================


class TestResolveports:
    """Test port resolution through entity mapping."""
    
    def test_resolve_ports_basic(self, simple_monomer):
        """Test resolving ports through entity map."""
        # Create entity mapping (simulating merge)
        old_atoms = simple_monomer.unwrap().atoms
        new_atoms = [Atom(symbol="C"), Atom(symbol="C")]
        emap = {old_atoms[0]: new_atoms[0], old_atoms[1]: new_atoms[1]}
        
        # Resolve ports
        resolved = simple_monomer.resolve_ports(emap)
        
        assert "in" in resolved
        assert "out" in resolved
        assert resolved["in"] is new_atoms[0]
        assert resolved["out"] is new_atoms[1]
    
    def test_resolve_ports_partial_mapping(self, simple_monomer):
        """Test resolving ports with partial entity map."""
        old_atoms = simple_monomer.unwrap().atoms
        new_atom = Atom(symbol="C")
        emap = {old_atoms[0]: new_atom}  # Only map first atom
        
        resolved = simple_monomer.resolve_ports(emap)
        
        # Only "in" port should be resolved
        assert "in" in resolved
        assert "out" not in resolved
        assert resolved["in"] is new_atom
