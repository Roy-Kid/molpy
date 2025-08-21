"""
Comprehensive tests for the MolPy wrapper system.

This module tests both explicit wrapping and multiple inheritance auto-composition
using the actual MolPy wrapper classes. Tests cover:
- Basic wrapper functionality
- Multiple inheritance composition
- Wrapper chains and delegation
- Post-initialization hooks
- Wrapper independence and state management
"""

import numpy as np
import pytest

from molpy.core.atomistic import Atomistic
from molpy.core.protocol import Struct
from molpy.core.wrapper import Spatial, Wrapper


class TestWrapperBasic:
    """Test basic wrapper functionality."""

    def test_explicit_wrapping_works(self):
        """Test that explicit wrapping still works as before."""
        # Create base struct
        base = Struct(name="test_struct")

        # Wrap with Spatial
        spatial_wrapper = Spatial(base)
        assert spatial_wrapper.name == "test_struct"

        # Test that wrapper chain works
        assert spatial_wrapper.wrapper_depth() == 1
        assert "Spatial" in spatial_wrapper.wrapper_types()

    def test_wrapped_passthrough(self):
        """Test that delegation methods remain functional."""
        # Test __getattr__
        wrapper = Atomistic(name="test")
        assert hasattr(wrapper, "name")
        assert wrapper.name == "test"

        # Test __getitem__
        assert wrapper["name"] == "test"

        # Test __setitem__
        wrapper["test_key"] = "test_value"
        assert wrapper["test_key"] == "test_value"

        # Test unwrap
        wrapped = wrapper.unwrap()
        assert isinstance(wrapped, Struct)

    def test_struct_creation_when_wrapped_is_none(self):
        """Test that Struct is created when wrapped is None."""
        wrapper = Atomistic(name="test")

        # Should have created a Struct internally
        wrapped = wrapper.unwrap()
        assert isinstance(wrapped, Struct)
        assert wrapped["name"] == "test"

        # Should have collections initialized
        assert "atoms" in wrapper
        assert "bonds" in wrapper
        assert "angles" in wrapper
        assert "dihedrals" in wrapper


class TestWrapperComposition:
    """Test multiple inheritance auto-composition."""

    def test_multiple_inheritance_auto_composition(self):
        """Test that multiple inheritance auto-composition works."""

        # Define class with multiple inheritance using actual MolPy classes
        class NewWrappedStruct(Atomistic, Spatial):
            def __init__(self, name="test"):
                super().__init__(name=name)

            def __post_init__(self):
                """Custom __post_init__ for the combined class."""
                self.combined_initialized = True

        # Create instance
        combined = NewWrappedStruct("combined_test")

        # Test that it has both Atomistic and Spatial capabilities
        assert hasattr(combined, "atoms")
        assert hasattr(combined, "bonds")
        assert hasattr(combined, "move")
        assert hasattr(combined, "rotate")

        # Test that all __post_init__ methods were called
        assert combined.combined_initialized

        # Test that collections are initialized
        assert "atoms" in combined
        assert "bonds" in combined
        assert "angles" in combined
        assert "dihedrals" in combined

    def test_post_init_runs_once(self):
        """Test that __post_init__ is called exactly once per instance."""
        call_count = 0

        class TestWrapper(Atomistic):
            def __post_init__(self):
                nonlocal call_count
                call_count += 1

        # Create instance
        wrapper = TestWrapper(name="test")
        assert call_count == 1

        # Create another instance
        wrapper2 = TestWrapper(name="test2")
        assert call_count == 2

    def test_multiple_inheritance_mro_safe(self):
        """Test that multiple inheritance respects MRO."""

        class BaseWrapper(Atomistic):
            def __post_init__(self):
                self.base_initialized = True

        class MiddleWrapper(BaseWrapper):
            def __post_init__(self):
                super().__post_init__()
                self.middle_initialized = True

        class TopWrapper(MiddleWrapper, Spatial):
            def __post_init__(self):
                super().__post_init__()
                self.top_initialized = True

        # Create instance
        wrapper = TopWrapper(name="test")

        # Check that all __post_init__ methods were called
        assert wrapper.base_initialized
        assert wrapper.middle_initialized
        assert wrapper.top_initialized

        # Check that it has both wrapper capabilities
        assert hasattr(wrapper, "atoms")
        assert hasattr(wrapper, "move")

    def test_wrapper_bases_recording(self):
        """Test that wrapper bases are correctly recorded in __new__."""

        class TestWrapper1(Atomistic):
            pass

        class TestWrapper2(TestWrapper1, Spatial):
            pass

        # Create instance
        wrapper = TestWrapper2(name="test")

        # Check that wrapper bases are recorded
        assert hasattr(wrapper, "_wrapper_bases")
        assert len(wrapper._wrapper_bases) == 2  # Atomistic and Spatial
        assert Atomistic in wrapper._wrapper_bases
        assert Spatial in wrapper._wrapper_bases

    def test_wrappers_initialized_flag(self):
        """Test that _wrappers_initialized flag prevents double initialization."""

        class TestWrapper(Atomistic):
            def __init__(self, **props):
                super().__init__(**props)
                # The flag should be True after initialization
                assert self._wrappers_initialized
                # Collections should be initialized
                assert hasattr(self, "atoms")

        wrapper = TestWrapper(name="test")
        # Flag should remain True
        assert wrapper._wrappers_initialized
        assert "atoms" in wrapper


class TestWrapperChains:
    """Test wrapper chain functionality and independence."""

    def test_wrapper_chain_functionality(self):
        """Test wrapper chain functionality."""
        # Create wrapper chain
        base = Atomistic(name="base")
        spatial_wrapper = Spatial(base)

        # Test that spatial_wrapper can access base through Atomistic
        assert spatial_wrapper.name == "base"
        assert hasattr(spatial_wrapper, "atoms")
        assert hasattr(spatial_wrapper, "move")

        # Test data modification through chain
        spatial_wrapper["chain_key"] = "chain_value"
        assert base["chain_key"] == "chain_value"

    def test_wrapper_independence(self):
        """Test that wrappers work independently."""
        base = Atomistic(name="base")

        # Create independent wrappers
        spatial1 = Spatial(base)
        spatial2 = Spatial(base)

        # Modify data through one wrapper
        spatial1["key"] = "value1"

        # Other wrapper should see the change
        assert spatial2["key"] == "value1"

        # But each wrapper maintains its own state
        # Note: depth is 2 because Atomistic is itself a wrapper, and Spatial wraps it
        assert spatial1.wrapper_depth() == 2
        assert spatial2.wrapper_depth() == 2

    def test_wrapper_with_custom_properties(self):
        """Test wrapper with custom properties."""

        class CustomWrapper(Atomistic):
            def __init__(self, name="custom", custom_prop="default"):
                super().__init__(name=name)
                self.custom_prop = custom_prop

            def __post_init__(self):
                super().__post_init__()
                self.custom_initialized = True

        wrapper = CustomWrapper("test", "custom_value")

        # Check custom property
        assert wrapper.custom_prop == "custom_value"
        assert wrapper.custom_initialized

        # Check wrapped functionality
        assert wrapper.name == "test"
        assert "atoms" in wrapper


class TestSpatialOperations:
    """Test spatial operations in wrapper composition."""

    def test_spatial_operations_in_multiple_inheritance(self):
        """Test that spatial operations work in multiple inheritance."""

        class H2O(Atomistic, Spatial):
            def __init__(self, name="water"):
                super().__init__(name=name)

                # Build molecule
                self.oxygen = self.def_atom(name="O", element="O", xyz=[0, 0, 0])
                self.hydrogen1 = self.def_atom(
                    name="H1", element="H", xyz=[0.9572, 0, 0]
                )
                self.hydrogen2 = self.def_atom(
                    name="H2", element="H", xyz=[-0.2400, 0, 0]
                )
                self.def_bond(self.oxygen, self.hydrogen1)
                self.def_bond(self.oxygen, self.hydrogen2)

        # Create instance
        water = H2O("water")

        # Test spatial operations
        assert hasattr(water, "move")
        assert hasattr(water, "rotate")
        assert hasattr(water, "scale")
        assert hasattr(water, "move_to")

        # Test that collections are initialized
        assert len(water.atoms) == 3
        assert len(water.bonds) == 2

    def test_spatial_properties_access(self):
        """Test that spatial properties are accessible through wrapper."""
        # Create a structure with spatial wrapper
        base = Atomistic(name="test")
        base.def_atom(name="C", element="C", xyz=[0, 0, 0])
        base.def_atom(name="H", element="H", xyz=[1, 0, 0])

        spatial_wrapper = Spatial(base)

        # Test spatial properties
        assert hasattr(spatial_wrapper, "xyz")
        assert hasattr(spatial_wrapper, "positions")
        assert hasattr(spatial_wrapper, "symbols")

        # Test that properties work correctly
        assert spatial_wrapper.xyz.shape == (2, 3)
        assert spatial_wrapper.positions.shape == (2, 3)
        assert spatial_wrapper.symbols == ["C", "H"]


class TestWrapperEdgeCases:
    """Test edge cases and error conditions."""

    def test_wrapper_with_empty_collections(self):
        """Test wrapper behavior with empty collections."""
        wrapper = Atomistic(name="empty")

        # Collections should exist but be empty
        assert len(wrapper.atoms) == 0
        assert len(wrapper.bonds) == 0
        assert len(wrapper.angles) == 0
        assert len(wrapper.dihedrals) == 0

    def test_wrapper_method_delegation(self):
        """Test that wrapper methods are properly delegated."""
        wrapper = Atomistic(name="test")

        # Test that wrapper methods work
        assert hasattr(wrapper, "add_atom")
        assert hasattr(wrapper, "def_atom")
        assert hasattr(wrapper, "add_bond")

        # Test that methods can be called
        atom = wrapper.def_atom(name="C", element="C", xyz=[0, 0, 0])
        assert atom.name == "C"
        assert len(wrapper.atoms) == 1

    def test_wrapper_repr_and_str(self):
        """Test wrapper string representations."""
        wrapper = Atomistic(name="test")

        # Test that repr works
        repr_str = repr(wrapper)
        assert "Atomistic" in repr_str
        assert "0 atoms" in repr_str

        # Test that str works
        str_str = str(wrapper)
        assert "Atomistic" in str_str
