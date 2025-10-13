"""
Comprehensive tests for the MolPy Wrapper base class.

Tests core Wrapper functionality including deep copy, attribute lookup,
assignment behavior, and post-init kwargs handling.
"""

import numpy as np
import pytest

from molpy.core.protocol import Struct
from molpy.core.wrapper import Wrapper


# Custom wrapper classes for testing
class CustomWrapperA(Wrapper):
    """Custom wrapper that consumes specific kwargs."""

    def __post_init__(self, **props):
        """Consume prop_a if present."""
        self.prop_a_value = props.pop("prop_a", None)
        return props


class CustomWrapperB(Wrapper):
    """Another custom wrapper for testing multi-layer."""

    def __post_init__(self, **props):
        """Consume prop_b if present."""
        self.prop_b_value = props.pop("prop_b", None)
        return props


class CombinedWrapper(CustomWrapperA, CustomWrapperB):
    """Combined wrapper for testing multiple inheritance."""

    def __post_init__(self, **props):
        """Consume prop_c if present."""
        self.prop_c_value = props.pop("prop_c", None)
        # Don't call super().__post_init__ - the framework handles calling all __post_init__ methods
        return props


class TestWrapperBasic:
    """Test basic wrapper functionality."""

    def test_explicit_wrapping_works(self):
        """Test that explicit wrapping works."""
        base = Struct(name="test_struct")
        wrapper = Wrapper(base)
        assert wrapper.name == "test_struct"
        assert wrapper.wrapper_depth() == 1
        assert "Wrapper" in wrapper.wrapper_types()

    def test_wrapped_passthrough(self):
        """Test that delegation methods remain functional."""
        base = Struct(name="test_struct", value=42)
        wrapper = Wrapper(base)
        assert hasattr(wrapper, "name")
        assert wrapper.name == "test_struct"
        assert wrapper["name"] == "test_struct"
        wrapper["test_key"] = "test_value"
        assert wrapper["test_key"] == "test_value"
        wrapped = wrapper.unwrap()
        assert isinstance(wrapped, Struct)

    def test_get_innermost(self):
        """Test _get_innermost method."""
        base = Struct(name="innermost")
        middle = Wrapper(base)
        outer = Wrapper(middle)
        innermost = outer._get_innermost()
        assert innermost is base
        assert isinstance(innermost, Struct)


class TestDeepCopy:
    """Test requirement 1: .copy() always returns a deep copy."""

    def test_copy_creates_deep_copy(self):
        """Test that copy method creates a deep copy."""
        base = Struct(name="test_struct", data=np.array([1, 2, 3]))
        wrapper = Wrapper(base)
        wrapper["key"] = "value"
        wrapper_copy = wrapper.copy()
        assert wrapper_copy is not wrapper
        assert wrapper_copy.unwrap() is not wrapper.unwrap()
        assert wrapper_copy["key"] == "value"
        wrapper_copy["key"] = "new_value"
        assert wrapper["key"] == "value"
        assert wrapper_copy["key"] == "new_value"
        wrapper_copy["data"][0] = 999
        assert wrapper["data"][0] == 1, "Original should be unaffected"

    def test_copy_with_nested_wrappers(self):
        """Test deep copy with nested wrapper chain."""
        base = Struct(name="innermost", inner_data=[1, 2, 3])
        middle = Wrapper(base)
        outer = Wrapper(middle)
        outer["middle_data"] = {"key": "value"}
        copied = outer.copy()
        base["inner_data"].append(4)
        outer["middle_data"]["key"] = "modified"
        assert len(copied["inner_data"]) == 3
        assert copied["middle_data"]["key"] == "value"


class TestAttributeLookup:
    """Test requirement 2: Attribute/key lookup searches from outer to inner."""

    def test_getattr_searches_outer_to_inner(self):
        """Test __getattr__ searches from outer layers to inner layers."""
        base = Struct(inner_attr="inner_value", shared_attr="from_inner")
        wrapper = Wrapper(base)
        wrapper["wrapper_attr"] = "wrapper_value"
        wrapper["shared_attr"] = "from_wrapper"
        assert wrapper.wrapper_attr == "wrapper_value"
        assert wrapper.inner_attr == "inner_value"
        assert wrapper.shared_attr == "from_wrapper"

    def test_getitem_searches_outer_to_inner(self):
        """Test __getitem__ searches from outer layers to inner layers."""
        base = Struct(inner_key="inner_value", shared_key="from_inner")
        wrapper = Wrapper(base)
        wrapper["wrapper_key"] = "wrapper_value"
        wrapper["shared_key"] = "from_wrapper"
        assert wrapper["wrapper_key"] == "wrapper_value"
        assert wrapper["inner_key"] == "inner_value"
        assert wrapper["shared_key"] == "from_wrapper"

    def test_getitem_raises_keyerror_when_not_found(self):
        """Test __getitem__ raises KeyError when key doesn't exist."""
        base = Struct(name="test")
        wrapper = Wrapper(base)
        with pytest.raises(KeyError):
            _ = wrapper["nonexistent_key"]

    def test_contains_searches_through_chain(self):
        """Test __contains__ searches through wrapper chain."""
        base = Struct(inner_key="value")
        wrapper = Wrapper(base)
        wrapper["outer_key"] = "value"
        assert "inner_key" in wrapper
        assert "outer_key" in wrapper
        assert "nonexistent" not in wrapper


class TestAttributeAssignment:
    """Test requirement 3: Assignment searches outer to inner, creates in innermost if not found."""

    def test_setattr_modifies_existing_in_place(self):
        """Test __setattr__ modifies existing attributes where they are."""
        base = Struct(existing_attr="original")
        wrapper = Wrapper(base)
        wrapper.existing_attr = "modified"
        assert base["existing_attr"] == "modified"
        assert wrapper.existing_attr == "modified"

    def test_setattr_creates_in_innermost_when_new(self):
        """Test __setattr__ creates new attributes in innermost Struct."""
        base = Struct(name="test")
        wrapper = Wrapper(base)
        wrapper.new_attr = "new_value"
        assert "new_attr" in base
        assert base["new_attr"] == "new_value"
        assert wrapper.new_attr == "new_value"

    def test_setitem_modifies_existing_in_place(self):
        """Test __setitem__ modifies existing keys where they are."""
        base = Struct(existing_key="original")
        wrapper = Wrapper(base)
        wrapper["existing_key"] = "modified"
        assert base["existing_key"] == "modified"
        assert wrapper["existing_key"] == "modified"

    def test_setitem_creates_in_innermost_when_new(self):
        """Test __setitem__ creates new keys in innermost Struct."""
        base = Struct(name="test")
        wrapper = Wrapper(base)
        wrapper["new_key"] = "new_value"
        assert "new_key" in base
        assert base["new_key"] == "new_value"
        assert wrapper["new_key"] == "new_value"


class TestPostInitKwargs:
    """Test requirement 4: Post-init kwargs handling."""

    def test_post_init_consumes_kwargs(self):
        """Test __post_init__ can consume specific kwargs."""
        wrapper = CustomWrapperA(
            wrapped=None, prop_a="value_a", prop_other="value_other"
        )
        assert wrapper.prop_a_value == "value_a"
        innermost = wrapper._get_innermost()
        assert isinstance(innermost, Struct)
        assert "prop_other" in innermost
        assert innermost["prop_other"] == "value_other"
        assert "prop_a" not in innermost

    def test_post_init_with_multiple_wrappers(self):
        """Test __post_init__ with multiple inheritance."""
        wrapper = CombinedWrapper(
            wrapped=None,
            prop_a="value_a",
            prop_b="value_b",
            prop_c="value_c",
            prop_struct="value_struct",
        )
        assert wrapper.prop_a_value == "value_a"
        assert wrapper.prop_b_value == "value_b"
        assert wrapper.prop_c_value == "value_c"
        innermost = wrapper._get_innermost()
        assert isinstance(innermost, Struct)
        assert "prop_struct" in innermost
        assert innermost["prop_struct"] == "value_struct"
        assert "prop_a" not in innermost
        assert "prop_b" not in innermost
        assert "prop_c" not in innermost
