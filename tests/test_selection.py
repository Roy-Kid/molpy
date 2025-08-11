import numpy as np
import pytest

from molpy.core.frame import Block
from molpy.core.selection import (
    AtomIndexSelection,
    AtomTypeSelection,
    CoordinateRangeSelection,
    DistanceSelection,
    ElementSelection,
    MaskPredicate,
    _And,
    _Not,
    _Or,
)


class TestAtomTypeSelection:
    """Test AtomTypeSelection class."""

    def test_init_with_int(self):
        """Test initialization with integer atom type."""
        selection = AtomTypeSelection(1)
        assert selection.atom_type == 1
        assert selection.field == "type"

    def test_init_with_string(self):
        """Test initialization with string atom type."""
        selection = AtomTypeSelection("protein")
        assert selection.atom_type == "protein"
        assert selection.field == "type"

    def test_init_with_custom_field(self):
        """Test initialization with custom field name."""
        selection = AtomTypeSelection(2, field="custom_type")
        assert selection.atom_type == 2
        assert selection.field == "custom_type"

    def test_mask_with_valid_data(self):
        """Test mask generation with valid data."""
        selection = AtomTypeSelection(1)
        block = Block({"type": [1, 2, 1, 3, 1], "x": [0.0, 1.0, 2.0, 3.0, 4.0]})

        mask = selection.mask(block)
        expected = np.array([True, False, True, False, True])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_missing_field(self):
        """Test mask generation when field is missing."""
        selection = AtomTypeSelection(1)
        block = Block({"x": [0.0, 1.0, 2.0]})

        mask = selection.mask(block)
        expected = np.zeros(3, dtype=bool)
        np.testing.assert_array_equal(mask, expected)


class TestAtomIndexSelection:
    """Test AtomIndexSelection class."""

    def test_init_with_list(self):
        """Test initialization with list of indices."""
        selection = AtomIndexSelection([0, 2, 4])
        np.testing.assert_array_equal(selection.indices, np.array([0, 2, 4]))
        assert selection.id_field == "id"

    def test_init_with_numpy_array(self):
        """Test initialization with numpy array of indices."""
        indices = np.array([1, 3, 5])
        selection = AtomIndexSelection(indices)
        np.testing.assert_array_equal(selection.indices, indices)

    def test_init_with_custom_id_field(self):
        """Test initialization with custom ID field name."""
        selection = AtomIndexSelection([0, 1], id_field="custom_id")
        assert selection.id_field == "custom_id"

    def test_init_with_invalid_type(self):
        """Test initialization with invalid type raises error."""
        with pytest.raises(
            TypeError, match="indices must be a list\\[int\\] or np.ndarray"
        ):
            AtomIndexSelection("invalid")

    def test_mask_with_valid_data(self):
        """Test mask generation with valid data."""
        selection = AtomIndexSelection([0, 2, 4])
        block = Block({"id": [0, 1, 2, 3, 4, 5], "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]})

        mask = selection.mask(block)
        expected = np.array([True, False, True, False, True, False])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_missing_field(self):
        """Test mask generation when ID field is missing."""
        selection = AtomIndexSelection([0, 1])
        block = Block({"x": [0.0, 1.0, 2.0]})

        mask = selection.mask(block)
        expected = np.zeros(3, dtype=bool)
        np.testing.assert_array_equal(mask, expected)


class TestElementSelection:
    """Test ElementSelection class."""

    def test_init_with_valid_element(self):
        """Test initialization with valid element symbol."""
        selection = ElementSelection("C")
        assert selection.element == "C"
        assert selection.field == "element"

    def test_init_with_custom_field(self):
        """Test initialization with custom field name."""
        selection = ElementSelection("H", field="atom_element")
        assert selection.element == "H"
        assert selection.field == "atom_element"

    def test_init_with_invalid_type(self):
        """Test initialization with invalid type raises error."""
        with pytest.raises(TypeError, match="element must be a string"):
            ElementSelection(1)

    def test_mask_with_valid_data(self):
        """Test mask generation with valid data."""
        selection = ElementSelection("C")
        block = Block(
            {"element": ["C", "H", "C", "O", "C"], "x": [0.0, 1.0, 2.0, 3.0, 4.0]}
        )

        mask = selection.mask(block)
        expected = np.array([True, False, True, False, True])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_missing_field(self):
        """Test mask generation when element field is missing."""
        selection = ElementSelection("C")
        block = Block({"x": [0.0, 1.0, 2.0]})

        mask = selection.mask(block)
        expected = np.zeros(3, dtype=bool)
        np.testing.assert_array_equal(mask, expected)


class TestCoordinateRangeSelection:
    """Test CoordinateRangeSelection class."""

    def test_init_with_valid_axis(self):
        """Test initialization with valid axis."""
        selection = CoordinateRangeSelection("x", min_value=0.0, max_value=5.0)
        assert selection.axis == "x"
        assert selection.min_value == 0.0
        assert selection.max_value == 5.0

    def test_init_with_invalid_axis(self):
        """Test initialization with invalid axis raises error."""
        with pytest.raises(ValueError, match="axis must be 'x', 'y', or 'z'"):
            CoordinateRangeSelection("w")

    def test_init_with_invalid_range(self):
        """Test initialization with invalid range raises error."""
        with pytest.raises(
            ValueError, match="min_value cannot be greater than max_value"
        ):
            CoordinateRangeSelection("x", min_value=5.0, max_value=0.0)

    def test_mask_with_both_bounds(self):
        """Test mask generation with both min and max bounds."""
        selection = CoordinateRangeSelection("x", min_value=1.0, max_value=4.0)
        block = Block({"x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0]})

        mask = selection.mask(block)
        expected = np.array([False, True, True, True, True, False])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_only_min(self):
        """Test mask generation with only min bound."""
        selection = CoordinateRangeSelection("y", min_value=2.0)
        block = Block({"y": [0.0, 1.0, 2.0, 3.0, 4.0]})

        mask = selection.mask(block)
        expected = np.array([False, False, True, True, True])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_only_max(self):
        """Test mask generation with only max bound."""
        selection = CoordinateRangeSelection("z", max_value=2.0)
        block = Block({"z": [0.0, 1.0, 2.0, 3.0, 4.0]})

        mask = selection.mask(block)
        expected = np.array([True, True, True, False, False])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_missing_field(self):
        """Test mask generation when coordinate field is missing."""
        selection = CoordinateRangeSelection("x", min_value=0.0)
        block = Block({"y": [0.0, 1.0, 2.0]})

        mask = selection.mask(block)
        expected = np.zeros(3, dtype=bool)
        np.testing.assert_array_equal(mask, expected)


class TestDistanceSelection:
    """Test DistanceSelection class."""

    def test_init_with_list_center(self):
        """Test initialization with list center."""
        selection = DistanceSelection([0.0, 0.0, 0.0], max_distance=2.0)
        np.testing.assert_array_equal(selection.center, np.array([0.0, 0.0, 0.0]))
        assert selection.max_distance == 2.0
        assert selection.min_distance is None

    def test_init_with_numpy_center(self):
        """Test initialization with numpy array center."""
        center = np.array([1.0, 2.0, 3.0])
        selection = DistanceSelection(center, max_distance=3.0)
        np.testing.assert_array_equal(selection.center, center)

    def test_init_with_min_distance(self):
        """Test initialization with min distance."""
        selection = DistanceSelection(
            [0.0, 0.0, 0.0], max_distance=5.0, min_distance=1.0
        )
        assert selection.min_distance == 1.0

    def test_init_with_invalid_center_type(self):
        """Test initialization with invalid center type raises error."""
        with pytest.raises(
            TypeError, match="center must be a list\\[float\\] or np.ndarray"
        ):
            DistanceSelection("invalid", max_distance=1.0)

    def test_init_with_invalid_center_length(self):
        """Test initialization with invalid center length raises error."""
        with pytest.raises(ValueError, match="center must have exactly 3 coordinates"):
            DistanceSelection([0.0, 0.0], max_distance=1.0)

    def test_init_with_negative_max_distance(self):
        """Test initialization with negative max distance raises error."""
        with pytest.raises(ValueError, match="max_distance must be non-negative"):
            DistanceSelection([0.0, 0.0, 0.0], max_distance=-1.0)

    def test_init_with_negative_min_distance(self):
        """Test initialization with negative min distance raises error."""
        with pytest.raises(ValueError, match="min_distance must be non-negative"):
            DistanceSelection([0.0, 0.0, 0.0], max_distance=5.0, min_distance=-1.0)

    def test_init_with_invalid_distance_range(self):
        """Test initialization with invalid distance range raises error."""
        with pytest.raises(
            ValueError, match="min_distance cannot be greater than max_distance"
        ):
            DistanceSelection([0.0, 0.0, 0.0], max_distance=1.0, min_distance=2.0)

    def test_mask_with_max_distance_only(self):
        """Test mask generation with only max distance."""
        selection = DistanceSelection([0.0, 0.0, 0.0], max_distance=2.0)
        block = Block(
            {
                "x": [0.0, 1.0, 2.0, 3.0, 4.0],
                "y": [0.0, 0.0, 0.0, 0.0, 0.0],
                "z": [0.0, 0.0, 0.0, 0.0, 0.0],
            }
        )

        mask = selection.mask(block)
        expected = np.array([True, True, True, False, False])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_both_distances(self):
        """Test mask generation with both min and max distances."""
        selection = DistanceSelection(
            [0.0, 0.0, 0.0], max_distance=4.0, min_distance=1.0
        )
        block = Block(
            {
                "x": [0.0, 1.0, 2.0, 3.0, 4.0],
                "y": [0.0, 0.0, 0.0, 0.0, 0.0],
                "z": [0.0, 0.0, 0.0, 0.0, 0.0],
            }
        )

        mask = selection.mask(block)
        expected = np.array([False, True, True, True, True])
        np.testing.assert_array_equal(mask, expected)

    def test_mask_with_missing_fields(self):
        """Test mask generation when coordinate fields are missing."""
        selection = DistanceSelection([0.0, 0.0, 0.0], max_distance=1.0)
        block = Block({"x": [0.0, 1.0, 2.0]})

        mask = selection.mask(block)
        expected = np.zeros(3, dtype=bool)
        np.testing.assert_array_equal(mask, expected)


class TestCombinators:
    """Test logical combinator classes."""

    def test_and_combination(self):
        """Test AND combination of two selections."""
        selection1 = AtomTypeSelection(1)
        selection2 = ElementSelection("C")

        combined = _And(selection1, selection2)
        assert isinstance(combined, _And)
        assert combined.a is selection1
        assert combined.b is selection2

    def test_or_combination(self):
        """Test OR combination of two selections."""
        selection1 = AtomTypeSelection(1)
        selection2 = AtomTypeSelection(2)

        combined = _Or(selection1, selection2)
        assert isinstance(combined, _Or)
        assert combined.a is selection1
        assert combined.b is selection2

    def test_not_combination(self):
        """Test NOT combination of a selection."""
        selection = AtomTypeSelection(1)

        negated = _Not(selection)
        assert isinstance(negated, _Not)
        assert negated.a is selection


class TestMaskPredicate:
    """Test base MaskPredicate class and operations."""

    def test_and_operator(self):
        """Test & operator between two selections."""
        selection1 = AtomTypeSelection(1)
        selection2 = ElementSelection("C")

        combined = selection1 & selection2
        assert isinstance(combined, _And)

    def test_or_operator(self):
        """Test | operator between two selections."""
        selection1 = AtomTypeSelection(1)
        selection2 = AtomTypeSelection(2)

        combined = selection1 | selection2
        assert isinstance(combined, _Or)

    def test_not_operator(self):
        """Test ~ operator on a selection."""
        selection = AtomTypeSelection(1)

        negated = ~selection
        assert isinstance(negated, _Not)

    def test_call_method(self):
        """Test calling selection on a block."""
        selection = AtomTypeSelection(1)
        block = Block({"type": [1, 2, 1, 3], "x": [0.0, 1.0, 2.0, 3.0]})

        result = selection(block)
        expected_mask = np.array([True, False, True, False])
        np.testing.assert_array_equal(result["type"], block["type"][expected_mask])


class TestComplexSelections:
    """Test complex selection combinations."""

    def test_nested_combinations(self):
        """Test nested logical combinations."""
        carbon = ElementSelection("C")
        hydrogen = ElementSelection("H")
        oxygen = ElementSelection("O")

        # (C OR H) AND NOT O
        complex_selection = (carbon | hydrogen) & ~oxygen

        block = Block(
            {
                "element": ["C", "H", "O", "C", "H", "N"],
                "x": [0.0, 1.0, 2.0, 3.0, 4.0, 5.0],
            }
        )

        mask = complex_selection.mask(block)
        expected = np.array([True, True, False, True, True, False])
        np.testing.assert_array_equal(mask, expected)

    def test_spatial_and_property_combination(self):
        """Test combining spatial and property-based selections."""
        carbon = ElementSelection("C")
        in_sphere = DistanceSelection([0.0, 0.0, 0.0], max_distance=2.0)

        # Carbon atoms within sphere
        carbon_in_sphere = carbon & in_sphere

        block = Block(
            {
                "element": ["C", "H", "C", "O", "C"],
                "x": [0.0, 1.0, 2.0, 3.0, 4.0],
                "y": [0.0, 0.0, 0.0, 0.0, 0.0],
                "z": [0.0, 0.0, 0.0, 0.0, 0.0],
            }
        )

        mask = carbon_in_sphere.mask(block)
        expected = np.array([True, False, True, False, False])
        np.testing.assert_array_equal(mask, expected)
