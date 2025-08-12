import numpy as np
import pytest

from molpy.analysis.base import ComputeContext
from molpy.analysis.neighbor.compute import (
    NeighborList,
    NeighborResult,
    UpdatableNeighborList,
)
from molpy.core.frame import Frame


class TestUpdatableNeighborList:
    """Test cases for UpdatableNeighborList wrapper class."""

    def setup_method(self):
        """Set up test fixtures."""
        # Create a simple test frame with 4 atoms in a square
        positions = np.array(
            [
                [0.0, 0.0, 0.0],  # atom 0
                [1.0, 0.0, 0.0],  # atom 1
                [0.0, 1.0, 0.0],  # atom 2
                [1.0, 1.0, 0.0],  # atom 3
            ]
        )

        # Create box (orthogonal)
        box = np.array([[5.0, 0.0, 0.0], [0.0, 5.0, 0.0], [0.0, 0.0, 5.0]])

        # Create frame
        self.frame = Frame()
        self.frame.box = box
        self.frame["atoms"] = {
            "id": np.arange(4),
            "position": positions,
            "type": np.array([1, 1, 1, 1]),
        }

        # Create context
        self.context = ComputeContext.attach_frame(self.frame)

    def test_updatable_wrapper_creation(self):
        """Test creating UpdatableNeighborList wrapper."""
        updatable = UpdatableNeighborList(
            name="updatable_base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        assert isinstance(updatable, UpdatableNeighborList)
        assert updatable.name == "updatable_base"

    def test_incremental_add_atom(self):
        """Test adding atoms incrementally."""
        updatable = UpdatableNeighborList(
            name="updatable_add", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        # Initial computation
        result = updatable(self.context)
        initial_count = len(result.indices)

        # Add a new atom
        result.add_atom(atom_id=4, position=np.array([2.0, 2.0, 0.0]), atom_type=1)

        # Check that neighbor data was updated
        assert len(result.frame["atoms"]["id"]) == 5
        assert not result.is_dirty()  # Should be cleared after update

        # Verify that we have more neighbors now
        new_count = len(result.indices)
        assert new_count >= initial_count

    def test_incremental_remove_atom(self):
        """Test removing atoms incrementally."""
        updatable = UpdatableNeighborList(
            name="base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        # Initial computation
        result = updatable(self.context)
        initial_count = len(result.indices)

        # Remove an atom
        result.remove_atom(atom_id=1)

        # Check that atom was removed
        assert len(result.frame["atoms"]["id"]) == 3
        assert 1 not in result.frame["atoms"]["id"]

        # Check that neighbor data was updated
        assert not result.is_dirty()

    def test_incremental_update_position(self):
        """Test updating atom positions incrementally."""
        updatable = UpdatableNeighborList(
            name="base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        # Initial computation
        result = updatable(self.context)
        initial_neighbors = result.indices.copy()

        # Update atom position
        new_position = np.array([3.0, 0.0, 0.0])  # Move atom 1 further away
        result.update_atom_position(atom_id=1, new_position=new_position)

        # Check that position was updated
        atom_idx = np.where(result.frame["atoms"]["id"] == 1)[0][0]
        assert np.allclose(result.frame["atoms"]["position"][atom_idx], new_position)

        # Check that neighbor data was updated
        assert not result.is_dirty()

    def test_dirty_tracking(self):
        """Test dirty flag tracking."""
        updatable = UpdatableNeighborList(
            name="base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        result = updatable(self.context)

        # Initially should not be dirty
        assert not result.is_dirty()

        # Mark atom as modified
        result.mark_atom_modified(0)
        assert result.is_dirty()
        assert 0 in result.dirty_atoms
        assert 0 in result.modified_atoms

        # Clear dirty flags
        result.clear_dirty_flags()
        assert not result.is_dirty()

    def test_context_chaining_with_incremental(self):
        """Test context chaining with incremental updates."""
        updatable = UpdatableNeighborList(
            name="base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        # Initial computation
        result1 = updatable(self.context)

        # Add atom to result1
        result1.add_atom(atom_id=4, position=np.array([2.0, 2.0, 0.0]), atom_type=1)

        # Create another computation that uses result1
        result2 = updatable(result1)

        # Check context chain
        assert result2.get_depth() == 2
        assert result2.unwrap() == result1

        # Both should have the new atom
        assert len(result1.frame["atoms"]["id"]) == 5
        assert len(result2.frame["atoms"]["id"]) == 5

    def test_incremental_disabled(self):
        """Test behavior when incremental updates are disabled."""
        updatable = UpdatableNeighborList(
            name="updatable_disabled",
            method="linkedcell",
            r_cut=2.0,
            cell_size=1.5,
            enable_incremental=False,
        )

        # Initial computation
        result = updatable(self.context)

        # Add atom (should not trigger automatic update)
        result.add_atom(atom_id=4, position=np.array([2.0, 2.0, 0.0]), atom_type=1)

        # Should still be dirty since no automatic update
        assert result.is_dirty()

        # Manual computation should work
        result2 = updatable(result)
        assert not result2.is_dirty()

    def test_method_passthrough(self):
        """Test that method calls are properly passed through to base compute."""
        updatable = UpdatableNeighborList(
            name="base", method="aabb", r_cut=2.0, cell_size=1.5
        )

        # Should use the same method as base
        assert updatable.method == "aabb"
        assert updatable.r_cut == 2.0

        # Computation should work the same
        result = updatable(self.context)
        assert len(result.indices) > 0

    def test_multiple_updates(self):
        """Test multiple incremental updates in sequence."""
        updatable = UpdatableNeighborList(
            name="base", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        result = updatable(self.context)

        # Multiple updates
        result.add_atom(atom_id=4, position=np.array([2.0, 2.0, 0.0]), atom_type=1)
        result.add_atom(atom_id=5, position=np.array([3.0, 3.0, 0.0]), atom_type=1)
        result.update_atom_position(atom_id=0, new_position=np.array([0.1, 0.0, 0.0]))

        # Should have 6 atoms and clean neighbor data
        assert len(result.frame["atoms"]["id"]) == 6
        assert not result.is_dirty()

        # Verify neighbor data is valid
        assert len(result.indices) > 0
        assert len(result.distances) > 0
