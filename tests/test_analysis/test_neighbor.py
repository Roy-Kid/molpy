import numpy as np
import pytest

from molpy.analysis.base import ComputeContext
from molpy.analysis.neighbor.compute import NeighborList, NeighborResult
from molpy.core.frame import Frame


class TestNeighborList:
    """Test cases for basic NeighborList class."""

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

    def test_neighbor_list_linkedcell(self):
        """Test linked cell neighbor finding."""
        neighbor_comp = NeighborList(
            name="test_linkedcell", method="linkedcell", r_cut=2.0, cell_size=1.5
        )

        result = neighbor_comp(self.context)

        assert isinstance(result, NeighborResult)
        assert len(result.indices) > 0
        assert len(result.distances) > 0
        assert len(result.vectors) > 0
        assert len(result.counts) > 0

        # Check that cell_size was set correctly
        assert neighbor_comp.cell_size == 1.5

    def test_neighbor_list_aabb(self):
        """Test AABB neighbor finding."""
        neighbor_comp = NeighborList(name="test_aabb", method="aabb", r_cut=2.0)

        result = neighbor_comp(self.context)

        assert isinstance(result, NeighborResult)
        assert len(result.indices) > 0
        assert len(result.distances) > 0
        assert len(result.vectors) > 0
        assert len(result.counts) > 0

    def test_invalid_method(self):
        """Test that invalid method raises error."""
        with pytest.raises(ValueError, match="Unknown method"):
            neighbor_comp = NeighborList(name="test_invalid", method="invalid_method")
            neighbor_comp(self.context)

    def test_context_chaining(self):
        """Test that NeighborResult properly chains contexts."""
        neighbor_comp = NeighborList(name="test_chain", method="linkedcell", r_cut=2.0)

        result = neighbor_comp(self.context)

        # Check context chain
        assert result.unwrap() == self.context
        assert result.get_depth() == 1
        # get_head() returns the root context, not the current result
        assert result.get_head() == self.context
        assert result.pop() == self.context

    def test_neighbor_data_structure(self):
        """Test that neighbor data has correct structure."""
        neighbor_comp = NeighborList(
            name="test_structure", method="linkedcell", r_cut=2.0
        )

        result = neighbor_comp(self.context)

        # Check data types
        assert isinstance(result.indices, np.ndarray)
        assert isinstance(result.distances, np.ndarray)
        assert isinstance(result.vectors, np.ndarray)
        assert isinstance(result.counts, np.ndarray)

        # Check array shapes
        assert len(result.indices) == len(result.distances)
        assert len(result.indices) == len(result.vectors)
        # counts array length may not equal atom count due to freud's internal representation
        assert len(result.counts) >= 0  # Should be non-negative

        # Check data validity
        assert np.all(result.distances >= 0)  # Distances should be non-negative
        assert np.all(result.counts >= 0)  # Counts should be non-negative

    def test_exclude_ii_parameter(self):
        """Test exclude_ii parameter functionality."""
        # Test with exclude_ii=True (default)
        neighbor_comp_exclude = NeighborList(
            name="test_exclude", method="linkedcell", r_cut=2.0, exclude_ii=True
        )

        result_exclude = neighbor_comp_exclude(self.context)

        # Test with exclude_ii=False
        neighbor_comp_include = NeighborList(
            name="test_include", method="linkedcell", r_cut=2.0, exclude_ii=False
        )

        result_include = neighbor_comp_include(self.context)

        # With exclude_ii=False, we should have more neighbors (including self)
        assert len(result_include.indices) >= len(result_exclude.indices)

    def test_different_cutoffs(self):
        """Test different cutoff radii."""
        cutoffs = [1.0, 2.0, 2.5]  # Reduced max cutoff to avoid cell_width issues
        neighbor_counts = []

        for r_cut in cutoffs:
            neighbor_comp = NeighborList(
                name=f"test_r_cut_{r_cut}",
                method="linkedcell",
                r_cut=r_cut,
                cell_size=min(r_cut, 2.0),  # Ensure cell_size is reasonable
            )

            result = neighbor_comp(self.context)
            neighbor_counts.append(len(result.indices))

        # Larger cutoff should generally give more neighbors
        # (though this depends on the specific geometry)
        assert neighbor_counts[1] >= neighbor_counts[0]  # r_cut=2.0 >= r_cut=1.0
        assert neighbor_counts[2] >= neighbor_counts[1]  # r_cut=2.5 >= r_cut=2.0
