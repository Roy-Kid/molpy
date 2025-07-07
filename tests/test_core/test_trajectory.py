#!/usr/bin/env python3
"""Unit tests for the refactored Trajectory and TrajectoryReader classes.

This module contains comprehensive tests for:
- Trajectory class functionality (caching, indexing, metadata)
- TrajectoryReader base class functionality
- LammpsTrajectoryReader implementation
- Integration between Trajectory and TrajectoryReader
- Edge cases and error handling

Uses pytest framework with modern Python 3.10+ type hints and Google-style docstrings.
"""
import numpy as np
import pytest

from molpy.core.box import Box
from molpy.core.frame import Block, Frame
from molpy.core.trajectory import Trajectory


class TestTrajectory:
    """Test suite for the refactored Trajectory class."""

    @pytest.fixture
    def trajectory(self) -> Trajectory:
        """Create an empty trajectory instance.
        
        Returns:
            Empty Trajectory instance for testing.
        """
        return Trajectory()

    @pytest.fixture
    def test_frames(self) -> list[Frame]:
        """Create test frame instances.
        
        Returns:
            List of 5 test Frame instances with random atom data.
        """
        frames = []
        for i in range(5):
            atoms_data = {
                'id': np.arange(1, 11),  # 10 atoms
                'type': np.ones(10, dtype=int),
                'x': np.random.random(10),
                'y': np.random.random(10),
                'z': np.random.random(10)
            }
            
            box = Box(
                matrix=np.diag([10.0, 10.0, 10.0]),
                origin=np.zeros(3)
            )
            
            frame = Frame(box=box, timestep=i)
            frame["atoms"] = Block(atoms_data)
            frames.append(frame)
        return frames
    
    def test_empty_trajectory_initialization(self, trajectory: Trajectory) -> None:
        """Test creating an empty trajectory.
        
        Args:
            trajectory: Empty trajectory fixture.
        """
        assert len(trajectory) == 0
        assert len(trajectory.get_loaded_indices()) == 0
        assert trajectory._total_frames is None

    def test_trajectory_with_frames_initialization(self, test_frames: list[Frame]) -> None:
        """Test creating trajectory with initial frames.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:3])
        assert len(traj.get_loaded_indices()) == 3
        assert traj.get_loaded_indices() == [0, 1, 2]

    def test_append_frame(self, trajectory: Trajectory, test_frames: list[Frame]) -> None:
        """Test appending frames to trajectory.
        
        Args:
            trajectory: Empty trajectory fixture.
            test_frames: List of test frames fixture.
        """
        trajectory.append(test_frames[0])
        assert len(trajectory.get_loaded_indices()) == 1
        assert 0 in trajectory.frames

    def test_frame_access_by_index(self, test_frames: list[Frame]) -> None:
        """Test accessing frames by index.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:3])
        
        # Test positive indexing
        frame = traj[1]
        assert isinstance(frame, Frame)
        
        # Test negative indexing (requires total frames to be set)
        traj.set_total_frames(3)
        frame = traj[-1]
        assert isinstance(frame, Frame)

    def test_frame_access_unloaded_raises_keyerror(self) -> None:
        """Test that accessing unloaded frames raises KeyError."""
        traj = Trajectory()
        traj.set_total_frames(5)
        
        with pytest.raises(KeyError):
            _ = traj[2]

    def test_slice_access(self, test_frames: list[Frame]) -> None:
        """Test slice access returns new trajectory.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:4])
        
        sliced = traj[1:3]
        assert isinstance(sliced, Trajectory)
        assert len(sliced.get_loaded_indices()) == 2

    def test_is_loaded(self, trajectory: Trajectory, test_frames: list[Frame]) -> None:
        """Test checking if frame is loaded.
        
        Args:
            trajectory: Empty trajectory fixture.
            test_frames: List of test frames fixture.
        """
        trajectory.append(test_frames[0])
        
        assert trajectory.is_loaded(0) is True
        assert trajectory.is_loaded(1) is False

    def test_need_more(self, test_frames: list[Frame]) -> None:
        """Test need_more method.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:2])
        
        # We have 2 frames loaded, need more if total > 2
        assert traj.need_more(5) is True  # This sets total frames to 5
        
        # Now we know there are 5 total frames and 2 loaded
        assert traj.need_more(5) is True
        
        # If total equals loaded, no need for more
        assert traj.need_more(2) is False

    def test_cache_size_limit(self, test_frames: list[Frame]) -> None:
        """Test LRU cache behavior with size limits.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(max_cache_size=2)
        
        # Add 3 frames - should evict the first one
        for i in range(3):
            traj.append(test_frames[i])
        
        # Should only have 2 frames cached
        assert len(traj.frames) == 2
        assert traj.is_loaded(0) is False  # First frame should be evicted

    def test_lru_access_order(self, test_frames: list[Frame]) -> None:
        """Test that LRU eviction works correctly.
        
        Args:
            test_frames: List of test frames fixture.
        """
        # Create trajectory with cache size 3, add 3 frames
        traj = Trajectory(max_cache_size=3)
        
        # Add frames manually to have more control
        for i in range(3):
            traj._add_frame(i, test_frames[i])
        
        # Access frame 0 to make it recently used
        _ = traj[0]
        
        # Add another frame - should evict frame 1, not frame 0
        traj._add_frame(3, test_frames[3])
        
        assert traj.is_loaded(0) is True   # Recently accessed
        assert traj.is_loaded(1) is False  # Should be evicted (LRU)
        assert traj.is_loaded(2) is True   # Still there
        assert traj.is_loaded(3) is True   # Just added

    def test_clear_cache(self, test_frames: list[Frame]) -> None:
        """Test clearing the cache.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:3])
        traj.clear_cache()
        
        assert len(traj.frames) == 0
        assert len(traj.get_loaded_indices()) == 0

    def test_copy(self, test_frames: list[Frame]) -> None:
        """Test copying trajectory.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:2], test_meta="value")
        copy_traj = traj.copy()
        
        assert len(copy_traj.get_loaded_indices()) == 2
        assert copy_traj.meta.get("test_meta") == "value"
        assert traj is not copy_traj

    def test_iteration(self, test_frames: list[Frame]) -> None:
        """Test iterating over trajectory.
        
        Args:
            test_frames: List of test frames fixture.
        """
        traj = Trajectory(test_frames[:3])
        frames = list(traj)
        
        assert len(frames) == 3
        for frame in frames:
            assert isinstance(frame, Frame)
