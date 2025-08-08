"""Comprehensive tests for the trajectory module."""

import pytest
import numpy as np
from typing import Generator, Iterator
from unittest.mock import Mock

from molpy.core.trajectory import (
    FrameProvider, 
    SizedFrameProvider,
    CachedFrames,
    FrameGenerator,
    Trajectory,
    SplitStrategy,
    FrameIntervalStrategy,
    TimeIntervalStrategy,
    CustomStrategy,
    TrajectorySplitter
)
from molpy.core.frame import Frame


@pytest.fixture
def mock_frame():
    """Create a mock frame for testing."""
    frame = Mock(spec=Frame)
    frame.metadata = {}
    return frame


@pytest.fixture
def mock_frames():
    """Create a list of mock frames for testing."""
    frames = []
    for i in range(10):
        frame = Mock(spec=Frame)
        frame.metadata = {'time': i * 0.5}  # Time: 0.0, 0.5, 1.0, ..., 4.5
        frames.append(frame)
    return frames


@pytest.fixture
def frame_list_provider(mock_frames):
    """Create a list-based frame provider."""
    class FrameList(list):
        def __getitem__(self, key):
            if isinstance(key, int):
                return super().__getitem__(key)
            elif isinstance(key, slice):
                return super().__getitem__(key)
            else:
                raise TypeError(f"Invalid key type: {type(key)}")
                
        def __iter__(self):
            return super().__iter__()
            
        def __len__(self):
            return super().__len__()
    
    return FrameList(mock_frames)


@pytest.fixture
def frame_generator(mock_frames):
    """Create a generator of frames for testing."""
    def gen():
        for frame in mock_frames:
            yield frame
    return gen()


class TestCachedFrames:
    """Test the CachedFrames dict subclass."""
    
    def test_cached_frames_is_dict(self, mock_frame):
        cache = CachedFrames()
        cache[0] = mock_frame
        assert cache[0] is mock_frame
        assert len(cache) == 1


class TestFrameGenerator:
    """Test the FrameGenerator class."""
    
    def test_init(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        assert isinstance(fg._cache, CachedFrames)
        assert not fg._exhausted
    
    def test_iteration(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        result = list(fg)
        assert len(result) == len(mock_frames)
        assert all(isinstance(f, Mock) for f in result)
    
    def test_iteration_caches_frames(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        
        # First iteration should cache frames
        list(fg)
        assert len(fg._cache) == len(mock_frames)
        
        # Second iteration should use cached frames
        result = list(fg)
        assert len(result) == len(mock_frames)
    
    def test_getitem_positive_index(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        
        # Access frame at index 3
        frame = fg[3]
        assert frame is mock_frames[3]
        
        # Should have cached frames 0-3
        assert len(fg._cache) == 4
    
    def test_getitem_negative_index_raises(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        
        with pytest.raises(IndexError, match="Negative indexing not supported"):
            fg[-1]
    
    def test_getitem_out_of_range(self, mock_frames):
        def gen():
            for frame in mock_frames[:3]:  # Only 3 frames
                yield frame
        
        fg = FrameGenerator(gen())
        
        with pytest.raises(IndexError, match="Index 5 out of range"):
            fg[5]
    
    def test_getitem_slice(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        
        # Test basic slice
        result = fg[2:5]
        assert isinstance(result, list)
        assert len(result) == 3
        assert result[0] is mock_frames[2]
        assert result[2] is mock_frames[4]
    
    def test_getitem_slice_with_step(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        
        # Test slice with step
        result = fg[1:8:2]
        assert isinstance(result, list)
        assert len(result) == 4
        assert result[0] is mock_frames[1]
        assert result[1] is mock_frames[3]
        assert result[2] is mock_frames[5]
        assert result[3] is mock_frames[7]
    
    def test_getitem_slice_negative_raises(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        
        with pytest.raises(IndexError, match="Reverse slicing not supported"):
            fg[-3:-1]
    
    def test_getitem_invalid_type(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        
        with pytest.raises(TypeError, match="Invalid key type"):
            fg[1.5]  # type: ignore[arg-type]


class TestTrajectory:
    """Test the Trajectory class."""
    
    def test_init_with_frame_provider(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        assert traj._frames is frame_list_provider
        assert traj._topology is None
    
    def test_init_with_topology(self, frame_list_provider):
        topology = Mock()
        traj = Trajectory(frame_list_provider, topology)
        assert traj._topology is topology
    
    def test_iteration(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        result = list(traj)
        assert len(result) == len(frame_list_provider)
    
    def test_len_with_sized_provider(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        assert len(traj) == len(frame_list_provider)
    
    def test_len_with_generator_raises(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        traj = Trajectory(fg)
        
        with pytest.raises(TypeError, match="Length not available for generator-based"):
            len(traj)
    
    def test_has_length_sized_provider(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        assert traj.has_length() is True
    
    def test_has_length_generator(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        traj = Trajectory(fg)
        assert traj.has_length() is False
    
    def test_getitem_int_returns_frame(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        frame = traj[3]
        assert frame is frame_list_provider[3]
    
    def test_getitem_slice_returns_trajectory(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        sub_traj = traj[2:5]
        
        assert isinstance(sub_traj, Trajectory)
        assert sub_traj._topology is traj._topology
    
    def test_getitem_invalid_type(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        
        with pytest.raises(TypeError, match="Invalid key type"):
            traj[1.5]  # type: ignore[arg-type]
    
    def test_map_function(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        traj = Trajectory(fg)
        
        def transform_frame(frame):
            new_frame = Mock(spec=Frame)
            new_frame.metadata = {'transformed': True}
            return new_frame
        
        mapped_traj = traj.map(transform_frame)
        
        assert isinstance(mapped_traj, Trajectory)
        assert isinstance(mapped_traj._frames, FrameGenerator)
        
        # Test that mapping actually works
        result = list(mapped_traj)
        assert len(result) == len(mock_frames)
        assert all(f.metadata.get('transformed') for f in result)


class TestSplitStrategy:
    """Test the abstract SplitStrategy class."""
    
    def test_abstract_method(self):
        # Cannot instantiate abstract class directly
        with pytest.raises(TypeError, match="Can't instantiate abstract class"):
            SplitStrategy()  # type: ignore[abstract]


class TestFrameIntervalStrategy:
    """Test the FrameIntervalStrategy class."""
    
    def test_init(self):
        strategy = FrameIntervalStrategy(5)
        assert strategy.interval == 5
    
    def test_get_split_indices_sized_trajectory(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        strategy = FrameIntervalStrategy(3)
        
        indices = strategy.get_split_indices(traj)
        expected = [0, 3, 6, 9, 10]  # For 10 frames with interval 3
        assert indices == expected
    
    def test_get_split_indices_exact_multiple(self, frame_list_provider):
        # Create a provider with exactly 9 frames
        provider = frame_list_provider[:9]
        traj = Trajectory(provider)
        strategy = FrameIntervalStrategy(3)
        
        indices = strategy.get_split_indices(traj)
        expected = [0, 3, 6, 9]
        assert indices == expected
    
    def test_get_split_indices_generator_raises(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        traj = Trajectory(fg)
        strategy = FrameIntervalStrategy(3)
        
        with pytest.raises(TypeError, match="Frame interval splitting requires trajectory with known length"):
            strategy.get_split_indices(traj)


class TestTimeIntervalStrategy:
    """Test the TimeIntervalStrategy class."""
    
    def test_init(self):
        strategy = TimeIntervalStrategy(1.0)
        assert strategy.interval == 1.0
    
    def test_get_split_indices_with_time_metadata(self, mock_frames):
        # Mock frames have time metadata: 0.0, 0.5, 1.0, 1.5, 2.0, ...
        provider = mock_frames
        traj = Trajectory(provider)
        strategy = TimeIntervalStrategy(1.0)  # Split every 1.0 time unit
        
        indices = strategy.get_split_indices(traj)
        # Should split at times: 0.0, 1.0, 2.0, 3.0, 4.0
        # Corresponding to frame indices: 0, 2, 4, 6, 8, (10)
        expected = [0, 2, 4, 6, 8, 10]
        assert indices == expected
    
    def test_get_split_indices_no_time_metadata(self, frame_list_provider):
        # Remove time metadata
        for frame in frame_list_provider:
            frame.metadata = {}
        
        traj = Trajectory(frame_list_provider)
        strategy = TimeIntervalStrategy(1.0)
        
        indices = strategy.get_split_indices(traj)
        # Should only have start and end indices
        expected = [0, 10]
        assert indices == expected
    
    def test_get_split_indices_generator_exhausts(self, mock_frames):
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        traj = Trajectory(fg)
        strategy = TimeIntervalStrategy(1.0)
        
        indices = strategy.get_split_indices(traj)
        # Should work but exhaust the generator
        assert len(indices) >= 2
        assert indices[0] == 0
        assert indices[-1] == len(mock_frames)


class TestCustomStrategy:
    """Test the CustomStrategy class."""
    
    def test_init_and_call(self, frame_list_provider):
        def custom_split_func(traj):
            return [0, 5, 10]
        
        strategy = CustomStrategy(custom_split_func)
        traj = Trajectory(frame_list_provider)
        
        indices = strategy.get_split_indices(traj)
        assert indices == [0, 5, 10]


class TestTrajectorySplitter:
    """Test the TrajectorySplitter class."""
    
    def test_init(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        splitter = TrajectorySplitter(traj)
        assert splitter.trajectory is traj
    
    def test_split_with_frame_interval(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        splitter = TrajectorySplitter(traj)
        strategy = FrameIntervalStrategy(3)
        
        segments = splitter.split(strategy)
        
        assert len(segments) == 4  # [0:3], [3:6], [6:9], [9:10]
        assert all(isinstance(seg, Trajectory) for seg in segments)
    
    def test_split_frames_convenience(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        splitter = TrajectorySplitter(traj)
        
        segments = splitter.split_frames(4)
        
        assert len(segments) == 3  # [0:4], [4:8], [8:10]
        assert all(isinstance(seg, Trajectory) for seg in segments)
    
    def test_split_time_convenience(self, mock_frames):
        provider = mock_frames
        traj = Trajectory(provider)
        splitter = TrajectorySplitter(traj)
        
        segments = splitter.split_time(1.0)
        
        assert len(segments) >= 1
        assert all(isinstance(seg, Trajectory) for seg in segments)
    
    def test_split_preserves_topology(self, frame_list_provider):
        topology = Mock()
        traj = Trajectory(frame_list_provider, topology)
        splitter = TrajectorySplitter(traj)
        
        segments = splitter.split_frames(3)
        
        assert all(seg._topology is topology for seg in segments)


class TestProtocolCompliance:
    """Test that classes properly implement the protocols."""
    
    def test_frame_generator_implements_frame_provider(self, frame_generator):
        fg = FrameGenerator(frame_generator)
        assert isinstance(fg, FrameProvider)
    
    def test_list_implements_sized_frame_provider(self, frame_list_provider):
        assert isinstance(frame_list_provider, SizedFrameProvider)
        assert isinstance(frame_list_provider, FrameProvider)


class TestErrorHandling:
    """Test error handling throughout the module."""
    
    def test_frame_generator_handles_empty_generator(self):
        def empty_gen():
            return
            yield  # This will never be reached
        
        fg = FrameGenerator(empty_gen())
        result = list(fg)
        assert result == []
        assert fg._exhausted
    
    def test_trajectory_slicing_edge_cases(self, frame_list_provider):
        traj = Trajectory(frame_list_provider)
        
        # Empty slice
        empty_traj = traj[5:5]
        assert isinstance(empty_traj, Trajectory)
        
        # Slice beyond bounds
        beyond_traj = traj[8:20]
        assert isinstance(beyond_traj, Trajectory)


class TestIntegration:
    """Integration tests combining multiple components."""
    
    def test_full_workflow_with_frame_generator(self, mock_frames):
        # Create trajectory from generator
        def gen():
            for frame in mock_frames:
                yield frame
        
        fg = FrameGenerator(gen())
        traj = Trajectory(fg)
        
        # Apply transformation
        def add_id(frame):
            new_frame = Mock(spec=Frame)
            new_frame.metadata = frame.metadata.copy()
            new_frame.metadata['id'] = id(frame)
            return new_frame
        
        transformed_traj = traj.map(add_id)
        
        # Collect results
        result = list(transformed_traj)
        assert len(result) == len(mock_frames)
        assert all('id' in f.metadata for f in result)
    
    def test_full_workflow_with_splitting(self, frame_list_provider):
        # Create trajectory
        topology = Mock()
        traj = Trajectory(frame_list_provider, topology)
        
        # Split trajectory
        splitter = TrajectorySplitter(traj)
        segments = splitter.split_frames(3)
        
        # Verify segments
        assert len(segments) == 4
        
        # Test accessing frames in segments
        first_segment = segments[0]
        frames_in_first = list(first_segment)
        assert len(frames_in_first) == 3
        
        # Verify topology is preserved
        assert all(seg._topology is topology for seg in segments)
