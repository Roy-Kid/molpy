from abc import ABC, abstractmethod
from typing import (
    Callable,
    Generator,
    Iterator,
    Protocol,
    Sequence,
    Union,
    overload,
    runtime_checkable,
)

import numpy as np

from .frame import Frame


@runtime_checkable
class FrameProvider(Protocol):
    """Protocol for frame providers that can be indexed."""

    def __getitem__(self, key: int | slice) -> Frame | list[Frame]: ...
    def __iter__(self) -> Iterator[Frame]: ...


@runtime_checkable
class SizedFrameProvider(FrameProvider, Protocol):
    """Protocol for frame providers that have a known length."""

    def __len__(self) -> int: ...


class CachedFrames(dict[int, Frame]): ...


class FrameGenerator(FrameProvider):
    """Generator-based frame provider with caching. No predefined length."""

    def __init__(self, frames: Generator[Frame, None, None]):
        """Initialize with a generator of frames."""
        self._frames = frames
        self._cache = CachedFrames()
        self._exhausted = False

    def __iter__(self) -> Iterator[Frame]:
        idx = 0
        # If we're exhausted, just iterate through cached frames
        if self._exhausted:
            for i in sorted(self._cache.keys()):
                yield self._cache[i]
            return

        # Otherwise, continue from where we left off
        while not self._exhausted:
            if idx in self._cache:
                yield self._cache[idx]
            else:
                try:
                    frame = next(self._frames)
                    self._cache[idx] = frame
                    yield frame
                except StopIteration:
                    self._exhausted = True
                    break
            idx += 1

    def __getitem__(self, key: int | slice) -> Frame | list[Frame]:
        if isinstance(key, int):
            if key < 0:
                raise IndexError(
                    "Negative indexing not supported for generator-based frames. "
                    "Consider using Framelist for full indexing support."
                )

            if key in self._cache:
                return self._cache[key]

            # Generate frames up to the requested index
            for idx, frame in enumerate(self):
                if idx == key:
                    return frame
            raise IndexError(f"Index {key} out of range (generator exhausted)")
        elif isinstance(key, slice):
            if key.start < 0 or key.stop < 0:
                raise IndexError("Reverse slicing not supported for generators")
            start = key.start or 0
            stop = key.stop
            step = key.step or 1

            result = []
            for i in range(start, stop, step):
                try:
                    result.append(self[i])
                except IndexError:
                    break
            return result
        else:
            raise TypeError(f"Invalid key type: {type(key)}")


class Trajectory:
    """A sequence of molecular frames with optional topology."""

    def __init__(self, frames: FrameProvider, topology=None):
        """Initialize trajectory with a frame provider."""
        self._frames = frames
        self._topology = topology

    def __iter__(self) -> Iterator[Frame]:
        return iter(self._frames)

    def __len__(self) -> int:
        """Return the number of frames in the trajectory."""
        if isinstance(self._frames, SizedFrameProvider):
            return len(self._frames)
        else:
            raise TypeError(
                "Length not available for generator-based trajectories. "
                "Use count_frames() to exhaust and count if needed."
            )

    def has_length(self) -> bool:
        """Check if this trajectory has a known length without computing it."""
        return isinstance(self._frames, SizedFrameProvider)

    def __getitem__(self, key: int | slice) -> "Frame | Trajectory":
        if isinstance(key, int):
            frame = self._frames[key]
            # Protocol guarantees Frame for int index
            return frame  # type: ignore[return-value]
        elif isinstance(key, slice):
            sliced_frames = self._frames[key]
            return Trajectory(sliced_frames, self._topology)  # type: ignore[arg-type]
        else:
            raise TypeError(f"Invalid key type: {type(key)}")

    def map(self, func: Callable[[Frame], Frame]) -> "Trajectory":
        """Apply a function to each frame, returning a new trajectory."""

        def mapped_generator() -> Generator[Frame, None, None]:
            for frame in self._frames:
                yield func(frame)

        return Trajectory(FrameGenerator(mapped_generator()), self._topology)


# ====================== Trajectory Splitters ====================


class SplitStrategy(ABC):
    """Abstract splitting strategy."""

    @abstractmethod
    def get_split_indices(self, trajectory: Trajectory) -> list[int]:
        """Get split point indices."""
        raise NotImplementedError


class FrameIntervalStrategy(SplitStrategy):
    """Split every N frames."""

    def __init__(self, interval: int):
        self.interval = interval

    def get_split_indices(self, trajectory: Trajectory) -> list[int]:
        if not trajectory.has_length():
            raise TypeError(
                "Frame interval splitting requires trajectory with known length"
            )

        length = len(trajectory)
        indices = list(range(0, length, self.interval))
        if indices[-1] != length:
            indices.append(length)
        return indices


class TimeIntervalStrategy(SplitStrategy):
    """Split by time intervals."""

    def __init__(self, interval: float):
        self.interval = interval

    def get_split_indices(self, trajectory: Trajectory) -> list[int]:
        indices = [0]
        start_time = None
        frame_count = 0

        for i, frame in enumerate(trajectory):
            frame_count = i + 1  # Keep track of total frames seen
            # Check if frame has time information in metadata
            frame_time = frame.metadata.get("time", None)
            if frame_time is not None:
                if start_time is None:
                    start_time = frame_time

                # Check if we've reached the next interval
                if frame_time >= start_time + len(indices) * self.interval:
                    indices.append(i)

        # Add final index using the count from iteration
        final_index = frame_count
        if indices[-1] != final_index:
            indices.append(final_index)

        return indices


class CustomStrategy(SplitStrategy):
    """Split using custom function."""

    def __init__(self, split_func: Callable[[Trajectory], list[int]]):
        self.split_func = split_func

    def get_split_indices(self, trajectory: Trajectory) -> list[int]:
        return self.split_func(trajectory)


class TrajectorySplitter:
    """Splits trajectories into lazy segments."""

    def __init__(self, trajectory: Trajectory):
        self.trajectory = trajectory

    def split(self, strategy: SplitStrategy) -> list[Trajectory]:
        """Split trajectory using strategy, returning lazy segments."""
        indices = strategy.get_split_indices(self.trajectory)

        segments = []
        for i in range(len(indices) - 1):
            start, end = indices[i], indices[i + 1]
            # Use trajectory slicing instead of accessing private members
            segment = self.trajectory[start:end]
            segments.append(segment)

        return segments

    def split_frames(self, interval: int) -> list[Trajectory]:
        """Convenience method for frame-based splitting."""
        return self.split(FrameIntervalStrategy(interval))

    def split_time(self, interval: float) -> list[Trajectory]:
        """Convenience method for time-based splitting."""
        return self.split(TimeIntervalStrategy(interval))
