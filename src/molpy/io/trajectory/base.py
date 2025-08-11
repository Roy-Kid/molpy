import mmap
from abc import ABC, abstractmethod
from pathlib import Path
from typing import TYPE_CHECKING, Iterator, List, NamedTuple, Optional, Union

if TYPE_CHECKING:
    from ...core.frame import Frame

PathLike = Union[str, bytes]  # type_check_only


class FrameLocation(NamedTuple):
    """Location information for a frame."""

    file_index: int
    byte_offset: int
    file_path: Path


class TrajectoryReader(ABC):
    """
    Base class for trajectory file readers that act as providers.

    This class provides memory-mapped file reading and directly returns Frame objects
    without needing to interact with Trajectory objects. Supports reading from multiple files.
    """

    def __init__(self, fpath: Union[Path, str, List[Path], List[str]]):
        """
        Initialize the trajectory reader.

        Args:
            fpath: Path to trajectory file or list of paths to multiple trajectory files
        """
        # Handle both single file and multiple files
        if isinstance(fpath, (str, Path)):
            self.fpaths = [Path(fpath)]
        else:
            self.fpaths = [Path(p) for p in fpath]

        # Validate all files exist
        for path in self.fpaths:
            if not path.exists():
                raise FileNotFoundError(f"File not found: {path}")

        self._frame_locations: List[FrameLocation] = []  # location info for each frame
        self._mms: List[mmap.mmap] = []  # memory-mapped file objects for each file
        self._total_frames = 0

        self._open_files()

    @property
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self._total_frames

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        for mm in self._mms:
            if mm is not None:
                mm.close()

    def read_frame(self, index: int) -> "Frame":
        """
        Read a specific frame from the trajectory file(s).

        Args:
            index: Global frame index to read

        Returns:
            The Frame object
        """
        if index < 0:
            index = self._total_frames + index

        if index < 0 or index >= self._total_frames:
            raise IndexError(
                f"Frame index {index} out of range [0, {self._total_frames})"
            )

        # Get location info for this frame
        location = self._get_frame_location(index)

        # Calculate frame end position
        if index + 1 < len(self._frame_locations):
            next_location = self._frame_locations[index + 1]
            if next_location.file_index == location.file_index:
                frame_end = next_location.byte_offset
            else:
                frame_end = None  # End of file
        else:
            frame_end = None  # Last frame

        # Get the memory-mapped file and read frame data
        mm = self._get_mmap(location.file_index)
        frame_bytes = mm[location.byte_offset : frame_end]
        frame_lines = frame_bytes.decode().splitlines()

        # Parse the frame lines using the derived class implementation
        return self._parse_frame(frame_lines)

    def read_frames(self, indices: List[int]) -> List["Frame"]:
        """
        Read multiple frames from the trajectory file.

        Args:
            indices: List of frame indices to read

        Returns:
            List of Frame objects
        """
        return [self.read_frame(i) for i in indices]

    def read_range(self, start: int, stop: int, step: int = 1) -> List["Frame"]:
        """
        Read a range of frames from the trajectory file.

        Args:
            start: Starting frame index
            stop: Stopping frame index (exclusive)
            step: Step size

        Returns:
            List of Frame objects
        """
        indices = list(range(start, stop, step))
        return self.read_frames(indices)

    def read_all(self) -> List["Frame"]:
        """Read all frames from the trajectory file."""
        return [self.read_frame(i) for i in range(self._total_frames)]

    @abstractmethod
    def _parse_frame(self, frame_lines: List[str]) -> "Frame":
        """
        Parse frame lines into a Frame object.

        Args:
            frame_lines: List of strings representing the frame data

        Returns:
            Frame object
        """
        pass

    @abstractmethod
    def _parse_trajectory(self, file_index: int):
        """Parse trajectory file at given index, storing frame locations."""
        pass

    def __len__(self) -> int:
        return self._total_frames

    def __iter__(self) -> Iterator["Frame"]:
        """Iterate over all frames."""
        for i in range(self._total_frames):
            yield self.read_frame(i)

    def _open_files(self):
        """Open trajectory files with memory mapping and build global index."""
        self._mms = []

        for file_index, fpath in enumerate(self.fpaths):
            # Open file
            fp = open(fpath, "rb")
            # Check if empty
            fp.seek(0, 2)
            if fp.tell() == 0:
                raise ValueError(f"File is empty: {fpath}")
            fp.seek(0)  # Seek back to beginning

            mm = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
            self._mms.append(mm)

            # Parse this file to get frame locations
            self._parse_trajectory(file_index)

    def _get_frame_location(self, index: int) -> FrameLocation:
        """Get location information for a frame."""
        if index >= len(self._frame_locations):
            raise IndexError(f"Frame index {index} out of range")
        return self._frame_locations[index]

    def _get_mmap(self, file_index: int) -> mmap.mmap:
        """Get the memory-mapped file object for a specific file."""
        if file_index >= len(self._mms) or self._mms[file_index] is None:
            raise ValueError(f"File {file_index} is not properly opened")
        return self._mms[file_index]

    @property
    def fpath(self) -> Path:
        """For backward compatibility - returns the first file path."""
        if not self.fpaths:
            raise ValueError("No files available")
        return self.fpaths[0]


class TrajectoryWriter(ABC):
    """Base class for all chemical file writers."""

    def __init__(self, fpath: Union[str, Path]):
        self.fpath = Path(fpath)
        self._fp = open(self.fpath, "w+b")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    @abstractmethod
    def write_frame(self, frame: "Frame"):
        """Write a single frame to the file."""
        pass

    def close(self):
        self._fp.close()
