from abc import ABC, abstractmethod
from typing import Iterator, Union, List, Optional, TYPE_CHECKING
from pathlib import Path
import mmap

if TYPE_CHECKING:
    from ...core.frame import Frame

PathLike = Union[str, bytes]  # type_check_only

class TrajectoryReader(ABC):
    """
    Base class for trajectory file readers that act as providers.
    
    This class provides memory-mapped file reading and directly returns Frame objects
    without needing to interact with Trajectory objects.
    """

    def __init__(self, fpath: Union[Path, str]):
        """
        Initialize the trajectory reader.
        
        Args:
            fpath: Path to trajectory file
        """
        self.fpath = Path(fpath)
        if not self.fpath.exists():
            raise FileNotFoundError(f"File not found: {self.fpath}")

        self._byte_offsets: List[int] = []  # list of byte offsets for each frame
        self._mm = None  # memory-mapped file object
        self._total_frames = 0

        self._open_file()

    @property
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self._total_frames

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._mm is not None:
            self._mm.close()

    def read_frame(self, index: int) -> "Frame":
        """
        Read a specific frame from the trajectory file.
        
        Args:
            index: Frame index to read
            
        Returns:
            The Frame object
        """
        if index < 0:
            index = self._total_frames + index
            
        if index < 0 or index >= self._total_frames:
            raise IndexError(f"Frame index {index} out of range [0, {self._total_frames})")
            
        # Read the frame directly
        frame = self._read_frame_data(index)
        return frame

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
    def _read_frame_data(self, index: int) -> "Frame":
        """
        Read frame data from file at the given index.
        
        Args:
            index: Frame index to read
            
        Returns:
            Frame object
        """
        pass

    def __len__(self) -> int:
        return self._total_frames

    def __iter__(self) -> Iterator["Frame"]:
        """Iterate over all frames."""
        for i in range(self._total_frames):
            yield self.read_frame(i)

    @abstractmethod
    def _parse_trajectory(self):
        """Parse trajectory file, storing frame offsets."""
        pass

    def _open_file(self):
        """Open trajectory file with memory mapping."""
        fp = open(self.fpath, "rb")
        # if empty, raise error
        fp.seek(0, 2)
        if fp.tell() == 0:
            raise ValueError("File is empty")
        fp.seek(0)  # Seek back to beginning
        self._mm = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
        self._parse_trajectory()

    def get_offset(self, index: int) -> int:
        """Get byte offset for a given frame index."""
        if index >= len(self._byte_offsets):
            raise IndexError(f"Frame index {index} out of range")
        return self._byte_offsets[index]

    def get_mmap(self) -> mmap.mmap:
        """Get the memory-mapped file object."""
        if self._mm is None:
            raise ValueError("File is empty or not properly opened")
        return self._mm


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

