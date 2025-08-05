from abc import ABC, abstractmethod
from typing import Iterator, Union, List, Optional, TYPE_CHECKING
from pathlib import Path
import mmap

if TYPE_CHECKING:
    from ...core.trajectory import Trajectory
    from ...core.frame import Frame

PathLike = Union[str, bytes]  # type_check_only

class TrajectoryReader(ABC):
    """
    Base class for trajectory file readers with lazy loading and caching support.
    
    This class provides memory-mapped file reading and works with Trajectory objects
    to enable on-demand frame loading and caching.
    """

    def __init__(self, trajectory: "Trajectory", fpath: Union[Path, str]):
        """
        Initialize the trajectory reader.
        
        Args:
            trajectory: Trajectory object to populate with frames
            fpath: Path to trajectory file
        """
        self.trajectory = trajectory
        self.fpath = Path(fpath)
        if not self.fpath.exists():
            raise FileNotFoundError(f"File not found: {self.fpath}")

        self._byte_offsets: List[int] = []  # list of byte offsets for each frame
        self._mm = None  # memory-mapped file object
        self._total_frames = 0

        self._open_file()
        self.trajectory.set_total_frames(self._total_frames)

    @property
    def n_frames(self) -> int:
        """Number of frames in the trajectory."""
        return self._total_frames

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._mm is not None:
            self._mm.close()

    def load_frame(self, index: int) -> "Frame":
        """
        Load a specific frame into the trajectory.
        
        Args:
            index: Frame index to load
            
        Returns:
            The loaded Frame object
        """
        if index < 0:
            index = self._total_frames + index
            
        if index < 0 or index >= self._total_frames:
            raise IndexError(f"Frame index {index} out of range [0, {self._total_frames})")
            
        # Check if frame is already loaded
        if self.trajectory.is_loaded(index):
            frame = self.trajectory.frames[index]  # Direct access to frame dict
            return frame
            
        # Load the frame
        frame = self.read_frame(index)
        self.trajectory._add_frame(index, frame)
        return frame

    def load_frames(self, indices: List[int]) -> List["Frame"]:
        """
        Load multiple frames into the trajectory.
        
        Args:
            indices: List of frame indices to load
            
        Returns:
            List of loaded Frame objects
        """
        return [self.load_frame(i) for i in indices]

    def load_range(self, start: int, stop: int, step: int = 1) -> List["Frame"]:
        """
        Load a range of frames into the trajectory.
        
        Args:
            start: Starting frame index
            stop: Stopping frame index (exclusive)
            step: Step size
            
        Returns:
            List of loaded Frame objects
        """
        indices = list(range(start, stop, step))
        return self.load_frames(indices)

    def preload_all(self) -> None:
        """Load all frames into the trajectory."""
        for i in range(self._total_frames):
            if not self.trajectory.is_loaded(i):
                self.load_frame(i)

    @abstractmethod
    def read_frame(self, index: int) -> "Frame":
        """
        Read a frame from file without adding it to the trajectory.
        
        Args:
            index: Frame index to read
            
        Returns:
            Frame object
        """
        pass

    def read_frames(self, indices: List[int]) -> List["Frame"]:
        """Read multiple frames from file without adding them to the trajectory."""
        return [self.read_frame(i) for i in indices]

    def __len__(self) -> int:
        return self._total_frames

    def __iter__(self) -> Iterator["Frame"]:
        """Iterate over all frames, loading them as needed."""
        for i in range(self._total_frames):
            yield self.load_frame(i)

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

