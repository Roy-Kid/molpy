from abc import ABC, abstractmethod
from collections.abc import Iterator
from io import BytesIO, StringIO, TextIOWrapper
from pathlib import Path
from typing import IO

from molpy.core.frame import Frame

PathLike = str | Path
FileLike = PathLike | IO | BytesIO | StringIO


# ─────────────────────────────────────────────────────────────────────
# Shared context-manager boilerplate
# ─────────────────────────────────────────────────────────────────────
class FileBase(ABC):
    """Common logic for Context-manager + lazy file handle."""

    def __init__(self, path: PathLike, mode: str, **open_kwargs):
        self._path = Path(path)
        self._mode = mode
        self._open_kwargs = open_kwargs
        self._fh: IO[str] | None = None

    # ---------- context-manager hooks ---------------------------------
    def __enter__(self):
        self._fh = self._path.open(self._mode, **self._open_kwargs)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fh:
            self._fh.close()
            self._fh = None

    # ---------- lazy accessor (works with or without `with`) ----------
    @property
    def fh(self) -> IO[str]:
        if self._fh is None:
            self._fh = self._path.open(self._mode, **self._open_kwargs)
        return self._fh


# ─────────────────────────────────────────────────────────────────────
# DataReader
# ─────────────────────────────────────────────────────────────────────
class DataReader(FileBase, ABC):
    """Base class for data file readers.

    Supports reading from both file paths and file-like objects (BytesIO, StringIO, etc.).
    """

    def __init__(self, source: FileLike, **open_kwargs):
        """Initialize the DataReader.

        Args:
            source: Either a file path (str/Path) or a file-like object (BytesIO, StringIO, etc.)
            **open_kwargs: Additional keyword arguments for file opening (only used for paths)
        """
        # Check if source is a file-like object
        if hasattr(source, "read"):
            # It's a file-like object
            self._is_file_object = True
            self._path = None
            self._mode = "r"
            self._open_kwargs = {}

            # Wrap BytesIO in TextIOWrapper if needed
            if isinstance(source, BytesIO):
                self._fh = TextIOWrapper(source, encoding="utf-8")
                self._owns_wrapper = True
            else:
                self._fh = source
                self._owns_wrapper = False
        else:
            # It's a path
            self._is_file_object = False
            self._path = Path(source)
            self._mode = "r"
            self._open_kwargs = open_kwargs
            self._fh = None
            self._owns_wrapper = False

    # Override parent methods to handle file-like objects
    def __enter__(self):
        if self._is_file_object:
            return self
        else:
            self._fh = self._path.open(self._mode, **self._open_kwargs)
            return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fh and not self._is_file_object:
            self._fh.close()
            self._fh = None
        elif self._owns_wrapper and self._fh:
            # Close the wrapper but not the underlying BytesIO
            self._fh.detach()
            self._fh = None

    @property
    def fh(self) -> IO[str]:
        if self._fh is None and not self._is_file_object:
            self._fh = self._path.open(self._mode, **self._open_kwargs)
        return self._fh

    # -- line helpers --------------------------------------------------
    def _iter_nonblank(self) -> Iterator[str]:
        """Iterate over non-blank, stripped lines."""
        self.fh.seek(0)
        for raw in self.fh:
            line = raw.strip()
            if line:
                yield line

    def __iter__(self) -> Iterator[str]:
        """`for line in reader:` yields non-blank, stripped lines."""
        return self._iter_nonblank()

    def read_lines(self) -> list[str]:
        """Return all lines at once."""
        return list(self.fh.readlines())

    # -- high-level parse ---------------------------------------------
    @abstractmethod
    def read(self, frame: Frame | None = None) -> Frame:
        """
        Populate / update a Frame from the underlying file.

        Args:
            frame: Optional existing Frame to populate. If None, creates a new one.

        Returns:
            The populated Frame object
        """
        ...


# ─────────────────────────────────────────────────────────────────────
# DataWriter
# ─────────────────────────────────────────────────────────────────────
class DataWriter(FileBase, ABC):
    """Base class for data file writers."""

    def __init__(self, path: PathLike, **open_kwargs):
        super().__init__(path, mode="w", **open_kwargs)

    @abstractmethod
    def write(self, frame: Frame) -> None:
        """
        Serialize frame into the underlying file.

        Args:
            frame: Frame object to write
        """
        ...
