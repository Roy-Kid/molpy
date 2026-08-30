"""Scientific-record I/O for ``*.mrec`` stores.

Path doors for in-memory :class:`~molpy.Record` and :class:`~molpy.Trajectory`:

* :func:`read_record` / :func:`write_record` — whole record
* :func:`write_trajectory` — trajectory-only record
* :class:`TrajectoryReader` — lazy frame cursor (one frame per
  :meth:`TrajectoryReader.read_frame`)

There is no :class:`FrameReader` and no :mod:`molpy.io.zarr`.
"""

from __future__ import annotations

from pathlib import Path

from molrs import Frame, Record, Trajectory
from molrs.io.mrec import TrajectoryReader as _MolrsTrajectoryReader
from molrs.io.mrec import read_record as _read_record
from molrs.io.mrec import write_record as _write_record
from molrs.io.mrec import write_trajectory as _write_trajectory


def _store_path(path: str | Path) -> str:
    return str(Path(path).expanduser())


class TrajectoryReader(_MolrsTrajectoryReader):
    """Open a lazy one-frame cursor over a ``*.mrec`` trajectory.

    Construction opens the store index. :meth:`read_frame` decodes exactly
    the asked-for frame.

    Args:
        path: Filesystem path of the record store.
    """

    def __init__(self, path: str | Path) -> None:
        super().__init__(_store_path(path))

    def read_frame(self, index: int) -> Frame:
        """Decode one committed frame.

        Args:
            index: Zero-based frame index.

        Returns:
            The frame at ``index``.

        Raises:
            IndexError: If ``index`` is past the commit marker.
            ValueError: If a section fails to decode.
        """
        return super().read_frame(index)


def read_record(path: str | Path) -> Record:
    """Read a scientific record from a ``*.mrec`` store.

    Args:
        path: Filesystem path of the record store.

    Returns:
        The in-memory :class:`~molpy.Record`.

    Raises:
        ValueError: If ``path`` uses a retired ``.zarr`` suffix, the store is
            not a readable record, or a section fails to decode.
    """
    return _read_record(_store_path(path))


def write_record(path: str | Path, record: Record) -> None:
    """Write a scientific record to a ``*.mrec`` store.

    Args:
        path: Destination filesystem path.
        record: In-memory :class:`~molpy.Record` to persist.

    Raises:
        ValueError: If ``path`` uses a retired ``.zarr`` suffix, the record
            has no state section, or a section fails to encode.
    """
    _write_record(_store_path(path), record)


def write_trajectory(path: str | Path, trajectory: Trajectory) -> None:
    """Write a trajectory as a record whose only state section is ``trajectory``.

    Args:
        path: Destination filesystem path.
        trajectory: In-memory :class:`~molpy.Trajectory` to persist.

    Raises:
        ValueError: If ``path`` uses a retired ``.zarr`` suffix, or a frame
            fails to encode.
    """
    _write_trajectory(_store_path(path), trajectory)


__all__ = ["TrajectoryReader", "read_record", "write_record", "write_trajectory"]
