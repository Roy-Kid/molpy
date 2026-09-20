"""Scientific-record I/O for ``*.mrec`` stores.

Path doors for in-memory :class:`~molpy.Frame` and :class:`~molpy.Trajectory`:

* :func:`read_frame` / :func:`write_frame` — Structure (``meta`` + ``frame/``)
* :func:`read_system` / :func:`write_system` — System-def (``meta`` + ``system/``)
* :func:`read_trajectory` / :func:`write_trajectory` — Trajectory shape
* :class:`TrajectoryReader` — lazy frame cursor (``len``, ``reader[i]``,
  iteration, ``.step`` / ``.time`` labels, ``has_block``)
* :class:`SequenceSchema` / :class:`TrajectoryWriter` — pin a schema and write
  a run frame by frame, without holding it all in memory
* the native ``mrec.schema`` re-exported as :mod:`molpy.io.mrec.schema`

There is no :class:`molpy.Record` and no :mod:`molpy.io.zarr`.
"""

from __future__ import annotations

from pathlib import Path

from molrs import Frame, Trajectory
from molrs.io.mrec import SequenceSchema, TrajectoryWriter
from molrs.io.mrec import TrajectoryReader as _MolrsTrajectoryReader
from molrs.io.mrec import read_frame as _read_frame
from molrs.io.mrec import read_meta as _read_meta
from molrs.io.mrec import read_system as _read_system
from molrs.io.mrec import read_trajectory as _read_trajectory
from molrs.io.mrec import schema
from molrs.io.mrec import sections as _sections
from molrs.io.mrec import write_frame as _write_frame
from molrs.io.mrec import write_system as _write_system
from molrs.io.mrec import write_trajectory as _write_trajectory

# `SequenceSchema` / `TrajectoryWriter` need no molpy-typed conversion, so they
# are re-exported from molrs verbatim (the path handling already lives there).
# Only `TrajectoryReader` earns a subclass — for the `_store_path` door below.


def _store_path(path: str | Path) -> str:
    return str(Path(path).expanduser())


class TrajectoryReader(_MolrsTrajectoryReader):
    """Open a lazy one-frame cursor over a ``*.mrec`` trajectory.

    Construction opens the store index only; each read decodes exactly the
    frame asked for. Beyond :meth:`read_frame` the cursor supports ``len()``,
    ``reader[i]`` (negative indices included), iteration, the ``.step`` /
    ``.time`` frame labels, and :meth:`has_block` — all inherited from the
    the native cursor.

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


def read_frame(path: str | Path) -> Frame:
    """Read the ``frame`` section of a ``*.mrec`` store."""
    return _read_frame(_store_path(path))


def write_frame(
    path: str | Path,
    frame: Frame,
    system: Frame | None = None,
    meta: dict | None = None,
) -> None:
    """Write a snapshot as a record whose only state section is ``frame``."""
    _write_frame(_store_path(path), frame, system, meta)


def read_system(path: str | Path) -> Frame:
    """Read the ``system`` section of a ``*.mrec`` store."""
    return _read_system(_store_path(path))


def write_system(path: str | Path, system: Frame, meta: dict | None = None) -> None:
    """Write a topology as a record whose only state section is ``system``."""
    _write_system(_store_path(path), system, meta)


def read_trajectory(path: str | Path) -> Trajectory:
    """Read the ``trajectory`` section of a ``*.mrec`` store."""
    return _read_trajectory(_store_path(path))


def write_trajectory(path: str | Path, trajectory: Trajectory) -> None:
    """Write a trajectory as a record whose only state section is ``trajectory``."""
    _write_trajectory(_store_path(path), trajectory)


def read_meta(path: str | Path) -> dict:
    """Read the mandatory ``meta`` document of a ``*.mrec`` store."""
    return _read_meta(_store_path(path))


def sections(path: str | Path) -> frozenset[str]:
    """Child group names at the record root."""
    return _sections(_store_path(path))


__all__ = [
    "SequenceSchema",
    "TrajectoryReader",
    "TrajectoryWriter",
    "read_frame",
    "read_meta",
    "read_system",
    "read_trajectory",
    "schema",
    "sections",
    "write_frame",
    "write_system",
    "write_trajectory",
]
