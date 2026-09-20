"""GROMACS .gro file I/O — thin native wrappers.

Parse/serialize live in the native ``io`` module, which also owns
:class:`GroFieldFormatter` (the .gro → canonical name map); there is no separate
Python parser and no molpy copy of that formatter.
"""

from __future__ import annotations

from pathlib import Path

import molrs.io
from molrs import Frame

from molpy.core.fields import GroFieldFormatter

from .base import DataReader, DataWriter


class GroReader(DataReader):
    """Read GRO natively (first frame if multi-frame)."""

    _formatter = GroFieldFormatter()

    def __init__(self, path: str | Path, **kwargs: object) -> None:
        super().__init__(Path(path), **kwargs)

    def read(self, frame: Frame | None = None) -> Frame:
        del frame
        frames = molrs.io.read_gro(str(self._path))
        if not frames:
            raise OSError(f"no frames parsed from GRO file: {self._path}")
        return frames[0]


class GroWriter(DataWriter):
    """Write GRO natively."""

    _formatter = GroFieldFormatter()

    def __init__(self, path: str | Path, **kwargs: object) -> None:
        super().__init__(Path(path), **kwargs)

    def write(self, frame: Frame) -> None:
        molrs.io.write_gro(str(self._path), frame)
