"""LAMMPS dump trajectory write — native-backed.

Incremental :meth:`write_frame` buffers frames; :meth:`close` flushes via
the native ``write_lammps_traj`` / the native ``write_lammps_dump_local``
(requires each frame to carry ``box``).
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import molrs.io
from molrs import Frame

from .base import TrajectoryWriter


class LammpsTrajectoryWriter(TrajectoryWriter):
    """Write a LAMMPS dump custom / atom trajectory (native)."""

    def __init__(self, fpath: str | Path, atom_style: str = "full") -> None:
        super().__init__(fpath)
        self.atom_style = atom_style
        self._frames: list[Frame] = []
        if self._fp is not None:
            self._fp.close()
            self._fp = None

    def write_frame(self, frame: Frame, timestep: int | None = None) -> None:
        if frame.box is None:
            raise ValueError(
                "LAMMPS trajectory write requires frame.box (the native core needs a simbox)"
            )
        if timestep is not None:
            frame.meta["timestep"] = int(timestep)
        elif "timestep" not in frame.meta:
            frame.meta["timestep"] = len(self._frames)
        self._frames.append(frame)

    def close(self) -> None:
        if self._frames:
            molrs.io.write_lammps_traj(str(self.fpath), self._frames)
            self._frames = []
        if self._fp is not None:
            self._fp.close()
            self._fp = None

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()


class LammpsDumpLocalWriter(TrajectoryWriter):
    """Write LAMMPS dump local (OVITO Load trajectory bonds) natively."""

    def __init__(self, fpath: str | Path) -> None:
        super().__init__(fpath)
        self._frames: list[Frame] = []
        if self._fp is not None:
            self._fp.close()
            self._fp = None

    def write_frame(self, frame: Frame, timestep: int | None = None) -> None:
        if frame.box is None:
            raise ValueError(
                "LAMMPS dump local write requires frame.box (the native core needs a simbox)"
            )
        if timestep is not None:
            frame.meta["timestep"] = int(timestep)
        elif "timestep" not in frame.meta:
            frame.meta["timestep"] = len(self._frames)
        self._frames.append(frame)

    def close(self) -> None:
        if self._frames:
            molrs.io.write_lammps_dump_local(str(self.fpath), self._frames)
            self._frames = []
        if self._fp is not None:
            self._fp.close()
            self._fp = None

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        self.close()
