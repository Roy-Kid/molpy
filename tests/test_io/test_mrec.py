"""Unit tests for :mod:`molpy.io.mrec` scientific-record I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from molpy import Block, Frame, Record, Trajectory
from molpy.io.mrec import (
    TrajectoryReader,
    read_record,
    write_record,
    write_trajectory,
)

# Å; dyadic so a bit-exact f64 round-trip is the golden, not a tolerance.
_N_ATOMS = 3
_ATOM_X = (0.0, 1.0, 0.5)
_ATOM_Y = (0.25, 0.0, 2.0)
_ATOM_Z = (0.0, 4.0, 0.125)


def _coords_frame() -> Frame:
    atoms = Block()
    atoms["x"] = np.array(_ATOM_X, dtype=np.float64)
    atoms["y"] = np.array(_ATOM_Y, dtype=np.float64)
    atoms["z"] = np.array(_ATOM_Z, dtype=np.float64)
    frame = Frame()
    frame["atoms"] = atoms
    return frame


def _assert_coords(frame: Frame) -> None:
    atoms = frame["atoms"]
    assert atoms.nrows == _N_ATOMS
    np.testing.assert_array_equal(
        np.asarray(atoms["x"]), np.array(_ATOM_X, dtype=np.float64)
    )
    np.testing.assert_array_equal(
        np.asarray(atoms["y"]), np.array(_ATOM_Y, dtype=np.float64)
    )
    np.testing.assert_array_equal(
        np.asarray(atoms["z"]), np.array(_ATOM_Z, dtype=np.float64)
    )


class TestTrajectoryReader:
    def test_read_frame(self, tmp_path: Path) -> None:
        path = tmp_path / "traj.mrec"
        write_trajectory(str(path), Trajectory([_coords_frame()]))
        reader = TrajectoryReader(str(path))
        _assert_coords(reader.read_frame(0))


class TestWriteRecord:
    def test_round_trips_system_coordinates(self, tmp_path: Path) -> None:
        path = tmp_path / "record.mrec"
        record = Record()
        record.set_system(_coords_frame())
        write_record(str(path), record)

        loaded = read_record(str(path))
        assert loaded.system is not None
        _assert_coords(loaded.system)


class TestReadRecord:
    def test_stamps_format_name_mrec(self, tmp_path: Path) -> None:
        path = tmp_path / "record.mrec"
        record = Record()
        record.set_system(_coords_frame())
        write_record(str(path), record)

        loaded = read_record(str(path))
        assert loaded.meta["format_name"] == "mrec"


class TestWriteTrajectory:
    def test_round_trips_coordinates(self, tmp_path: Path) -> None:
        path = tmp_path / "traj.mrec"
        write_trajectory(str(path), Trajectory([_coords_frame()]))

        loaded = read_record(str(path))
        assert loaded.trajectory is not None
        assert len(loaded.trajectory) == 1
        _assert_coords(loaded.trajectory.frames[0])


class TestMrecSurface:
    def test_import_zarr_raises_import_error(self) -> None:
        with pytest.raises(ImportError):
            import molpy.io.zarr  # noqa: F401

    def test_has_no_frame_reader(self) -> None:
        import molpy.io.mrec as mrec

        assert not hasattr(mrec, "FrameReader")
