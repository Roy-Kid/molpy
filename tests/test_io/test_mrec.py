"""Unit tests for :mod:`molpy.io.mrec` scientific-record I/O."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from molpy import Block, Frame, Trajectory
from molpy.io.mrec import (
    TrajectoryReader,
    read_frame,
    read_meta,
    read_system,
    read_trajectory,
    schema,
    sections,
    write_frame,
    write_system,
    write_trajectory,
)


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
        write_trajectory(path, Trajectory([_coords_frame()]))
        reader = TrajectoryReader(path)
        _assert_coords(reader.read_frame(0))


class TestWriteFrame:
    def test_round_trips_coordinates(self, tmp_path: Path) -> None:
        path = tmp_path / "snapshot.mrec"
        write_frame(path, _coords_frame())
        _assert_coords(read_frame(path))
        assert sections(path) == frozenset({"meta", "frame"})
        meta = read_meta(path)
        schema.validate_meta(meta)
        assert meta["molrec_version"] == schema.MOLREC_VERSION
        assert "format_name" not in meta


class TestWriteSystem:
    def test_round_trips_coordinates(self, tmp_path: Path) -> None:
        path = tmp_path / "system.mrec"
        write_system(path, _coords_frame())
        _assert_coords(read_system(path))
        assert "frame" not in sections(path)


class TestWriteTrajectory:
    def test_round_trips_coordinates(self, tmp_path: Path) -> None:
        path = tmp_path / "traj.mrec"
        write_trajectory(path, Trajectory([_coords_frame()]))
        loaded = read_trajectory(path)
        assert len(loaded) == 1
        _assert_coords(loaded[0])


class TestSchema:
    def test_sole_version_key_is_molrec_version(self) -> None:
        assert schema.MOLREC_VERSION == 1
        with pytest.raises(Exception, match="molrec_version"):
            schema.validate_meta({"record_schema_version": 1, "format_name": "mrec"})


class TestMrecSurface:
    def test_import_zarr_raises_import_error(self) -> None:
        with pytest.raises(ImportError):
            import molpy.io.zarr  # noqa: F401

    def test_has_no_frame_reader(self) -> None:
        import molpy.io.mrec as mrec

        assert not hasattr(mrec, "FrameReader")

    def test_has_no_record(self) -> None:
        import molpy
        import molpy.io.mrec as mrec

        assert not hasattr(molpy, "Record")
        assert not hasattr(mrec, "read_record")
        assert not hasattr(mrec, "write_record")
