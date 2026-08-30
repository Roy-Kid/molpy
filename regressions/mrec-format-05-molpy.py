r"""Write and read a ``*.mrec`` record through ``molpy.io.mrec``.

Scientific-record I/O belongs on ``molpy.io.mrec`` (spec mrec-format-05-molpy).
This script writes a one-system record to a directory whose name ends in
``.mrec``, reads it back with ``read_record``, and asserts ``format_name`` as
the literal ``"mrec"``. ``Trajectory.read`` / ``Trajectory.write`` are called
only to show they fail loud naming ``molpy.io.mrec``.

Provenance of the goldens: hand-written literals, no external oracle and no
third-party scientific package at run time (``molpy`` + ``numpy`` only).
Runner:

    python regressions/mrec-format-05-molpy.py

(2026-08-30).
"""

from __future__ import annotations

import os
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore", category=FutureWarning)

import molpy as mp
from molpy.io import mrec

# Å; dyadic so a bit-exact f64 round-trip is the golden, not a tolerance.
ATOM_X = (0.0, 1.0, 0.5)
ATOM_Y = (0.25, 0.0, 2.0)
ATOM_Z = (0.0, 4.0, 0.125)
N_ATOMS = 3

with tempfile.TemporaryDirectory() as tmp:
    store = os.path.join(tmp, "record.mrec")
    assert store.endswith(".mrec"), store

    atoms = mp.Block()
    atoms["x"] = np.array(ATOM_X, dtype=np.float64)
    atoms["y"] = np.array(ATOM_Y, dtype=np.float64)
    atoms["z"] = np.array(ATOM_Z, dtype=np.float64)
    system = mp.Frame()
    system["atoms"] = atoms

    record = mp.Record()
    record.set_system(system)
    record.meta = {"creator": {"name": "mrec-format-05-molpy"}}
    mrec.write_record(store, record)

    loaded = mrec.read_record(store)
    meta = loaded.meta
    assert meta["format_name"] == "mrec", meta.get("format_name")
    assert loaded.system is not None, "system section missing after read_record"
    got = loaded.system["atoms"]
    assert got.nrows == N_ATOMS, got.nrows
    np.testing.assert_array_equal(
        np.asarray(got["x"]), np.array(ATOM_X, dtype=np.float64)
    )
    np.testing.assert_array_equal(
        np.asarray(got["y"]), np.array(ATOM_Y, dtype=np.float64)
    )
    np.testing.assert_array_equal(
        np.asarray(got["z"]), np.array(ATOM_Z, dtype=np.float64)
    )
    assert loaded.meta["creator"]["name"] == "mrec-format-05-molpy"
    assert os.path.isdir(store), store

    loud = os.path.join(tmp, "loud.mrec")
    traj = mp.Trajectory([system])
    try:
        traj.read(loud)
    except (TypeError, RuntimeError) as exc:
        assert "molpy.io.mrec" in str(exc), exc
    else:
        raise AssertionError("Trajectory.read must fail loud naming molpy.io.mrec")
    try:
        traj.write(loud)
    except (TypeError, RuntimeError) as exc:
        assert "molpy.io.mrec" in str(exc), exc
    else:
        raise AssertionError("Trajectory.write must fail loud naming molpy.io.mrec")
    assert not os.path.exists(loud), loud

print("mrec-format-05-molpy ok: write_record/read_record record.mrec format_name=mrec")
