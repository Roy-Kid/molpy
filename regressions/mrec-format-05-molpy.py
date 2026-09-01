r"""Write and read a ``*.mrec`` store through ``molpy.io.mrec``.

Scientific-record I/O belongs on ``molpy.io.mrec``. This script writes a
one-system store to a directory whose name ends in ``.mrec``, reads it back
with ``read_system``, and asserts ``molrec_version`` as the literal ``1``.
There is no ``molpy.Record``.

Provenance of the goldens: hand-written literals, no external oracle and no
third-party scientific package at run time (``molpy`` + ``numpy`` only).
Runner:

    python regressions/mrec-format-05-molpy.py

(2026-08-31).
"""

from __future__ import annotations

import tempfile
import warnings
from pathlib import Path

import numpy as np

warnings.filterwarnings("ignore", category=FutureWarning)

import molpy as mp
from molpy.io import mrec

ATOM_X = (0.0, 1.0, 0.5)
ATOM_Y = (0.25, 0.0, 2.0)
ATOM_Z = (0.0, 4.0, 0.125)
N_ATOMS = 3

with tempfile.TemporaryDirectory() as tmp:
    store = Path(tmp) / "record.mrec"
    assert store.suffix == ".mrec", store

    atoms = mp.Block()
    atoms["x"] = np.array(ATOM_X, dtype=np.float64)
    atoms["y"] = np.array(ATOM_Y, dtype=np.float64)
    atoms["z"] = np.array(ATOM_Z, dtype=np.float64)
    system = mp.Frame()
    system["atoms"] = atoms

    mrec.write_system(store, system)
    loaded = mrec.read_system(store)
    got = loaded["atoms"]
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

    meta = mrec.read_meta(store)
    mrec.schema.validate_meta(meta)
    assert meta["molrec_version"] == mrec.schema.MOLREC_VERSION
    assert "format_name" not in meta, meta
    assert store.is_dir(), store
    assert not hasattr(mp, "Record")

    loud = Path(tmp) / "loud.mrec"
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
    assert not loud.exists(), loud

print("mrec-format-05-molpy ok: write_system/read_system record.mrec molrec_version=1")
