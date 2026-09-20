"""LAMMPSEngine helpers: style lines from a force field, coordinate splicing."""

import numpy as np
import pytest

import molrs

import molpy as mp
from molpy.engine.lammps import _splice_coords, _style_lines


def test_style_lines_name_the_bonded_styles_only():
    ff = mp.ForceField(name="t", units="real")
    atoms = ff.def_atomstyle("full")
    ct = atoms.def_type("CT", mass=12.011)
    ff.def_bondstyle("harmonic").def_type(ct, ct, k=1.0, r0=1.5)
    ff.def_pairstyle("lj/cut", cutoff=10.0).def_type(ct, ct, epsilon=0.1, sigma=3.0)
    assert _style_lines(ff) == ["bond_style harmonic"]


def _frame(ids, x):
    frame = molrs.Frame()
    frame["atoms"] = {
        "id": np.array(ids),
        "x": np.array(x, dtype=float),
        "y": np.zeros(len(ids)),
        "z": np.zeros(len(ids)),
        "type": np.array(["A"] * len(ids)),
    }
    frame.box = mp.Box.cubic(10.0)
    return frame


def test_splice_matches_relaxed_coordinates_by_id():
    original = _frame([1, 2], [0.0, 1.0])
    relaxed = _frame([2, 1], [7.0, 5.0])  # reversed order, moved atoms
    out = _splice_coords(original, relaxed)
    assert out["atoms"]["x"].tolist() == [5.0, 7.0]
    assert list(out["atoms"]["type"]) == ["A", "A"]
    assert original["atoms"]["x"].tolist() == [0.0, 1.0], "input is not mutated"
    np.testing.assert_allclose(out.box.matrix, original.box.matrix)


def test_splice_rejects_a_changed_atom_count():
    with pytest.raises(RuntimeError, match="atom count changed"):
        _splice_coords(_frame([1, 2], [0.0, 1.0]), _frame([1], [0.0]))
