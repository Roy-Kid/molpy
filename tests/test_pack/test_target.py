"""Target: a frame replicated `number` times inside a constraint."""

import numpy as np

import molrs

from molpy.pack.constraint import MinDistanceConstraint
from molpy.pack.target import Target


def _two_atoms() -> molrs.Frame:
    frame = molrs.Frame()
    frame["atoms"] = {
        "id": np.array([1, 2]),
        "x": np.array([0.0, 1.0]),
        "y": np.array([0.0, 0.0]),
        "z": np.array([0.0, 0.5]),
    }
    return frame


def test_n_points_counts_every_copy():
    target = Target(_two_atoms(), number=3, constraint=MinDistanceConstraint(2.0))
    assert target.n_points == 6


def test_points_tile_the_coordinates_per_copy():
    target = Target(_two_atoms(), number=2, constraint=MinDistanceConstraint(2.0))
    points = target.points
    assert points.shape == (4, 3)
    np.testing.assert_array_equal(points[:2], points[2:])
    np.testing.assert_array_equal(points[1], [1.0, 0.0, 0.5])


def test_empty_frame_has_no_points():
    frame = molrs.Frame()
    frame["atoms"] = molrs.Block()
    target = Target(frame, number=5, constraint=MinDistanceConstraint(2.0), name="void")
    assert target.n_points == 0
    assert target.points.shape == (0, 3)
    assert "void" in repr(target)
