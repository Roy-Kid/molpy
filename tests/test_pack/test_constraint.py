"""Packing constraints: the min-distance penalty and its gradient."""

import numpy as np
import pytest

from molpy.pack.constraint import MinDistanceConstraint


def _finite_difference(constraint, points, h=1e-6):
    grad = np.zeros_like(points)
    for index in np.ndindex(points.shape):
        plus = points.copy()
        plus[index] += h
        minus = points.copy()
        minus[index] -= h
        grad[index] = (constraint.penalty(plus) - constraint.penalty(minus)) / (2 * h)
    return grad


class TestMinDistanceConstraint:
    def test_penalty_is_the_squared_shortfall_of_each_close_pair(self):
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [10.0, 0.0, 0.0]])
        # One pair at distance 1 under dmin = 2.5: (2.5 - 1)^2; the far point is free.
        assert MinDistanceConstraint(2.5).penalty(points) == pytest.approx(2.25)

    def test_penalty_is_zero_when_every_pair_is_far_enough(self):
        points = np.array([[0.0, 0.0, 0.0], [3.0, 0.0, 0.0], [0.0, 3.0, 0.0]])
        assert MinDistanceConstraint(2.5).penalty(points) == 0.0
        assert MinDistanceConstraint(2.5).penalty(points[:1]) == 0.0

    def test_gradient_matches_finite_differences(self):
        rng = np.random.default_rng(3)
        points = rng.uniform(0.0, 4.0, size=(12, 3))
        constraint = MinDistanceConstraint(2.0)
        np.testing.assert_allclose(
            constraint.dpenalty(points),
            _finite_difference(constraint, points),
            atol=1e-5,
        )

    def test_gradient_points_toward_the_neighbour_and_sums_to_zero(self):
        # Descending the gradient separates the pair: d(penalty)/dx0 > 0 for the
        # point on the left, so a step against it moves it further left.
        points = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        grad = MinDistanceConstraint(2.5).dpenalty(points)
        assert grad[0, 0] > 0 > grad[1, 0]
        np.testing.assert_allclose(grad.sum(axis=0), 0.0, atol=1e-12)
