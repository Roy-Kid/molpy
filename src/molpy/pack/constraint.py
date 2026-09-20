import numpy as np

import molrs

from molpy.core.region import BoxRegion, SphereRegion


# === Base Constraint class ===
class Constraint:
    """Base class for all packing constraints."""

    def penalty(self, points: np.ndarray) -> float:
        """Calculate penalty for given points. Lower is better."""
        raise NotImplementedError

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        """Calculate gradient of penalty with respect to points."""
        raise NotImplementedError

    def __and__(self, other: "Constraint") -> "AndConstraint":
        """Combine constraints with AND (both must be satisfied)."""
        return AndConstraint(self, other)

    def __or__(self, other: "Constraint") -> "OrConstraint":
        """Combine constraints with OR (either can be satisfied)."""
        return OrConstraint(self, other)


class AndConstraint(Constraint):
    def __init__(self, a: Constraint, b: Constraint):
        self.a = a
        self.b = b

    def penalty(self, points):
        return self.a.penalty(points) + self.b.penalty(points)

    def dpenalty(self, points):
        return self.a.dpenalty(points) + self.b.dpenalty(points)


class OrConstraint(Constraint):
    def __init__(self, a: Constraint, b: Constraint):
        self.a = a
        self.b = b

    def penalty(self, points):
        pa = self.a.penalty(points)
        pb = self.b.penalty(points)
        return min(pa, pb)

    def dpenalty(self, points):
        pa = self.a.penalty(points)
        pb = self.b.penalty(points)
        if pa < pb:
            return self.a.dpenalty(points)
        elif pb < pa:
            return self.b.dpenalty(points)
        else:
            # When penalties are equal, use average gradient for stability
            return (self.a.dpenalty(points) + self.b.dpenalty(points)) / 2.0


class InsideBoxConstraint(Constraint):
    """Points must lie inside an axis-aligned box: ``penalty = Σ |r - box|²``.

    A point outside pays the squared distance to the nearest face(s); the
    penalty and its gradient are smooth, so a minimiser walks a stray point
    back in.
    """

    def __init__(self, length, origin=np.array([0, 0, 0])):
        length = np.asarray(length, dtype=float)
        if length.ndim == 0:  # scalar edge -> cube
            length = np.full(3, float(length))
        self.region = BoxRegion(length, origin)
        self.lengths = np.array(length)
        self.origin = np.array(origin, dtype=float)
        self.upper = self.origin + self.lengths

    def _excess(self, points: np.ndarray) -> np.ndarray:
        """Per-axis signed overshoot past the nearest face (0 inside)."""
        points = np.asarray(points, dtype=float)
        return np.minimum(points - self.origin, 0.0) + np.maximum(
            points - self.upper, 0.0
        )

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(self._excess(points) ** 2))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        return 2.0 * self._excess(points)

    def __invert__(self):
        return OutsideBoxConstraint(self.origin, self.upper - self.origin)


class OutsideBoxConstraint(Constraint):
    """Points must lie outside an axis-aligned box: ``penalty = Σ d_in²``.

    A point inside pays the squared distance to the nearest face (the shortest
    way out); the gradient points inward, so descending it leaves the box.
    """

    def __init__(self, origin, lengths):
        self.region = BoxRegion(lengths, origin)
        self.origin = np.array(origin, dtype=float)
        self.upper = self.origin + np.array(lengths, dtype=float)

    def _depth(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """``(inside, depth, axis)``: the shortest exit per interior point."""
        points = np.asarray(points, dtype=float)
        to_lower = points - self.origin
        to_upper = self.upper - points
        inside = np.all(to_lower > 0.0, axis=1) & np.all(to_upper > 0.0, axis=1)
        gap = np.minimum(to_lower, to_upper)  # (n, 3): distance to each face pair
        axis = np.argmin(gap, axis=1)
        depth = gap[np.arange(len(points)), axis]
        return inside, depth, axis

    def penalty(self, points: np.ndarray) -> float:
        inside, depth, _ = self._depth(points)
        return float(np.sum(depth[inside] ** 2))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        points = np.asarray(points, dtype=float)
        inside, depth, axis = self._depth(points)
        grad = np.zeros_like(points)
        rows = np.flatnonzero(inside)
        # d(depth²)/dr along the exit axis: +2·depth if the nearest face is the
        # lower one (moving up deepens), −2·depth if it is the upper one.
        centre = 0.5 * (self.origin + self.upper)
        sign = np.where(points[rows, axis[rows]] < centre[axis[rows]], 1.0, -1.0)
        grad[rows, axis[rows]] = 2.0 * depth[rows] * sign
        return grad

    def __invert__(self):
        return InsideBoxConstraint(self.upper - self.origin, self.origin)


class InsideSphereConstraint(Constraint):
    """Points must lie inside a sphere: ``penalty = Σ max(0, |r - c| - R)²``."""

    def __init__(self, radius, center):
        self.region = SphereRegion(radius, center)
        self.radius = float(radius)
        self.center = np.array(center, dtype=np.float64)

    def _overshoot(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        diff = np.asarray(points, dtype=float) - self.center
        dist = np.linalg.norm(diff, axis=1)
        return diff, dist, np.maximum(0.0, dist - self.radius)

    def penalty(self, points: np.ndarray) -> float:
        _, _, over = self._overshoot(points)
        return float(np.sum(over**2))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff, dist, over = self._overshoot(points)
        grad = np.zeros_like(diff)
        out = over > 0.0
        grad[out] = (2.0 * over[out] / dist[out])[:, None] * diff[out]
        return grad

    def __invert__(self):
        return OutsideSphereConstraint(self.radius, self.center)


class OutsideSphereConstraint(Constraint):
    """Points must lie outside a sphere: ``penalty = Σ max(0, R - |r - c|)²``."""

    def __init__(self, radius, center):
        self.region = SphereRegion(radius, center)
        self.radius = float(radius)
        self.center = np.array(center, dtype=np.float64)

    def _penetration(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        diff = np.asarray(points, dtype=float) - self.center
        dist = np.linalg.norm(diff, axis=1)
        return diff, dist, np.maximum(0.0, self.radius - dist)

    def penalty(self, points: np.ndarray) -> float:
        _, _, depth = self._penetration(points)
        return float(np.sum(depth**2))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff, dist, depth = self._penetration(points)
        grad = np.zeros_like(diff)
        inside = (depth > 0.0) & (dist > 1e-12)
        grad[inside] = (-2.0 * depth[inside] / dist[inside])[:, None] * diff[inside]
        return grad

    def __invert__(self):
        return InsideSphereConstraint(self.radius, self.center)


# === Min-distance constraint (pairwise distances) ===
class MinDistanceConstraint(Constraint):
    def __init__(self, dmin: float):
        self.dmin = dmin

    def _close_pairs(
        self, points: np.ndarray
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        """Pairs closer than ``dmin``: ``(i, j, displacement_ij, distance)``.

        The molrs neighbour query is O(N) in the number of points; only the
        pairs that can carry a penalty are ever materialised.
        """
        points = np.ascontiguousarray(points, dtype=np.float64)
        if len(points) < 2:
            empty = np.zeros(0, dtype=np.int64)
            return empty, empty, np.zeros((0, 3)), np.zeros(0)
        found = molrs.NeighborQuery.free(points, self.dmin).query_self()
        i = np.asarray(found.query_point_indices(), dtype=np.int64)
        j = np.asarray(found.point_indices(), dtype=np.int64)
        disp = np.asarray(found.disp(), dtype=np.float64)  # points[j] - points[i]
        return i, j, disp, np.sqrt(np.asarray(found.dist_sq(), dtype=np.float64))

    def penalty(self, points: np.ndarray) -> float:
        _, _, _, dist = self._close_pairs(points)
        violations = np.maximum(0.0, self.dmin - dist)
        return float(np.sum(violations**2))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        i, j, disp, dist = self._close_pairs(points)
        grad = np.zeros_like(points, dtype=np.float64)
        keep = (dist < self.dmin) & (dist > 1e-8)
        i, j, disp, dist = i[keep], j[keep], disp[keep], dist[keep]
        # d/dr_i of (dmin - |r_i - r_j|)^2 = -2 (dmin - d) (r_i - r_j) / d
        #                                  =  2 (dmin - d) (r_j - r_i) / d.
        # This is the gradient (uphill); descending it moves i away from j.
        grad_i = (2.0 * (self.dmin - dist) / dist)[:, None] * disp
        np.add.at(grad, i, grad_i)
        np.add.at(grad, j, -grad_i)
        return grad
