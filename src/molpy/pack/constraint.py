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
    def __init__(self, length, origin=np.array([0, 0, 0])):
        length = np.asarray(length, dtype=float)
        if length.ndim == 0:  # scalar edge -> cube
            length = np.full(3, float(length))
        self.region = BoxRegion(length, origin)
        self.lengths = np.array(length)
        self.origin = np.array(origin)
        self.upper = self.origin + self.lengths

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(~self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        not_in = ~self.region.isin(points)
        grad = np.zeros_like(points)
        lower_mask = points < self.origin
        upper_mask = points > self.origin + self.lengths
        grad[lower_mask & not_in[:, None]] = 1
        grad[upper_mask & not_in[:, None]] = -1
        return grad

    def __invert__(self):
        return OutsideBoxConstraint(self.origin, self.upper - self.origin)


class OutsideBoxConstraint(Constraint):
    def __init__(self, origin, lengths):
        self.region = BoxRegion(lengths, origin)
        self.origin = np.array(origin)
        self.upper = self.origin + np.array(lengths)

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        is_in = self.region.isin(points)
        grad = np.zeros_like(points)
        lower_mask = points > self.origin
        upper_mask = points < self.upper
        grad[lower_mask & is_in[:, None]] = -1
        grad[upper_mask & is_in[:, None]] = 1

        return grad

    def __invert__(self):
        return InsideBoxConstraint(self.origin, self.upper - self.origin)


class InsideSphereConstraint(Constraint):
    def __init__(self, radius, center):
        self.region = SphereRegion(radius, center)
        self.radius = radius
        self.center = np.array(center, dtype=np.float64)

    def penalty(self, points: np.ndarray) -> float:
        # Check how many points are outside the sphere (should be inside)
        return float(np.sum(~self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff = points - self.center
        dist = np.linalg.norm(diff, axis=1)
        not_in = dist > self.radius
        grad = np.zeros_like(points)
        # For points outside, gradient should point toward center (negative direction)
        grad[not_in] = -diff[not_in] / (dist[not_in, np.newaxis] + 1e-8)
        return grad

    def __invert__(self):
        return OutsideSphereConstraint(self.radius, self.center)


class OutsideSphereConstraint(Constraint):
    def __init__(self, radius, center):
        self.region = SphereRegion(radius, center)
        self.radius = radius
        self.center = np.array(center, dtype=np.float64)

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff = points - self.center
        dist = np.linalg.norm(diff, axis=1)
        is_in = dist <= self.radius
        grad = np.zeros_like(points)
        # For points inside the sphere, push them outward
        mask = is_in & (dist > 1e-8)  # Avoid division by zero
        grad[mask] = -diff[mask] / dist[mask, np.newaxis]
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
