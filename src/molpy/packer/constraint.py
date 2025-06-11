import numpy as np
import molpy as mp

# === Base Constraint class ===
class Constraint:
    def penalty(self, points: np.ndarray) -> float:
        raise NotImplementedError
    
    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        raise NotImplementedError
    
    def __and__(self, other):
        return AndConstraint(self, other)

    def __or__(self, other):
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
        else:
            return self.b.dpenalty(points)
    
class InsideBoxConstraint(Constraint):
    def __init__(self, length, origin=np.array([0, 0, 0])):
        self.region = mp.BoxRegion(length, origin)
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
        grad[lower_mask&not_in[:, None]] = 1
        grad[upper_mask&not_in[:, None]] = -1
        return grad

    def __invert__(self):
        return OutsideBoxConstraint(self.origin, self.upper - self.origin)

class OutsideBoxConstraint(Constraint):
    def __init__(self, origin, lengths):
        self.region = mp.BoxRegion(lengths, origin)
        self.origin = np.array(origin)
        self.upper = self.origin + np.array(lengths)

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        is_in = self.region.isin(points)
        grad = np.zeros_like(points)
        lower_mask = points > self.origin
        upper_mask = points < self.upper
        grad[lower_mask&is_in[:, None]] = -1
        grad[upper_mask&is_in[:, None]] = 1

        return grad

    def __invert__(self):
        return InsideBoxConstraint(self.origin, self.upper - self.origin)
    
class InsideSphereConstraint(Constraint):
    def __init__(self, radius, center):
        self.region = mp.SphereRegion(radius, center)
        self.radius = radius
        self.center = np.array(center, dtype=np.float64)

    def penalty(self, points: np.ndarray) -> float:
        # dist = np.linalg.norm(points - self.center, axis=1)
        # return float(np.sum(dist > self.radius))
        return float(np.sum(self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff = points - self.center
        dist = np.linalg.norm(diff, axis=1)
        not_in = dist > self.radius
        grad = np.zeros_like(points)
        grad[not_in] = diff[not_in] / (dist[not_in] + 1e-8)

        return grad

    def __invert__(self):
        return OutsideSphereConstraint(self.radius, self.center)

class OutsideSphereConstraint(Constraint):

    def __init__(self, radius, center):
        self.region = mp.SphereRegion(radius, center)
        self.radius = radius
        self.center = np.array(center, dtype=np.float64)

    def penalty(self, points: np.ndarray) -> float:
        return float(np.sum(self.region.isin(points)))

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff = points - self.center
        dist = np.linalg.norm(diff, axis=1)
        is_in = dist <= self.radius
        grad = np.zeros_like(points)
        grad[is_in] = 1/diff[is_in]
        return grad      

    def __invert__(self):
        return InsideSphereConstraint(self.radius, self.center)

    
# === Min-distance constraint (pairwise distances) ===
class MinDistanceConstraint(Constraint):
    def __init__(self, dmin: float):
        self.dmin = dmin

    def penalty(self, points: np.ndarray) -> float:
        dist = np.linalg.norm(points[:, None, :] - points[None, :, :], axis=-1)
        i, j = np.triu_indices(len(points), k=1)
        d_ij = dist[i, j]
        violations = np.maximum(0, self.dmin - d_ij)
        penalty = np.sum(violations ** 2)
        return penalty

    def dpenalty(self, points: np.ndarray) -> np.ndarray:
        diff = points[:, None, :] - points[None, :, :]
        dist = np.linalg.norm(diff, axis=-1)
        i, j = np.triu_indices(len(points), k=1)
        too_close = np.zeros(len(points), dtype=bool)
        too_close = np.any(dist < self.dmin, axis=-1)

        grad = np.zeros_like(points)
        np.add.at(
            grad, i, -diff[i, j]
        )
        np.add.at(
            grad, j, diff[i, j]
        )
        grad[~too_close] = 0 
        return grad