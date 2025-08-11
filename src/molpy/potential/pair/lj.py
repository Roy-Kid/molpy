from typing import Dict, Tuple, Union, overload

import numpy as np
from numpy.typing import NDArray

import molpy as mp

from .base import PairPotential


class LJ126(PairPotential):
    """
    Lennard-Jones 12-6 pair potential with cutoff.

    The potential is defined as:
    V(r) = 4 * ε * ((σ/r)^12 - (σ/r)^6)

    The force is:
    F(r) = 24 * ε * (2 * (σ/r)^12 - (σ/r)^6) * dr / r^2

    Neighbor list and periodicity corrections are assumed done externally;
    dr and dr_norm are provided by the caller.

    Attributes:
        epsilon: Depth of potential well for each atom type [energy]
        sigma: Finite distance at which potential is zero [length]
    """

    name = "lj126/cut"
    type = "pair"

    def __init__(
        self,
        epsilon: Union[float, NDArray[np.float64]],
        sigma: Union[float, NDArray[np.float64]],
    ):
        """
        Initialize LJ126 potential.

        Args:
            epsilon: Depth of potential well, can be scalar or array for multiple types
            sigma: Finite distance at which potential is zero, can be scalar or array for multiple types
        """
        self.epsilon = np.array(epsilon, dtype=np.float64).reshape(-1, 1)
        self.sigma = np.array(sigma, dtype=np.float64).reshape(-1, 1)

        if self.epsilon.shape != self.sigma.shape:
            raise ValueError("epsilon and sigma must have the same shape")

    def _extract_data_from_frame(self, frame):
        """Extract necessary data from a Frame object."""
        # Extract atom coordinates
        atoms = frame["atoms"]
        coords = (
            atoms["xyz"]
            if "xyz" in atoms
            else np.column_stack([atoms["x"], atoms["y"], atoms["z"]])
        )
        atom_types = (
            atoms["type"] if "type" in atoms else np.zeros(len(coords), dtype=np.int32)
        )

        # If pairs information is pre-computed in frame, use it
        if "pairs" in frame:
            pairs = frame["pairs"]
            dr = pairs["dr"]
            dr_norm = np.linalg.norm(dr, axis=1, keepdims=True)
            pair_types = pairs["type"]
        else:
            # Simple all-pairs calculation (for testing/small systems)
            # In production, you'd use a neighbor list
            n_atoms = len(coords)
            dr_list = []
            pair_types_list = []

            for i in range(n_atoms):
                for j in range(i + 1, n_atoms):
                    dr_vec = coords[j] - coords[i]
                    dr_list.append(dr_vec)
                    # Use combination of atom types for pair type
                    pair_types_list.append(max(atom_types[i], atom_types[j]))

            if dr_list:
                dr = np.array(dr_list)
                dr_norm = np.linalg.norm(dr, axis=1, keepdims=True)
                pair_types = np.array(pair_types_list, dtype=np.int32)
            else:
                # No pairs found
                dr = np.empty((0, 3))
                dr_norm = np.empty((0, 1))
                pair_types = np.empty(0, dtype=np.int32)

        return dr, dr_norm, pair_types

    @staticmethod
    def energy(
        epsilon: np.ndarray, sigma: np.ndarray, dr_norm: np.ndarray
    ) -> np.ndarray:
        """
        Lennard-Jones pair energy: 4*epsilon*((sigma/dr)^12 - (sigma/dr)^6)
        Returns per-pair energy, shape (n_pairs, 1).
        """
        return 4 * epsilon * ((sigma / dr_norm) ** 12 - (sigma / dr_norm) ** 6)

    @staticmethod
    def force(
        epsilon: np.ndarray, sigma: np.ndarray, dr: np.ndarray, dr_norm: np.ndarray
    ) -> np.ndarray:
        """
        Lennard-Jones pair force vector:
        24*epsilon*(2*(sigma/dr)^12 - (sigma/dr)^6) * dr / dr_norm^2
        Returns per-pair force vectors, shape (n_pairs, 3).
        """
        return (
            24
            * epsilon
            * (2 * (sigma / dr_norm) ** 12 - (sigma / dr_norm) ** 6)
            * dr
            / (dr_norm**2)
        )

    @overload
    def calc_energy(self, frame: "mp.Frame") -> NDArray[np.float64]: ...

    @overload
    def calc_energy(
        self, dr_norm: NDArray[np.float64], pair_types: NDArray[np.int32]
    ) -> NDArray[np.float64]: ...

    def calc_energy(self, *args, **kwargs):
        """
        Compute pair energies.

        Can be called with either:
        1. A Frame object: calc_energy(frame)
        2. Direct data: calc_energy(dr_norm, pair_types)

        Returns:
            energies: Array of shape (n_pairs, 1), pair energies
        """
        if len(args) == 1 and hasattr(args[0], "__getitem__") and "atoms" in args[0]:
            # Frame input
            frame = args[0]
            dr, dr_norm, pair_types = self._extract_data_from_frame(frame)
        elif len(args) == 2:
            # Direct input
            dr_norm, pair_types = args
        else:
            raise ValueError(
                "calc_energy expects either a Frame or (dr_norm, pair_types)"
            )

        if len(pair_types) == 0:
            return np.array([]).reshape(0, 1)

        if np.any(pair_types >= len(self.epsilon)) or np.any(pair_types < 0):
            raise ValueError(
                f"pair_types contains invalid indices. Must be in range [0, {len(self.epsilon)})"
            )

        eps = self.epsilon[pair_types]
        sig = self.sigma[pair_types]
        return self.energy(eps, sig, dr_norm)

    @overload
    def calc_force(self, frame: "mp.Frame") -> NDArray[np.float64]: ...

    @overload
    def calc_force(
        self,
        dr: NDArray[np.float64],
        dr_norm: NDArray[np.float64],
        pair_types: NDArray[np.int32],
    ) -> NDArray[np.float64]: ...

    def calc_force(self, *args, **kwargs):
        """
        Compute pair force vectors.

        Can be called with either:
        1. A Frame object: calc_force(frame)
        2. Direct data: calc_force(dr, dr_norm, pair_types)

        Returns:
            forces: Array of shape (n_pairs, 3), force vectors
        """
        if len(args) == 1 and hasattr(args[0], "__getitem__") and "atoms" in args[0]:
            # Frame input
            frame = args[0]
            dr, dr_norm, pair_types = self._extract_data_from_frame(frame)
        elif len(args) == 3:
            # Direct input
            dr, dr_norm, pair_types = args
        else:
            raise ValueError(
                "calc_force expects either a Frame or (dr, dr_norm, pair_types)"
            )

        if len(pair_types) == 0:
            return np.array([]).reshape(0, 3)

        if dr.shape[0] != dr_norm.shape[0] or dr.shape[0] != pair_types.shape[0]:
            raise ValueError(
                "dr, dr_norm, and pair_types must have the same number of pairs"
            )

        if dr.shape[1] != 3:
            raise ValueError("dr must have shape (n_pairs, 3)")

        if dr_norm.shape[1] != 1:
            raise ValueError("dr_norm must have shape (n_pairs, 1)")

        if np.any(pair_types >= len(self.epsilon)) or np.any(pair_types < 0):
            raise ValueError(
                f"pair_types contains invalid indices. Must be in range [0, {len(self.epsilon)})"
            )

        eps = self.epsilon[pair_types]
        sig = self.sigma[pair_types]
        return self.force(eps, sig, dr, dr_norm)
