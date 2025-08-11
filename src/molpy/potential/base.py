from collections import UserList

import numpy as np

from molpy.core.forcefield import KernelMeta
from molpy.core.frame import Frame


class Potential(metaclass=KernelMeta):
    """
    Base class for all potential functions in MolPy.

    This class provides a template for defining potential functions that can be used in molecular simulations.
    It includes methods for evaluating the potential and its derivatives, as well as a method to check if the potential is periodic.
    """

    def __call__(self, *args, **kwargs):
        """Evaluate the potential."""
        raise NotImplementedError("Subclasses must implement this method.")

    def calc_energy(self, frame: Frame) -> float: ...

    def calc_forces(self, frame: Frame) -> np.ndarray:
        """
        Calculate the forces acting on the particles in the given frame.

        Parameters
        ----------
        frame : Frame
            The frame containing the particle positions and other relevant data.

        Returns
        -------
        np.ndarray
            An array of forces acting on each particle.
        """
        raise NotImplementedError("Subclasses must implement this method.")


class Potentials(UserList[Potential]):

    def calc_energy(self, frame: Frame) -> float:
        """
        Calculate the total energy of the system by summing the energies from all potentials.

        Parameters
        ----------
        frame : Frame
            The frame containing the particle positions and other relevant data.

        Returns
        -------
        float
            The total energy of the system.
        """
        return sum(pot.calc_energy(frame) for pot in self)

    def calc_forces(self, frame: Frame) -> np.ndarray:
        """
        Calculate the total forces acting on the particles by summing the forces from all potentials.

        Parameters
        ----------
        frame : Frame
            The frame containing the particle positions and other relevant data.

        Returns
        -------
        np.ndarray
            An array of total forces acting on each particle.
        """
        return sum(
            (pot.calc_forces(frame) for pot in self),
            start=np.zeros_like((frame["atoms", "xyz"])),
        )
