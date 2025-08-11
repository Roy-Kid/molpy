"""Unit tests for LJ126 potential implementation."""

import numpy as np
import pytest

from molpy.potential.pair.lj import LJ126


class MockFrame:
    """Mock Frame object for testing."""

    def __init__(self, atoms_data):
        self.data = {"atoms": atoms_data}

    def __getitem__(self, key):
        return self.data[key]

    def __contains__(self, key):
        return key in self.data


class MockAtoms:
    """Mock Atoms object for testing."""

    def __init__(self, xyz, atom_types):
        self.xyz = xyz
        self.type = atom_types

    def __getitem__(self, key):
        return getattr(self, key)

    def __contains__(self, key):
        return hasattr(self, key)


class TestLJ126:
    """Test suite for LJ126 potential."""

    def test_init_scalar_params(self):
        """Test initialization with scalar parameters."""
        lj = LJ126(epsilon=1.0, sigma=2.0)
        assert lj.epsilon.shape == (1, 1)
        assert lj.sigma.shape == (1, 1)
        assert lj.epsilon[0, 0] == 1.0
        assert lj.sigma[0, 0] == 2.0

    def test_init_array_params(self):
        """Test initialization with array parameters."""
        epsilon = np.array([1.0, 2.0, 3.0])
        sigma = np.array([1.0, 1.5, 2.0])
        lj = LJ126(epsilon=epsilon, sigma=sigma)
        assert lj.epsilon.shape == (3, 1)
        assert lj.sigma.shape == (3, 1)
        np.testing.assert_array_equal(lj.epsilon.flatten(), epsilon)
        np.testing.assert_array_equal(lj.sigma.flatten(), sigma)

    def test_energy_static_method(self):
        """Test static energy calculation method."""
        epsilon = np.array([[1.0], [2.0]])
        sigma = np.array([[1.0], [1.5]])
        dr_norm = np.array([[1.0], [1.5]])

        energy = LJ126.energy(epsilon, sigma, dr_norm)
        assert energy.shape == (2, 1)
        assert not np.isnan(energy).any()

    def test_force_static_method(self):
        """Test static force calculation method."""
        epsilon = np.array([[1.0], [2.0]])
        sigma = np.array([[1.0], [1.5]])
        dr = np.array([[1.0, 0.0, 0.0], [0.0, 1.5, 0.0]])
        dr_norm = np.array([[1.0], [1.5]])

        force = LJ126.force(epsilon, sigma, dr, dr_norm)
        assert force.shape == (2, 3)
        assert not np.isnan(force).any()

    def test_calc_energy_direct(self):
        """Test calc_energy with direct arguments."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        dr_norm = np.array([[2.0], [1.0], [0.5]])
        pair_types = np.array([0, 0, 0])

        energies = lj.calc_energy(dr_norm, pair_types)
        assert energies.shape == (3, 1)
        assert not np.isnan(energies).any()

    def test_calc_force_direct(self):
        """Test calc_force with direct arguments."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        dr = np.array([[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        dr_norm = np.array([[1.0], [1.0]])
        pair_types = np.array([0, 0])

        forces = lj.calc_force(dr, dr_norm, pair_types)
        assert forces.shape == (2, 3)
        assert not np.isnan(forces).any()

    def test_calc_energy_frame(self):
        """Test calc_energy with Frame object."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Create a mock frame with 2 atoms
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        atom_types = np.array([0, 0])

        atoms = MockAtoms(coords, atom_types)
        frame = MockFrame(atoms)

        energies = lj.calc_energy(frame)  # type: ignore
        assert energies.shape[1] == 1
        assert not np.isnan(energies).any()

    def test_calc_force_frame(self):
        """Test calc_force with Frame object."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Create a mock frame with 2 atoms
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
        atom_types = np.array([0, 0])

        atoms = MockAtoms(coords, atom_types)
        frame = MockFrame(atoms)

        forces = lj.calc_force(frame)  # type: ignore
        assert forces.shape[1] == 3
        assert not np.isnan(forces).any()

    def test_empty_system(self):
        """Test with empty system."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        dr_norm = np.array([]).reshape(0, 1)
        pair_types = np.array([], dtype=np.int32)

        energies = lj.calc_energy(dr_norm, pair_types)
        assert energies.shape == (0, 1)

        dr = np.array([]).reshape(0, 3)
        forces = lj.calc_force(dr, dr_norm, pair_types)
        assert forces.shape == (0, 3)

    def test_multiple_atom_types(self):
        """Test with multiple atom types."""
        lj = LJ126(epsilon=np.array([1.0, 2.0]), sigma=np.array([1.0, 1.5]))

        dr_norm = np.array([[1.0], [1.5]])
        pair_types = np.array([0, 1])

        energies = lj.calc_energy(dr_norm, pair_types)
        assert energies.shape == (2, 1)

        dr = np.array([[1.0, 0.0, 0.0], [0.0, 1.5, 0.0]])
        forces = lj.calc_force(dr, dr_norm, pair_types)
        assert forces.shape == (2, 3)

    def test_validation_errors(self):
        """Test input validation errors."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Test invalid pair type indices
        with pytest.raises(ValueError, match="pair_types contains invalid indices"):
            dr_norm = np.array([[1.0]])
            pair_types = np.array([1])  # Out of range
            lj.calc_energy(dr_norm, pair_types)

        # Test mismatched array shapes
        with pytest.raises(ValueError, match="must have the same number of pairs"):
            dr = np.array([[1.0, 0.0, 0.0]])
            dr_norm = np.array([[1.0], [2.0]])  # Different length
            pair_types = np.array([0])
            lj.calc_force(dr, dr_norm, pair_types)

    def test_extract_data_from_frame(self):
        """Test _extract_data_from_frame method."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Create mock frame with 3 atoms
        coords = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        atom_types = np.array([0, 0, 0])

        atoms = MockAtoms(coords, atom_types)
        frame = MockFrame(atoms)

        dr, dr_norm, pair_types = lj._extract_data_from_frame(frame)  # type: ignore

        # Should have 3 pairs from 3 atoms
        expected_pairs = 3
        assert dr.shape == (expected_pairs, 3)
        assert dr_norm.shape == (expected_pairs, 1)
        assert pair_types.shape == (expected_pairs,)

    def test_lj_physics(self):
        """Test that LJ potential behaves physically correctly."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Test at sigma distance (should be near zero energy)
        dr_norm = np.array([[1.0]])
        pair_types = np.array([0])
        energy = lj.calc_energy(dr_norm, pair_types)
        assert abs(energy[0, 0]) < 1e-10  # Should be very close to zero

        # Test at very close distance (should be high positive energy)
        dr_norm = np.array([[0.5]])
        energy = lj.calc_energy(dr_norm, pair_types)
        assert energy[0, 0] > 0  # Repulsive

        # Test at far distance (should be negative energy)
        dr_norm = np.array([[1.5]])
        energy = lj.calc_energy(dr_norm, pair_types)
        assert energy[0, 0] < 0  # Attractive

    def test_cutoff_effect(self):
        """Test behavior at cutoff distance."""
        lj = LJ126(epsilon=1.0, sigma=1.0)

        # Test just below cutoff distance
        dr_norm = np.array([[2.4]])
        pair_types = np.array([0])
        energy_below = lj.calc_energy(dr_norm, pair_types)

        # Test at very far distance (should be very small)
        dr_norm = np.array([[10.0]])
        energy_far = lj.calc_energy(dr_norm, pair_types)
        assert abs(energy_far[0, 0]) < abs(energy_below[0, 0])  # Should be smaller

        # Test forces also decrease with distance
        dr = np.array([[10.0, 0.0, 0.0]])
        dr_norm = np.array([[10.0]])
        force_far = lj.calc_force(dr, dr_norm, pair_types)
        assert np.linalg.norm(force_far[0]) < 1e-5  # Should be very small
