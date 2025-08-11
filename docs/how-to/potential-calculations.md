# How to Calculate Potential Energies and Forces

This guide shows practical examples of working with potential energy functions in MolPy. You'll learn how to calculate energies, forces, and analyze molecular interactions.

## Understanding Potential Functions

Potential functions describe the energy landscape of molecular systems. They include bonded interactions (bonds, angles, dihedrals) and non-bonded interactions (van der Waals, electrostatic).

## Basic Potential Calculations

### Creating Simple Potentials

Start with basic potential functions for common interactions:

```python
import molpy as mp
import numpy as np

class HarmonicBondPotential(mp.Potential):
    """Simple harmonic bond potential."""

    def __init__(self, k=100.0, r0=1.0):
        self.k = k  # Force constant
        self.r0 = r0  # Equilibrium distance

    def calc_energy(self, frame):
        """Calculate bond energy."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        # Simple example: calculate energy for first two atoms
        if len(positions) >= 2:
            r = np.linalg.norm(positions[1] - positions[0])
            energy = 0.5 * self.k * (r - self.r0)**2
            return energy
        return 0.0

    def calc_forces(self, frame):
        """Calculate forces from bond potential."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        forces = np.zeros_like(positions)

        if len(positions) >= 2:
            r_vec = positions[1] - positions[0]
            r = np.linalg.norm(r_vec)

            if r > 0:
                # Force magnitude
                force_mag = -self.k * (r - self.r0)
                # Force vector
                force_vec = force_mag * r_vec / r

                forces[0] = -force_vec  # Force on atom 0
                forces[1] = force_vec   # Force on atom 1

        return forces

# Create potential
bond_potential = HarmonicBondPotential(k=100.0, r0=1.0)
print(f"Created harmonic bond potential with k={bond_potential.k}, r0={bond_potential.r0}")
```

### Using Potentials with Frames

Apply potentials to molecular frames to calculate energies and forces:

```python
# Create a simple molecular frame
water_atoms = mp.Block({
    'x': [0.0, 0.9572, -0.2400],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['O', 'H', 'H']
})

water_frame = mp.Frame({'atoms': water_atoms})

# Calculate energy and forces
energy = bond_potential.calc_energy(water_frame)
forces = bond_potential.calc_forces(water_frame)

print(f"Bond energy: {energy:.3f} kJ/mol")
print(f"Forces shape: {forces.shape}")
print(f"Force on oxygen: {forces[0]}")
```

## Multiple Potential Functions

### Combining Different Potentials

Use the Potentials container to combine multiple interaction types:

```python
class LennardJonesPotential(mp.Potential):
    """Lennard-Jones potential for non-bonded interactions."""

    def __init__(self, epsilon=0.066, sigma=3.5):
        self.epsilon = epsilon
        self.sigma = sigma

    def calc_energy(self, frame):
        """Calculate LJ energy."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        energy = 0.0
        n_atoms = len(positions)

        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                r_vec = positions[j] - positions[i]
                r = np.linalg.norm(r_vec)

                if r > 0:
                    r6 = (self.sigma / r)**6
                    r12 = r6**2
                    energy += 4 * self.epsilon * (r12 - r6)

        return energy

    def calc_forces(self, frame):
        """Calculate LJ forces."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        forces = np.zeros_like(positions)
        n_atoms = len(positions)

        for i in range(n_atoms):
            for j in range(i+1, n_atoms):
                r_vec = positions[j] - positions[i]
                r = np.linalg.norm(r_vec)

                if r > 0:
                    r6 = (self.sigma / r)**6
                    r12 = r6**2

                    # Force magnitude
                    force_mag = 24 * self.epsilon * (2*r12 - r6) / r
                    # Force vector
                    force_vec = force_mag * r_vec / r

                    forces[i] += force_vec
                    forces[j] -= force_vec

        return forces

# Create LJ potential
lj_potential = LennardJonesPotential(epsilon=0.066, sigma=3.5)

# Combine potentials
all_potentials = mp.Potentials([bond_potential, lj_potential])

# Calculate total energy and forces
total_energy = all_potentials.calc_energy(water_frame)
total_forces = all_potentials.calc_forces(water_frame)

print(f"Total energy: {total_energy:.3f} kJ/mol")
print(f"Total forces shape: {total_forces.shape}")
```

## Advanced Potential Applications

### Custom Potential Functions

Create specialized potentials for specific applications:

```python
class MorsePotential(mp.Potential):
    """Morse potential for anharmonic bonds."""

    def __init__(self, D=100.0, alpha=2.0, r0=1.0):
        self.D = D      # Well depth
        self.alpha = alpha  # Width parameter
        self.r0 = r0    # Equilibrium distance

    def calc_energy(self, frame):
        """Calculate Morse potential energy."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        energy = 0.0
        if len(positions) >= 2:
            r = np.linalg.norm(positions[1] - positions[0])
            exp_term = np.exp(-self.alpha * (r - self.r0))
            energy = self.D * (1 - exp_term)**2

        return energy

    def calc_forces(self, frame):
        """Calculate Morse potential forces."""
        atoms = frame['atoms']
        positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

        forces = np.zeros_like(positions)

        if len(positions) >= 2:
            r_vec = positions[1] - positions[0]
            r = np.linalg.norm(r_vec)

            if r > 0:
                exp_term = np.exp(-self.alpha * (r - self.r0))
                force_mag = -2 * self.D * self.alpha * exp_term * (1 - exp_term)
                force_vec = force_mag * r_vec / r

                forces[0] = -force_vec
                forces[1] = force_vec

        return forces

# Create Morse potential
morse_potential = MorsePotential(D=100.0, alpha=2.0, r0=1.0)

# Test with water frame
morse_energy = morse_potential.calc_energy(water_frame)
morse_forces = morse_potential.calc_forces(water_frame)

print(f"Morse energy: {morse_energy:.3f} kJ/mol")
```

### Potential Energy Landscapes

Analyze potential energy landscapes for molecular systems:

```python
def analyze_potential_landscape(potential, frame, scan_range=(-2.0, 2.0), n_points=50):
    """Analyze potential energy landscape by scanning coordinates."""

    original_positions = frame['atoms'][['x', 'y', 'z']].copy()
    energies = []
    distances = []

    # Scan along x-axis for first atom
    for i, x_offset in enumerate(np.linspace(scan_range[0], scan_range[1], n_points)):
        # Modify position
        frame['atoms']['x'][0] = original_positions['x'][0] + x_offset

        # Calculate energy
        energy = potential.calc_energy(frame)
        energies.append(energy)
        distances.append(x_offset)

    # Restore original positions
    frame['atoms']['x'][0] = original_positions['x'][0]

    return distances, energies

# Analyze bond potential landscape
distances, energies = analyze_potential_landscape(bond_potential, water_frame)

print("Potential energy landscape:")
for d, e in zip(distances[::5], energies[::5]):  # Show every 5th point
    print(f"  Distance: {d:6.2f} Å, Energy: {e:8.3f} kJ/mol")
```

## Practical Applications

### Energy Minimization

Use potentials for simple energy minimization:

```python
def simple_energy_minimization(potential, frame, max_steps=100, step_size=0.01):
    """Simple gradient descent energy minimization."""

    positions = frame['atoms'][['x', 'y', 'z']].copy()
    energies = []

    for step in range(max_steps):
        # Calculate current energy
        energy = potential.calc_energy(frame)
        energies.append(energy)

        # Calculate forces
        forces = potential.calc_forces(frame)

        # Update positions (simple gradient descent)
        for i in range(len(positions)):
            frame['atoms']['x'][i] += step_size * forces[i][0]
            frame['atoms']['y'][i] += step_size * forces[i][1]
            frame['atoms']['z'][i] += step_size * forces[i][2]

        # Check convergence
        if step > 0 and abs(energies[-1] - energies[-2]) < 1e-6:
            print(f"Converged at step {step}")
            break

    return energies

# Minimize water structure
initial_energy = bond_potential.calc_energy(water_frame)
print(f"Initial energy: {initial_energy:.3f} kJ/mol")

energies = simple_energy_minimization(bond_potential, water_frame)

final_energy = bond_potential.calc_energy(water_frame)
print(f"Final energy: {final_energy:.3f} kJ/mol")
print(f"Energy change: {final_energy - initial_energy:.3f} kJ/mol")
```

### Force Analysis

Analyze forces to understand molecular behavior:

```python
def analyze_forces(potential, frame):
    """Analyze forces in a molecular system."""

    forces = potential.calc_forces(frame)
    atoms = frame['atoms']

    print("Force analysis:")
    print(f"Total force magnitude: {np.linalg.norm(forces.sum(axis=0)):.6f}")

    for i, element in enumerate(atoms['element']):
        force_mag = np.linalg.norm(forces[i])
        print(f"  {element}{i}: force magnitude = {force_mag:.6f}")

    # Check if forces sum to zero (conservation of momentum)
    total_force = forces.sum(axis=0)
    print(f"Net force: {total_force}")

    return forces

# Analyze forces
forces = analyze_forces(all_potentials, water_frame)
```

## Summary

This guide covered practical potential energy calculations:

- Create custom potential functions
- Combine multiple interaction types
- Calculate energies and forces
- Analyze potential landscapes
- Perform energy minimization
- Analyze molecular forces

MolPy's potential system provides a flexible framework for defining and using various interaction potentials. The key is understanding how to implement the required methods and combine potentials effectively.

### Next Steps

Continue exploring MolPy by learning about packing algorithms, analysis techniques, and advanced simulation workflows.
