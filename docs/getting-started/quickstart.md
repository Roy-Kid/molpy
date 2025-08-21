# Quickstart Guide

Get up and running with MolPy in minutes. This guide walks you through installation verification, creating your first molecular system, manipulating it, and running basic analysis.

## Installation Verification

### 1. Check Your Installation

First, verify that MolPy is properly installed:

```python
import molpy as mp
import numpy as np

print(f"MolPy version: {mp.version.version}")
print(f"NumPy version: {np.version.version}")

# Test basic functionality
print("✓ MolPy imported successfully")
```

## Your First MolPy Session

### 1. Create a Simulation Box

MolPy uses `Box` objects to define periodic boundary conditions:

```python
# Cubic box with 10 Å sides
box = mp.Box.cubic(10.0)

# Orthogonal box with different dimensions
box_orth = mp.Box.orth([8.0, 10.0, 12.0])

# Check box properties
print(f"Box volume: {box.volume:.1f} Å³")
print(f"Box lengths: {box.lengths}")
print(f"Periodic boundaries: {box.pbc}")
```

### 2. Create Atomic Data

MolPy uses `Block` objects to store atomic properties:

```python
# Create atomic coordinates and properties
atoms_data = {
    'x': [0.0, -0.06, 0.0],
    'y': [0.75, 0.52, 0.0],
    'z': [-0.75, 0.52, 0.0],
    'element': ['O', 'H', 'H'],
    'type': [1, 2, 2],
    'mass': [15.999, 1.008, 1.008]
}

# Create a Block from the data
atoms_block = mp.Block(atoms_data)
print(f"Created {len(atoms_block)} atoms")
```

### 3. Build Your First Frame

A `Frame` is the core container that holds atomic data and simulation box:

```python
# Create a frame with atoms and box
frame = mp.Frame(data={'atoms': atoms_data}, box=box)

# Access data
print(f"Frame has {len(frame['atoms'])} atoms")
print(f"Box volume: {frame.box.volume:.1f} Å³")

# Check what's in the frame
print(f"Available blocks: {list(frame.blocks())}")
print(f"Available variables: {list(frame.variables('atoms'))}")
```

### 4. Define Force Field

```python

ff = mp.ForceField()
atomstyle = ff.def_atomstyle("full")
atomtype1 = atomstyle.def_type("C", mass=12.01)
atomtype2 = atomstyle.def_type("H", mass=1.008)

bondstyle = ff.def_bondstyle("harmonic")
bondtype1 = bondstyle.def_type(atomtype1, atomtype1, k=1000.0, r0=1.0)

anglestyle = ff.def_anglestyle("harmonic")

```

### 5. Basic Atom Selection

Select atoms blocks using various criteria:

```python
# Select by element type
carbons = frame["atoms"][ frame['atoms']['element'] == 'C' ]
# Select by atom type
select_type1 = frame["atoms"][mp.AtomTypeSelection(1)]
```

### 8. Save and Load

MolPy provides unified serialization and io functions:

```python
# simple frame
mp.io.write_pdb("frame.pdb", frame)
frame = mp.io.read_pdb("frame.pdb")

# read lammps trajectory
for frame in mp.io.read_lammps_trajectory("trajectory.lammpstrj"):
    print(frame)

# forcefield
mp.io.write_lammps_forcefield("forcefield.lammps", ff)
```

## What You've Learned

✅ **Installation verification** - Ensure everything is working
✅ **Box creation** - Define simulation boundaries
✅ **Atomic data** - Store coordinates and properties
✅ **Frame assembly** - Combine data and geometry
✅ **Atom selection** - Filter atoms by criteria
✅ **Basic analysis** - Calculate molecular properties
✅ **Data persistence** - Save and load your work
✅ **Complete workflow** - Build a real molecular system

## Next Steps

- **[Block & Frame Basics](../tutorials/block-frame.md)** - Deep dive into data structures
- **[Selection & Filtering](../tutorials/selection.md)** - Advanced atom selection
- **[Analysis Tools](../tutorials/trajectory.md)** - Built-in analysis functions
- **[API Reference](../reference/index.md)** - Complete function reference

Ready to explore more? Check out the [tutorials](../tutorials/block-frame.md) for advanced features!
