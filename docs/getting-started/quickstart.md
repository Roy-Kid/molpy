# Quickstart Guide

Get up and running with MolPy in minutes. This guide walks you through installation verification, creating your first molecular system, manipulating it, and running basic analysis.

## Installation Verification

### 1. Check Your Installation

First, verify that MolPy is properly installed:

```python
import molpy as mp
import numpy as np

print(f"MolPy version: {mp.__version__}")
print(f"NumPy version: {np.__version__}")

# Test basic functionality
print("✓ MolPy imported successfully")
```

### 2. Verify Core Components

Test that the main components are working:

```python
# Test Block creation
test_block = mp.Block({'x': [1.0, 2.0], 'y': [0.0, 0.0]})
print(f"✓ Block creation: {test_block.nrows} rows")

# Test Frame creation
test_frame = mp.Frame({'atoms': test_block})
print(f"✓ Frame creation: {len(test_frame['atoms'])} atoms")

# Test Box creation
test_box = mp.Box.cubic(10.0)
print(f"✓ Box creation: {test_box.volume:.1f} Å³")

print("All core components working correctly!")
```

## Your First MolPy Session

### 3. Create a Simulation Box

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

### 4. Create Atomic Data

MolPy uses `Block` objects to store atomic properties:

```python
# Create atomic coordinates and properties
atoms_data = {
    'x': [0.0, 1.0, 2.0, 3.0],
    'y': [0.0, 0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C', 'H'],
    'type': [1, 1, 1, 2],
    'mass': [12.01, 12.01, 12.01, 1.008]
}

# Create a Block from the data
atoms_block = mp.Block(atoms_data)
print(f"Created {len(atoms_block)} atoms")
```

### 5. Build Your First Frame

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

### 6. Basic Atom Selection

Select atoms using various criteria:

```python
# Select by element type
carbons = frame['atoms']['element'] == 'C'
print(f"Found {np.sum(carbons)} carbon atoms")

# Select by position (atoms in first half of box)
first_half = frame['atoms']['x'] < 5.0
print(f"Atoms in first half: {np.sum(first_half)}")

# Select by atom type
type_1 = frame['atoms']['type'] == 1
print(f"Type 1 atoms: {np.sum(type_1)}")
```

### 7. Simple Analysis

Calculate basic properties:

```python
# Calculate center of mass
masses = frame['atoms']['mass']
positions = np.column_stack([
    frame['atoms']['x'],
    frame['atoms']['y'],
    frame['atoms']['z']
])

com = np.average(positions, weights=masses, axis=0)
print(f"Center of mass: {com}")

# Calculate total mass
total_mass = np.sum(masses)
print(f"Total mass: {total_mass:.2f} amu")

# Calculate radius of gyration
distances = positions - com
rg_squared = np.average(np.sum(distances**2, axis=1), weights=masses)
rg = np.sqrt(rg_squared)
print(f"Radius of gyration: {rg:.2f} Å")
```

### 8. Save and Load

MolPy provides unified serialization:

```python
# Save to dictionary (JSON-serializable)
frame_dict = frame.to_dict()

# Save to file (HDF5 format)
frame.save('my_system.h5')

# Load from file
loaded_frame = mp.Frame.load('my_system.h5')

# Verify data integrity
print(f"Loaded frame has {len(loaded_frame['atoms'])} atoms")
print(f"Data identical: {np.allclose(frame['atoms']['x'], loaded_frame['atoms']['x'])}")
```

## Complete First Workflow

### Building a Water Dimer System

Let's create a complete molecular system from scratch:

```python
def create_water_dimer():
    """Create a complete water dimer system."""

    # Create water molecule 1
    water1_atoms = mp.Block({
        'x': [0.0, 0.9572, -0.2400],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H'],
        'type': [1, 2, 2],
        'mass': [15.999, 1.008, 1.008]
    })

    # Create water molecule 2 (offset)
    water2_atoms = mp.Block({
        'x': [3.0, 3.9572, 2.7600],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H'],
        'type': [1, 2, 2],
        'mass': [15.999, 1.008, 1.008]
    })

    # Combine into dimer
    dimer_atoms = mp.Block({
        'x': np.concatenate([water1_atoms['x'], water2_atoms['x']]),
        'y': np.concatenate([water1_atoms['y'], water2_atoms['y']]),
        'z': np.concatenate([water1_atoms['z'], water2_atoms['z']]),
        'element': np.concatenate([water1_atoms['element'], water2_atoms['element']),
        'type': np.concatenate([water1_atoms['type'], water2_atoms['type']),
        'mass': np.concatenate([water1_atoms['mass'], water2_atoms['mass'])
    })

    # Create frame with box
    box = mp.Box.cubic(20.0)
    dimer_frame = mp.Frame(data={'atoms': dimer_atoms}, box=box)

    return dimer_frame

# Create the water dimer
water_dimer = create_water_dimer()
print(f"✓ Created water dimer with {len(water_dimer['atoms'])} atoms")

# Analyze the dimer
dimer_masses = water_dimer['atoms']['mass']
dimer_positions = np.column_stack([
    water_dimer['atoms']['x'],
    water_dimer['atoms']['y'],
    water_dimer['atoms']['z']
])

dimer_com = np.average(dimer_positions, weights=dimer_masses, axis=0)
print(f"✓ Dimer center of mass: {dimer_com}")

# Save the system
water_dimer.save('water_dimer.h5')
print("✓ Water dimer saved to 'water_dimer.h5'")
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

## Common Patterns

### Working with Multiple Frames

```python
# Create trajectory-like data
frames = []
for i in range(5):
    # Simulate some movement
    new_positions = positions + np.random.normal(0, 0.1, positions.shape)
    new_data = atoms_data.copy()
    new_data['x'] = new_positions[:, 0]
    new_data['y'] = new_positions[:, 1]
    new_data['z'] = new_positions[:, 2]

    frame = mp.Frame(data={'atoms': new_data}, box=box)
    frames.append(frame)

print(f"Created {len(frames)} frames")
```

### Batch Processing

```python
# Process multiple frames
results = []
for frame in frames:
    # Calculate property for each frame
    com = np.average(positions, weights=masses, axis=0)
    results.append(com)

results = np.array(results)
print(f"Center of mass trajectory shape: {results.shape}")
```

## Troubleshooting Common Issues

### Import Errors
```python
# If you get import errors, check your Python environment
import sys
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version}")

# Make sure you're in the right environment
import molpy
print(f"MolPy location: {molpy.__file__}")
```

### Data Type Issues
```python
# MolPy automatically converts data to NumPy arrays
# If you have issues, ensure your data is numeric
atoms_data = {
    'x': [0.0, 1.0, 2.0],  # Must be numeric
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C']  # Can be strings
}
```

Ready to explore more? Check out the [tutorials](../tutorials/block-frame.md) for advanced features!
