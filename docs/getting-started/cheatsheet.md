# MolPy Cheat Sheet

Quick reference for the most common MolPy operations.

## Import and Setup

```python
import molpy as mp
import numpy as np
```

## Core Data Structures

### Creating Frames

```python
# Basic frame
atoms_data = {
    'x': [0.0, 1.0, 2.0],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C']
}
frame = mp.Frame(data={'atoms': atoms_data})

# Frame with box and metadata
box = mp.Box.cubic(10.0)
frame = mp.Frame(
    data={'atoms': atoms_data},
    box=box,
    metadata={'time': 0.0, 'temperature': 300.0}
)
```

### Creating Boxes

```python
# Different box types
cubic = mp.Box.cubic(10.0)                    # 10×10×10 Å³
orth = mp.Box.orth([8.0, 10.0, 12.0])        # 8×10×12 Å³
tric = mp.Box.tric([10.0, 10.0, 10.0],       # Triclinic with tilts
                   [0.1, 0.0, 0.0])

# Box with specific PBC
box = mp.Box.cubic(10.0, pbc=[True, True, False])  # Periodic in X,Y only
```

### Creating Blocks

```python
# From dictionary
atoms = mp.Block({
    'x': [0.0, 1.0, 2.0],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C']
})

# From existing data
positions = np.array([[0, 0, 0], [1, 0, 0], [2, 0, 0]])
atoms = mp.Block({
    'x': positions[:, 0],
    'y': positions[:, 1],
    'z': positions[:, 2]
})
```

## Data Access

### Frame Access

```python
# Access blocks
atoms = frame['atoms']
bonds = frame['bonds']

# Access variables within blocks
x_coords = frame['atoms', 'x']           # Shorthand
x_coords = frame['atoms']['x']           # Explicit
positions = frame['atoms', ['x', 'y', 'z']]  # Multiple variables

# Access metadata
time = frame.metadata['time']
box = frame.box
```

### Block Access

```python
# Get variables
x = atoms['x']
positions = atoms[['x', 'y', 'z']]

# Get by index
first_atom = atoms[0]
first_two = atoms[:2]

# Get by selection
carbons = atoms[atoms['element'] == 'C']
```

## Atom Selection

### Basic Selection

```python
# Select by element
carbons = frame['atoms']['element'] == 'C'
oxygens = frame['atoms']['element'] == 'O'

# Select by position
in_box = (frame['atoms']['x'] > 0) & (frame['atoms']['x'] < 5)
in_sphere = np.sqrt(frame['atoms']['x']**2 + frame['atoms']['y']**2) < 3.0

# Select by type
type_1 = frame['atoms']['type'] == 1
```

### Advanced Selection

```python
# Combine selections
carbons_in_box = carbons & in_box

# Select specific atoms
specific_atoms = frame['atoms'][[0, 2, 5]]  # Atoms 0, 2, 5

# Select by multiple criteria
heavy_atoms = (frame['atoms']['element'] != 'H') & (frame['atoms']['mass'] > 10)
```

## Analysis Operations

### Basic Calculations

```python
# Center of mass
masses = frame['atoms']['mass']
positions = frame['atoms'][['x', 'y', 'z']]
com = np.average(positions, weights=masses, axis=0)

# Total mass
total_mass = np.sum(masses)

# Radius of gyration
distances = positions - com
rg_squared = np.average(np.sum(distances**2, axis=1), weights=masses)
rg = np.sqrt(rg_squared)
```

### Distance Calculations

```python
# Pairwise distances
positions = frame['atoms'][['x', 'y', 'z']]
distances = np.linalg.norm(positions[1:] - positions[:-1], axis=1)

# Distance from point
center = np.array([5.0, 5.0, 5.0])
distances_from_center = np.linalg.norm(positions - center, axis=1)
```

## File I/O

### Reading Files

```python
# Read different formats
pdb_frame = mp.read_pdb("protein.pdb")
lammps_frame = mp.read_lammps_data("data.lammps", "full")
amber_frame, amber_ff = mp.read_amber("protein.prmtop", "protein.inpcrd")

# Read trajectories
traj_reader = mp.read_lammps_trajectory("trajectory.lammpstrj")
```

### Writing Files

```python
# Write to different formats
mp.write_pdb("output.pdb", frame)
mp.write_lammps_data("output.data", frame)

# Save to HDF5
frame.save("system.h5")

# Load from HDF5
loaded_frame = mp.Frame.load("system.h5")
```

## System Operations

### Creating Systems

```python
# Basic system
system = mp.FrameSystem(frame=frame, box=box)

# System with force field
system = mp.FrameSystem(
    frame=frame,
    box=box,
    forcefield=forcefield
)
```

### Box Operations

```python
# Wrap coordinates
wrapped_positions = box.wrap(positions)

# Calculate distances with PBC
distances = box.dist_all(pos1, pos2)

# Transform coordinates
transformed = box.transform(transformation_matrix)
```

## Common Patterns

### Working with Trajectories

```python
# Process multiple frames
for frame in trajectory:
    # Your analysis here
    com = calculate_center_of_mass(frame)
    results.append(com)

# Process in chunks
for i in range(0, len(trajectory), 100):
    chunk = trajectory[i:i+100]
    process_chunk(chunk)
```

### Batch Processing

```python
# Process multiple files
files = ["file1.pdb", "file2.pdb", "file3.pdb"]
results = []

for file in files:
    frame = mp.read_pdb(file)
    result = analyze_frame(frame)
    results.append(result)
```

## Quick Tips

- **Use vectorized operations** instead of Python loops
- **Check data types** with `frame['atoms'].dtypes`
- **Access metadata** with `frame.metadata`
- **Combine selections** with `&` (and) and `|` (or)
- **Save intermediate results** to avoid recomputation
- **Use lazy loading** for large trajectories

## Common Errors

```python
# ❌ Wrong: mixing data types
atoms_data = {'x': [0, 1, 2], 'element': ['C', 'C', 'C']}

# ✅ Correct: consistent data types
atoms_data = {'x': [0.0, 1.0, 2.0], 'element': ['C', 'C', 'C']}

# ❌ Wrong: accessing non-existent variables
mass = frame['atoms']['mass']  # Error if 'mass' doesn't exist

# ✅ Correct: check first
if 'mass' in frame['atoms'].variables():
    mass = frame['atoms']['mass']
```
