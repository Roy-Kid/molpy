# How to Pack Molecular Systems

This guide shows practical examples of using MolPy's packing tools to create molecular assemblies, mixtures, and complex molecular systems. You'll learn how to pack molecules into simulation boxes with various constraints.

## Understanding Molecular Packing

Molecular packing involves placing molecules in space to create realistic molecular assemblies. This is essential for creating initial configurations for simulations, studying molecular mixtures, and building complex molecular systems.

## Basic Packing Workflows

### Creating a Packing Session

Start by setting up a packing session with your working directory:

```python
import molpy as mp
from pathlib import Path
import numpy as np

# Create working directory
workdir = Path("packing_work")
workdir.mkdir(exist_ok=True)

# Create packing session
session = mp.pack.Session(workdir, packer="packmol")
print(f"Created packing session in {workdir}")
```

### Adding Molecules to Pack

Add different molecular species with their target numbers and constraints:

```python
# Create water molecule
water_atoms = mp.Block({
    'x': [0.0, 0.9572, -0.2400],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['O', 'H', 'H']
})

water_frame = mp.Frame({'atoms': water_atoms})

# Create methane molecule
methane_atoms = mp.Block({
    'x': [0.0, 0.0, 0.0, 0.0, 0.0],
    'y': [0.0, 0.0, 0.0, 0.0, 0.0],
    'z': [0.0, 1.09, -1.09, 0.0, 0.0],
    'element': ['C', 'H', 'H', 'H', 'H']
})

methane_frame = mp.Frame({'atoms': methane_atoms})

# Add molecules to packing session
water_target = session.add_target(water_frame, number=100, constraint="inside box 0 0 0 20 20 20")
methane_target = session.add_target(methane_frame, number=50, constraint="inside box 0 0 0 20 20 20")

print(f"Added {water_target.number} water molecules")
print(f"Added {methane_target.number} methane molecules")
```

## Packing Constraints

### Box Constraints

Define spatial regions where molecules can be placed:

```python
# Pack molecules in a cubic box
cubic_constraint = "inside box 0 0 0 15 15 15"
water_in_cube = session.add_target(water_frame, number=50, constraint=cubic_constraint)

# Pack in a rectangular box
rect_constraint = "inside box 0 0 0 20 15 10"
methane_in_rect = session.add_target(methane_frame, number=30, constraint=rect_constraint)

print("Added molecules with box constraints")
```

### Distance Constraints

Control how close molecules can be to each other:

```python
# Pack with minimum distance between molecules
distance_constraint = "inside box 0 0 0 20 20 20 outside sphere 0 0 0 5"
water_with_distance = session.add_target(water_frame, number=25, constraint=distance_constraint)

print("Added water molecules with distance constraint")
```

### Complex Constraints

Combine multiple constraints for sophisticated packing:

```python
# Pack water in outer shell, methane in center
outer_constraint = "inside box 0 0 0 20 20 20 outside sphere 10 10 10 8"
inner_constraint = "inside sphere 10 10 10 6"

water_outer = session.add_target(water_frame, number=75, constraint=outer_constraint)
methane_inner = session.add_target(methane_frame, number=20, constraint=inner_constraint)

print("Added molecules with complex spatial constraints")
```

## Running Packing Simulations

### Basic Packing

Execute the packing simulation with basic parameters:

```python
# Run packing simulation
result = session.optimize(max_steps=1000, seed=42)

if result.success:
    print("Packing successful!")
    print(f"Final system has {len(result.frame['atoms'])} atoms")
else:
    print("Packing failed")
    print(f"Error: {result.error}")
```

### Advanced Packing

Use more sophisticated packing parameters:

```python
# Run with more steps and different seed
result = session.optimize(max_steps=5000, seed=12345)

if result.success:
    print("Advanced packing successful!")

    # Analyze the packed system
    packed_frame = result.frame
    print(f"Packed system: {packed_frame['atoms'].nrows} atoms")

    # Check density
    box_volume = 20 * 20 * 20  # Assuming 20x20x20 box
    total_molecules = 100 + 50  # water + methane
    density = total_molecules / box_volume
    print(f"Molecular density: {density:.3f} molecules/Å³")
else:
    print("Advanced packing failed")
```

## Packing Strategies

### Layered Packing

Create layered molecular systems:

```python
def create_layered_system():
    """Create a layered molecular system."""

    # Create new session for layered packing
    layered_session = mp.pack.Session(Path("layered_packing"), packer="packmol")

    # Bottom layer (water)
    bottom_constraint = "inside box 0 0 0 20 20 5"
    bottom_water = layered_session.add_target(water_frame, number=40, constraint=bottom_constraint)

    # Middle layer (methane)
    middle_constraint = "inside box 0 0 5 20 20 10"
    middle_methane = layered_session.add_target(methane_frame, number=25, constraint=middle_constraint)

    # Top layer (water)
    top_constraint = "inside box 0 0 10 20 20 20"
    top_water = layered_session.add_target(water_frame, number=40, constraint=top_constraint)

    return layered_session

# Create and run layered packing
layered_session = create_layered_system()
layered_result = layered_session.optimize(max_steps=2000)

if layered_result.success:
    print("Layered packing successful!")
    print(f"Layered system: {layered_result.frame['atoms'].nrows} atoms")
```

### Concentration Gradients

Create systems with concentration gradients:

```python
def create_concentration_gradient():
    """Create a system with concentration gradient."""

    gradient_session = mp.pack.Session(Path("gradient_packing"), packer="packmol")

    # Create gradient by varying molecule density along x-axis
    for i in range(5):
        x_start = i * 4
        x_end = (i + 1) * 4
        n_molecules = 20 - i * 3  # Decreasing concentration

        constraint = f"inside box {x_start} 0 0 {x_end} 20 20"
        gradient_session.add_target(water_frame, number=n_molecules, constraint=constraint)

    return gradient_session

# Create and run gradient packing
gradient_session = create_concentration_gradient()
gradient_result = gradient_session.optimize(max_steps=1500)

if gradient_result.success:
    print("Gradient packing successful!")
```

## Packing Analysis

### Analyzing Packed Systems

Examine the results of your packing simulations:

```python
def analyze_packed_system(packed_frame):
    """Analyze a packed molecular system."""

    atoms = packed_frame['atoms']
    elements = atoms['element']

    # Count different molecule types
    element_counts = {}
    for element in elements:
        element_counts[element] = element_counts.get(element, 0) + 1

    print("Packed system analysis:")
    print(f"Total atoms: {atoms.nrows}")
    print(f"Element composition: {element_counts}")

    # Calculate molecular density
    positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

    # Estimate box dimensions
    box_x = max(positions[:, 0]) - min(positions[:, 0])
    box_y = max(positions[:, 1]) - min(positions[:, 1])
    box_z = max(positions[:, 2]) - min(positions[:, 2])

    box_volume = box_x * box_y * box_z
    print(f"Estimated box dimensions: {box_x:.1f} x {box_y:.1f} x {box_z:.1f} Å")
    print(f"Estimated box volume: {box_volume:.1f} Å³")

    # Calculate density
    if 'O' in element_counts:
        n_water = element_counts['O']  # Each water has one O
        water_density = n_water / box_volume
        print(f"Water density: {water_density:.3f} molecules/Å³")

    return element_counts, box_volume

# Analyze packed system
if result.success:
    element_counts, box_volume = analyze_packed_system(result.frame)
```

### Quality Assessment

Assess the quality of your packed systems:

```python
def assess_packing_quality(packed_frame):
    """Assess the quality of a packed molecular system."""

    atoms = packed_frame['atoms']
    positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])

    # Check for overlapping atoms (simplified)
    min_distance = 1.0  # Minimum allowed distance in Å
    overlaps = 0

    for i in range(len(positions)):
        for j in range(i+1, len(positions)):
            distance = np.linalg.norm(positions[i] - positions[j])
            if distance < min_distance:
                overlaps += 1

    print("Packing quality assessment:")
    print(f"Total atom pairs: {len(positions) * (len(positions) - 1) // 2}")
    print(f"Overlapping pairs: {overlaps}")

    if overlaps == 0:
        print("✓ No overlapping atoms detected")
    else:
        print(f"⚠ {overlaps} overlapping atom pairs detected")

    # Check spatial distribution
    x_coords = positions[:, 0]
    y_coords = positions[:, 1]
    z_coords = positions[:, 2]

    print(f"X range: {x_coords.min():.2f} to {x_coords.max():.2f} Å")
    print(f"Y range: {y_coords.min():.2f} to {y_coords.max():.2f} Å")
    print(f"Z range: {z_coords.min():.2f} to {z_coords.max():.2f} Å")

    return overlaps == 0

# Assess packing quality
if result.success:
    is_good = assess_packing_quality(result.frame)
    if is_good:
        print("Packing quality is acceptable")
    else:
        print("Packing quality needs improvement")
```

## Practical Applications

### Creating Simulation Boxes

Use packing to create initial configurations for molecular dynamics:

```python
def create_simulation_box():
    """Create a simulation box with packed molecules."""

    # Create packing session
    sim_session = mp.pack.Session(Path("simulation_box"), packer="packmol")

    # Add molecules for a realistic liquid simulation
    water_constraint = "inside box 0 0 0 30 30 30"
    sim_session.add_target(water_frame, number=200, constraint=water_constraint)

    # Add some ions if needed
    # ion_constraint = "inside box 0 0 0 30 30 30"
    # sim_session.add_target(ion_frame, number=10, constraint=ion_constraint)

    # Run packing
    sim_result = sim_session.optimize(max_steps=3000)

    if sim_result.success:
        # Create simulation system
        sim_system = mp.FrameSystem(
            frame=sim_result.frame,
            box=mp.Box.cubic(30.0)
        )

        # Save for simulation
        mp.write_lammps(Path("simulation_box"), sim_system)
        print("Simulation box created and saved")

        return sim_system

    return None

# Create simulation box
sim_box = create_simulation_box()
```

## Summary

This guide covered practical molecular packing workflows:

- Set up packing sessions with different packers
- Add molecules with various constraints
- Create complex spatial arrangements
- Run packing simulations
- Analyze packed systems
- Assess packing quality
- Create simulation boxes

MolPy's packing tools provide flexible ways to create molecular assemblies and mixtures. The key is understanding how to define appropriate constraints and analyze the results effectively.

### Next Steps

Continue exploring MolPy by learning about analysis techniques, advanced simulation workflows, and integration with other molecular modeling tools.
