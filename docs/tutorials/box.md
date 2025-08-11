# Box & Simulation Boundaries

Simulation boxes define the geometric boundaries and periodic conditions for molecular systems. This tutorial teaches you how to work with MolPy's Box class to create, manipulate, and understand simulation boundaries.

## Understanding Simulation Boxes

Simulation boxes define the spatial domain where your molecular system exists. They can be free (no boundaries), orthogonal (rectangular with right angles), or triclinic (parallelepiped with arbitrary angles). The box also controls periodic boundary conditions, which determine how atoms interact across boundaries.

MolPy's Box class provides a comprehensive interface for managing simulation boundaries. It handles different box types automatically, manages periodic boundary conditions, and provides methods for coordinate wrapping and distance calculations. The box is essential for molecular dynamics simulations, crystal structures, and any system where spatial boundaries matter.

Boxes are fundamental to understanding how your molecular system behaves in space. They affect everything from atom interactions to system properties like density and pressure.

## Creating Simulation Boxes

### Basic Box Types

Box provides several convenient constructors for common box types. These make it easy to create appropriate simulation boundaries for different systems.

```python
import molpy as mp
import numpy as np

# Create a cubic box (equal dimensions in all directions)
cubic_box = mp.Box.cubic(
    length=20.0,           # 20 Å in each direction
    pbc=[True, True, True] # Periodic in all directions
)

# Create an orthogonal box (different lengths, right angles)
orth_box = mp.Box.orth(
    lengths=[15.0, 20.0, 25.0],  # Different dimensions
    pbc=[True, True, False]       # Periodic in x, y; free in z
)

# Create a triclinic box (arbitrary angles)
tric_box = mp.Box.tric(
    lengths=[20.0, 20.0, 20.0],  # Box dimensions
    tilts=[0.0, 0.0, 5.0],       # Tilt angles
    pbc=[True, True, True]        # Periodic in all directions
)

# Check box properties
print(f"Cubic box: {cubic_box}")
print(f"Orthogonal box: {orth_box}")
print(f"Triclinic box: {tric_box}")
```

Each box type serves different purposes. Cubic boxes are common for simple liquids and gases, orthogonal boxes are useful for interfaces and membranes, and triclinic boxes are essential for crystal structures with non-orthogonal unit cells.

### Box Properties and Information

Boxes provide various properties that help you understand their geometry and behavior.

```python
def analyze_box(box, name=""):
    """Analyze box properties and characteristics."""
    print(f"=== {name} Analysis ===")

    # Basic properties
    print(f"Style: {box.style}")
    print(f"Volume: {box.volume:.1f} Å³")
    print(f"Lengths: {box.lengths}")
    print(f"Origin: {box.origin}")

    # Periodic boundary conditions
    print(f"PBC: {box.pbc}")
    print(f"Periodic: {box.is_periodic}")

    # Box boundaries
    print(f"Bounds: {box.bounds}")
    print(f"X range: [{box.xlo:.1f}, {box.xhi:.1f}]")
    print(f"Y range: [{box.ylo:.1f}, {box.yhi:.1f}]")
    print(f"Z range: [{box.zlo:.1f}, {box.zhi:.1f}]")

    # For triclinic boxes
    if box.style == mp.Box.Style.TRICLINIC:
        print(f"Tilts: {box.tilts}")
        print(f"Angles: {box.angles}")

# Analyze different box types
analyze_box(cubic_box, "Cubic")
analyze_box(orth_box, "Orthogonal")
analyze_box(tric_box, "Triclinic")
```

This analysis reveals the key characteristics of each box type. Understanding these properties is essential for setting up simulations and interpreting results.

## Working with Box Geometry

### Box Transformations

Boxes can be transformed and modified to suit different simulation needs.

```python
# Scale a box
scaled_box = cubic_box * 2.0
print(f"Original volume: {cubic_box.volume:.1f} Å³")
print(f"Scaled volume: {scaled_box.volume:.1f} Å³")

# Transform box with a matrix
transformation = np.array([
    [2.0, 0.0, 0.0],
    [0.0, 1.5, 0.0],
    [0.0, 0.0, 1.0]
])

transformed_box = cubic_box.transform(transformation)
print(f"Transformed box: {transformed_box}")
print(f"New volume: {transformed_box.volume:.1f} Å³")

# Merge boxes (useful for combining systems)
merged_box = cubic_box.merge(orth_box)
print(f"Merged box: {merged_box}")
```

These transformations are useful for adapting boxes to different simulation conditions, combining systems, or creating supercells.

### Coordinate Operations

Boxes provide methods for working with atomic coordinates, including wrapping and distance calculations.

```python
# Create some atomic coordinates
positions = np.array([
    [0.0, 0.0, 0.0],      # Center of box
    [25.0, 0.0, 0.0],     # Outside box (x direction)
    [0.0, 25.0, 0.0],     # Outside box (y direction)
    [0.0, 0.0, 25.0]      # Outside box (z direction)
])

# Wrap coordinates to box (periodic boundary conditions)
wrapped_positions = cubic_box.wrap(positions)
print("Original positions:")
print(positions)
print("Wrapped positions:")
print(wrapped_positions)

# Calculate distances with periodic boundary conditions
distances = cubic_box.dist_all(positions[0], positions[1:])
print(f"Distances from center: {distances}")

# Convert to fractional coordinates
fractional = cubic_box.make_fractional(positions)
print("Fractional coordinates:")
print(fractional)
```

Coordinate wrapping is essential for periodic systems, ensuring that atoms that move outside the box reappear on the opposite side. Distance calculations with periodic boundary conditions give you the shortest distance between atoms, accounting for box periodicity.

## Periodic Boundary Conditions

### Understanding PBC

Periodic boundary conditions determine how your system behaves at boundaries. They're essential for simulating bulk materials and avoiding surface effects.

```python
# Create boxes with different PBC settings
all_periodic = mp.Box.cubic(20.0, pbc=[True, True, True])
xy_periodic = mp.Box.cubic(20.0, pbc=[True, True, False])
free_box = mp.Box.cubic(20.0, pbc=[False, False, False])

print(f"All periodic: {all_periodic.is_periodic}")
print(f"XY periodic: {xy_periodic.is_periodic}")
print(f"Free box: {free_box.is_periodic}")

# Check individual directions
print(f"X periodic: {all_periodic.periodic_x}")
print(f"Y periodic: {all_periodic.periodic_y}")
print(f"Z periodic: {all_periodic.periodic_z}")

# Modify PBC settings
all_periodic.periodic_z = False
print(f"Modified PBC: {all_periodic.pbc}")
```

Different PBC configurations are useful for different simulation types. All periodic is common for bulk materials, XY periodic is useful for membranes and interfaces, and free boundaries are used for isolated molecules.

### PBC-Aware Calculations

When working with periodic systems, you need to account for boundary conditions in your calculations.

```python
# Create a periodic box
periodic_box = mp.Box.cubic(10.0, pbc=[True, True, True])

# Test coordinates that cross boundaries
test_positions = np.array([
    [1.0, 1.0, 1.0],      # Inside box
    [9.5, 1.0, 1.0],      # Near boundary
    [0.5, 1.0, 1.0]       # Near opposite boundary
])

# Calculate distances with PBC
for i in range(1, len(test_positions)):
    distance = periodic_box.dist(test_positions[0], test_positions[i])
    print(f"Distance 0 to {i}: {distance:.3f} Å")

# Check if positions are in box
in_box = periodic_box.isin(test_positions)
print(f"Positions in box: {in_box}")
```

PBC-aware calculations are crucial for accurate simulations. They ensure that interactions across boundaries are handled correctly and that system properties are calculated properly.

## Box Applications

### Crystal Structure Boxes

Boxes are essential for crystal structures, where the unit cell defines the repeating pattern.

```python
# Create a box for a simple cubic crystal
crystal_box = mp.Box.cubic(5.0, pbc=[True, True, True])

# Create atomic positions for a simple cubic lattice
lattice_positions = []
for i in range(2):
    for j in range(2):
        for k in range(2):
            pos = [i * 2.5, j * 2.5, k * 2.5]
            lattice_positions.append(pos)

lattice_positions = np.array(lattice_positions)

# Wrap positions to box
wrapped_lattice = crystal_box.wrap(lattice_positions)

print(f"Crystal box: {crystal_box}")
print(f"Lattice positions: {len(lattice_positions)}")
print(f"Box volume: {crystal_box.volume:.1f} Å³")
```

This example shows how boxes define the unit cell for crystal structures. The periodic boundary conditions ensure that the crystal pattern repeats infinitely in all directions.

### Interface and Surface Boxes

For systems with interfaces or surfaces, you need boxes that are periodic in some directions but free in others.

```python
# Create a box for a water-vapor interface
interface_box = mp.Box.orth(
    lengths=[20.0, 20.0, 40.0],  # Taller in z for interface
    pbc=[True, True, False]       # Periodic in x,y; free in z
)

# Add vacuum regions
interface_box.origin = [0.0, 0.0, -20.0]  # Center the box

print(f"Interface box: {interface_box}")
print(f"Z range: [{interface_box.zlo:.1f}, {interface_box.zhi:.1f}]")
print(f"PBC: {interface_box.pbc}")
```

Interface boxes are essential for studying surfaces, membranes, and other systems where you need to maintain periodicity in some directions while allowing free boundaries in others.

## Box Best Practices

### Box Sizing

Always ensure your box is appropriately sized for your molecular system. Too small boxes can cause artifacts, while too large boxes waste computational resources.

```python
def recommend_box_size(molecular_size, buffer=5.0):
    """Recommend appropriate box size for a molecular system."""
    recommended_size = molecular_size + 2 * buffer
    return mp.Box.cubic(recommended_size)

# Example: water molecule with ~3 Å diameter
water_size = 3.0
water_box = recommend_box_size(water_size, buffer=8.0)
print(f"Recommended water box: {water_box}")
print(f"Box volume: {water_box.volume:.1f} Å³")
```

### PBC Selection

Choose periodic boundary conditions based on your simulation goals. Consider the physical meaning of your system and what you're trying to simulate.

```python
def create_appropriate_box(system_type, dimensions):
    """Create appropriate box for different system types."""
    if system_type == "bulk":
        # Bulk materials need full periodicity
        return mp.Box.cubic(dimensions, pbc=[True, True, True])
    elif system_type == "interface":
        # Interfaces need periodicity in plane only
        return mp.Box.orth(dimensions, pbc=[True, True, False])
    elif system_type == "isolated":
        # Isolated molecules need no periodicity
        return mp.Box.cubic(dimensions, pbc=[False, False, False])
    else:
        raise ValueError(f"Unknown system type: {system_type}")

# Create boxes for different systems
bulk_box = create_appropriate_box("bulk", 20.0)
interface_box = create_appropriate_box("interface", [20.0, 20.0, 40.0])
isolated_box = create_appropriate_box("isolated", 30.0)

print(f"Bulk box PBC: {bulk_box.pbc}")
print(f"Interface box PBC: {interface_box.pbc}")
print(f"Isolated box PBC: {isolated_box.pbc}")
```

## Summary

This tutorial covered the fundamental concepts of MolPy's Box class. You learned how to create different types of simulation boxes, manage periodic boundary conditions, and perform coordinate operations with proper boundary handling.

Simulation boxes are essential for defining the spatial domain of your molecular system. They control how atoms interact across boundaries, affect system properties, and determine the computational approach needed for your simulation. Understanding box behavior is crucial for setting up accurate and efficient molecular simulations.

### Next Steps

Continue your MolPy journey by learning about trajectory analysis, understanding region-based operations, exploring selection and filtering, and mastering system organization for complex assemblies.

Understanding simulation boxes is essential for molecular dynamics and structural analysis in MolPy. The Box class provides the foundation for proper spatial boundary management and periodic system handling.
