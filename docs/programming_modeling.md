# Programming-based Material Modeling

This guide shows how to build complex molecular systems using MolPy's core capabilities in a programmatic way.

## Core Concepts

MolPy provides fundamental building blocks:
- **`copy()`** - Create an independent copy of a molecule
- **`move(delta)`** - Translate coordinates (returns self for chaining)
- **`rotate(axis, angle, about=None)`** - Rotate around an axis (returns self for chaining)
- **`scale(factor, about=None)`** - Scale coordinates (returns self for chaining)
- **`merge(other)`** - Deep copy another system into this one
- **`+=`** - Merge operator (in-place)
- **`+`** - Combination operator (new object)
- **`replicate(n, transform=None)`** - Create n copies with optional transformation

These primitives allow you to **programmatically construct any structure**.

## Basic Patterns

### 1. Simple Molecular Templates

Define reusable molecular structures:

```python
from molpy import Atomistic

class Water(Atomistic):
    def __post_init__(self, **props):
        o = self.add_atom(symbol="O", xyz=[0, 0, 0])
        h1 = self.add_atom(symbol="H", xyz=[0.757, 0.586, 0])
        h2 = self.add_atom(symbol="H", xyz=[-0.757, 0.586, 0])
        self.add_bond(o, h1)
        self.add_bond(o, h2)

class Methane(Atomistic):
    def __post_init__(self, **props):
        c = self.add_atom(symbol="C", xyz=[0, 0, 0])
        h1 = self.add_atom(symbol="H", xyz=[0.629, 0.629, 0.629])
        h2 = self.add_atom(symbol="H", xyz=[-0.629, -0.629, 0.629])
        h3 = self.add_atom(symbol="H", xyz=[-0.629, 0.629, -0.629])
        h4 = self.add_atom(symbol="H", xyz=[0.629, -0.629, -0.629])
        for h in [h1, h2, h3, h4]:
            self.add_bond(c, h)
```

### 2. Building Systems with `+=`

Merge molecules into a system:

```python
# Create an empty system
system = Atomistic()

# Add molecules at different positions
system += Water()
system += Water().move([5, 0, 0])
system += Water().move([10, 0, 0])
system += Methane().move([2.5, 5, 0])

print(system)  # <Atomistic, 17 atoms (C:1 H:7 O:3), 9 bonds, with coords>
```

### 3. Chaining Transformations

Operations return `self` for fluent interface:

```python
# Move and rotate in one line
system += Water().move([5, 0, 0]).rotate([0, 0, 1], 45)

# Complex transformations
system += Methane() \
    .move([3, 0, 0]) \
    .rotate([1, 0, 0], 90) \
    .scale(1.5)
```

## Recipes for Common Patterns

### Linear Array

Create molecules in a line:

```python
system = Atomistic()
for i in range(10):
    system += Water().move([i * 5, 0, 0])
```

Or more concise with `replicate()`:

```python
system = Water().replicate(10, lambda mol, i: mol.move([i * 5, 0, 0]))
```

### 2D Grid

```python
system = Atomistic()
nx, ny = 5, 5
spacing = 5.0

for i in range(nx):
    for j in range(ny):
        system += Water().move([i * spacing, j * spacing, 0])
```

With `replicate()`:

```python
n = nx * ny
system = Water().replicate(
    n, 
    lambda mol, idx: mol.move([
        (idx % nx) * spacing,
        (idx // nx) * spacing,
        0
    ])
)
```

### 3D Lattice

```python
system = Atomistic()
nx, ny, nz = 3, 3, 3
spacing = 5.0

for i in range(nx):
    for j in range(ny):
        for k in range(nz):
            system += Water().move([i * spacing, j * spacing, k * spacing])
```

### Random Selection

Mix different molecules with controlled ratios:

```python
import random

molecules = [Water, Methane, Benzene]
weights = [0.5, 0.3, 0.2]  # 50% water, 30% methane, 20% benzene

system = Atomistic()
for i in range(100):
    MolClass = random.choices(molecules, weights=weights)[0]
    system += MolClass().move([
        random.uniform(0, 50),
        random.uniform(0, 50),
        random.uniform(0, 50)
    ])
```

### Cumulative Array (MolTemplate-style)

Each molecule transformed relative to the previous one:

```python
system = Atomistic()
current = Water()
system += current

for i in range(9):  # Add 9 more
    current = current.copy().move([5, 0, 0])  # Cumulative
    system += current
```

### Linear Polymer Chain

Connect monomers end-to-end:

```python
class Monomer(Atomistic):
    def __post_init__(self, **props):
        # Simple CH2 unit
        c1 = self.add_atom(symbol="C", xyz=[0, 0, 0], name="head")
        c2 = self.add_atom(symbol="C", xyz=[1.54, 0, 0], name="tail")
        h1 = self.add_atom(symbol="H", xyz=[0, 1, 0])
        h2 = self.add_atom(symbol="H", xyz=[0, -1, 0])
        h3 = self.add_atom(symbol="H", xyz=[1.54, 1, 0])
        h4 = self.add_atom(symbol="H", xyz=[1.54, -1, 0])
        self.add_bond(c1, c2)
        for h in [h1, h2]:
            self.add_bond(c1, h)
        for h in [h3, h4]:
            self.add_bond(c2, h)

# Build chain
chain = Atomistic()
for i in range(10):
    chain += Monomer().move([i * 1.54, 0, 0])
```

### Helical Polymer

```python
import math

system = Atomistic()
n_units = 20
radius = 10.0
pitch = 34.0  # Rise per turn

for i in range(n_units):
    angle = (i / n_units) * 2 * math.pi * 3  # 3 turns
    z = (i / n_units) * pitch * 3
    
    system += Monomer().move([
        radius * math.cos(angle),
        radius * math.sin(angle),
        z
    ]).rotate([0, 0, 1], math.degrees(angle))
```

### Ring Structure

```python
import math

n_molecules = 8
radius = 10.0

system = Atomistic()
for i in range(n_molecules):
    angle = (i / n_molecules) * 2 * math.pi
    
    system += Water().move([
        radius * math.cos(angle),
        radius * math.sin(angle),
        0
    ]).rotate([0, 0, 1], math.degrees(angle))
```

### Sphere Surface

```python
import math
import random

def fibonacci_sphere(n):
    """Generate ~uniformly distributed points on sphere"""
    points = []
    phi = math.pi * (3 - math.sqrt(5))  # Golden angle
    
    for i in range(n):
        y = 1 - (i / (n - 1)) * 2  # y from 1 to -1
        r = math.sqrt(1 - y * y)
        theta = phi * i
        x = math.cos(theta) * r
        z = math.sin(theta) * r
        points.append([x, y, z])
    
    return points

# Place molecules on sphere
radius = 20.0
n = 100

system = Atomistic()
for point in fibonacci_sphere(n):
    system += Water().move([
        point[0] * radius,
        point[1] * radius,
        point[2] * radius
    ])
```

### Bilayer Structure

```python
# Top layer
top = Atomistic()
for i in range(10):
    for j in range(10):
        top += Molecule().move([i * 5, j * 5, 5])

# Bottom layer (inverted)
bottom = Atomistic()
for i in range(10):
    for j in range(10):
        bottom += Molecule().move([i * 5, j * 5, -5]).rotate([1, 0, 0], 180)

# Combine
bilayer = top + bottom
```

### Random Perturbations

Add noise to perfect structures:

```python
import random

# Start with perfect grid
system = Water().replicate(25, lambda mol, i: mol.move([
    (i % 5) * 5,
    (i // 5) * 5,
    0
]))

# Add random noise
for atom in system.atoms:
    pos = atom.get('xyz', [0, 0, 0])
    atom['xyz'] = [
        pos[0] + random.uniform(-0.5, 0.5),
        pos[1] + random.uniform(-0.5, 0.5),
        pos[2] + random.uniform(-0.5, 0.5)
    ]
```

### Random Rotations

```python
import random

system = Atomistic()
for i in range(20):
    # Random axis
    axis = [random.random(), random.random(), random.random()]
    angle = random.uniform(0, 360)
    
    system += Water() \
        .rotate(axis, angle) \
        .move([i * 5, 0, 0])
```

## Advanced Patterns

### Hierarchical Assembly

Build complex structures from sub-assemblies:

```python
# Create a cluster
def water_cluster(n=4):
    """Create a cluster of n water molecules"""
    cluster = Atomistic()
    for i in range(n):
        angle = (i / n) * 360
        cluster += Water().move([3, 0, 0]).rotate([0, 0, 1], angle)
    return cluster

# Use clusters as building blocks
system = Atomistic()
for i in range(5):
    system += water_cluster(4).move([i * 10, 0, 0])
```

### Parameterized Structures

```python
def create_tube(molecule_class, radius, height, n_rings, n_per_ring):
    """Create a tubular structure"""
    tube = Atomistic()
    
    for ring_idx in range(n_rings):
        z = (ring_idx / (n_rings - 1)) * height
        
        for i in range(n_per_ring):
            angle = (i / n_per_ring) * 360
            x = radius * math.cos(math.radians(angle))
            y = radius * math.sin(math.radians(angle))
            
            tube += molecule_class() \
                .move([x, y, z]) \
                .rotate([0, 0, 1], angle)
    
    return tube

# Use it
nanotube = create_tube(Methane, radius=10, height=50, n_rings=10, n_per_ring=8)
```

### Conditional Assembly

```python
def is_valid_position(pos, existing_positions, min_dist=3.0):
    """Check if position is far enough from existing ones"""
    for ex_pos in existing_positions:
        dist = sum((a - b)**2 for a, b in zip(pos, ex_pos)) ** 0.5
        if dist < min_dist:
            return False
    return True

# Add molecules with minimum distance constraint
import random

system = Atomistic()
positions = []

while len(positions) < 100:
    pos = [random.uniform(0, 50) for _ in range(3)]
    if is_valid_position(pos, positions, min_dist=4.0):
        system += Water().move(pos)
        positions.append(pos)
```

### Combining Operations

```python
# Create a multi-component system
system = Atomistic()

# Water layer at bottom
water_layer = Water().replicate(25, lambda mol, i: mol.move([
    (i % 5) * 5,
    (i // 5) * 5,
    0
]))

# Methane layer in middle
methane_layer = Methane().replicate(16, lambda mol, i: mol.move([
    (i % 4) * 6,
    (i // 4) * 6,
    10
]))

# Ring of benzene on top
benzene_ring = Atomistic()
for i in range(8):
    angle = (i / 8) * 360
    benzene_ring += Benzene() \
        .move([15, 0, 20]) \
        .rotate([0, 0, 1], angle) \
        .move([10 * math.cos(math.radians(angle)), 10 * math.sin(math.radians(angle)), 0])

# Combine all
system = water_layer + methane_layer + benzene_ring
```

## Key Principles

1. **Molecules are templates**: Define once, instantiate many times
2. **Transformations return self**: Enable method chaining
3. **`copy()` before modifying**: If you need to preserve the original
4. **Use `+=` to build systems**: Clean and readable
5. **`replicate()` for repetitive patterns**: More concise than loops
6. **Python is your friend**: Use loops, functions, conditionals, random, math

## Philosophy

MolPy doesn't provide pre-built functions for every possible structure. Instead, it gives you **primitive operations** that you can **compose programmatically** to build anything.

This approach is:
- **More flexible**: Not limited to predefined patterns
- **More powerful**: Combine operations in creative ways
- **More maintainable**: Clear, readable code
- **More Pythonic**: Use standard Python constructs

If you find yourself repeating a pattern, **write a function**. That's the MolPy way.
