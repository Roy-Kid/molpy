# How to Use Packmol API

This guide explains how to use MolPy's Packmol API to automate molecular packing workflows.

## Overview

The Packmol API provides a clean interface for molecular packing using the Packmol software:
- **Packmol**: Main packing component for creating molecular systems
- **Target**: Defines molecules and their packing constraints
- **Constraints**: Spatial constraints for packing (boxes, spheres, etc.)

All components inherit from appropriate base classes, providing consistent working directory management and executable configuration.

## Why This API Design?

The new API addresses several key problems with traditional Packmol workflows:

1. **Complex Input Generation**: Instead of manually writing Packmol input files, you define targets and constraints
2. **File Management**: Automatic working directory creation and file organization
3. **Constraint Handling**: Easy definition of complex spatial constraints
4. **Integration**: Seamless integration with MolPy data structures
5. **Workflow Control**: Generate input files, run packing, or do both in one call

## Basic Setup

```python
import molpy as mp
import molpy.pack as mpk
from molpy.pack.packer.packmol import Packmol

# Create Packmol instance
packmol = Packmol(executable="packmol")
```

## Understanding the Workflow

The Packmol workflow follows this sequence:

```
Molecular Targets → Spatial Constraints → Packmol Input → Packing → Final System
```

Each step transforms the input:

1. **Targets** define what molecules to pack and how many
2. **Constraints** define where molecules can be placed
3. **Packmol Input** is automatically generated from targets and constraints
4. **Packing** creates the final molecular system
5. **Final System** is returned as a MolPy Frame

## Component-by-Component Guide

### Packmol Component

**What it does**: Orchestrates the entire molecular packing process using Packmol.

**When to use**: When you need to pack multiple molecules into a defined space.

```python
packmol = Packmol(executable="packmol")

# Basic usage with stored targets
result = packmol(max_steps=1000, seed=4628)

# With custom targets
result = packmol(targets=my_targets, max_steps=1000, seed=4628)

# With custom working directory
result = packmol(
    targets=my_targets,
    max_steps=1000,
    seed=4628,
    workdir="/scratch/packing_work"
)
```

**Key Parameters Explained**:
- `targets`: List of packing targets (molecules + constraints)
- `max_steps`: Maximum optimization steps for Packmol
- `seed`: Random seed for reproducible results
- `workdir`: Optional custom working directory

**What you get**:
- Packed molecular system as a MolPy Frame
- Optimized coordinates from Packmol
- Proper atom indexing and molecular IDs

**Additional Methods**:
```python
# Generate input file without running packing
input_file = packmol.generate_input_only(
    targets=my_targets,
    max_steps=1000,
    seed=4628
)

# Read existing Packmol output
frame = packmol.read_packmol_output("/path/to/packmol/workdir")
```

### Target Component

**What it does**: Defines a molecule and how many copies to pack.

**When to use**: Always needed to specify what molecules to pack.

```python
from molpy.pack.target import Target

# Create a target
target = Target(
    frame=molecule_frame,      # MolPy Frame of the molecule
    number=100,                # Number of copies to pack
    constraint=spatial_constraint,  # Where to pack them
    name="water_molecules"     # Descriptive name
)

# Or use the convenience method
target = packmol.def_target(
    frame=molecule_frame,
    number=100,
    constraint=spatial_constraint,
    name="water_molecules"
)
```

**Key Parameters Explained**:
- `frame`: MolPy Frame containing the molecule structure
- `number`: How many copies of this molecule to pack
- `constraint`: Spatial constraint defining where molecules can be placed
- `name`: Descriptive name for the target

### Constraint System

**What it does**: Defines spatial regions where molecules can be packed.

**When to use**: To control the spatial distribution of packed molecules.

```python
import molpy.pack.constraint as mpc

# Box constraint - pack inside a rectangular box
box_constraint = mpc.InsideBoxConstraint(
    region=mp.Box.orth([20.0, 20.0, 20.0])  # 20Å cubic box
)

# Sphere constraint - pack inside a sphere
sphere_constraint = mpc.InsideSphereConstraint(
    region=mp.Box.sphere(origin=[0, 0, 0], radius=10.0)  # 10Å radius sphere
)

# Outside constraint - pack outside a region
outside_constraint = mpc.OutsideBoxConstraint(
    region=mp.Box.orth([5.0, 5.0, 5.0])  # Outside 5Å box
)

# Combined constraints
combined = mpc.AndConstraint(
    a=mpc.InsideBoxConstraint(box_region),
    b=mpc.OutsideSphereConstraint(sphere_region)
)
```

## Complete Workflow Example

Let's walk through packing water molecules into a box:

```python
import molpy as mp
import molpy.pack as mpk
from molpy.pack.packer.packmol import Packmol
from molpy.pack.constraint import InsideBoxConstraint

# Step 1: Create water molecule
water = mp.Atomistic()
water["name"] = "water"

# Add water atoms (simplified)
water.def_atom(xyz=[0.0, 0.0, 0.0], element='O', name='O1')
water.def_atom(xyz=[0.96, 0.0, 0.0], element='H', name='H1')
water.def_atom(xyz=[-0.24, 0.93, 0.0], element='H', name='H2')

# Add bonds
water.def_bond(i=0, j=1)
water.def_bond(i=0, j=2)

# Convert to Frame
water_frame = water.to_frame()

# Step 2: Define packing constraint
box_region = mp.Box.orth([20.0, 20.0, 20.0])  # 20Å cubic box
constraint = InsideBoxConstraint(region=box_region)

# Step 3: Create packing target
target = mpk.Target(
    frame=water_frame,
    number=100,  # Pack 100 water molecules
    constraint=constraint,
    name="water_system"
)

# Step 4: Create Packmol instance
packmol = Packmol(executable="packmol")

# Step 5: Run packing
packed_system = packmol(
    targets=[target],
    max_steps=1000,
    seed=4628
)

print(f"✓ Packing completed: {len(packed_system['atoms'])} atoms")
print(f"✓ System dimensions: {packed_system.box.lengths}")
```

## Advanced Usage Patterns

### Multiple Molecule Types

```python
# Create different molecule targets
water_target = mpk.Target(
    frame=water_frame,
    number=50,
    constraint=water_constraint,
    name="water"
)

ion_target = mpk.Target(
    frame=ion_frame,
    number=10,
    constraint=ion_constraint,
    name="ions"
)

# Pack multiple types
packed_system = packmol(
    targets=[water_target, ion_target],
    max_steps=1000,
    seed=4628
)
```

### Complex Spatial Constraints

```python
# Create a complex constraint: inside box but outside central sphere
box_constraint = mpc.InsideBoxConstraint(
    region=mp.Box.orth([30.0, 30.0, 30.0])
)

sphere_constraint = mpc.OutsideSphereConstraint(
    region=mp.Box.sphere(origin=[15, 15, 15], radius=5.0)
)

# Combine constraints
complex_constraint = mpc.AndConstraint(
    a=box_constraint,
    b=sphere_constraint
)

# Use in target
target = mpk.Target(
    frame=molecule_frame,
    number=200,
    constraint=complex_constraint,
    name="shell_molecules"
)
```

### Custom Working Directories

```python
# Specify custom working directory
packmol = Packmol(
    executable="packmol",
    workdir="/scratch/packmol_work"
)

# Or override per call
result = packmol(
    targets=my_targets,
    workdir="/scratch/specific_packing"
)
```

### Generate Input Only

```python
# Generate Packmol input file without running
input_file = packmol.generate_input_only(
    targets=my_targets,
    max_steps=1000,
    seed=4628,
    workdir="/scratch/input_only"
)

print(f"Generated input file: {input_file}")

# You can then run packmol manually if needed
# packmol < /scratch/input_only/.packmol.inp > /scratch/input_only/.packmol.out
```

## Troubleshooting Common Issues

### Packmol Not Found

```bash
# Check if packmol is in PATH
which packmol

# Install packmol if needed
conda install -c conda-forge packmol

# Or specify full path
packmol = Packmol(executable="/usr/local/bin/packmol")
```

### Packing Fails

```python
# Check constraint definitions
print(f"Constraint: {target.constraint}")
print(f"Target region: {target.constraint.region}")

# Verify molecule structures
print(f"Molecule atoms: {len(target.frame['atoms'])}")
print(f"Molecule box: {target.frame.box}")

# Try with fewer molecules first
target.number = 10  # Start with small number
```

### Memory Issues

```python
# Large systems may need more memory
# Consider packing in stages

# Stage 1: Pack core molecules
core_target = mpk.Target(frame=core_frame, number=100, constraint=core_constraint)
core_system = packmol(targets=[core_target])

# Stage 2: Pack around core
shell_target = mpk.Target(frame=shell_frame, number=200, constraint=shell_constraint)
shell_system = packmol(targets=[shell_target])
```

## Performance Optimization

### Parallel Processing

```python
import concurrent.futures

def pack_molecules(targets, workdir):
    local_packmol = Packmol(workdir=workdir)
    return local_packmol(targets=targets)

# Pack different regions in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    futures = []
    for i, region_targets in enumerate(regions):
        workdir = f"/scratch/region_{i}"
        future = executor.submit(pack_molecules, region_targets, workdir)
        futures.append(future)

    results = [future.result() for future in futures]
```

### Memory Management

```python
# Clean up intermediate files after successful packing
import shutil

# After successful completion
shutil.rmtree(packmol.workdir)

# Or keep specific files
for file in packmol.intermediate_files:
    if file != ".optimized.pdb":  # Keep the result
        (packmol.workdir / file).unlink()
```

## Integration with Other MolPy Components

### Working with Polymers

```python
from molpy.builder.polymer import Monomer

# Create polymer monomer
monomer = Monomer()
# ... add atoms and bonds ...

# Pack polymer chains
polymer_target = mpk.Target(
    frame=monomer.to_frame(),
    number=50,
    constraint=polymer_constraint,
    name="polymer_chains"
)
```

### Reading Existing Structures

```python
# Read a PDB file and pack it
frame = mp.io.read_pdb("molecule.pdb")
struct = mp.Atomistic.from_frame(frame)

# Create target
target = mpk.Target(
    frame=frame,
    number=100,
    constraint=packing_constraint
)
```

## Best Practices

1. **Start Small**: Begin with few molecules to test constraints
2. **Validate Constraints**: Ensure constraints are physically reasonable
3. **Use Descriptive Names**: Name targets and constraints clearly
4. **Monitor Progress**: Check Packmol output for convergence
5. **Clean Up**: Remove intermediate files after successful packing

## Related Resources

- [Packmol Official Documentation](http://leandro.iqm.unicamp.br/m3g/packmol/)
- [MolPy Core API](https://molcrafts.github.io/molpy/)
- [MolQ Job Management](https://github.com/molcrafts/molq)

---

This API design follows MolPy's core principles: **LLM-friendly, modular, and composable**. The Packmol component can be used independently or combined with other MolPy components to create complex molecular systems.
