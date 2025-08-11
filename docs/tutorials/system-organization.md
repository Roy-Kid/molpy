# System Organization

Molecular systems can range from simple molecules to complex assemblies with multiple components. This tutorial teaches you how to work with MolPy's system classes to organize, manage, and manipulate complex molecular systems.

## Understanding System Organization

Molecular systems often contain multiple components that need to be organized and managed together. A single molecule might have multiple conformations, a mixture might contain different molecular species, or a crystal might need periodic replication.

MolPy's system classes provide different approaches to organizing these complex assemblies. `FrameSystem` uses the Frame/Block API for columnar data storage, `StructSystem` uses the Struct/Wrapper API for object graph storage, and `PeriodicSystem` adds periodic boundary condition operations. All systems inherit from `Systemic`, which provides common functionality like box management and force field association.

This organization makes it easy to work with systems of different complexity levels. Simple systems can use basic organization, while complex systems can leverage advanced features like periodic replication and vacuum addition.

## Basic System Types

### FrameSystem for Columnar Data

`FrameSystem` is designed for systems that use the Frame/Block API for data storage. It's ideal for molecular dynamics trajectories, analysis workflows, and systems where you need efficient access to atomic properties.

```python
import molpy as mp
import numpy as np

# Create a simple water molecule frame
water_atoms = mp.Block({
    'x': [0.0, 0.9572, -0.2400],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['O', 'H', 'H'],
    'type': [1, 2, 2]
})

water_frame = mp.Frame({'atoms': water_atoms})
water_frame.metadata = {'molecule': 'water', 'charge': 0.0}

# Create a system with the water frame
water_system = mp.FrameSystem(
    frame=water_frame,
    box=mp.Box.cubic(20.0),
    forcefield=mp.ForceField("water_ff")
)

# Access system components
print(f"System frame: {water_system.frame}")
print(f"System box: {water_system.box}")
print(f"System force field: {water_system.forcefield}")
```

FrameSystem provides direct access to the underlying frame, making it easy to work with atomic data and perform analysis operations. The system maintains the frame, box, and force field as separate components that can be accessed and modified independently.

### StructSystem for Object Graphs

`StructSystem` is designed for systems that use the Struct/Wrapper API for object graph storage. It's ideal for complex molecular assemblies, systems with hierarchical structure, and cases where you need object-oriented access to molecular components.

```python
# Create a struct system
struct_system = mp.StructSystem(
    box=mp.Box.cubic(30.0),
    forcefield=mp.ForceField("organic_ff")
)

# Add structures to the system
water_struct = mp.Struct()  # Simplified - in practice this would be more complex
methane_struct = mp.Struct()

struct_system.add_struct(water_struct)
struct_system.add_struct(methane_struct)

# Access system information
print(f"System box: {struct_system.box}")
print(f"System force field: {struct_system.forcefield}")
print(f"Number of structures: {len(struct_system.structs)}")
```

StructSystem maintains a list of structures and provides methods for adding and managing them. This approach is useful when you need to work with complex molecular assemblies that have multiple distinct components.

## System Properties and Management

### Box and Force Field Management

All system types inherit from `Systemic`, which provides common functionality for managing simulation boxes and force fields.

```python
# Create a system with custom box and force field
custom_box = mp.Box.orth(10.0, 15.0, 20.0)
custom_ff = mp.ForceField("custom_ff", units="real")

# Create system
system = mp.FrameSystem(
    frame=water_frame,
    box=custom_box,
    forcefield=custom_ff
)

# Modify system properties
new_box = mp.Box.cubic(25.0)
system.set_box(new_box)

new_ff = mp.ForceField("new_ff")
system.set_forcefield(new_ff)

# Access current properties
print(f"Current box: {system.box}")
print(f"Current force field: {system.forcefield}")
```

The system provides a unified interface for managing these fundamental properties. You can change the box dimensions, switch force fields, or modify other system parameters as needed for your simulation or analysis.

### System Information and Status

Systems provide various ways to access information about their current state and organization.

```python
def analyze_system(system):
    """Analyze system properties and organization."""
    print("=== System Analysis ===")

    # Basic information
    print(f"System type: {type(system).__name__}")
    print(f"Box: {system.box}")
    print(f"Force field: {system.forcefield.name}")

    # Frame-specific information
    if hasattr(system, 'frame') and system.frame is not None:
        frame = system.frame
        print(f"Frame blocks: {list(frame.blocks())}")
        if 'atoms' in frame:
            atoms = frame['atoms']
            print(f"Atoms: {atoms.nrows}")

    # Struct-specific information
    if hasattr(system, 'structs'):
        print(f"Structures: {len(system.structs)}")

# Analyze the water system
analyze_system(water_system)
```

This analysis provides a comprehensive view of your system's organization and current state. It shows what type of system you have, what components it contains, and how they're organized.

## Advanced System Operations

### Periodic System Operations

`PeriodicSystem` wraps other system types and adds periodic boundary condition operations. This is essential for crystal structures, bulk materials, and systems that need periodic replication.

```python
# Create a periodic system from the water system
# First, ensure the box is periodic
water_system.box = mp.Box.cubic(20.0)  # This creates a periodic box

# Create periodic system
periodic_system = mp.PeriodicSystem(water_system)

# Create a supercell (2x2x2 replication)
transformation_matrix = [[2, 0, 0], [0, 2, 0], [0, 0, 2]]
supercell = periodic_system.make_supercell(transformation_matrix)

print(f"Original system atoms: {water_system.frame['atoms'].nrows}")
print(f"Supercell atoms: {supercell.frame['atoms'].nrows}")
print(f"Expected atoms: {water_system.frame['atoms'].nrows * 8}")  # 2³ = 8
```

Periodic systems enable advanced operations like supercell creation, which is essential for studying bulk properties, defects, and crystal structures. The system automatically handles the replication of atomic coordinates and maintains proper periodic boundary conditions.

### Vacuum Addition

Periodic systems can also add vacuum regions, which is useful for studying surfaces, interfaces, and systems that need to transition from periodic to non-periodic conditions.

```python
# Add vacuum in the z-direction
vacuum_system = periodic_system.add_vacuum(
    vacuum_thickness=10.0,
    direction='z'
)

print(f"Original box z-length: {water_system.box.lengths[2]}")
print(f"Vacuum box z-length: {vacuum_system.box.lengths[2]}")
```

Vacuum addition is useful for studying surface phenomena, creating interfaces between different materials, or transitioning from bulk to surface calculations.

## System Organization Patterns

### Multi-Component Systems

Complex systems often contain multiple molecular components that need to be organized and managed together.

```python
# Create a system for a molecular mixture
mixture_system = mp.FrameSystem(
    box=mp.Box.cubic(50.0),
    forcefield=mp.ForceField("mixture_ff")
)

# Create frames for different components
water_frame = mp.Frame({
    'atoms': mp.Block({
        'x': [0.0, 0.9572, -0.2400],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H']
    })
})

methane_frame = mp.Frame({
    'atoms': mp.Block({
        'x': [10.0, 10.0, 10.0, 10.0, 10.0],
        'y': [10.0, 10.0, 10.0, 10.0, 10.0],
        'z': [10.0, 11.09, 9.91, 10.0, 10.0],
        'element': ['C', 'H', 'H', 'H', 'H']
    })
})

# Combine components (simplified - in practice you'd merge the frames)
mixture_system.frame = water_frame  # For this example, just use water
mixture_system.metadata = {
    'components': ['water', 'methane'],
    'composition': 'binary mixture'
}

print(f"Mixture system: {mixture_system.frame['atoms'].nrows} total atoms")
```

This pattern allows you to organize complex mixtures, multi-component systems, and assemblies with different molecular species.

### System Evolution and Trajectories

Systems can represent different states of the same molecular assembly, enabling trajectory analysis and system evolution studies.

```python
# Create multiple system states for trajectory analysis
trajectory_systems = []

for step in range(5):
    # Create a frame for this step (simplified)
    step_frame = mp.Frame({
        'atoms': mp.Block({
            'x': [0.0 + step * 0.1, 0.9572 + step * 0.1, -0.2400 + step * 0.1],
            'y': [0.0, 0.0, 0.0],
            'z': [0.0, 0.0, 0.0],
            'element': ['O', 'H', 'H']
        })
    })

    step_frame.metadata = {'step': step, 'time': step * 0.001}

    # Create system for this step
    step_system = mp.FrameSystem(
        frame=step_frame,
        box=mp.Box.cubic(20.0),
        forcefield=mp.ForceField("trajectory_ff")
    )

    trajectory_systems.append(step_system)

# Analyze trajectory
for i, system in enumerate(trajectory_systems):
    frame = system.frame
    metadata = frame.metadata
    print(f"Step {metadata['step']}: {frame['atoms']['x'][0]:.3f}")
```

This pattern enables trajectory analysis, where you can study how systems evolve over time or across different conditions.

## Data Organization Principles

### System Hierarchy

Organize your systems with a clear hierarchy: System → Components → Properties. Each system should have a clear purpose and contain only the components needed for that purpose.

Use consistent naming conventions for systems, components, and properties. This makes your code more maintainable and easier for others to understand.

### Component Management

When working with multi-component systems, organize components logically. Group related molecules together, use consistent indexing, and maintain clear relationships between components.

Consider the computational requirements of your system organization. Frame-based systems are efficient for large numbers of atoms, while struct-based systems are better for complex molecular assemblies.

### System Validation

Always validate your system organization to ensure it makes sense for your intended use. Check that boxes are appropriate for your molecular content, that force fields are compatible with your molecular types, and that periodic conditions are set correctly when needed.

## Summary

This tutorial covered the fundamental concepts of MolPy's system organization. You learned how to create and manage different types of systems, how to organize complex molecular assemblies, and how to perform advanced operations like supercell creation and vacuum addition.

The system classes provide flexible organization for molecular systems of varying complexity. FrameSystem handles columnar data efficiently, StructSystem manages object graphs, and PeriodicSystem adds periodic boundary operations. This organization enables complex molecular modeling workflows while maintaining clear structure and relationships.

### Next Steps

Continue your MolPy journey by exploring advanced analysis techniques for complex molecular systems, understanding how to integrate different system types, and learning about performance optimization for large-scale simulations.

Understanding system organization is essential for working with complex molecular assemblies in MolPy. The organized approach provides the foundation for sophisticated molecular modeling and simulation workflows.
