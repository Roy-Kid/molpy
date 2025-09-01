# How to Use Programming-Based Molecular Modeling

This comprehensive guide demonstrates how to use MolPy's programming-based approach for building complex molecular systems. This approach emphasizes the separation of molecular structure definition from spatial operations, making it easier to create, manipulate, and combine molecular components.

## Overview

Programming-based molecular modeling in MolPy follows a fundamental principle: **separate structure from space**. This approach allows you to:

- Define molecular structures once and reuse them multiple times
- Perform spatial operations (translation, rotation, scaling) independently
- Combine molecular components in a modular, maintainable way
- Build complex systems from simple, well-defined building blocks

The key insight is that a molecule's chemical structure (atoms, bonds, angles, dihedrals) is conceptually different from its spatial arrangement (position, orientation, scale). By keeping these concerns separate, you can build more complex systems with cleaner, more maintainable code.

## Core Concepts

### Atomistic vs Spatial: Separation of Concerns

In MolPy's programming-based approach, we distinguish between two fundamental aspects of molecular modeling:

**Atomistic (Structure)**: This represents the intrinsic properties of a molecule - its chemical composition, bonding patterns, and internal geometry. An `Atomistic` object contains:
- Atomic coordinates relative to the molecule's internal reference frame
- Chemical bonds, angles, and dihedrals
- Force field parameters and atom types
- Molecular topology and connectivity

**Spatial (Position)**: This represents where and how the molecule is positioned in 3D space. The `Spatial` wrapper provides:
- Translation operations to move molecules to specific positions
- Rotation operations to orient molecules in desired directions
- Scaling operations to adjust molecular size
- Spatial queries and distance calculations

This separation allows you to create a molecular template once and then place multiple copies of it at different locations and orientations throughout your system.

### The Wrapper Pattern

The `Spatial` wrapper is a key component of programming-based molecular modeling. It acts as a "spatial interface" for any molecular structure:

```python
from molpy.core.atomistic import Atomistic
from molpy.core.wrapper import Spatial

# Create a molecular structure
water = Atomistic(name="water")
# ... add atoms and bonds ...

# Wrap it with spatial capabilities
spatial_water = Spatial(water)

# Now you can perform spatial operations
spatial_water.move([10.0, 5.0, 2.0])  # Translate to position
spatial_water.rotate([0, 0, 1], np.pi/4)  # Rotate around Z-axis
spatial_water.scale(1.5)  # Scale by 50%
```

The wrapper pattern ensures that:
- The original molecular structure remains unchanged
- Spatial operations are applied consistently across all atoms
- You can chain multiple operations together
- The interface is clean and intuitive

### Modular Design Philosophy

Programming-based molecular modeling encourages a modular approach where you build complex systems from simple, reusable components:

1. **Define Templates**: Create well-defined molecular building blocks using class inheritance
2. **Spatial Operations**: Position and orient templates as needed
3. **Combine Components**: Assemble templates into larger systems
4. **Iterate and Refine**: Modify individual components without affecting others

This modularity makes your code more maintainable, testable, and reusable across different projects.

## Basic Operations

### Creating Molecular Structures

The foundation of MolTemplate-style modeling is creating well-defined molecular structures. For one-time use, create `Atomistic` instances directly:

```python
from molpy.core.atomistic import Atomistic
import numpy as np

# Create the molecular structure directly
water = Atomistic(name="water")

# Add atoms with internal coordinates
# Oxygen at origin, hydrogens at typical H-O-H angles
oxygen = water.def_atom(
    name="O",
    element="O",
    xyz=np.array([0.0, 0.0, 0.0])
)

hydrogen1 = water.def_atom(
    name="H1",
    element="H",
    xyz=np.array([0.9572, 0.0, 0.0])  # ~1 Å O-H bond
)

hydrogen2 = water.def_atom(
    name="H2",
    element="H",
    xyz=np.array([-0.2400, 0.0, 0.0])  # ~1 Å O-H bond
)

# Define chemical bonds
water.def_bond(oxygen, hydrogen1)
water.def_bond(oxygen, hydrogen2)

print(f"Water molecule created with {len(water.atoms)} atoms")
```

For reusable molecular templates, use the new wrapper system with multiple inheritance from `Atomistic` and `Spatial`:

```python
from molpy.core.atomistic import Atomistic
from molpy.core.wrapper import Spatial
import numpy as np

class H2O(Atomistic, Spatial):
    """Water molecule template class with built-in spatial capabilities."""

    def __post_init__(self):
        """Build the water molecule structure after initialization."""
        # Oxygen at origin
        self.oxygen = self.def_atom(
            name="O",
            element="O",
            xyz=np.array([0.0, 0.0, 0.0])
        )

        # Hydrogen atoms at typical H-O-H angles
        self.hydrogen1 = self.def_atom(
            name="H1",
            element="H",
            xyz=np.array([0.9572, 0.0, 0.0])
        )

        self.hydrogen2 = self.def_atom(
            name="H2",
            element="H",
            xyz=np.array([-0.2400, 0.0, 0.0])
        )

        # Define chemical bonds
        self.def_bond(self.oxygen, self.hydrogen1)
        self.def_bond(self.oxygen, self.hydrogen2)

# Create water template instances with full spatial capabilities
water1 = H2O("water1")
water2 = H2O("water2")

# Now you can use spatial operations directly on the template instances!
water1.move([5.0, 0.0, 0.0])  # Move first water molecule
water2.rotate([0, 0, 1], np.pi/2)  # Rotate second water molecule

print(f"Created {len(water1.atoms)} atoms in water template")
print(f"Water1 position: {water1.positions[0]}")  # Access positions directly
print(f"Water2 symbols: {water2.symbols}")  # Access symbols directly
```

**Key points about structure creation with the new wrapper system**:
- **Direct Creation**: For one-time use, create `Atomistic` instances directly
- **Multiple Inheritance**: For reusable templates, inherit from both `Atomistic` and `Spatial` to get full functionality
- **Built-in Spatial Operations**: Templates automatically have spatial capabilities without explicit wrapping
- **Internal Coordinates**: Define atoms relative to a molecular reference frame
- **Chemical Bonds**: Establish the molecular topology
- **Naming Convention**: Use descriptive names for atoms and molecules
- **Unified Interface**: Access both structural and spatial properties through the same object
- **Simplified Class Design**: No need for `__init__` or `_build_molecule` methods - use `__post_init__()` instead

### Spatial Transformations

With the new wrapper system, you have two approaches for spatial operations:

#### Approach 1: Direct Spatial Operations on Template Instances

When you inherit from both `Atomistic` and `Spatial`, your template instances automatically have spatial capabilities:

```python
# Create water template with built-in spatial operations
water = H2O("water")

# Translation: Move molecule to specific position
water.move([5.0, 3.0, 1.0])
print(f"Water moved to position: {water.positions[0]}")

# Rotation: Orient molecule around specific axis
rotation_axis = np.array([0.0, 0.0, 1.0])  # Z-axis
rotation_angle = np.pi / 2  # 90 degrees
water.rotate(rotation_axis, rotation_angle)

# Scaling: Adjust molecular size
water.scale(1.2)  # Scale by 20%

# Access spatial properties directly
print(f"Water center of mass: {water.xyz.mean(axis=0)}")
print(f"Water symbols: {water.symbols}")
```

#### Approach 2: Explicit Wrapping for Existing Structures

For existing structures or when you need to add spatial capabilities dynamically:

```python
from molpy.core.wrapper import Spatial

# Wrap existing structure for spatial operations
existing_water = Atomistic(name="existing_water")
# ... add atoms and bonds ...

spatial_water = Spatial(existing_water)

# Now you can perform spatial operations
spatial_water.move([10.0, 5.0, 2.0])
spatial_water.rotate([0, 0, 1], np.pi/4)
spatial_water.scale(1.5)
```

**Understanding Spatial Operations**:
- **Translation (`move`)**: Shifts all atoms by a vector while preserving internal geometry
- **Rotation (`rotate`)**: Rotates the molecule around an axis through a specified origin
- **Scaling (`scale`)**: Expands or contracts the molecule around a center point
- **Chaining**: Multiple operations can be combined for complex transformations
- **Direct Property Access**: Access `xyz`, `positions`, and `symbols` directly from template instances
- **Unified Interface**: No need to unwrap or access wrapped objects - everything is available directly

### The New Wrapper System: Benefits and Features

The refactored wrapper system provides several key improvements over the old approach:

#### 0. **Simplified Class Design**
With the new wrapper system, you can create molecular templates with minimal boilerplate code:

```python
# Minimal class definition - no __init__ needed!
class SimpleWater(Atomistic, Spatial):
    def __post_init__(self):
        """Called automatically after initialization."""
        # Build molecule structure
        self.oxygen = self.def_atom(name="O", element="O", xyz=[0, 0, 0])
        self.hydrogen1 = self.def_atom(name="H1", element="H", xyz=[0.9572, 0, 0])
        self.hydrogen2 = self.def_atom(name="H2", element="H", xyz=[-0.2400, 0, 0])

        # Create bonds
        self.def_bond(self.oxygen, self.hydrogen1)
        self.def_bond(self.oxygen, self.hydrogen2)

# Usage - just create an instance!
water = SimpleWater("water")  # No need to call __init__ or _build_molecule
print(f"Water has {len(water.atoms)} atoms")  # Works immediately
```

**Why this works:**
- The wrapper system automatically calls `__post_init__()` after all wrappers are initialized
- No need for `__init__` methods unless you need custom initialization logic
- The `__post_init__()` method is the perfect place to build molecular structures
- This pattern eliminates boilerplate code and makes templates cleaner

#### 1. **Multiple Inheritance Auto-Composition**
```python
class Methane(Atomistic, Spatial):
    """Methane molecule with built-in spatial capabilities."""

    def __post_init__(self):
        """Build methane molecule (CH4)."""
        # Carbon at center
        self.carbon = self.def_atom(name="C", element="C", xyz=[0, 0, 0])

        # Hydrogens at tetrahedral positions
        self.h1 = self.def_atom(name="H1", element="H", xyz=[1.09, 1.09, 1.09])
        self.h2 = self.def_atom(name="H2", element="H", xyz=[1.09, -1.09, -1.09])
        self.h3 = self.def_atom(name="H3", element="H", xyz=[-1.09, 1.09, -1.09])
        self.h4 = self.def_atom(name="H4", element="H", xyz=[-1.09, -1.09, 1.09])

        # Define bonds
        for h in [self.h1, self.h2, self.h3, self.h4]:
            self.def_bond(self.carbon, h)

# Create and use methane with full capabilities
methane = Methane("CH4")
methane.move([10, 0, 0])  # Move to position
methane.rotate([0, 1, 0], np.pi/4)  # Rotate around Y-axis
print(f"Methane symbols: {methane.symbols}")  # Access symbols directly
print(f"Methane positions: {methane.positions}")  # Access positions directly
```

#### 2. **Automatic Post-Initialization**
The wrapper system automatically calls `__post_init__()` methods in the correct order:
```python
class CustomMolecule(Atomistic, Spatial):
    def __post_init__(self):
        """Called automatically after all wrappers are initialized."""
        self.custom_property = "initialized"
        print("Custom molecule fully initialized!")

# The __post_init__ is called automatically
mol = CustomMolecule("custom")
assert mol.custom_property == "initialized"
```

#### 3. **Wrapper Chain Management**
```python
# Create wrapper chain
base = Atomistic(name="base")
spatial_wrapper = Spatial(base)

# Query wrapper information
print(f"Wrapper depth: {spatial_wrapper.wrapper_depth()}")  # 2 (Atomistic + Spatial)
print(f"Wrapper types: {spatial_wrapper.wrapper_types()}")  # ['Spatial', 'Atomistic', 'Struct']
print(f"Unwrapped object: {type(spatial_wrapper.unwrap())}")  # <class 'molpy.core.atomistic.Atomistic'>
```

#### 4. **Seamless Property Access**
```python
class Protein(Atomistic, Spatial):
    def __post_init__(self):
        """Add protein atoms after initialization."""
        # Add protein atoms...
        pass

protein = Protein("my_protein")

# Access both structural and spatial properties seamlessly
print(f"Protein has {len(protein.atoms)} atoms")  # Structural
print(f"Protein center: {protein.xyz.mean(axis=0)}")  # Spatial
print(f"Protein symbols: {protein.symbols}")  # Combined
print(f"Protein positions: {protein.positions}")  # Combined

# Perform spatial operations
protein.move([100, 0, 0])  # Move protein
protein.scale(1.5)  # Scale protein
```

### Molecular Assembly

The power of programming-based molecular modeling comes from combining multiple molecular components:

```python
def create_water_system(n_molecules=10):
    """Create a system with multiple water molecules."""

    # Create the base system
    system = Atomistic(name="water_system")

    # Add multiple water molecules at different positions
    for i in range(n_molecules):
        # Create a copy of the water template
        water_copy = H2O(f"water_{i}")

        # Wrap for spatial operations
        spatial_water = Spatial(water_copy)

        # Position each water molecule
        x_pos = i * 3.0  # Space molecules 3 Å apart
        y_pos = np.random.uniform(-5.0, 5.0)  # Random Y position
        z_pos = np.random.uniform(-5.0, 5.0)  # Random Z position

        spatial_water.move([x_pos, y_pos, z_pos])

        # Add to system
        system.add_struct(water_copy)

    return system

# Create the multi-molecule system
water_system = create_water_system(15)
print(f"Created system with {len(water_system.atoms)} total atoms")
```

**Assembly Principles**:
- **Template Instances**: Use class instances for reusable molecular templates
- **Spatial Independence**: Each molecule can be positioned independently
- **System Integration**: Use `add_struct()` to combine components
- **Scalability**: This approach scales well to hundreds or thousands of molecules

## Practical Patterns

### Template Creation and Reuse

Creating reusable templates using class inheritance is fundamental to programming-based molecular modeling:

```python
class Benzene(Atomistic, Spatial):
    """Benzene ring template class with proper geometry."""

    def __post_init__(self):
        """Build the benzene ring structure."""
        # Carbon atoms in a regular hexagon
        self.carbon_atoms = []
        radius = 1.40  # C-C bond length in Å
        n_carbons = 6

        for i in range(n_carbons):
            angle = i * 2 * np.pi / n_carbons
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)

            carbon = self.def_atom(
                name=f"C{i+1}",
                element="C",
                xyz=np.array([x, y, 0.0])
            )
            self.carbon_atoms.append(carbon)

        # Hydrogen atoms attached to each carbon
        self.hydrogen_atoms = []
        h_radius = radius + 1.09  # C-H bond length

        for i, carbon in enumerate(self.carbon_atoms):
            angle = i * 2 * np.pi / n_carbons
            h_x = h_radius * np.cos(angle)
            h_y = h_radius * np.sin(angle)

            hydrogen = self.def_atom(
                name=f"H{i+1}",
                element="H",
                xyz=np.array([h_x, h_y, 0.0])
            )
            self.hydrogen_atoms.append(hydrogen)

        # Create bonds: C-C ring and C-H bonds
        for i in range(n_carbons):
            j = (i + 1) % n_carbons
            self.def_bond(self.carbon_atoms[i], self.carbon_atoms[j])  # C-C ring
            self.def_bond(self.carbon_atoms[i], self.hydrogen_atoms[i])  # C-H

# Create and verify the template
benzene1 = Benzene("benzene1")
benzene2 = Benzene("benzene2")
print(f"Benzene template: {len(benzene1.atoms)} atoms, {len(benzene1.bonds)} bonds")
```

**Template Design Guidelines**:
- **Class Inheritance**: Inherit from `(Atomistic, Spatial)` for reusable templates with built-in spatial capabilities
- **Use __post_init__()**: Build molecular structures in the `__post_init__()` method, not in `__init__`
- **Minimal Boilerplate**: No need for `__init__` methods unless you need custom initialization logic
- **Geometric Accuracy**: Use realistic bond lengths and angles
- **Naming Convention**: Use systematic names for easy identification
- **Modularity**: Design templates that can be easily connected
- **Documentation**: Include clear descriptions of connection points

### Batch Operations and Spatial Arrangement

Programming-based molecular modeling excels at creating multiple instances and arranging them systematically:

```python
def create_benzene_monolayer(n_x=5, n_y=5, spacing=5.0):
    """Create a 2D monolayer of benzene molecules."""

    monolayer = Atomistic(name="benzene_monolayer")

    # Create a grid of benzene molecules
    for i in range(n_x):
        for j in range(n_y):
            # Create a new instance of the benzene template
            benzene = Benzene(f"benzene_{i}_{j}")

            # Wrap for spatial operations
            spatial_benzene = Spatial(benzene)

            # Calculate position in the grid
            x_pos = i * spacing
            y_pos = j * spacing

            # Position the molecule
            spatial_benzene.move([x_pos, y_pos, 0.0])

            # Optional: Add some random rotation for realism
            if np.random.random() > 0.5:
                rotation_axis = np.array([0.0, 0.0, 1.0])
                rotation_angle = np.random.uniform(0, 2 * np.pi)
                spatial_benzene.rotate(rotation_axis, rotation_angle)

            # Add to monolayer
            monolayer.add_struct(benzene)

    return monolayer

# Create the monolayer
monolayer = create_benzene_monolayer(6, 6, 6.0)
print(f"Created monolayer with {len(monolayer.atoms)} total atoms")
```

**Batch Operation Strategies**:
- **Template Instances**: Create new instances for each molecule
- **Grid Patterns**: Use nested loops for regular arrangements
- **Random Variations**: Add randomness for realistic systems
- **Spacing Control**: Maintain appropriate intermolecular distances
- **Performance**: This approach is efficient for hundreds of molecules

### Connection and Bonding Operations

One of the most powerful features of programming-based molecular modeling is the ability to connect molecular components:

```python
class SimpleMonomer(Atomistic, Spatial):
    """Simple 2-atom monomer template for demonstration."""

    def __post_init__(self):
        """Build the monomer structure."""
        self.atom1 = self.def_atom(
            name="A1",
            element="C",
            xyz=np.array([0.0, 0.0, 0.0])
        )
        self.atom2 = self.def_atom(
            name="A2",
            element="C",
            xyz=np.array([1.54, 0.0, 0.0])
        )

        self.def_bond(self.atom1, self.atom2)

def create_polymer_chain(n_monomers=10, bond_length=1.54):
    """Create a linear polymer chain from monomer templates."""

    polymer = Atomistic(name="polymer_chain")

    # Add first monomer
    first_monomer = SimpleMonomer("monomer_0")
    polymer.add_struct(first_monomer)

    # Add subsequent monomers
    for i in range(1, n_monomers):
        # Create new monomer instance
        monomer = SimpleMonomer(f"monomer_{i}")

        # Wrap for spatial operations
        spatial_monomer = Spatial(monomer)

        # Position monomer for connection
        # This assumes the monomer has connection points at specific atoms
        connection_point = polymer.atoms[-1]["xyz"]  # Last atom of previous monomer
        monomer_connection = monomer.atom1["xyz"]  # First atom of new monomer

        # Calculate translation to connect monomers
        translation_vector = connection_point - monomer_connection + np.array([bond_length, 0.0, 0.0])
        spatial_monomer.move(translation_vector)

        # Add to polymer
        polymer.add_struct(monomer)

        # Create connecting bond (this would need to be implemented)
        # polymer.def_bond(previous_atom, current_atom)

    return polymer

# Create polymer chain
polymer = create_polymer_chain(8)
print(f"Created polymer with {len(polymer.atoms)} atoms")
```

**Connection Strategies**:
- **Template Classes**: Use class instances for consistent molecular structures
- **Anchor Points**: Define clear connection points in your templates
- **Bond Lengths**: Use realistic chemical bond lengths
- **Orientation Control**: Ensure proper molecular alignment
- **Topology Management**: Maintain correct bonding patterns

### Constraint Handling and Optimization

Real molecular systems often have spatial constraints that need to be satisfied:

```python
def create_constrained_system(molecule_class, surface_positions, min_distance=3.0):
    """Create a system with molecules constrained to surface positions."""

    constrained_system = Atomistic(name="constrained_system")

    for i, surface_pos in enumerate(surface_positions):
        # Create new molecule instance
        molecule = molecule_class(f"molecule_{i}")

        # Wrap for spatial operations
        spatial_molecule = Spatial(molecule)

        # Check distance constraints
        valid_position = True
        for existing_molecule in constrained_system.atoms:
            distance = np.linalg.norm(surface_pos - existing_molecule["xyz"])
            if distance < min_distance:
                valid_position = False
                break

        if valid_position:
            # Position molecule at surface
            spatial_molecule.move_to(surface_pos)

            # Add to system
            constrained_system.add_struct(molecule)
        else:
            print(f"Position {i} too close to existing molecules, skipping")

    return constrained_system

# Example surface positions
surface_positions = [
    np.array([0.0, 0.0, 0.0]),
    np.array([5.0, 0.0, 0.0]),
    np.array([0.0, 5.0, 0.0]),
    np.array([5.0, 5.0, 0.0]),
    np.array([2.5, 2.5, 0.0])
]

# Create constrained system using H2O class
constrained = create_constrained_system(H2O, surface_positions, min_distance=4.0)
print(f"Constrained system: {len(constrained.atoms)} molecules placed")
```

**Constraint Handling Approaches**:
- **Distance Checking**: Verify minimum intermolecular distances
- **Position Validation**: Ensure molecules fit within boundaries
- **Iterative Placement**: Try alternative positions if constraints fail
- **Optimization**: Use algorithms to find optimal arrangements

## Advanced Techniques

### Layered Construction

Complex molecular systems can be built in layers, with each layer adding new functionality:

```python
class BaseSubstrate(Atomistic, Spatial):
    """Base substrate layer template."""

    def __post_init__(self):
        """Build the base substrate structure."""
        # Implementation would create a 2D sheet structure
        pass

class FunctionalLayer(Atomistic, Spatial):
    """Functional molecules layer template."""

    def __post_init__(self):
        """Build functional molecules on the base."""
        # Implementation would place molecules on the substrate
        pass

class ProtectiveCoating(Atomistic, Spatial):
    """Protective coating layer template."""

    def __post_init__(self):
        """Build a protective coating layer."""
        # Implementation would add protective molecules
        pass

def create_layered_system():
    """Create a multi-layered molecular system."""

    system = Atomistic(name="layered_system")

    # Layer 1: Base substrate
    base_layer = BaseSubstrate()
    system.add_struct(base_layer)

    # Layer 2: Functional molecules
    functional_layer = FunctionalLayer(base_layer)
    system.add_struct(functional_layer)

    # Layer 3: Protective coating
    protective_layer = ProtectiveCoating(functional_layer)
    system.add_struct(protective_layer)

    return system
```

**Layered Construction Benefits**:
- **Modularity**: Each layer can be developed and tested independently
- **Complexity Management**: Break complex systems into manageable components
- **Reusability**: Layers can be reused in different combinations
- **Maintainability**: Easier to modify individual layers

### Dynamic Structure Modification

Programming-based molecular modeling allows for runtime modification of molecular structures:

```python
def create_adaptive_system(base_class, environment_conditions):
    """Create a system that adapts to environmental conditions."""

    # Start with base template instance
    adaptive_system = base_class("adaptive_system")

    # Modify based on conditions
    if environment_conditions.get('temperature', 300) > 350:
        # High temperature: add stabilizing molecules
        stabilizers = create_stabilizing_molecules()
        spatial_stabilizers = Spatial(stabilizers)
        spatial_stabilizers.move([0.0, 0.0, 5.0])
        adaptive_system.add_struct(stabilizers)

    if environment_conditions.get('pressure', 1.0) > 2.0:
        # High pressure: compress the system
        spatial_system = Spatial(adaptive_system)
        spatial_system.scale(0.9)  # Compress by 10%

    return adaptive_system

# Example usage
conditions = {'temperature': 400, 'pressure': 2.5}
adaptive = create_adaptive_system(H2O, conditions)
```

**Dynamic Modification Capabilities**:
- **Conditional Logic**: Modify structures based on external parameters
- **Runtime Adaptation**: Adjust systems during simulation
- **Performance Optimization**: Only apply modifications when needed
- **Flexibility**: Handle changing requirements dynamically

### Optimization Strategies

For complex systems, optimization algorithms can help find optimal molecular arrangements:

```python
def optimize_molecular_arrangement(molecule_class, n_molecules, target_function, constraints):
    """Optimize molecular arrangement using iterative improvement."""

    best_arrangement = None
    best_score = float('-inf')

    for iteration in range(100):  # Max 100 iterations
        # Create candidate arrangement
        candidate = create_candidate_arrangement(molecule_class, n_molecules, constraints)

        # Evaluate candidate
        candidate_score = target_function(candidate)

        # Accept if better
        if candidate_score > best_score:
            best_arrangement = candidate
            best_score = candidate_score

        # Optional: Add some randomness to avoid local optima
        if np.random.random() < 0.1:
            candidate = randomize_positions(candidate)

    return best_arrangement

def target_function(arrangement):
    """Example target function: maximize intermolecular interactions."""
    # This would calculate some measure of system quality
    # e.g., total interaction energy, packing efficiency, etc.
    return 0.0  # Placeholder

def create_candidate_arrangement(molecule_class, n_molecules, constraints):
    """Create a candidate arrangement by perturbing current positions."""
    # Implementation would create variations of the current arrangement
    pass
```

**Optimization Approaches**:
- **Iterative Improvement**: Gradually improve arrangements
- **Random Sampling**: Explore different configurations
- **Constraint Satisfaction**: Ensure all constraints are met
- **Multi-objective**: Balance multiple optimization goals

## Best Practices

### Performance Considerations

When building large molecular systems, performance becomes important:

```python
def create_efficient_system(molecule_class, n_molecules=1000):
    """Create a large system efficiently."""

    # Pre-allocate arrays for better performance
    system = Atomistic(name="large_system")

    # Batch operations when possible
    batch_size = 100
    for batch_start in range(0, n_molecules, batch_size):
        batch_end = min(batch_start + batch_size, n_molecules)

        # Create batch of molecules
        batch_molecules = []
        for i in range(batch_start, batch_end):
            molecule = molecule_class(f"molecule_{i}")
            spatial_molecule = Spatial(molecule)

            # Calculate position efficiently
            x_pos = (i % 50) * 3.0
            y_pos = (i // 50) * 3.0
            spatial_molecule.move([x_pos, y_pos, 0.0])

            batch_molecules.append(molecule)

        # Add batch to system
        for molecule in batch_molecules:
            system.add_struct(molecule)

    return system
```

**Performance Tips**:
- **Batch Operations**: Process multiple molecules together
- **Vectorization**: Use numpy operations when possible
- **Memory Management**: Avoid unnecessary object creation
- **Algorithm Choice**: Select appropriate algorithms for your system size

### Code Organization

Well-organized code makes MolTemplate-style modeling more maintainable:

```python
# molecular_templates.py
class MolecularTemplates:
    """Collection of molecular template classes for reuse."""

    @staticmethod
    def water():
        return H2O()

    @staticmethod
    def benzene():
        return Benzene()

    @staticmethod
    def monomer():
        return SimpleMonomer()

# spatial_operations.py
class SpatialOperations:
    """Common spatial operations for molecular systems."""

    @staticmethod
    def arrange_in_grid(template_class, n_x, n_y, spacing):
        return create_grid_arrangement(template_class, n_x, n_y, spacing)

    @staticmethod
    def arrange_in_sphere(template_class, n_molecules, radius):
        return create_spherical_arrangement(template_class, n_molecules, radius)

# system_builder.py
class SystemBuilder:
    """High-level system construction interface."""

    def __init__(self, template_classes, operations):
        self.template_classes = template_classes
        self.operations = operations

    def build_monolayer(self, template_name, dimensions):
        template_class = self.template_classes.get(template_name)
        return self.operations.arrange_in_grid(template_class, *dimensions)
```

**Code Organization Principles**:
- **Separation of Concerns**: Keep templates, operations, and building separate
- **Class-based Templates**: Use class inheritance for molecular templates with `(Atomistic, Spatial)`
- **Simplified Class Design**: Use `__post_init__()` instead of `__init__` for structure building
- **Reusability**: Design components that can be reused across projects
- **Configuration**: Use configuration files for system parameters
- **Documentation**: Document interfaces and usage patterns

### Debugging and Validation

Robust debugging and validation are essential for complex molecular systems:

```python
def validate_molecular_system(system):
    """Validate a molecular system for common issues."""

    issues = []

    # Check for overlapping atoms
    for i, atom1 in enumerate(system.atoms):
        for j, atom2 in enumerate(system.atoms[i+1:], i+1):
            distance = np.linalg.norm(atom1["xyz"] - atom2["xyz"])
            if distance < 0.5:  # Atoms too close
                issues.append(f"Atoms {i} and {j} too close: {distance:.3f} Å")

    # Check for disconnected components
    if len(system.bonds) == 0:
        issues.append("System has no bonds - may be disconnected")

    # Check for unrealistic bond lengths
    for bond in system.bonds:
        bond_length = np.linalg.norm(
            bond.itom["xyz"] - bond.jtom["xyz"]
        )
        if bond_length > 3.0:  # Unrealistically long bond
            issues.append(f"Bond {bond.itom['name']}-{bond.jtom['name']} too long: {bond_length:.3f} Å")

    return issues

# Usage
system = create_complex_system()
validation_issues = validate_molecular_system(system)

if validation_issues:
    print("Validation issues found:")
    for issue in validation_issues:
        print(f"  - {issue}")
else:
    print("System validation passed!")
```

**Validation Strategies**:
- **Geometric Checks**: Verify realistic bond lengths and angles
- **Topology Validation**: Ensure proper molecular connectivity
- **Spatial Constraints**: Check for overlapping or misplaced molecules
- **Performance Monitoring**: Track system creation time and memory usage

## Comparison with Other Approaches

### Programming-Based vs Traditional Methods

Programming-based molecular modeling offers several advantages over traditional approaches:

**Traditional Approach**:
```python
# Old way: Everything mixed together
def create_water_system_old_way():
    system = []
    for i in range(10):
        # Create atoms directly
        o_atom = {'element': 'O', 'x': i*3.0, 'y': 0.0, 'z': 0.0}
        h1_atom = {'element': 'H', 'x': i*3.0 + 0.9572, 'y': 0.0, 'z': 0.0}
        h2_atom = {'element': 'H', 'x': i*3.0 - 0.2400, 'y': 0.0, 'z': 0.0}

        # Add to system
        system.extend([o_atom, h1_atom, h2_atom])

    return system
```

**Programming-Based Approach with New Wrapper System**:
```python
# New way: Clean separation of concerns using the new wrapper system
def create_water_system_moltemplate():
    # Create system
    system = Atomistic(name="water_system")

    # Position multiple instances of the H2O class (now with built-in spatial capabilities)
    for i in range(10):
        water = H2O(f"water_{i}")  # H2O inherits from (Atomistic, Spatial)
        water.move([i*3.0, 0.0, 0.0])  # Direct spatial operations - no wrapping needed!
        system.add_struct(water)

    return system

# Even simpler: Create a water cluster with different orientations
def create_water_cluster():
    system = Atomistic(name="water_cluster")

    # Create water molecules at different positions and orientations
    positions = [
        ([0, 0, 0], 0),      # Center, no rotation
        ([3, 0, 0], np.pi/4), # Right, 45° rotation
        ([0, 3, 0], np.pi/2), # Up, 90° rotation
        ([-3, 0, 0], -np.pi/4) # Left, -45° rotation
    ]

    for i, (pos, angle) in enumerate(positions):
        water = H2O(f"water_{i}")
        water.move(pos)  # Move to position
        water.rotate([0, 0, 1], angle)  # Rotate around Z-axis
        system.add_struct(water)

    return system
```

**Key Advantages**:
- **Maintainability**: Changes to molecular structure only need to be made in the class definition
- **Reusability**: Template classes can be used across multiple projects
- **Readability**: Code clearly separates structure definition from spatial arrangement
- **Scalability**: Easy to add hundreds or thousands of molecules
- **Testing**: Individual components can be tested independently
- **OOP Principles**: Follows proper object-oriented design patterns
- **Unified Interface**: New wrapper system provides seamless access to both structural and spatial properties
- **Auto-Composition**: Multiple inheritance automatically combines functionality without manual wrapping

### Integration with Other Tools

Programming-based molecular modeling integrates well with other molecular modeling tools:

```python
def export_to_lammps(system, filename):
    """Export system to LAMMPS format."""
    # Convert MolPy system to LAMMPS format
    lammps_data = convert_to_lammps_format(system)

    with open(filename, 'w') as f:
        f.write(lammps_data)

def export_to_pdb(system, filename):
    """Export system to PDB format."""
    # Convert MolPy system to PDB format
    pdb_data = convert_to_pdb_format(system)

    with open(filename, 'w') as f:
        f.write(pdb_data)

def import_from_amber(amber_file):
    """Import system from AMBER format."""
    # Import AMBER structure and convert to MolPy format
    amber_system = read_amber_file(amber_file)
    molpy_system = convert_from_amber_format(amber_system)

    return molpy_system
```

**Integration Benefits**:
- **Format Flexibility**: Export to various simulation and visualization formats
- **Tool Interoperability**: Work with existing molecular modeling workflows
- **Data Exchange**: Share systems with collaborators using different tools
- **Pipeline Integration**: Fit into larger computational workflows

---

This guide provides a comprehensive introduction to programming-based molecular modeling in MolPy. By following these principles and patterns, you can build complex molecular systems that are maintainable, scalable, and easy to understand. The key is to always separate structure from space, use class inheritance for reusable templates, leverage the new wrapper system for seamless functionality, and use the simplified `__post_init__()` pattern for cleaner code.

## Summary: The New Wrapper System

The refactored wrapper system represents a significant improvement in MolPy's molecular modeling capabilities:

### **What Changed**
- **Before**: Required explicit wrapping with `Spatial(atomistic)` to access spatial operations
- **After**: Use multiple inheritance `class MyMolecule(Atomistic, Spatial)` for automatic composition

### **Key Benefits**
1. **Simplified API**: No more manual wrapping - spatial capabilities are built-in
2. **Better Performance**: Eliminates wrapper object creation overhead
3. **Cleaner Code**: Single class definition provides all functionality
4. **Type Safety**: Better IDE support and type checking
5. **Consistent Interface**: All properties and methods accessible from the same object

### **Migration Guide**
```python
# Old way (still supported)
class OldWater(Atomistic):
    def __init__(self, name="water"):
        super().__init__(name=name)
        # ... build molecule

water = OldWater("water")
spatial_water = Spatial(water)  # Manual wrapping
spatial_water.move([1, 0, 0])

# New way (recommended)
class NewWater(Atomistic, Spatial):
    def __init__(self, name="water"):
        super().__init__(name=name)
        # ... build molecule

water = NewWater("water")
water.move([1, 0, 0])  # Direct access - no wrapping needed
```

### **Best Practices**
- Use multiple inheritance `(Atomistic, Spatial)` for new molecular templates
- Keep explicit wrapping for dynamic addition of spatial capabilities
- Leverage the unified interface for cleaner, more maintainable code
- Take advantage of automatic post-initialization for complex setup logic

The new wrapper system maintains full backward compatibility while providing a more intuitive and efficient interface for molecular modeling.
