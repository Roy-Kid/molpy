# Block & Frame Basics

Blocks and Frames are the fundamental data structures in MolPy that organize molecular simulation data. This tutorial teaches you how to understand and work with these core components that form the foundation of molecular data management.

## Understanding the Data Model

Molecular simulation data has a natural hierarchical structure. Atoms have properties like coordinates, mass, charge, and element type. Molecules contain multiple atoms with connectivity information. Systems include multiple molecules, simulation boxes, and metadata. Trajectories track how systems evolve over time.

MolPy's Block and Frame design mirrors this hierarchy, making it intuitive to work with molecular data. Instead of scattered arrays and dictionaries, you get organized containers that understand the relationships between different types of molecular information.

A `Block` is a container that maps variable names to NumPy arrays. Think of it as a table where rows represent individual entities like atoms or bonds, and columns represent different properties like x, y, z coordinates, mass, or charge. Data types are automatically optimized for scientific computing, with float64 for coordinates and int32 for indices.

A `Frame` is a hierarchical container that represents a complete molecular system at a point in time. It contains multiple blocks for different data types like atoms, bonds, angles, and dihedrals. It includes metadata like time, temperature, pressure, and step number. This structure allows you to organize related data logically, access data efficiently through consistent interfaces, and maintain data integrity across operations.

## Working with Blocks

### Creating and Accessing Data

Blocks behave like dictionaries but automatically convert values to NumPy arrays. This ensures type consistency and enables efficient numerical operations.

```python
import molpy as mp
import numpy as np

# Create atomic data for a water molecule
water_atoms = mp.Block({
    'x': [0.0, 0.9572, -0.2400],    # X coordinates in Å
    'y': [0.0, 0.0, 0.0],           # Y coordinates in Å
    'z': [0.0, 0.0, 0.0],           # Z coordinates in Å
    'element': ['O', 'H', 'H'],      # Chemical elements
    'type': [1, 2, 2],               # Atom types (O=1, H=2)
    'mass': [15.999, 1.008, 1.008]  # Atomic masses in amu
})

# Access individual variables
oxygen_x = water_atoms['x'][0]  # Oxygen X coordinate
all_positions = water_atoms[['x', 'y', 'z']]  # All atomic coordinates

# Check block properties
print(f"Number of atoms: {water_atoms.nrows}")
print(f"Available variables: {list(water_atoms.keys())}")
print(f"X coordinates shape: {water_atoms['x'].shape}")
```

When you create a block, MolPy validates that all arrays have the same length and converts them to NumPy arrays. The `nrows` property tells you how many entities (atoms, bonds, etc.) are in the block, and you can access variables individually or combine multiple variables into coordinate matrices.

### Modifying and Organizing Data

Blocks provide flexible ways to modify and organize your data while maintaining consistency.

```python
# Add new properties
water_atoms['charge'] = [-0.8476, 0.4238, 0.4238]  # Partial charges in e
water_atoms['radius'] = [1.52, 1.20, 1.20]         # Van der Waals radii in Å

# Modify existing properties
water_atoms['x'] = water_atoms['x'] + 5.0  # Shift all X coordinates by 5 Å

# Remove properties
del water_atoms['radius']

# Create a new block for bonds
water_bonds = mp.Block({
    'i': [0, 0],      # First atom index in each bond
    'j': [1, 2],      # Second atom index in each bond
    'type': [1, 1]    # Bond type for each bond
})

print(f"Water molecule: {water_atoms.nrows} atoms, {water_bonds.nrows} bonds")
```

Blocks maintain data consistency by ensuring all variables have the same number of rows. When you modify coordinates, the changes affect the original data. You can organize different types of data into separate blocks, making it clear what each block represents.

### Data Access Patterns

Blocks support various access patterns that are useful for different types of analysis.

```python
# Row-based access for individual atoms
first_atom = water_atoms[0]  # Get first atom's properties
print(f"First atom: {first_atom}")

# Slicing for subsets
first_two = water_atoms[:2]  # Get first two atoms
print(f"First two atoms: {first_two.nrows}")

# Boolean indexing for selection
oxygens = water_atoms[water_atoms['element'] == 'O']
hydrogens = water_atoms[water_atoms['element'] == 'H']
print(f"Oxygens: {oxygens.nrows}, Hydrogens: {hydrogens.nrows}")

# Multiple variable access
coordinates = water_atoms[['x', 'y', 'z']]  # 3x3 coordinate matrix
print(f"Coordinate matrix shape: {coordinates.shape}")
```

These access patterns let you work with your data efficiently. Row-based access is useful for examining individual entities, slicing creates subsets for analysis, boolean indexing selects entities based on properties, and multiple variable access creates matrices for vector operations.

## Working with Frames

### Creating and Organizing Frames

Frames organize multiple blocks into a hierarchical structure that represents a complete molecular system.

```python
# Create a frame with water molecule data
water_frame = mp.Frame({
    'atoms': water_atoms,
    'bonds': water_bonds
})

# Add metadata
water_frame.metadata = {
    'molecule': 'water',
    'charge': 0.0,
    'description': 'Single water molecule in vacuum'
}

# Access different parts of the frame
atoms_block = water_frame['atoms']
bonds_block = water_frame['bonds']
metadata = water_frame.metadata

print(f"Frame contains {len(list(water_frame.blocks()))} blocks")
print(f"Available variables in atoms: {list(water_frame.variables('atoms'))}")
```

Frames use a dictionary-like interface where block names map to Block objects. The metadata stores system-level information that describes the overall state of your molecular system. This organization makes it easy to work with different types of data while maintaining clear relationships.

### Frame Access Patterns

Frames support both block-level and variable-level access, giving you flexibility in how you work with your data.

```python
# Block-level access
atoms = water_frame['atoms']
bonds = water_frame['bonds']

# Variable-level access using tuples
oxygen_x = water_frame['atoms', 'x']  # Same as water_frame['atoms']['x']
hydrogen_positions = water_frame['atoms', ['x', 'y', 'z']]

# Check what's available
print(f"Blocks: {list(water_frame.blocks())}")
print(f"Variables in atoms: {list(water_frame.variables('atoms'))}")
print(f"Variables in bonds: {list(water_frame.variables('bonds'))}")

# Add new blocks
water_frame['angles'] = mp.Block({
    'i': [1],      # First atom in angle
    'j': [0],      # Central atom
    'k': [2],      # Third atom in angle
    'value': [104.5]  # Angle in degrees
})
```

The tuple-based access `frame['block', 'variable']` provides a convenient shorthand for accessing specific variables. This pattern is especially useful when you need to work with variables from specific blocks in your analysis.

## Data Organization Principles

### Block Naming and Structure

Use descriptive names for your blocks that clearly indicate what they contain. Common names include 'atoms', 'bonds', 'angles', 'dihedrals', and 'impropers'. This naming convention makes your code self-documenting and easier for others to understand.

Organize related variables within each block. Keep all atomic properties together in the atoms block, all connectivity information in the bonds block, and so on. This logical grouping makes it easier to work with specific types of data and maintains clear relationships between different properties.

### Frame Organization

Use frames to organize complete molecular systems. Each frame should represent a snapshot of a system at a specific point in time or a specific configuration. Include relevant metadata that describes the system state, simulation parameters, or other important information.

When working with trajectories or multiple systems, use consistent frame structures. This consistency makes it easier to write analysis code that works across different frames and enables efficient batch processing.

### Data Consistency

Maintain consistency within your blocks by ensuring all variables have the same number of rows. This consistency is enforced automatically by MolPy, but it's important to understand why it's necessary. When you add or remove variables, make sure they align with your existing data structure.

Use consistent data types for similar properties. Coordinates should typically be float64 for precision, while indices and types can be integers. Element symbols and names are typically strings. This consistency makes your data easier to work with and prevents type-related errors.

## (De)serialization

### Dictionary Conversion

Blocks and frames can be converted to and from dictionaries, making them easy to save, load, and share.

```python
# Convert frame to dictionary
water_dict = water_frame.to_dict()
print(f"Serialized frame has {len(water_dict)} keys")

# Convert back to frame
restored_frame = mp.Frame.from_dict(water_dict)
print(f"Restored frame has {len(list(restored_frame.blocks()))} blocks")

# Verify data integrity
original_atoms = water_frame['atoms']
restored_atoms = restored_frame['atoms']
print(f"Data identical: {np.allclose(original_atoms['x'], restored_atoms['x'])}")
```

This serialization capability is essential for data persistence, sharing between different programs, and building complex workflows that involve multiple steps of data processing.

### CSV Export

Blocks can be exported to CSV format for analysis in external tools.

```python
# Export atomic data to CSV
water_atoms.to_csv('water_atoms.csv')

# Import from CSV
imported_atoms = mp.Block.from_csv('water_atoms.csv')
print(f"Imported {imported_atoms.nrows} atoms from CSV")
```

CSV export is useful for data analysis in spreadsheet applications, statistical software, or other tools that don't directly support MolPy's data structures.

## Summary

This tutorial covered the fundamental concepts of MolPy's data model. You learned how blocks organize data by mapping variable names to arrays, providing efficient storage and fast access to molecular properties. You explored how frames organize complete molecular systems by combining multiple blocks with metadata, creating a hierarchical structure that mirrors the natural organization of molecular data.

The organized structure of blocks and frames makes molecular data easier to work with, analyze, and share. Blocks contain related variables with consistent row counts, while frames combine multiple blocks with metadata to provide a unified view of molecular systems. This organization enables efficient analysis workflows and makes your code more maintainable and understandable.

### Next Steps

Continue your MolPy journey by exploring atomistic structures to understand how molecular connectivity works, learning about element and chemistry handling for chemical calculations, mastering force field basics for molecular mechanics simulations, and understanding topology management for complex molecular graphs.

Understanding Blocks and Frames is essential for working with MolPy. These data structures provide the foundation for all molecular simulation and analysis tasks, giving you the tools to organize, manipulate, and analyze molecular data efficiently and reliably.
