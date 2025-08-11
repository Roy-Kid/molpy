# Force Field Basics

Force fields define the mathematical functions and parameters that describe molecular interactions in simulations. This tutorial teaches you how to work with MolPy's force field system to define, organize, and apply molecular mechanics parameters.

## Understanding Force Fields

Force fields are collections of mathematical functions and parameters that describe how atoms interact with each other. They include bonded interactions like bonds, angles, and dihedrals, as well as non-bonded interactions like van der Waals forces and electrostatic interactions.

MolPy's force field system organizes these parameters hierarchically. At the top level, a `ForceField` contains multiple styles for different interaction types. Each style manages multiple types of that interaction, and each type stores the specific parameters needed for calculations.

This organization makes it easy to define complex force fields, combine parameters from different sources, and apply them consistently across molecular systems. Instead of scattered parameter files, you get a unified system that understands the relationships between different interaction types.

## Working with Force Field Types

### Basic Type Structure

Each interaction type in a force field stores specific parameters and metadata. The base `Type` class provides a flexible container for these parameters.

```python
import molpy as mp

# Create an atom type with parameters
carbon_type = mp.AtomType(
    name="C",
    parms=[12.01, 1.7],  # mass, radius
    charge=0.0,
    epsilon=0.066,        # LJ parameter
    sigma=3.5             # LJ parameter
)

# Create a bond type with parameters
ch_bond = mp.BondType(
    itype=carbon_type,
    jtype=carbon_type,
    name="C-C",
    parms=[1.53, 322.0],  # equilibrium length, force constant
    bond_order=1
)

# Access type properties
print(f"Carbon type: {carbon_type.name}")
print(f"Bond type: {ch_bond.name}")
print(f"Bond parameters: {ch_bond.parms}")
print(f"Bond atom types: {ch_bond.atomtypes}")
```

Types store parameters in two ways: positional parameters in the `parms` list and keyword parameters as attributes. This design accommodates different parameter formats while maintaining a consistent interface.

### Type Collections and Management

Types are organized into collections that make it easy to manage related parameters and find specific types when needed.

```python
# Create a type container
atom_types = mp.TypeContainer()

# Add types to the container
atom_types.add(carbon_type)
atom_types.add(mp.AtomType("H", parms=[1.008, 1.2], charge=0.0))
atom_types.add(mp.AtomType("O", parms=[15.999, 1.52], charge=-0.5))

# Access types by index or name
first_type = atom_types[0]
hydrogen_type = atom_types.get("H")
oxygen_type = atom_types.get("O")

# Find types by condition
charged_types = atom_types.get_all_by(lambda t: abs(t.charge) > 0.1)
print(f"Charged types: {[t.name for t in charged_types]}")

# Iterate through all types
for atom_type in atom_types:
    print(f"{atom_type.name}: mass={atom_type.parms[0]}, charge={atom_type.charge}")
```

Type containers provide efficient access to types and support filtering operations. This makes it easy to find specific types or groups of types that meet certain criteria.

## Working with Force Field Styles

### Style Organization

Styles organize related types and provide methods for creating and managing them. Each style corresponds to a specific interaction category.

```python
# Create an atom style
atom_style = mp.AtomStyle("organic")

# Define atom types within the style
atom_style.def_type("C", parms=[12.01, 1.7], charge=0.0)
atom_style.def_type("H", parms=[1.008, 1.2], charge=0.0)
atom_style.def_type("O", parms=[15.999, 1.52], charge=-0.5)

# Create a bond style
bond_style = mp.BondStyle("harmonic")

# Define bond types
bond_style.def_type(
    itype=atom_style.get("C"),
    jtype=atom_style.get("H"),
    name="C-H",
    parms=[1.09, 340.0]  # length, force constant
)

# Access style information
print(f"Atom style: {atom_style.name}")
print(f"Number of types: {atom_style.n_types}")
print(f"Available types: {[t.name for t in atom_style.get_types()]}")
```

Styles provide a logical grouping for related types and ensure consistency in parameter definitions. They also support type creation with automatic validation and organization.

### Style Containers

Style containers manage multiple styles of the same type, allowing you to organize force field parameters into logical groups.

```python
# Create a style container
atom_styles = mp.StyleContainer()

# Add different atom styles
atom_styles.add(atom_style)
atom_styles.add(mp.AtomStyle("inorganic"))

# Access styles
organic_style = atom_styles.get("organic")
inorganic_style = atom_styles.get("inorganic")

# Find styles by condition
large_styles = atom_styles.get_all_by(lambda s: s.n_types > 5)
print(f"Large styles: {[s.name for s in large_styles]}")
```

Style containers enable you to organize force field parameters into logical groups. This is useful when working with complex force fields that combine parameters from different sources or when you need to maintain separate parameter sets for different types of systems.

## Building Complete Force Fields

### Force Field Assembly

A `ForceField` brings together all the styles and types into a unified system that can be applied to molecular structures.

```python
# Create a complete force field
ff = mp.ForceField(name="organic_ff", units="real")

# Add atom styles
ff.def_atomstyle("organic", parms=[], charge=0.0)
ff.def_bondstyle("harmonic", parms=[], bond_order=1)
ff.def_anglestyle("harmonic", parms=[], angle_type="tetrahedral")
ff.def_dihedralstyle("harmonic", parms=[], dihedral_type="proper")

# Access force field information
print(f"Force field: {ff.name}")
print(f"Units: {ff.units}")
print(f"Atom styles: {ff.n_atomstyles}")
print(f"Bond styles: {ff.n_bondstyles}")
print(f"Angle styles: {ff.n_anglestyles}")
print(f"Dihedral styles: {ff.n_dihedralstyles}")
```

The force field provides a unified interface for all interaction types. It maintains counts of styles and types, making it easy to understand the scope and complexity of your force field.

### Type Management and Access

Once you have a complete force field, you can access types and parameters through the organized structure.

```python
# Get all atom types from the force field
atom_types = ff.get_atomtypes()
print(f"Total atom types: {len(atom_types)}")

# Get specific styles
atom_style = ff.get_atomstyle("organic")
bond_style = ff.get_bondstyle("harmonic")

# Check if specific types exist
if "C" in atom_style:
    carbon_type = atom_style["C"]
    print(f"Carbon type found: {carbon_type.name}")

# Access style information
print(f"Bond style types: {bond_style.n_types}")
```

This organized access makes it easy to find and use specific parameters. You can check what types are available, access their parameters, and understand how they're organized within the force field.

## Force Field Operations

### Merging and Combining

Force fields can be combined to create more comprehensive parameter sets or to merge parameters from different sources.

```python
# Create a second force field
ff2 = mp.ForceField(name="inorganic_ff", units="real")
ff2.def_atomstyle("metal", parms=[], charge=0.0)

# Merge force fields
combined_ff = mp.ForceField.from_forcefields("combined_ff", ff, ff2)
print(f"Combined force field: {combined_ff.name}")
print(f"Total atom styles: {combined_ff.n_atomstyles}")

# Alternative merge method
ff.merge(ff2)
print(f"Merged force field atom styles: {ff.n_atomstyles}")
```

Merging is useful when you need to combine parameters from different sources or when building complex force fields that cover multiple chemical systems.

### Serialization and Persistence

Force fields can be saved and loaded, making it easy to share parameter sets and maintain consistency across different projects.

```python
# Convert force field to dictionary
ff_dict = ff.to_dict()
print(f"Serialized force field has {len(ff_dict)} keys")

# Convert back to force field
restored_ff = mp.ForceField.from_dict(ff_dict)
print(f"Restored force field: {restored_ff.name}")
print(f"Atom styles: {restored_ff.n_atomstyles}")

# Verify data integrity
print(f"Force field identical: {ff.name == restored_ff.name}")
```

This serialization capability is essential for sharing force field parameters, maintaining parameter databases, and building reproducible simulation workflows.

## Data Organization Principles

### Type and Style Naming

Use descriptive names for your types and styles that clearly indicate their purpose and content. Names like "organic", "harmonic", and "tetrahedral" make your force field self-documenting.

Organize related types within styles logically. Group atom types by chemical character, bond types by interaction model, and so on. This organization makes it easier to find specific parameters and understand the structure of your force field.

### Parameter Consistency

Maintain consistency in parameter units and formats across your force field. The `units` attribute helps track what units your parameters use, but you should also ensure that parameters within the same style use consistent formats.

Use consistent naming conventions for types and styles. This makes your force field easier to work with and reduces the chance of errors when accessing parameters.

### Force Field Structure

Organize your force field with a clear hierarchy: ForceField → Styles → Types → Parameters. Each level should have a clear purpose and consistent interface.

When working with multiple force fields, use consistent structures and naming conventions. This consistency makes it easier to combine force fields and write analysis code that works across different parameter sets.

## Summary

This tutorial covered the fundamental concepts of MolPy's force field system. You learned how types store individual interaction parameters, how styles organize related types, and how force fields bring everything together into unified systems.

The hierarchical organization of force fields makes parameter management efficient and intuitive. Types store specific parameters, styles group related types, and force fields provide unified access to all interaction parameters. This organization enables complex force field construction, parameter sharing, and consistent application across molecular systems.

### Next Steps

Continue your MolPy journey by understanding topology management for complex molecular graphs, learning about system organization for working with large assemblies, and exploring advanced analysis techniques for complex molecular systems.

Understanding force fields is essential for molecular mechanics simulations in MolPy. The organized parameter management system provides the foundation for accurate and efficient molecular dynamics calculations.
