# Atomistic Structures

Molecular systems are more than just collections of atoms with coordinates. They have connectivity, topology, and chemical relationships that define their behavior and properties. This tutorial teaches you how to work with MolPy's atomistic classes to build, manipulate, and analyze molecular structures with full connectivity information.

## Understanding Molecular Connectivity

### Why Connectivity Matters

In molecular simulations and analysis, knowing how atoms are connected is essential. Bonds define which atoms are chemically linked, angles describe the geometry around central atoms, and dihedrals control rotational freedom. This connectivity information determines molecular flexibility, chemical reactivity, and physical properties.

MolPy's atomistic classes provide a systematic way to represent these relationships. Instead of working with disconnected coordinate arrays, you can create structures that understand the chemical bonds and geometric relationships between atoms. This makes it possible to analyze molecular properties, perform structural modifications, and generate force field parameters automatically.

### The Hierarchy of Molecular Entities

Molecular connectivity follows a natural hierarchy. Atoms are the fundamental units, each with specific properties like element type, mass, and charge. Bonds connect pairs of atoms, defining the primary structure. Angles involve three atoms and describe local geometry. Dihedrals involve four atoms and control rotational degrees of freedom. Improper dihedrals maintain planarity in groups like aromatic rings.

MolPy represents each level of this hierarchy with dedicated classes. The `Atom` class handles individual atomic properties, `Bond` manages pairwise connections, `Angle` describes three-atom geometries, and `Dihedral` controls four-atom rotations. The `Atomistic` class brings all these together into a complete molecular structure.

## Working with Individual Atoms

### Creating and Understanding Atoms

The `Atom` class represents individual atoms with all their properties. Each atom can have coordinates, chemical properties, physical properties, and custom attributes. This flexibility allows you to represent everything from simple atomic systems to complex biomolecules with specialized properties.

```python
import molpy as mp
import numpy as np

# Create individual atoms with different properties
carbon_atom = mp.Atom(
    name="C",
    x=0.0,
    y=0.0,
    z=0.0,
    element="C",
    mass=12.01,
    charge=0.0,
    type=1
)

oxygen_atom = mp.Atom(
    name="O",
    x=1.4,
    y=0.0,
    z=0.0,
    element="O",
    mass=15.999,
    charge=-0.5,
    type=2
)

# Access atom properties
print(f"Carbon atom: {carbon_atom.name}, mass: {carbon_atom.mass} amu")
print(f"Oxygen atom: {oxygen_atom.name}, charge: {oxygen_atom.charge} e")
print(f"Distance between atoms: {np.sqrt((oxygen_atom.x - carbon_atom.x)**2)} Å")
```

Atoms behave like dictionaries, so you can add, modify, or remove properties as needed. This makes them flexible for different types of molecular systems and analysis requirements. You can store experimental data, simulation parameters, or custom properties alongside the standard atomic information.

### Atom Collections and Operations

When working with multiple atoms, you often need to perform operations on collections. MolPy provides efficient ways to work with groups of atoms, from simple property extraction to complex geometric calculations.

```python
# Create a collection of atoms
atoms = [
    mp.Atom(name="H1", x=0.0, y=1.0, z=0.0, element="H", mass=1.008, type=3),
    mp.Atom(name="H2", x=0.0, y=-1.0, z=0.0, element="H", mass=1.008, type=3),
    mp.Atom(name="H3", x=1.0, y=0.0, z=0.0, element="H", mass=1.008, type=3),
    mp.Atom(name="H4", x=-1.0, y=0.0, z=0.0, element="H", mass=1.008, type=3)
]

# Extract properties across all atoms
positions = np.array([[atom.x, atom.y, atom.z] for atom in atoms])
masses = np.array([atom.mass for atom in atoms])
elements = [atom.element for atom in atoms]

# Calculate center of mass
center_of_mass = np.average(positions, weights=masses, axis=0)
print(f"Center of mass: {center_of_mass}")

# Find atoms by property
hydrogens = [atom for atom in atoms if atom.element == "H"]
print(f"Number of hydrogen atoms: {len(hydrogens)}")
```

This approach gives you the flexibility to work with atoms individually or as groups, depending on your analysis needs. You can extract specific properties for vectorized calculations, filter atoms by chemical properties, or perform geometric operations on selected subsets.

## Building Molecular Bonds

### Understanding Bond Connectivity

Bonds represent the primary chemical connections between atoms. In MolPy, a `Bond` object connects exactly two atoms and can carry additional information like bond type, order, or force field parameters. Bonds are fundamental for understanding molecular structure and for many types of analysis.

```python
# Create bonds between atoms
c_o_bond = mp.Bond(
    itom=carbon_atom,
    jtom=oxygen_atom,
    type=1,
    order=2,  # Double bond
    length=1.4
)

# Bonds automatically sort atoms by ID for consistency
print(f"Bond connects: {c_o_bond.itom.name} and {c_o_bond.jtom.name}")
print(f"Bond type: {c_o_bond.type}, order: {c_o_bond.order}")
print(f"Bond length: {c_o_bond.length} Å")
```

Bonds provide several important features. They automatically sort the connected atoms by ID, ensuring consistent ordering regardless of how you create them. They can store bond-specific properties like type, order, or force field parameters. And they maintain references to the connected atoms, allowing you to traverse molecular connectivity.

### Working with Bond Collections

When building complex molecules, you'll work with many bonds simultaneously. MolPy provides efficient ways to create, manage, and analyze bond collections.

```python
# Create multiple bonds for a molecule
bonds = []

# Carbon-hydrogen bonds
for i, h_atom in enumerate(hydrogens):
    bond = mp.Bond(
        itom=carbon_atom,
        jtom=h_atom,
        type=2,  # C-H bond type
        order=1,  # Single bond
        length=1.09
    )
    bonds.append(bond)

# Add the C-O bond
bonds.append(c_o_bond)

# Analyze bond collection
print(f"Total bonds: {len(bonds)}")
print(f"Bond types: {set(bond.type for bond in bonds)}")

# Find bonds by type
ch_bonds = [bond for bond in bonds if bond.type == 2]
print(f"Number of C-H bonds: {len(ch_bonds)}")
```

This approach lets you build complex molecular structures systematically. You can create bonds based on geometric proximity, chemical rules, or predefined connectivity patterns. The bond collection provides a foundation for more complex analysis like ring detection, path finding, or force field parameter assignment.

## Managing Molecular Angles

### Understanding Angle Geometry

Angles describe the geometry around central atoms and are crucial for understanding molecular shape and flexibility. In MolPy, an `Angle` object involves three atoms and can store geometric and force field information. Angles are essential for maintaining proper molecular geometry during simulations.

```python
# Create angles for the molecular structure
angles = []

# H-C-H angles (around carbon)
for i in range(len(hydrogens)):
    for j in range(i + 1, len(hydrogens)):
        angle = mp.Angle(
            itom=hydrogens[i],
            jtom=carbon_atom,  # Central atom
            ktom=hydrogens[j],
            type=1,  # H-C-H angle type
            value=109.5  # Tetrahedral angle in degrees
        )
        angles.append(angle)

# H-C-O angles
for h_atom in hydrogens:
    angle = mp.Angle(
        itom=h_atom,
        jtom=carbon_atom,
        ktom=oxygen_atom,
        type=2,  # H-C-O angle type
        value=109.5
    )
    angles.append(angle)

print(f"Total angles: {len(angles)}")
print(f"Angle types: {set(angle.type for angle in angles)}")
```

Angles provide important geometric constraints. They help maintain proper molecular geometry during energy minimization and molecular dynamics. They can store equilibrium values and force constants for force field calculations. And they provide a way to analyze molecular flexibility and conformational preferences.

## Controlling Molecular Dihedrals

### Understanding Dihedral Rotations

Dihedrals control the rotational freedom around bonds and are essential for describing molecular conformations. In MolPy, a `Dihedral` object involves four atoms and defines the rotation around the central bond. Dihedrals are crucial for understanding molecular flexibility and for generating diverse conformations.

```python
# Create dihedrals for rotational degrees of freedom
dihedrals = []

# H-C-C-O dihedral (rotation around C-C bond)
for h_atom in hydrogens:
    dihedral = mp.Dihedral(
        itom=h_atom,
        jtom=carbon_atom,
        ktom=carbon_atom,  # Same as jtom for this example
        ltom=oxygen_atom,
        type=1,  # H-C-C-O dihedral type
        value=0.0  # Initial dihedral angle
    )
    dihedrals.append(dihedral)

print(f"Total dihedrals: {len(dihedrals)}")
print(f"Dihedral types: {set(dihedral.type for dihedral in dihedrals)}")
```

Dihedrals provide several important capabilities. They define the rotational degrees of freedom in your molecule, allowing you to explore different conformations. They can store force field parameters for proper energy calculations. And they provide a framework for conformational analysis and structure generation.

## Building Complete Molecular Structures

### The Atomistic Class

The `Atomistic` class brings together all the molecular components into a unified structure. It manages atoms, bonds, angles, and dihedrals as a cohesive system, providing methods for analysis, modification, and conversion to other formats.

```python
# Create a complete molecular structure
molecule = mp.Atomistic()

# Add atoms
molecule.add_atom(carbon_atom)
molecule.add_atom(oxygen_atom)
for h_atom in hydrogens:
    molecule.add_atom(h_atom)

# Add bonds
for bond in bonds:
    molecule.add_bond(bond)

# Add angles
for angle in angles:
    molecule.add_angle(angle)

# Add dihedrals
for dihedral in dihedrals:
    molecule.add_dihedral(dihedral)

print(f"Molecule contains:")
print(f"  {len(molecule.atoms)} atoms")
print(f"  {len(molecule.bonds)} bonds")
print(f"  {len(molecule.angles)} angles")
print(f"  {len(molecule.dihedrals)} dihedrals")
```

The `Atomistic` class provides a comprehensive interface for working with molecular structures. It maintains the relationships between all components, ensures consistency, and provides methods for common operations like adding or removing atoms, analyzing connectivity, and converting to other formats.

### Molecular Analysis and Properties

Once you have a complete molecular structure, you can perform comprehensive analysis to understand its properties and behavior.

```python
def analyze_molecular_structure(molecule):
    """Analyze the complete molecular structure."""
    print("=== Molecular Structure Analysis ===")

    # Basic composition
    elements = [atom.element for atom in molecule.atoms]
    element_counts = {}
    for element in elements:
        element_counts[element] = element_counts.get(element, 0) + 1

    print(f"Molecular formula: {''.join(f'{element}{count}' for element, count in element_counts.items())}")
    print(f"Total mass: {sum(atom.mass for atom in molecule.atoms):.3f} amu")
    print(f"Total charge: {sum(atom.charge for atom in molecule.atoms):.3f} e")

    # Connectivity analysis
    print(f"\nConnectivity:")
    print(f"  Bonds per atom: {len(molecule.bonds) / len(molecule.atoms):.1f}")

    # Check for rings (simplified)
    if len(molecule.bonds) >= len(molecule.atoms):
        print("  Structure may contain rings")
    else:
        print("  Structure is acyclic")

# Analyze the molecule
analyze_molecular_structure(molecule)
```

This analysis provides a comprehensive view of your molecular system. It shows the chemical composition, physical properties, and structural characteristics. This information is essential for understanding how your molecule will behave in simulations, what properties it will have, and how it compares to experimental data.

## Converting Between Formats

### Frame Conversion

The `Atomistic` class can convert to and from other MolPy formats, making it easy to integrate with the rest of the framework.

```python
# Convert to Frame format
frame = molecule.to_frame()
print(f"Converted to frame with {len(frame.blocks())} blocks")

# Convert back to Atomistic
restored_molecule = mp.Atomistic.from_frame(frame)
print(f"Restored molecule has {len(restored_molecule.atoms)} atoms")

# Verify conversion integrity
print(f"Conversion successful: {len(molecule.atoms) == len(restored_molecule.atoms)}")
```

This conversion capability is essential for integrating atomistic structures with other MolPy components. You can use the rich connectivity information of `Atomistic` for analysis and modification, then convert to `Frame` format for trajectory analysis or I/O operations.

## Best Practices

### Structure Building Guidelines

When building molecular structures, start with the core atoms and work outward. Begin with the central atoms that have the most connections, then add peripheral atoms. This approach makes it easier to maintain proper geometry and connectivity.

Always validate your structures after creation and modification. Check that bond counts are reasonable for each element type, that angles are physically possible, and that the overall structure makes chemical sense. This validation prevents problems in subsequent analysis steps.

Use consistent naming conventions for atoms, bonds, and other components. Descriptive names make your structures easier to understand and debug. Consider using systematic names that indicate the position and type of each component.

### Integration with Other Components

The `Atomistic` class is designed to work seamlessly with other MolPy components. You can convert to `Frame` format for trajectory analysis, use with force field classes for parameter assignment, and integrate with topology classes for graph-based analysis.

When converting between formats, always verify that the conversion preserves the essential information you need. Some formats may not preserve all the detailed connectivity information, so choose the appropriate format for your specific use case.

## Summary

This tutorial covered the comprehensive atomistic modeling capabilities of MolPy. You learned how to create and work with individual atoms, build molecular connectivity through bonds, angles, and dihedrals, and assemble complete molecular structures using the `Atomistic` class.

You discovered how to analyze molecular properties, validate structural integrity, and perform modifications while maintaining consistency. You explored the conversion capabilities that allow integration with other MolPy components, and you learned best practices for building and working with complex molecular systems.

### Key Takeaways

Molecular connectivity is essential for understanding chemical behavior and physical properties. The atomistic classes provide a systematic way to represent and work with this connectivity, from simple bonds to complex molecular assemblies.

The `Atomistic` class brings together all molecular components into a unified structure, providing methods for analysis, modification, and conversion. This integration makes it possible to build complex workflows that combine structure building, analysis, and simulation.

### Next Steps

Continue your MolPy journey by exploring element and chemistry handling for chemical calculations, mastering force field basics for molecular mechanics simulations, understanding topology management for complex molecular graphs, and learning about system organization for working with large assemblies.

The atomistic classes provide the foundation for building and analyzing molecular structures in MolPy. Understanding how to work with these classes enables you to create complex molecular systems, analyze their properties, and integrate them into larger simulation and analysis workflows.
