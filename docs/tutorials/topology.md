# Topology Management

Molecular topology describes how atoms are connected to form molecules. This tutorial teaches you how to work with MolPy's Topology class to create, analyze, and manipulate molecular connectivity graphs.

## Understanding Molecular Topology

Molecular topology represents the connectivity between atoms as a mathematical graph. Atoms are vertices, bonds are edges, and the overall structure describes how molecules are assembled from their constituent parts.

MolPy's Topology class extends the igraph library to provide specialized functionality for molecular systems. It automatically detects higher-order connectivity like angles, dihedrals, and impropers from the basic bond network. This makes it easy to analyze molecular structure, identify functional groups, and understand chemical relationships.

The graph-based approach provides powerful tools for topology analysis. You can find paths between atoms, identify rings, detect functional groups, and perform complex structural queries. Instead of manually tracking connectivity, the topology system automatically maintains these relationships.

## Creating Molecular Topologies

### Basic Topology Construction

Topologies are built by adding atoms and bonds. The system automatically maintains the graph structure and provides access to various connectivity measures.

```python
import molpy as mp
import numpy as np

# Create a topology for a water molecule
water_topo = mp.Topology()

# Add atoms (vertices)
water_topo.add_atoms(3)  # 3 atoms: O, H, H

# Add bonds (edges)
water_topo.add_bond(0, 1)  # O-H bond
water_topo.add_bond(0, 2)  # O-H bond

# Check topology properties
print(f"Water topology: {water_topo.n_atoms} atoms, {water_topo.n_bonds} bonds")
print(f"Atoms: {water_topo.atoms}")
print(f"Bonds: {water_topo.bonds}")
```

The topology automatically tracks the number of atoms and bonds. The `atoms` property gives you the vertex indices, and the `bonds` property provides the edge list as pairs of atom indices.

### Higher-Order Connectivity

Topology automatically detects higher-order connectivity patterns like angles and dihedrals from the bond network.

```python
# Check for angles and dihedrals
print(f"Angles: {water_topo.n_angles}")
print(f"Dihedrals: {water_topo.n_dihedrals}")

# Get specific connectivity patterns
if water_topo.n_angles > 0:
    angles = water_topo.angles
    print(f"Angle patterns: {angles}")

# Create a more complex topology (propane)
propane_topo = mp.Topology()
propane_topo.add_atoms(11)  # C3H8

# Add C-C bonds
propane_topo.add_bond(0, 1)  # C1-C2
propane_topo.add_bond(1, 2)  # C2-C3

# Add C-H bonds
for i in range(3):  # For each carbon
    for j in range(3, 11):  # Add hydrogens
        if j < 6:  # First carbon gets 3 H
            if i == 0:
                propane_topo.add_bond(i, j)
        elif j < 9:  # Second carbon gets 2 H
            if i == 1:
                propane_topo.add_bond(i, j)
        else:  # Third carbon gets 3 H
            if i == 2:
                propane_topo.add_bond(i, j)

print(f"Propane topology: {propane_topo.n_atoms} atoms, {propane_topo.n_bonds} bonds")
print(f"Angles: {propane_topo.n_angles}")
print(f"Dihedrals: {propane_topo.n_dihedrals}")
```

The topology automatically identifies angles (three-atom patterns) and dihedrals (four-atom patterns) from the bond network. This provides a complete picture of molecular connectivity without requiring manual specification.

## Topology Analysis

### Connectivity Patterns

Topology provides various ways to analyze molecular connectivity and identify structural patterns.

```python
def analyze_topology(topo):
    """Analyze topology connectivity patterns."""
    print("=== Topology Analysis ===")

    print(f"Atoms: {topo.n_atoms}")
    print(f"Bonds: {topo.n_bonds}")
    print(f"Angles: {topo.n_angles}")
    print(f"Dihedrals: {topo.n_dihedrals}")

    # Analyze bond connectivity
    bonds = topo.bonds
    if len(bonds) > 0:
        # Count bonds per atom
        bond_counts = {}
        for bond in bonds:
            for atom in bond:
                bond_counts[atom] = bond_counts.get(atom, 0) + 1

        print(f"Bond counts per atom: {bond_counts}")

        # Find terminal atoms (1 bond)
        terminal_atoms = [atom for atom, count in bond_counts.items() if count == 1]
        print(f"Terminal atoms: {terminal_atoms}")

        # Find branching atoms (3+ bonds)
        branching_atoms = [atom for atom, count in bond_counts.items() if count >= 3]
        print(f"Branching atoms: {branching_atoms}")

# Analyze the propane topology
analyze_topology(propane_topo)
```

This analysis reveals important structural information about your molecule. Terminal atoms are often hydrogens or functional groups, while branching atoms indicate structural complexity.

### Structural Features

Topology can identify various structural features that are important for understanding molecular properties.

```python
def identify_structural_features(topo):
    """Identify important structural features in the topology."""
    print("=== Structural Features ===")

    bonds = topo.bonds
    if len(bonds) == 0:
        return

    # Check for rings (simplified)
    n_atoms = topo.n_atoms
    n_bonds = topo.n_bonds

    if n_bonds >= n_atoms:
        print("Structure may contain rings")
    else:
        print("Structure is acyclic")

    # Check for linear chains
    bond_counts = {}
    for bond in bonds:
        for atom in bond:
            bond_counts[atom] = bond_counts.get(atom, 0) + 1

    linear_atoms = [atom for atom, count in bond_counts.items() if count == 2]
    if len(linear_atoms) >= 3:
        print(f"Linear chain detected with {len(linear_atoms)} atoms")

    # Check for branching
    branching_atoms = [atom for atom, count in bond_counts.items() if count >= 3]
    if branching_atoms:
        print(f"Branching points at atoms: {branching_atoms}")

# Identify features in propane
identify_structural_features(propane_topo)
```

These structural features are crucial for understanding molecular behavior. Ring systems have different properties than linear chains, and branching affects molecular flexibility and reactivity.

## Topology Manipulation

### Adding and Removing Connectivity

Topology provides methods for modifying molecular connectivity while maintaining graph consistency.

```python
# Add new atoms and bonds
water_topo.add_atoms(1)  # Add a fourth atom
water_topo.add_bond(0, 3)  # Connect to oxygen

# Remove bonds
water_topo.delete_bond((0, 3))  # Remove the bond we just added

# Add multiple bonds at once
new_bonds = [(1, 3), (2, 3)]
water_topo.add_bonds(new_bonds)

print(f"Modified water topology: {water_topo.n_atoms} atoms, {water_topo.n_bonds} bonds")
```

The topology automatically updates all connectivity measures when you modify the structure. This ensures that angles, dihedrals, and other patterns remain consistent with the current bond network.

### Topology Combination

Topologies can be combined to create larger molecular systems or to merge separate molecular fragments.

```python
# Create a second topology
methane_topo = mp.Topology()
methane_topo.add_atoms(5)  # CH4
methane_topo.add_bond(0, 1)  # C-H bonds
methane_topo.add_bond(0, 2)
methane_topo.add_bond(0, 3)
methane_topo.add_bond(0, 4)

# Combine topologies
combined_topo = water_topo.union(methane_topo)
print(f"Combined topology: {combined_topo.n_atoms} atoms, {combined_topo.n_bonds} bonds")
```

Combining topologies is useful when building complex molecular systems from simpler components or when analyzing interactions between different molecular fragments.

## Working with Topology Data

### Accessing Connectivity Information

Topology provides various ways to access connectivity information for analysis and manipulation.

```python
def explore_connectivity(topo):
    """Explore the connectivity information in a topology."""
    print("=== Connectivity Exploration ===")

    # Get all atoms and their properties
    atoms = topo.atoms
    print(f"Atom indices: {atoms}")

    # Get all bonds
    bonds = topo.bonds
    print(f"Bond pairs: {bonds}")

    # Get angles if they exist
    if topo.n_angles > 0:
        angles = topo.angles
        print(f"Angle patterns: {angles}")

    # Get dihedrals if they exist
    if topo.n_dihedrals > 0:
        dihedrals = topo.dihedrals
        print(f"Dihedral patterns: {dihedrals}")

    # Get impropers if they exist
    impropers = topo.improper
    if len(impropers) > 0:
        print(f"Improper patterns: {impropers}")

# Explore the propane topology
explore_connectivity(propane_topo)
```

This exploration reveals the complete connectivity structure of your molecule. Each connectivity pattern provides information about molecular geometry and flexibility.

### Topology Validation

Topology provides methods to validate the consistency and reasonableness of molecular connectivity.

```python
def validate_topology(topo):
    """Validate topology consistency and reasonableness."""
    print("=== Topology Validation ===")

    # Check basic properties
    n_atoms = topo.n_atoms
    n_bonds = topo.n_bonds

    print(f"Atoms: {n_atoms}, Bonds: {n_bonds}")

    # Check for reasonable bond counts
    if n_bonds > n_atoms * 2:
        print("Warning: High bond count - may indicate errors")
    elif n_bonds < n_atoms - 1:
        print("Warning: Low bond count - structure may be disconnected")

    # Check for isolated atoms
    bonds = topo.bonds
    if len(bonds) > 0:
        connected_atoms = set()
        for bond in bonds:
            connected_atoms.update(bond)

        isolated_atoms = set(range(n_atoms)) - connected_atoms
        if isolated_atoms:
            print(f"Warning: Isolated atoms: {isolated_atoms}")
        else:
            print("✓ All atoms are connected")

    # Check for duplicate bonds
    bond_set = set()
    for bond in bonds:
        sorted_bond = tuple(sorted(bond))
        if sorted_bond in bond_set:
            print(f"Warning: Duplicate bond: {bond}")
        bond_set.add(sorted_bond)

# Validate the propane topology
validate_topology(propane_topo)
```

This validation helps ensure that your topology makes chemical sense and identifies potential errors in connectivity specification.

## Data Organization Principles

### Topology Structure

Organize your topology with a clear hierarchy: atoms → bonds → higher-order patterns. Start with the basic atomic structure, then add bonds to create connectivity, and let the system automatically identify higher-order patterns.

Use consistent atom indexing throughout your topology. Atom indices should correspond to the order in which atoms were added, and bond specifications should use these consistent indices.

### Connectivity Management

When building complex topologies, consider the logical structure of your molecule. Group related atoms together and add bonds systematically to avoid errors and maintain clarity.

Use the topology's automatic pattern detection rather than manually specifying angles and dihedrals. This ensures consistency and reduces the chance of errors.

### Topology Validation

Always validate your topology after construction to ensure it makes chemical sense. Check for isolated atoms, unreasonable bond counts, and other potential issues.

Use the topology's built-in validation methods and add custom checks for specific chemical requirements in your system.

## Summary

This tutorial covered the fundamental concepts of MolPy's topology management system. You learned how to create molecular topologies by adding atoms and bonds, how the system automatically detects higher-order connectivity patterns, and how to analyze and manipulate molecular structure.

The graph-based approach to topology provides powerful tools for understanding molecular connectivity. Automatic pattern detection, structural analysis, and topology manipulation enable complex molecular modeling workflows. This system makes it easy to build, analyze, and modify molecular structures while maintaining consistency.

### Next Steps

Continue your MolPy journey by learning about system organization for working with large assemblies, exploring advanced analysis techniques for complex molecular systems, and understanding how to integrate topology with other MolPy components.

Understanding topology management is essential for working with molecular connectivity in MolPy. The graph-based system provides the foundation for structural analysis and molecular modeling.
