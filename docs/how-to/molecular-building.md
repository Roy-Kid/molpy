# How to Build Molecular Systems

This guide shows practical examples of building molecular systems using MolPy's builder tools. You'll learn how to construct polymers, assemble complex molecules, and create custom building blocks.

## Building Polymer Chains

### Creating Monomer Templates

Start by defining monomer templates that can be reused to build polymers:

```python
import molpy as mp
import numpy as np

# Create a simple monomer (ethylene-like)
def create_ethylene_monomer():
    """Create an ethylene monomer template."""
    atoms = mp.Block({
        'x': [0.0, 1.34, 0.0, 0.0],
        'y': [0.0, 0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0, 0.0],
        'element': ['C', 'C', 'H', 'H'],
        'type': [1, 1, 2, 2]
    })

    # Create anchor rules for connection points
    anchor_rules = [
        mp.AnchorRule(init=0, end=1, deletes=[2, 3]),  # Left anchor
        mp.AnchorRule(init=1, end=0, deletes=[2, 3])   # Right anchor
    ]

    monomer = mp.Monomer(atoms, anchors=anchor_rules)
    return monomer

# Create monomer template
ethylene = create_ethylene_monomer()
print(f"Ethylene monomer created with {len(ethylene.anchors)} anchors")
```

### Building Polymer Chains

Use the PolymerBuilder to assemble monomers into chains:

```python
# Create polymer builder with monomer templates
monomers = {
    'ethylene': ethylene
}

builder = mp.PolymerBuilder(monomers)

# Build a simple chain
def build_ethylene_chain(n_monomers=5):
    """Build a chain of ethylene monomers."""

    # Start with first monomer
    chain = monomers['ethylene']

    # Add monomers one by one
    for i in range(1, n_monomers):
        # Connect new monomer
        new_monomer = monomers['ethylene']
        # In practice, you'd use the builder's connection logic
        # This is a simplified example

        # Combine atomic coordinates
        combined_x = np.concatenate([chain['atoms']['x'], new_monomer['atoms']['x']])
        combined_y = np.concatenate([chain['atoms']['y'], new_monomer['atoms']['y']])
        combined_z = np.concatenate([chain['atoms']['z'], new_monomer['atoms']['z']])
        combined_elements = np.concatenate([chain['atoms']['element'], new_monomer['atoms']['element']])

        # Create new chain
        chain = mp.Block({
            'x': combined_x,
            'y': combined_y,
            'z': combined_z,
            'element': combined_elements
        })

    return chain

# Build the chain
polymer_chain = build_ethylene_chain(5)
print(f"Built polymer chain with {polymer_chain.nrows} atoms")
```

## Building Complex Molecules

### Assembling Multi-Component Systems

Build complex molecules by combining different building blocks:

```python
def build_water_dimer():
    """Build a water dimer system."""

    # Create first water molecule
    water1 = mp.Block({
        'x': [0.0, 0.9572, -0.2400],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H']
    })

    # Create second water molecule (offset)
    water2 = mp.Block({
        'x': [3.0, 3.9572, 2.7600],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H']
    })

    # Combine into dimer
    dimer = mp.Block({
        'x': np.concatenate([water1['x'], water2['x']]),
        'y': np.concatenate([water1['y'], water2['y']]),
        'z': np.concatenate([water1['z'], water2['z']]),
        'element': np.concatenate([water1['element'], water2['element'])
    })

    return dimer

# Build water dimer
water_dimer = build_water_dimer()
print(f"Water dimer built with {water_dimer.nrows} atoms")

# Create frame with the dimer
dimer_frame = mp.Frame({'atoms': water_dimer})
print(f"Created frame with {dimer_frame['atoms'].nrows} atoms")
```

### Building with Custom Geometries

Create molecules with specific geometric arrangements:

```python
def build_tetrahedral_molecule(center_atom='C', terminal_atoms=['H', 'H', 'H', 'H']):
    """Build a tetrahedral molecule (e.g., methane)."""

    # Tetrahedral angles
    angles = [109.5, 109.5, 109.5, 109.5]  # degrees
    bond_length = 1.09  # Angstroms

    # Calculate positions
    positions = []
    elements = [center_atom]

    # Center atom
    positions.append([0.0, 0.0, 0.0])

    # Terminal atoms in tetrahedral arrangement
    for i, (angle, element_symbol) in enumerate(zip(angles, terminal_atoms)):
        # Simplified tetrahedral positioning
        phi = i * 90  # azimuthal angle
        theta = 109.5  # polar angle

        x = bond_length * np.sin(np.radians(theta)) * np.cos(np.radians(phi))
        y = bond_length * np.sin(np.radians(theta)) * np.sin(np.radians(phi))
        z = bond_length * np.cos(np.radians(theta))

        positions.append([x, y, z])
        elements.append(element_symbol)

    # Create block
    molecule = mp.Block({
        'x': [pos[0] for pos in positions],
        'y': [pos[1] for pos in positions],
        'z': [pos[2] for pos in positions],
        'element': elements
    })

    return molecule

# Build methane
methane = build_tetrahedral_molecule('C', ['H', 'H', 'H', 'H'])
print(f"Methane built with {methane.nrows} atoms")

# Build carbon tetrachloride
ccl4 = build_tetrahedral_molecule('C', ['Cl', 'Cl', 'Cl', 'Cl'])
print(f"CCl4 built with {ccl4.nrows} atoms")
```

## Building Crystal Structures

### Creating Unit Cells

Build crystal structures by creating and replicating unit cells:

```python
def build_simple_cubic_crystal(unit_cell, repetitions=[2, 2, 2]):
    """Build a crystal by repeating a unit cell."""

    # Get unit cell dimensions
    cell_x = max(unit_cell['atoms']['x']) - min(unit_cell['atoms']['x'])
    cell_y = max(unit_cell['atoms']['y']) - min(unit_cell['atoms']['y'])
    cell_z = max(unit_cell['atoms']['z']) - min(unit_cell['atoms']['z'])

    # Build crystal
    all_positions = []
    all_elements = []

    for i in range(repetitions[0]):
        for j in range(repetitions[1]):
            for k in range(repetitions[2]):
                # Offset for this unit cell
                offset_x = i * cell_x
                offset_y = j * cell_y
                offset_z = k * cell_z

                # Add atoms from unit cell
                for atom_idx in range(unit_cell.nrows):
                    x = unit_cell['atoms']['x'][atom_idx] + offset_x
                    y = unit_cell['atoms']['y'][atom_idx] + offset_y
                    z = unit_cell['atoms']['z'][atom_idx] + offset_z

                    all_positions.append([x, y, z])
                    all_elements.append(unit_cell['atoms']['element'][atom_idx])

    # Create crystal block
    crystal = mp.Block({
        'x': [pos[0] for pos in all_positions],
        'y': [pos[1] for pos in all_positions],
        'z': [pos[2] for pos in all_positions],
        'element': all_elements
    })

    return crystal

# Build crystal from methane unit cell
methane_crystal = build_simple_cubic_crystal(methane, [3, 3, 3])
print(f"Methane crystal built with {methane_crystal.nrows} atoms")

# Create frame with crystal
crystal_frame = mp.Frame({'atoms': methane_crystal})
print(f"Crystal frame created")
```

## Building with Constraints

### Applying Geometric Constraints

Build molecules that satisfy specific geometric constraints:

```python
def build_linear_chain(elements, bond_lengths):
    """Build a linear chain with specified bond lengths."""

    if len(elements) != len(bond_lengths) + 1:
        raise ValueError("Number of elements must be number of bonds + 1")

    positions = [[0.0, 0.0, 0.0]]  # Start at origin

    # Build chain along x-axis
    current_x = 0.0
    for i, bond_length in enumerate(bond_lengths):
        current_x += bond_length
        positions.append([current_x, 0.0, 0.0])

    # Create block
    chain = mp.Block({
        'x': [pos[0] for pos in positions],
        'y': [pos[1] for pos in positions],
        'z': [pos[2] for pos in positions],
        'element': elements
    })

    return chain

# Build a simple linear molecule (e.g., H-C-C-H)
linear_molecule = build_linear_chain(['H', 'C', 'C', 'H'], [1.09, 1.54, 1.09])
print(f"Linear molecule built with {linear_molecule.nrows} atoms")

# Build a longer chain
long_chain = build_linear_chain(['H'] + ['C'] * 5 + ['H'], [1.09] + [1.54] * 5 + [1.09])
print(f"Long chain built with {long_chain.nrows} atoms")
```

## Building Workflows

### Automated Building Process

Create automated workflows for building complex systems:

```python
def build_protein_system(sequence, template_monomers):
    """Build a protein system from amino acid sequence."""

    # This is a simplified example - real protein building is more complex
    protein_atoms = []
    protein_elements = []

    current_pos = 0.0

    for i, residue in enumerate(sequence):
        if residue in template_monomers:
            monomer = template_monomers[residue]

            # Get monomer coordinates
            monomer_x = monomer['atoms']['x']
            monomer_y = monomer['atoms']['y']
            monomer_z = monomer['atoms']['z']
            monomer_elements = monomer['atoms']['element']

            # Offset monomer position
            offset_x = [x + current_pos for x in monomer_x]

            # Add to protein
            protein_atoms.extend(zip(offset_x, monomer_y, monomer_z))
            protein_elements.extend(monomer_elements)

            # Move to next position
            current_pos += max(monomer_x) - min(monomer_x) + 3.8  # Add spacing

    # Create protein block
    protein = mp.Block({
        'x': [pos[0] for pos in protein_atoms],
        'y': [pos[1] for pos in protein_atoms],
        'z': [pos[2] for pos in protein_atoms],
        'element': protein_elements
    })

    return protein

# Example usage (simplified)
# template_monomers = {'A': alanine_monomer, 'G': glycine_monomer, ...}
# protein = build_protein_system("AGAGAG", template_monomers)
```

## Summary

This guide covered practical molecular building workflows:

- Create monomer templates for polymers
- Build polymer chains and complex molecules
- Assemble crystal structures
- Apply geometric constraints
- Automate building processes

MolPy's builder tools provide flexible ways to construct molecular systems from simple building blocks. The key is designing reusable templates and understanding how to connect them properly.

### Next Steps

Continue exploring MolPy by learning about potential energy calculations, packing algorithms, and analysis techniques.
