# Selection & Filtering

Molecular analysis often requires selecting specific atoms or groups based on their properties. This tutorial teaches you how to work with MolPy's selection system to create complex queries and filter molecular data efficiently.

## Understanding Selection

Selection in MolPy is based on boolean masks that identify atoms or other entities meeting specific criteria. The system uses a predicate-based approach where each selection criterion produces a boolean array, and these can be combined using logical operators to create complex queries.

The selection system is built around the `MaskPredicate` class, which provides a clean interface for creating and combining selection criteria. Each predicate produces a boolean mask when applied to a block, and predicates can be combined using `&` (AND), `|` (OR), and `~` (NOT) operators.

This approach makes it easy to build complex selections from simple building blocks. Instead of writing complex conditional logic, you can combine predicates in intuitive ways to create sophisticated queries.

## Basic Selection Predicates

### Atom Type Selection

`AtomType` predicate selects atoms based on their type identifier. This is useful for selecting specific chemical elements or atom types in your system.

```python
import molpy as mp
import numpy as np

# Create a water molecule block
water_atoms = mp.Block({
    'x': [0.0, 0.9572, -0.2400],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['O', 'H', 'H'],
    'type': [1, 2, 2]  # O=1, H=2
})

# Create selection predicates
oxygen_sel = mp.AtomType(atom_type=1, field="type")
hydrogen_sel = mp.AtomType(atom_type=2, field="type")

# Apply selections
oxygen_atoms = oxygen_sel(water_atoms)
hydrogen_atoms = hydrogen_sel(water_atoms)

print(f"Oxygen atoms: {oxygen_atoms.nrows}")
print(f"Hydrogen atoms: {hydrogen_atoms.nrows}")
print(f"Oxygen positions: {oxygen_atoms['x']}")
```

The `AtomType` predicate takes an atom type value and the field name to match against. When applied to a block, it returns a new block containing only the atoms that match the criterion.

### Atom Index Selection

`AtomIndex` predicate selects atoms based on their index values. This is useful when you need to select specific atoms by their position in the system.

```python
# Create index-based selections
first_atom = mp.AtomIndex(indices=[0])
last_two_atoms = mp.AtomIndex(indices=[1, 2])

# Apply selections
first_atom_block = first_atom(water_atoms)
last_two_block = last_two_atoms(water_atoms)

print(f"First atom: {first_atom_block['element']}")
print(f"Last two atoms: {last_two_block['element']}")
```

Index selection is particularly useful when you need to work with specific atoms identified by their position or when you want to create custom groupings based on atom order.

## Combining Selection Predicates

### Logical Operations

Selection predicates can be combined using logical operators to create complex queries. This enables sophisticated atom selection based on multiple criteria.

```python
# Create a more complex system (methane)
methane_atoms = mp.Block({
    'x': [0.0, 0.0, 0.0, 0.0, 0.0],
    'y': [0.0, 0.0, 0.0, 0.0, 0.0],
    'z': [0.0, 1.09, -1.09, 0.0, 0.0],
    'element': ['C', 'H', 'H', 'H', 'H'],
    'type': [1, 2, 2, 2, 2],  # C=1, H=2
    'charge': [0.0, 0.0, 0.0, 0.0, 0.0]
})

# Create selection predicates
carbon_sel = mp.AtomTypeSelection(atom_type=1, field="type")
hydrogen_sel = mp.AtomTypeSelection(atom_type=2, field="type")
top_hydrogens = mp.AtomIndexSelection(indices=[1, 2])  # Top two hydrogens

# Combine selections using logical operators
carbon_and_top_h = carbon_sel | top_hydrogens
non_carbon = ~carbon_sel
carbon_or_bottom_h = carbon_sel | mp.AtomIndexSelection(indices=[3, 4])

# Apply combined selections
selected_atoms = carbon_and_top_h(methane_atoms)
non_carbon_atoms = non_carbon(methane_atoms)
carbon_or_bottom = carbon_or_bottom_h(methane_atoms)

print(f"Carbon + top H: {selected_atoms.nrows} atoms")
print(f"Non-carbon atoms: {non_carbon_atoms.nrows} atoms")
print(f"Carbon + bottom H: {carbon_or_bottom.nrows} atoms")
```

The logical operators work intuitively: `&` selects atoms that meet both criteria, `|` selects atoms that meet either criterion, and `~` inverts the selection.

### Complex Selection Patterns

You can build complex selection patterns by combining multiple predicates in various ways.

```python
def create_complex_selection():
    """Create a complex selection pattern."""
    # Define various selection criteria
    carbon = mp.AtomTypeSelection(atom_type=1, field="type")
    hydrogen = mp.AtomTypeSelection(atom_type=2, field="type")
    top_layer = mp.AtomIndexSelection(indices=[1, 2])
    bottom_layer = mp.AtomIndexSelection(indices=[3, 4])

    # Complex selection: carbon OR (hydrogen AND (top OR bottom))
    complex_sel = carbon | (hydrogen & (top_layer | bottom_layer))

    return complex_sel

# Apply complex selection
complex_sel = create_complex_selection()
result = complex_sel(methane_atoms)

print(f"Complex selection result: {result.nrows} atoms")
print(f"Selected elements: {result['element']}")
```

This pattern demonstrates how you can build sophisticated queries by combining simple predicates. The resulting selection can be as complex as needed while remaining readable and maintainable.

## Selection Applications

### Data Analysis

Selection predicates are powerful tools for analyzing specific subsets of your molecular data.

```python
def analyze_selected_atoms(block, selection, name=""):
    """Analyze properties of selected atoms."""
    selected = selection(block)

    print(f"=== {name} Analysis ===")
    print(f"Number of atoms: {selected.nrows}")

    if 'element' in selected:
        elements = selected['element']
        element_counts = {}
        for element in elements:
            element_counts[element] = element_counts.get(element, 0) + 1
        print(f"Element composition: {element_counts}")

    if 'charge' in selected:
        total_charge = np.sum(selected['charge'])
        print(f"Total charge: {total_charge:.3f}")

    if 'x' in selected and 'y' in selected and 'z' in selected:
        positions = np.column_stack([selected['x'], selected['y'], selected['z']])
        center = np.mean(positions, axis=0)
        print(f"Center of mass: {center}")

# Analyze different selections
analyze_selected_atoms(methane_atoms, carbon_sel, "Carbon")
analyze_selected_atoms(methane_atoms, hydrogen_sel, "Hydrogen")
analyze_selected_atoms(methane_atoms, top_hydrogens, "Top Hydrogens")
```

This analysis shows how selection can be used to examine specific parts of your molecular system. You can analyze properties, calculate centers of mass, or perform other calculations on selected subsets.

### Conditional Processing

Selection predicates enable conditional processing where different operations are applied to different parts of your system.

```python
def process_by_selection(block):
    """Process different atom types differently."""
    carbon = mp.AtomType(atom_type=1, field="type")
    hydrogen = mp.AtomType(atom_type=2, field="type")

    # Process carbon atoms
    carbon_atoms = carbon(block)
    if carbon_atoms.nrows > 0:
        print("Processing carbon atoms:")
        print(f"  Count: {carbon_atoms.nrows}")
        print(f"  Positions: {carbon_atoms[['x', 'y', 'z']]}")

    # Process hydrogen atoms
    hydrogen_atoms = hydrogen(block)
    if hydrogen_atoms.nrows > 0:
        print("Processing hydrogen atoms:")
        print(f"  Count: {hydrogen_atoms.nrows}")
        print(f"  Positions: {hydrogen_atoms[['x', 'y', 'z']]}")

# Process the methane system
process_by_selection(methane_atoms)
```

This pattern is useful when you need to apply different algorithms or parameters to different types of atoms in your system.

## Selection Best Practices

### Predicate Design

Design your selection predicates to be clear and reusable. Use descriptive names and organize predicates logically.

```python
# Good: Clear, descriptive predicates
carbon_atoms = mp.AtomType(atom_type=1, field="type")
hydrogen_atoms = mp.AtomType(atom_type=2, field="type")
terminal_atoms = mp.AtomIndex(indices=[1, 2, 3, 4])

# Good: Logical combinations
carbon_or_terminal = carbon_atoms | terminal_atoms
non_terminal_carbon = carbon_atoms & ~terminal_atoms
```

### Selection Efficiency

Consider the efficiency of your selections, especially for large systems. Combine predicates logically to minimize the number of operations.

```python
# Efficient: Combine predicates first, then apply
efficient_sel = (carbon_sel | hydrogen_sel) & top_hydrogens
result = efficient_sel(methane_atoms)

# Less efficient: Apply selections separately
carbon_result = carbon_sel(methane_atoms)
hydrogen_result = hydrogen_sel(methane_atoms)
combined = mp.Block({
    'x': np.concatenate([carbon_result['x'], hydrogen_result['x']]),
    'y': np.concatenate([carbon_result['y'], hydrogen_result['y']]),
    'z': np.concatenate([carbon_result['z'], hydrogen_result['z']])
})
```

### Selection Validation

Always validate your selections to ensure they produce the expected results.

```python
def validate_selection(block, selection, expected_count, name=""):
    """Validate that selection produces expected results."""
    result = selection(block)
    actual_count = result.nrows

    print(f"{name}: Expected {expected_count}, Got {actual_count}")
    if actual_count == expected_count:
        print("  ✓ Selection validated")
    else:
        print("  ✗ Selection validation failed")

    return actual_count == expected_count

# Validate selections
validate_selection(methane_atoms, carbon_sel, 1, "Carbon selection")
validate_selection(methane_atoms, hydrogen_sel, 4, "Hydrogen selection")
validate_selection(methane_atoms, top_hydrogens, 2, "Top hydrogen selection")
```

## Summary

This tutorial covered the fundamental concepts of MolPy's selection system. You learned how to create basic selection predicates, combine them using logical operators, and apply them to molecular data for analysis and processing.

The predicate-based selection system provides a powerful and intuitive way to filter molecular data. By combining simple predicates with logical operators, you can create complex queries that select exactly the atoms you need for your analysis. This system makes it easy to work with specific subsets of your molecular system while maintaining clean, readable code.

### Next Steps

Continue your MolPy journey by learning about trajectory analysis, understanding region-based operations, exploring box and boundary conditions, and mastering element and chemistry handling.

Understanding selection and filtering is essential for effective molecular analysis in MolPy. The predicate-based system provides the foundation for sophisticated data filtering and analysis workflows.
