# Custom Selection Predicates: Building Complex Queries

This chapter explains how to create custom selection predicates in MolPy. Selection predicates are the foundation of MolPy's filtering system, providing a functional and composable way to select atoms and other entities based on various criteria.

## Understanding Selection Predicates

### What Are Predicates?

In MolPy, a **selection predicate** is a function that takes a `Block` of data and returns a boolean mask. This mask can then be used to filter the data or combine with other predicates using logical operators.

```python
class MaskPredicate(ABC):
    """Boolean mask producer combinable with &, |, ~."""

    @abstractmethod
    def mask(self, block: "Block") -> np.ndarray: ...

    def __call__(self, block: "Block") -> "Block":
        return block[self.mask(block)]

    # Compositional logic
    def __and__(self, other: "MaskPredicate") -> "MaskPredicate":
        return _And(self, other)
    def __or__(self, other: "MaskPredicate") -> "MaskPredicate":
        return _Or(self, other)
    def __invert__(self) -> "MaskPredicate":
        return _Not(self)
```

**Key Insight**: Predicates are **functions that produce boolean masks**. They don't modify data; they create filters that can be applied to data.

### The Predicate Philosophy

MolPy's selection system is based on **functional composition** rather than object-oriented selection:

1. **Predicates are functions**: Each predicate is a callable that produces a mask
2. **Composable**: Predicates can be combined using logical operators (`&`, `|`, `~`)
3. **Immutable**: Predicates don't change data, they filter it
4. **Efficient**: Boolean masks enable vectorized operations

### How Predicates Work

```python
# Create a predicate
carbon_predicate = ElementSelection("C")

# Apply to data (produces a mask)
block = Frame(...)['atoms']
mask = carbon_predicate.mask(block)  # Returns: [True, False, True, False, ...]

# Use the mask to filter data
selected_atoms = block[mask]  # Only carbon atoms

# Or use the predicate directly
selected_atoms = carbon_predicate(block)  # Same result
```

## Built-in Selection Predicates

### AtomTypeSelection

Selects atoms by their type (integer or string):

```python
class AtomTypeSelection(MaskPredicate):
    """Select atoms by their type (integer or string)."""

    def __init__(self, atom_type: Union[int, str], field: str = "type"):
        self.atom_type = atom_type
        self.field = field

    def mask(self, block: "Block") -> np.ndarray:
        if self.field not in block:
            return np.zeros(block.nrows, dtype=bool)
        return (block[self.field] == self.atom_type)

# Usage
carbons = AtomTypeSelection(atom_type=1, field="type")
proteins = AtomTypeSelection(atom_type="protein", field="resname")
```

### AtomIndexSelection

Selects atoms by their indices:

```python
class AtomIndexSelection(MaskPredicate):
    """Select atoms by their indices."""

    def __init__(self, indices: Union[List[int], np.ndarray], id_field: str = "id"):
        # Convert to numpy array and validate
        if isinstance(indices, list):
            self.indices = np.array(indices, dtype=int)
        elif isinstance(indices, np.ndarray):
            self.indices = indices.astype(int)
        else:
            raise TypeError("indices must be a list[int] or np.ndarray")

        self.id_field = id_field

    def mask(self, block: "Block") -> np.ndarray:
        if self.id_field not in block:
            return np.zeros(block.nrows, dtype=bool)
        return np.isin(block[self.id_field], self.indices)

# Usage
specific_atoms = AtomIndexSelection(indices=[0, 5, 10])
chain_a = AtomIndexSelection(indices=range(100, 200), id_field="chain_id")
```

### ElementSelection

Selects atoms by their element symbol:

```python
class ElementSelection(MaskPredicate):
    """Select atoms by their element symbol."""

    def __init__(self, element: str, field: str = "element"):
        if not isinstance(element, str):
            raise TypeError("element must be a string")
        self.element = element
        self.field = field

    def mask(self, block: "Block") -> np.ndarray:
        if self.field not in block:
            return np.zeros(block.nrows, dtype=bool)
        return (block[self.field] == self.element)

# Usage
carbons = ElementSelection("C")
hydrogens = ElementSelection("H")
oxygens = ElementSelection("O")
```

### CoordinateRangeSelection

Selects atoms within a coordinate range:

```python
class CoordinateRangeSelection(MaskPredicate):
    """Select atoms within a coordinate range."""

    def __init__(self,
                 axis: str,
                 min_value: Optional[float] = None,
                 max_value: Optional[float] = None):
        if axis not in ["x", "y", "z"]:
            raise ValueError("axis must be 'x', 'y', or 'z'")

        if min_value is not None and max_value is not None and min_value > max_value:
            raise ValueError("min_value cannot be greater than max_value")

        self.axis = axis
        self.min_value = min_value
        self.max_value = max_value

    def mask(self, block: "Block") -> np.ndarray:
        if self.axis not in block:
            return np.zeros(block.nrows, dtype=bool)

        values = block[self.axis]
        mask = np.ones(block.nrows, dtype=bool)

        if self.min_value is not None:
            mask &= (values >= self.min_value)

        if self.max_value is not None:
            mask &= (values <= self.max_value)

        return mask

# Usage
top_layer = CoordinateRangeSelection("z", min_value=5.0)
in_box = CoordinateRangeSelection("x", min_value=0.0, max_value=10.0)
```

### DistanceSelection

Selects atoms within a distance from a reference point:

```python
class DistanceSelection(MaskPredicate):
    """Select atoms within a distance from a reference point."""

    def __init__(self,
                 center: Union[List[float], np.ndarray],
                 max_distance: float,
                 min_distance: Optional[float] = None):
        # Convert center to numpy array and validate
        if isinstance(center, list):
            self.center = np.array(center, dtype=float)
        elif isinstance(center, np.ndarray):
            self.center = center.astype(float)
        else:
            raise TypeError("center must be a list[float] or np.ndarray")

        if len(self.center) != 3:
            raise ValueError("center must have exactly 3 coordinates")

        if max_distance < 0:
            raise ValueError("max_distance must be non-negative")

        if min_distance is not None:
            if min_distance < 0:
                raise ValueError("min_distance must be non-negative")
            if min_distance > max_distance:
                raise ValueError("min_distance cannot be greater than max_distance")

        self.max_distance = max_distance
        self.min_distance = min_distance

    def mask(self, block: "Block") -> np.ndarray:
        required_fields = ["x", "y", "z"]
        if not all(field in block for field in required_fields):
            return np.zeros(block.nrows, dtype=bool)

        positions = np.column_stack([block["x"], block["y"], block["z"]])
        distances = np.linalg.norm(positions - self.center, axis=1)

        mask = distances <= self.max_distance

        if self.min_distance is not None:
            mask &= (distances >= self.min_distance)

        return mask

# Usage
near_origin = DistanceSelection([0.0, 0.0, 0.0], max_distance=5.0)
shell = DistanceSelection([0.0, 0.0, 0.0], max_distance=10.0, min_distance=5.0)
```

## Creating Custom Selection Predicates

### When to Create Custom Predicates

Create custom selection predicates when you need to:

1. **Select by custom criteria** not covered by built-in predicates
2. **Implement domain-specific selection logic** (chemical properties, spatial regions)
3. **Optimize selection performance** for specific use cases
4. **Combine multiple selection criteria** in complex ways
5. **Create reusable selection patterns** for your specific domain

### Basic Predicate Structure

Every custom predicate must inherit from `MaskPredicate` and implement the `mask` method:

```python
class CustomPredicate(MaskPredicate):
    """Template for custom selection predicates."""

    def __init__(self, param1, param2, **kwargs):
        # Store parameters
        self.param1 = param1
        self.param2 = param2
        # ... other initialization

    def mask(self, block: "Block") -> np.ndarray:
        """
        Generate boolean mask for the block.

        Args:
            block: Block of data to filter

        Returns:
            Boolean array with same length as block
        """
        # Your selection logic here
        # Must return np.ndarray with dtype=bool and length=block.nrows
        pass
```

### Predicate Design Patterns

#### Pattern 1: Property-Based Selection

Select atoms based on custom properties:

```python
class PropertyRangeSelection(MaskPredicate):
    """Select atoms with properties in a specific range."""

    def __init__(self, property_name: str, min_value: Optional[float] = None, max_value: Optional[float] = None):
        if not isinstance(property_name, str):
            raise TypeError("property_name must be a string")

        if min_value is not None and max_value is not None and min_value > max_value:
            raise ValueError("min_value cannot be greater than max_value")

        self.property_name = property_name
        self.min_value = min_value
        self.max_value = max_value

    def mask(self, block: "Block") -> np.ndarray:
        if self.property_name not in block:
            return np.zeros(block.nrows, dtype=bool)

        values = block[self.property_name]
        mask = np.ones(block.nrows, dtype=bool)

        if self.min_value is not None:
            mask &= (values >= self.min_value)

        if self.max_value is not None:
            mask &= (values <= self.max_value)

        return mask

# Usage
heavy_atoms = PropertyRangeSelection('mass', min_value=10.0)
charged_atoms = PropertyRangeSelection('charge', min_value=-1.0, max_value=1.0)
high_temp = PropertyRangeSelection('temperature', min_value=300.0)
```

**Key Design Principles**:
- **Validate inputs**: Check parameter types and logical consistency
- **Handle missing fields**: Return zero mask when required fields don't exist
- **Efficient operations**: Use vectorized operations when possible

#### Pattern 2: Spatial Selection

Select atoms based on spatial criteria:

```python
class WithinBoxSelection(MaskPredicate):
    """Select atoms within a 3D box."""

    def __init__(self,
                 min_corner: Union[List[float], np.ndarray],
                 max_corner: Union[List[float], np.ndarray]):
        # Convert to numpy arrays and validate
        if isinstance(min_corner, list):
            self.min_corner = np.array(min_corner, dtype=float)
        elif isinstance(min_corner, np.ndarray):
            self.min_corner = min_corner.astype(float)
        else:
            raise TypeError("min_corner must be a list[float] or np.ndarray")

        if isinstance(max_corner, list):
            self.max_corner = np.array(max_corner, dtype=float)
        elif isinstance(max_corner, np.ndarray):
            self.max_corner = max_corner.astype(float)
        else:
            raise TypeError("max_corner must be a list[float] or np.ndarray")

        if len(self.min_corner) != 3 or len(self.max_corner) != 3:
            raise ValueError("corners must have exactly 3 coordinates")

        if np.any(self.min_corner >= self.max_corner):
            raise ValueError("min_corner must be less than max_corner in all dimensions")

    def mask(self, block: "Block") -> np.ndarray:
        required_fields = ["x", "y", "z"]
        if not all(field in block for field in required_fields):
            return np.zeros(block.nrows, dtype=bool)

        positions = np.column_stack([block["x"], block["y"], block["z"]])

        # Check if positions are within box bounds
        within_x = (positions[:, 0] >= self.min_corner[0]) & (positions[:, 0] <= self.max_corner[0])
        within_y = (positions[:, 1] >= self.min_corner[1]) & (positions[:, 1] <= self.max_corner[1])
        within_z = (positions[:, 2] >= self.min_corner[2]) & (positions[:, 2] <= self.max_corner[2])

        return within_x & within_y & within_z

# Usage
unit_cell = WithinBoxSelection([0.0, 0.0, 0.0], [10.0, 10.0, 10.0])
active_region = WithinBoxSelection([-5.0, -5.0, -5.0], [5.0, 5.0, 5.0])
```

**Key Design Principles**:
- **Type conversion**: Accept both lists and numpy arrays for flexibility
- **Boundary validation**: Ensure logical consistency of parameters
- **Vectorized operations**: Use numpy operations for efficiency

#### Pattern 3: Chemical Environment Selection

Select atoms based on their chemical environment:

```python
class ChemicalEnvironmentSelection(MaskPredicate):
    """Select atoms based on their chemical environment."""

    def __init__(self,
                 central_atom_type: str,
                 neighbor_types: List[str],
                 max_distance: float,
                 min_neighbors: int = 1,
                 neighbor_field: str = "element"):
        if not isinstance(central_atom_type, str):
            raise TypeError("central_atom_type must be a string")

        if not isinstance(neighbor_types, list) or not all(isinstance(t, str) for t in neighbor_types):
            raise TypeError("neighbor_types must be a list of strings")

        if max_distance <= 0:
            raise ValueError("max_distance must be positive")

        if min_neighbors < 0:
            raise ValueError("min_neighbors must be non-negative")

        self.central_atom_type = central_atom_type
        self.neighbor_types = neighbor_types
        self.max_distance = max_distance
        self.min_neighbors = min_neighbors
        self.neighbor_field = neighbor_field

    def mask(self, block: "Block") -> np.ndarray:
        required_fields = ["x", "y", "z", self.neighbor_field]
        if not all(field in block for field in required_fields):
            return np.zeros(block.nrows, dtype=bool)

        # Find central atoms
        central_mask = block[self.neighbor_field] == self.central_atom_type
        central_indices = np.where(central_mask)[0]

        if len(central_indices) == 0:
            return np.zeros(block.nrows, dtype=bool)

        positions = np.column_stack([block["x"], block["y"], block["z"]])
        elements = block[self.neighbor_field]

        result_mask = np.zeros(block.nrows, dtype=bool)

        for central_idx in central_indices:
            central_pos = positions[central_idx]

            # Find neighbors within distance
            distances = np.linalg.norm(positions - central_pos, axis=1)
            neighbor_mask = distances <= self.max_distance

            # Check if neighbors are of required types
            neighbor_elements = elements[neighbor_mask]
            valid_neighbors = sum(1 for elem in neighbor_elements if elem in self.neighbor_types)

            if valid_neighbors >= self.min_neighbors:
                result_mask[central_idx] = True

        return result_mask

# Usage
# Select carbon atoms with at least 2 hydrogen neighbors within 1.5 Å
carbon_with_hydrogens = ChemicalEnvironmentSelection(
    central_atom_type='C',
    neighbor_types=['H'],
    max_distance=1.5,
    min_neighbors=2
)

# Select oxygen atoms with exactly 1 hydrogen neighbor within 1.0 Å
oxygen_with_hydrogen = ChemicalEnvironmentSelection(
    central_atom_type='O',
    neighbor_types=['H'],
    max_distance=1.0,
    min_neighbors=1
)
```

**Key Design Principles**:
- **Complex logic**: Handle multi-step selection processes
- **Performance consideration**: Use vectorized operations where possible
- **Flexible parameters**: Allow customization of selection criteria

#### Pattern 4: Composite Selection

Create complex selections by combining simple ones:

```python
class CompositeSelection(MaskPredicate):
    """Select atoms based on multiple criteria combined with custom logic."""

    def __init__(self,
                 predicates: List[MaskPredicate],
                 combination_logic: str = "all"):
        if not isinstance(predicates, list) or not all(isinstance(p, MaskPredicate) for p in predicates):
            raise TypeError("predicates must be a list of MaskPredicate objects")

        if combination_logic not in ["all", "any", "majority"]:
            raise ValueError("combination_logic must be 'all', 'any', or 'majority'")

        self.predicates = predicates
        self.combination_logic = combination_logic

    def mask(self, block: "Block") -> np.ndarray:
        if not self.predicates:
            return np.zeros(block.nrows, dtype=bool)

        # Generate masks for all predicates
        masks = [pred.mask(block) for pred in self.predicates]

        if self.combination_logic == "all":
            # All predicates must be True
            result = masks[0]
            for mask in masks[1:]:
                result &= mask
        elif self.combination_logic == "any":
            # Any predicate can be True
            result = masks[0]
            for mask in masks[1:]:
                result |= mask
        elif self.combination_logic == "majority":
            # Majority of predicates must be True
            result = np.zeros(block.nrows, dtype=bool)
            for mask in masks:
                result = result.astype(int) + mask.astype(int)
            threshold = len(masks) // 2 + 1
            result = result >= threshold

        return result

# Usage
# Select atoms that meet ALL criteria
strict_selection = CompositeSelection([
    ElementSelection("C"),
    PropertyRangeSelection("mass", min_value=10.0),
    DistanceSelection([0.0, 0.0, 0.0], max_distance=5.0)
], combination_logic="all")

# Select atoms that meet ANY criteria
loose_selection = CompositeSelection([
    ElementSelection("C"),
    ElementSelection("H"),
    ElementSelection("O")
], combination_logic="any")
```

**Key Design Principles**:
- **Flexible combination**: Support different logical combination strategies
- **Efficient evaluation**: Generate all masks once and combine efficiently
- **Clear interface**: Provide intuitive parameter names for combination logic

## Advanced Predicate Techniques

### Predicate Caching

For expensive predicates, implement caching to avoid recomputation:

```python
class CachedPredicate(MaskPredicate):
    """Predicate with caching for expensive computations."""

    def __init__(self, predicate: MaskPredicate, cache_key: Optional[str] = None):
        self.predicate = predicate
        self.cache_key = cache_key
        self._cached_mask = None
        self._cached_block_id = None

    def mask(self, block: "Block") -> np.ndarray:
        # Generate a unique identifier for the block
        block_id = id(block)

        # Check if we can use cached result
        if (self._cached_mask is not None and
            self._cached_block_id == block_id):
            return self._cached_mask

        # Compute the mask
        mask = self.predicate.mask(block)

        # Cache the result
        self._cached_mask = mask
        self._cached_block_id = block_id

        return mask

    def clear_cache(self):
        """Clear the cached mask."""
        self._cached_mask = None
        self._cached_block_id = None

# Usage
expensive_predicate = ChemicalEnvironmentSelection(
    central_atom_type='C',
    neighbor_types=['H', 'O', 'N'],
    max_distance=3.0,
    min_neighbors=3
)
cached_predicate = CachedPredicate(expensive_predicate)

# First call computes the mask
mask1 = cached_predicate.mask(block)

# Second call uses cached result
mask2 = cached_predicate.mask(block)  # Much faster!

# Clear cache when needed
cached_predicate.clear_cache()
```

### Conditional Predicates

Create predicates that only apply under certain conditions:

```python
class ConditionalPredicate(MaskPredicate):
    """Predicate that only applies when a condition is met."""

    def __init__(self,
                 predicate: MaskPredicate,
                 condition: Callable[["Block"], bool],
                 fallback_mask: Optional[np.ndarray] = None):
        self.predicate = predicate
        self.condition = condition
        self.fallback_mask = fallback_mask

    def mask(self, block: "Block") -> np.ndarray:
        if self.condition(block):
            return self.predicate.mask(block)
        elif self.fallback_mask is not None:
            # Return fallback mask (must match block length)
            if len(self.fallback_mask) == block.nrows:
                return self.fallback_mask
            else:
                return np.zeros(block.nrows, dtype=bool)
        else:
            return np.zeros(block.nrows, dtype=bool)

# Usage
def has_force_field(block):
    """Check if block has force field information."""
    return 'type' in block and 'charge' in block

# Only apply force field selection when force field data is available
force_field_selection = ConditionalPredicate(
    AtomTypeSelection(1),
    has_force_field,
    fallback_mask=None  # No selection when no force field
)
```

### Performance-Optimized Predicates

For high-performance applications, optimize predicate evaluation:

```python
class OptimizedPredicate(MaskPredicate):
    """Optimized predicate for high-performance applications."""

    def __init__(self, property_name: str, values: Union[List, np.ndarray]):
        self.property_name = property_name

        # Pre-compute lookup set for O(1) membership testing
        if isinstance(values, list):
            self.values = set(values)
        elif isinstance(values, np.ndarray):
            self.values = set(values.tolist())
        else:
            raise TypeError("values must be a list or numpy array")

    def mask(self, block: "Block") -> np.ndarray:
        if self.property_name not in block:
            return np.zeros(block.nrows, dtype=bool)

        # Use vectorized operations and pre-computed set
        return np.vectorize(lambda x: x in self.values)(block[self.property_name])

# Usage
# Pre-compute common element sets
common_elements = OptimizedPredicate("element", ["C", "H", "O", "N"])
rare_elements = OptimizedPredicate("element", ["Fe", "Cu", "Zn", "Mn"])

# These predicates will be very fast for repeated use
```

## Best Practices for Predicate Development

### 1. **Always Handle Missing Fields**

```python
def mask(self, block: "Block") -> np.ndarray:
    # Check if required fields exist
    if self.field_name not in block:
        return np.zeros(block.nrows, dtype=bool)

    # Proceed with selection logic
    return block[self.field_name] == self.value
```

### 2. **Validate Inputs Thoroughly**

```python
def __init__(self, param1, param2):
    # Type checking
    if not isinstance(param1, (int, float)):
        raise TypeError("param1 must be numeric")

    # Logical validation
    if param2 <= 0:
        raise ValueError("param2 must be positive")

    # Store validated parameters
    self.param1 = param1
    self.param2 = param2
```

### 3. **Use Vectorized Operations**

```python
# ✅ Good: Vectorized operations
def mask(self, block: "Block") -> np.ndarray:
    values = block[self.field]
    return (values >= self.min_value) & (values <= self.max_value)

# ❌ Avoid: Loops when possible
def mask(self, block: "Block") -> np.ndarray:
    result = np.zeros(block.nrows, dtype=bool)
    for i in range(block.nrows):
        value = block[self.field][i]
        result[i] = self.min_value <= value <= self.max_value
    return result
```

### 4. **Provide Clear Error Messages**

```python
def mask(self, block: "Block") -> np.ndarray:
    required_fields = ["x", "y", "z"]
    missing_fields = [field for field in required_fields if field not in block]

    if missing_fields:
        raise ValueError(f"Missing required fields: {missing_fields}")

    # Proceed with selection logic
    return self._compute_mask(block)
```

### 5. **Document Your Predicates**

```python
class CustomPredicate(MaskPredicate):
    """
    Select atoms based on custom criteria.

    This predicate selects atoms that meet specific requirements
    based on their properties and environment.

    Examples:
        >>> predicate = CustomPredicate(param1=1.0, param2=2.0)
        >>> mask = predicate.mask(block)
        >>> selected = block[mask]
    """

    def __init__(self, param1: float, param2: float):
        """
        Initialize the predicate.

        Args:
            param1: First parameter for selection
            param2: Second parameter for selection

        Raises:
            ValueError: If parameters are invalid
        """
        # Implementation...
```

## Summary

MolPy's selection predicate system provides a powerful and flexible way to filter molecular data:

- **Functional composition**: Predicates are functions that produce boolean masks
- **Logical combination**: Use `&`, `|`, and `~` operators to combine predicates
- **Efficient evaluation**: Vectorized operations for high performance
- **Extensible design**: Easy to create custom predicates for any purpose
- **Type safety**: Comprehensive input validation and error handling

The key insight is that predicates create **reusable, composable filters** that can be combined in complex ways to build sophisticated selection queries.

### Next Steps

Continue exploring MolPy's extension capabilities by learning about custom I/O systems, analysis pipelines, and integration with external tools.
