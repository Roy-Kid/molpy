# Custom Wrappers: Extending MolPy with Composition

This chapter explains how to create custom wrappers that extend MolPy's functionality using the composition pattern. Wrappers are the foundation of MolPy's extensibility, providing a clean way to add new capabilities without modifying existing classes.

## Understanding the Wrapper Stack

### Core Delegation Mechanism

MolPy's `Wrapper` class implements a sophisticated delegation system that creates a **stack of functionality**. When you access an attribute or method on a wrapped object, the system automatically searches through the wrapper stack until it finds what you're looking for.

```python
class Wrapper(Generic[T]):
    """Base class for all wrappers."""

    def __init__(self, wrapped: T):
        self._wrapped = wrapped

    def __getattr__(self, name: str) -> Any:
        """Delegate attribute access to the wrapped entity."""
        return getattr(self._wrapped, name)

    def __setattr__(self, name: str, value: Any) -> None:
        """Set attributes on wrapper or delegate to wrapped entity."""
        if name.startswith('_') or name in self.__dict__ or hasattr(type(self), name):
            super().__setattr__(name, value)
        else:
            setattr(self._wrapped, name, value)
```

**Key Insight**: The `__getattr__` method is only called when Python can't find the attribute in the current object. This creates a **fallback chain** that automatically searches through the wrapper stack.

### How the Stack Works

```python
# Create a molecule
molecule = Molecule(...)

# Wrap it with multiple capabilities
spatial_molecule = Spatial(molecule)           # Layer 1: Spatial operations
visual_molecule = Visual(spatial_molecule)     # Layer 2: Visual properties
hierarchical_molecule = HierarchyWrapper(visual_molecule)  # Layer 3: Hierarchy

# When you call a method, Python searches the stack:
hierarchical_molecule.move([1.0, 0.0, 0.0])
# 1. Does hierarchical_molecule have a 'move' method? No
# 2. Call __getattr__('move') → delegate to visual_molecule
# 3. Does visual_molecule have a 'move' method? No
# 4. Call __getattr__('move') → delegate to spatial_molecule
# 5. Does spatial_molecule have a 'move' method? Yes!
# 6. Execute the method
```

### Wrapper Chain Management

MolPy provides sophisticated tools to inspect and manage the wrapper stack:

```python
def wrapper_chain(self) -> list[Any]:
    """Return the list of objects from outermost wrapper to final unwrapped entity."""
    chain = []
    obj = self
    while isinstance(obj, Wrapper):
        chain.append(obj)
        obj = obj.unwrap()
    chain.append(obj)  # the final, non-wrapper object
    return chain

def wrapper_depth(self) -> int:
    """How many Wrapper layers are there?"""
    return len(self.wrapper_chain()) - 1

def wrapper_types(self) -> list[str]:
    """Return the list of class names for each layer in the wrapper chain."""
    return [type(o).__name__ for o in self.wrapper_chain()]
```

**Usage Example**:
```python
# Inspect the wrapper stack
print(f"Wrapper depth: {hierarchical_molecule.wrapper_depth()}")
print(f"Wrapper types: {hierarchical_molecule.wrapper_types()}")
print(f"Wrapper chain: {hierarchical_molecule.wrapper_chain()}")

# Output:
# Wrapper depth: 3
# Wrapper types: ['HierarchyWrapper', 'VisualWrapper', 'SpatialWrapper', 'Molecule']
# Wrapper chain: [<HierarchyWrapper>, <VisualWrapper>, <SpatialWrapper>, <Molecule>]
```

## Built-in Wrapper Classes

### Spatial Wrapper

The `Spatial` wrapper adds geometric operations to any entity with coordinates:

```python
class Spatial(Wrapper):
    """Wrapper providing spatial operations for entities."""

    @property
    def xyz(self) -> np.ndarray:
        """Get xyz coordinates as numpy array."""
        atoms = self.atoms
        return np.array([atom["xyz"] for atom in atoms], dtype=float)

    def move(self, vector: ArrayLike) -> None:
        """Translate the entity by a vector."""
        self.xyz = self.xyz + np.array(vector)

    def rotate(self, axis: ArrayLike, angle: float, origin: Optional[ArrayLike] = None) -> None:
        """Rotate around an axis by an angle."""
        if origin is None:
            origin = np.zeros(3)
        else:
            origin = np.array(origin)
        axis = np.array(axis, dtype=float)
        translated = self.xyz - origin
        rotated = rotate_by_rodrigues(translated, axis, angle)
        self.xyz = origin + rotated
```

**Key Features**:
- **Automatic coordinate extraction**: Assumes the wrapped entity has an `atoms` attribute
- **Vectorized operations**: All operations work on the entire coordinate set
- **Fallback delegation**: If `self.atoms` doesn't exist, it will delegate to the wrapped entity

### Hierarchy Wrapper

The `HierarchyWrapper` adds parent-child relationship management:

```python
class HierarchyWrapper(Wrapper):
    """Wrapper providing hierarchical structure functionality."""

    def __init__(self, wrapped):
        super().__init__(wrapped)
        if not hasattr(self._wrapped, '_parent'):
            self._wrapped._parent = None
        if not hasattr(self._wrapped, '_children'):
            self._wrapped._children = []

    def add_child(self, child) -> None:
        """Add a child entity."""
        child_entity = child.unwrap() if hasattr(child, 'unwrap') else child
        if child_entity not in self._wrapped._children:
            self._wrapped._children.append(child_entity)
            child_entity._parent = self._wrapped

    def get_descendants(self) -> list:
        """Get all descendants recursively."""
        descendants = []
        for child in self._wrapped._children:
            descendants.append(child)
            if hasattr(child, '_children'):
                child_wrapper = HierarchyWrapper(child)
                descendants.extend(child_wrapper.get_descendants())
        return descendants
```

**Key Features**:
- **Safe attribute access**: Checks if attributes exist before using them
- **Recursive operations**: Can traverse the entire hierarchy tree
- **Wrapper-aware**: Handles both wrapped and unwrapped children

### Visual Wrapper

The `VisualWrapper` adds visualization properties:

```python
class VisualWrapper(Wrapper):
    """Wrapper providing visualization functionality for entities."""

    def __init__(self, wrapped, color: Optional[str] = None, size: Optional[float] = None, **visual_props):
        super().__init__(wrapped)
        if not hasattr(self._wrapped, '_visual_props'):
            self._wrapped._visual_props = {}

        if color is not None:
            self._wrapped._visual_props['color'] = color
        if size is not None:
            self._wrapped._visual_props['size'] = size

        self._wrapped._visual_props.update(visual_props)

    @property
    def color(self) -> Optional[str]:
        """Get the color property."""
        return self._wrapped._visual_props.get('color')

    @property
    def visual_props(self) -> dict:
        """Get all visual properties."""
        return getattr(self._wrapped, '_visual_props', {}).copy()
```

**Key Features**:
- **Property management**: Provides getter/setter for visual properties
- **Safe defaults**: Uses `.get()` and `.copy()` to avoid errors
- **Extensible**: Accepts arbitrary visual properties via `**kwargs`

## Creating Custom Wrappers

### When to Create Wrappers

Create custom wrappers when you need to:

1. **Add new capabilities** to existing entities without inheritance
2. **Implement cross-cutting concerns** (logging, caching, validation)
3. **Provide domain-specific functionality** (analysis, visualization)
4. **Maintain backward compatibility** while adding features
5. **Create composable functionality** that can be mixed and matched

### Wrapper Design Patterns

#### Pattern 1: Capability Wrapper

Add specific capabilities to entities:

```python
class AnalysisWrapper(Wrapper):
    """Wrapper that adds analysis capabilities."""

    def __init__(self, wrapped):
        super().__init__(wrapped)
        # Validate that wrapped entity supports analysis
        if not hasattr(wrapped, 'atoms'):
            raise ValueError("Wrapped entity must have 'atoms' attribute")

    def calculate_center_of_mass(self) -> np.ndarray:
        """Calculate center of mass."""
        atoms = self.atoms
        if 'mass' in atoms:
            positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])
            masses = atoms['mass']
            return np.average(positions, weights=masses, axis=0)
        else:
            # Fallback to geometric center
            positions = np.column_stack([atoms['x'], atoms['y'], atoms['z']])
            return np.mean(positions, axis=0)

    def calculate_radius_of_gyration(self) -> float:
        """Calculate radius of gyration."""
        com = self.calculate_center_of_mass()
        positions = np.column_stack([self.atoms['x'], self.atoms['y'], self.atoms['z']])
        distances = np.linalg.norm(positions - com, axis=1)
        return np.sqrt(np.mean(distances**2))

# Usage
molecule = Molecule(...)
analysis_molecule = AnalysisWrapper(molecule)
com = analysis_molecule.calculate_center_of_mass()
rg = analysis_molecule.calculate_radius_of_gyration()
```

**Key Design Principles**:
- **Validate requirements**: Check that the wrapped entity has necessary attributes
- **Provide fallbacks**: Offer alternative implementations when possible
- **Delegate naturally**: Let the wrapper stack handle missing methods

#### Pattern 2: Validation Wrapper

Add validation and error checking:

```python
class ValidationWrapper(Wrapper):
    """Wrapper that adds validation capabilities."""

    def __init__(self, wrapped, validation_rules: dict = None):
        super().__init__(wrapped)
        self.validation_rules = validation_rules or {}
        self._validate_entity()

    def _validate_entity(self):
        """Validate the wrapped entity."""
        errors = []

        # Check required attributes
        if not hasattr(self._wrapped, 'atoms'):
            errors.append("Entity must have 'atoms' attribute")

        # Check data consistency
        if hasattr(self._wrapped, 'atoms'):
            atoms = self._wrapped.atoms
            if 'x' in atoms and 'y' in atoms and 'z' in atoms:
                if len(atoms['x']) != len(atoms['y']) or len(atoms['y']) != len(atoms['z']):
                    errors.append("Coordinate arrays must have same length")

        if errors:
            raise ValueError(f"Validation failed: {'; '.join(errors)}")

    def safe_operation(self, operation_name: str, *args, **kwargs):
        """Perform operation with validation."""
        try:
            if hasattr(self._wrapped, operation_name):
                return getattr(self._wrapped, operation_name)(*args, **kwargs)
            else:
                raise AttributeError(f"Operation '{operation_name}' not supported")
        except Exception as e:
            raise RuntimeError(f"Operation '{operation_name}' failed: {e}")

# Usage
molecule = Molecule(...)
validated_molecule = ValidationWrapper(molecule)
result = validated_molecule.safe_operation('calculate_energy')
```

**Key Design Principles**:
- **Comprehensive validation**: Check all requirements upfront
- **Clear error messages**: Provide specific information about what failed
- **Safe operations**: Wrap operations in try-catch blocks

#### Pattern 3: Composite Wrapper

Combine multiple capabilities:

```python
class CompositeWrapper(Wrapper):
    """Wrapper that combines multiple capabilities."""

    def __init__(self, wrapped, *capabilities):
        super().__init__(wrapped)
        self.capabilities = {}

        for capability in capabilities:
            if isinstance(capability, Wrapper):
                # Extract capability name and methods
                capability_name = capability.__class__.__name__.lower().replace('wrapper', '')
                self.capabilities[capability_name] = capability

                # Add capability methods to this wrapper
                for attr_name in dir(capability):
                    if not attr_name.startswith('_') and callable(getattr(capability, attr_name)):
                        setattr(self, f"{capability_name}_{attr_name}", getattr(capability, attr_name))

    def get_capability(self, name: str):
        """Get a specific capability wrapper."""
        return self.capabilities.get(name)

    def has_capability(self, name: str) -> bool:
        """Check if a capability is available."""
        return name in self.capabilities

# Usage
molecule = Molecule(...)
enhanced_molecule = CompositeWrapper(
    molecule,
    Spatial(molecule),
    AnalysisWrapper(molecule),
    ValidationWrapper(molecule)
)

# Use combined capabilities
enhanced_molecule.move([1.0, 0.0, 0.0])  # From Spatial
com = enhanced_molecule.calculate_center_of_mass()  # From AnalysisWrapper
enhanced_molecule.safe_operation('calculate_energy')  # From ValidationWrapper
```

**Key Design Principles**:
- **Method composition**: Combine methods from multiple wrappers
- **Namespace management**: Use prefixes to avoid method name conflicts
- **Capability discovery**: Provide methods to inspect available capabilities

## Advanced Wrapper Techniques

### Method Overriding

You can override methods from the wrapped entity to provide enhanced functionality:

```python
class EnhancedSpatialWrapper(Spatial):
    """Enhanced spatial wrapper with additional features."""

    def move(self, vector: ArrayLike) -> None:
        """Enhanced move with logging and validation."""
        # Pre-move validation
        if not self._validate_coordinates():
            raise ValueError("Invalid coordinates before move")

        # Log the move operation
        print(f"Moving entity by vector: {vector}")

        # Call the original method
        super().move(vector)

        # Post-move validation
        if not self._validate_coordinates():
            raise ValueError("Invalid coordinates after move")

    def _validate_coordinates(self) -> bool:
        """Validate that coordinates are reasonable."""
        coords = self.xyz
        return np.all(np.isfinite(coords)) and np.all(np.abs(coords) < 1e6)
```

**Key Insight**: Use `super().method()` to call the original implementation while adding your own logic.

### Conditional Wrapping

Create wrappers that only apply in certain conditions:

```python
class ConditionalWrapper(Wrapper):
    """Wrapper that only applies functionality under certain conditions."""

    def __init__(self, wrapped, condition: Callable[[Any], bool]):
        super().__init__(wrapped)
        self.condition = condition

    def __getattr__(self, name: str) -> Any:
        """Conditionally delegate attribute access."""
        if self.condition(self._wrapped):
            # Apply wrapper logic
            return getattr(self._wrapped, name)
        else:
            # Skip wrapper, go directly to wrapped entity
            return getattr(self._wrapped, name)

    def conditional_method(self, *args, **kwargs):
        """Method that only works when condition is met."""
        if not self.condition(self._wrapped):
            raise RuntimeError("Condition not met for this operation")
        # Perform the operation
        return self._wrapped.some_method(*args, **kwargs)

# Usage
def is_large_system(entity):
    return hasattr(entity, 'atoms') and len(entity.atoms) > 1000

large_system_wrapper = ConditionalWrapper(molecule, is_large_system)
```

### Wrapper Factories

Create functions that return appropriate wrappers:

```python
def create_analysis_wrapper(entity, analysis_type: str = "basic"):
    """Factory function for creating analysis wrappers."""
    if analysis_type == "basic":
        return AnalysisWrapper(entity)
    elif analysis_type == "advanced":
        return AdvancedAnalysisWrapper(entity)
    elif analysis_type == "custom":
        return CustomAnalysisWrapper(entity)
    else:
        raise ValueError(f"Unknown analysis type: {analysis_type}")

def wrap_for_visualization(entity, style: str = "default"):
    """Factory function for creating visualization wrappers."""
    base_wrapper = VisualWrapper(entity)

    if style == "scientific":
        base_wrapper.set_visual_prop("color_scheme", "element")
        base_wrapper.set_visual_prop("atom_size", "vdw")
    elif style == "cartoon":
        base_wrapper.set_visual_prop("color_scheme", "chain")
        base_wrapper.set_visual_prop("representation", "cartoon")

    return base_wrapper

# Usage
molecule = Molecule(...)
analysis_mol = create_analysis_wrapper(molecule, "advanced")
visual_mol = wrap_for_visualization(molecule, "scientific")
```

## Best Practices for Wrapper Development

### 1. **Follow the Delegation Pattern**

```python
class GoodWrapper(Wrapper):
    def __init__(self, wrapped):
        super().__init__(wrapped)

    def __getattr__(self, name):
        # Always delegate to wrapped object
        return getattr(self._wrapped, name)

# ❌ Avoid: Breaking the delegation chain
class BadWrapper(Wrapper):
    def __init__(self, wrapped):
        super().__init__(wrapped)

    def some_method(self):
        # This breaks delegation for 'some_method'
        return "local implementation"
```

### 2. **Maintain Transparency**

```python
class TransparentWrapper(Wrapper):
    def __init__(self, wrapped):
        super().__init__(wrapped)

    def __getattr__(self, name):
        # Always delegate to wrapped object
        return getattr(self._wrapped, name)

    def __getitem__(self, key):
        # Always delegate item access
        return self._wrapped[key]

    def __len__(self):
        # Always delegate length
        return len(self._wrapped)
```

### 3. **Validate Inputs Carefully**

```python
class SafeWrapper(Wrapper):
    def __init__(self, wrapped):
        super().__init__(wrapped)
        self._validate_wrapped()

    def _validate_wrapped(self):
        """Ensure wrapped object supports required operations."""
        required_attrs = ['atoms', 'box']
        for attr in required_attrs:
            if not hasattr(self._wrapped, attr):
                raise ValueError(f"Wrapped object must have '{attr}' attribute")
```

### 4. **Use Type Hints**

```python
from typing import TypeVar, Generic

T = TypeVar('T')

class TypedWrapper(Wrapper[T], Generic[T]):
    def __init__(self, wrapped: T):
        super().__init__(wrapped)

    def unwrap(self) -> T:
        return self._wrapped
```

### 5. **Document Your Extensions**

```python
class DocumentedWrapper(Wrapper):
    """
    Wrapper that adds specific functionality.

    This wrapper provides [describe functionality] for entities
    that have [describe requirements].

    Examples:
        >>> entity = SomeEntity(...)
        >>> wrapped = DocumentedWrapper(entity)
        >>> result = wrapped.some_method()
    """

    def __init__(self, wrapped):
        """
        Initialize the wrapper.

        Args:
            wrapped: Entity to wrap. Must have [describe requirements].

        Raises:
            ValueError: If wrapped entity doesn't meet requirements.
        """
        super().__init__(wrapped)
        self._validate_requirements()
```

## Summary

MolPy's wrapper system provides a powerful and flexible way to extend functionality:

- **Stack-based delegation**: Wrappers automatically fall back through the stack
- **Composition over inheritance**: Add capabilities without modifying existing classes
- **Transparent operation**: Wrapped objects behave like unwrapped ones
- **Chain management**: Built-in tools to inspect and manage wrapper stacks
- **Extensible design**: Easy to create new wrappers for any purpose

The key insight is that wrappers create a **stack of functionality** where each layer can add new capabilities while maintaining access to all underlying functionality through automatic delegation.

### Next Steps

Continue exploring MolPy's extension capabilities by learning about custom selection predicates, I/O systems, and analysis pipelines.
