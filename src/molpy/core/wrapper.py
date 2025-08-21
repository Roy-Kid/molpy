"""
Wrapper classes for the molpy framework.

This module provides wrapper classes that can be composed to add functionality
to molecular entities. Replaces the old mixin-based approach with a more
flexible composition pattern.
"""

from typing import Any, Callable, Generic, TypeVar, overload

import numpy as np
from numpy.typing import ArrayLike

from ..op import rotate_by_rodrigues
from .protocol import Struct

T = TypeVar("T")


class Wrapper(Generic[T]):
    """
    Base class for all wrappers.

    Wrappers provide a composable way to add functionality to entities
    without using inheritance mixins.

    Supports both explicit wrapping and multiple inheritance auto-composition.
    """

    def __new__(cls, *args, **kwargs):
        """Record other Wrapper bases from MRO for auto-composition."""
        instance = super().__new__(cls)

        # Find other Wrapper bases from MRO, but only MolPy wrapper classes
        wrapper_bases = []
        for base in cls.__mro__:
            if (
                issubclass(base, Wrapper)
                and base != Wrapper
                and base != cls
                and base.__module__.startswith("molpy.")
            ):
                wrapper_bases.append(base)

        # Store wrapper bases for later initialization
        instance._wrapper_bases = wrapper_bases
        instance._wrappers_initialized = False

        return instance

    def __init__(
        self,
        wrapped: T | None = None,
        key_mapping: dict[str, str | list[str]] | None = None,
        **props,
    ):
        """
        Initialize wrapper with an entity to wrap.

        Args:
            wrapped: The entity to wrap, or None to create a base Struct(**props)
            key_mapping: Optional key mapping for delegation
            **props: Properties to pass to Struct if wrapped is None
        """
        # Create/attach wrapped entity
        if wrapped is None:
            self._wrapped = Struct(**props)
        else:
            self._wrapped = wrapped

        self._key_mapping = key_mapping if key_mapping else {}

        # Run __post_init__ hooks for auto-composition
        if not self._wrappers_initialized:
            self._run_post_init_hooks()
            self._wrappers_initialized = True

    def _run_post_init_hooks(self):
        """Run __post_init__ hooks for auto-composition."""
        # Call self.__post_init__ if present (Python standard naming)
        if hasattr(self, "__post_init__"):
            self.__post_init__()

        # Call __post_init__ on other wrapper bases
        if hasattr(self, "_wrapper_bases"):
            for base in self._wrapper_bases:
                if hasattr(base, "__post_init__"):
                    base.__post_init__(self)

    def __post_init__(self):
        """
        Post-initialization hook for wrapper setup.

        Override this method in subclasses to perform one-time initialization
        that depends on the wrapped entity being fully constructed.
        """
        pass

    def unwrap(self) -> T:
        """Get the wrapped entity."""
        return self._wrapped

    def __getattr__(self, name: str) -> Any:
        """Delegate attribute access to the wrapped entity."""
        # For dict-like objects (like Struct), try to get as key first
        try:
            # Check if it's a dict-like object by trying to access as key
            return self._wrapped[name]  # type: ignore[attr-defined]
        except (KeyError, TypeError, AttributeError):
            # Not a dict-like object or key doesn't exist, try attribute access
            pass

        # Try direct attribute access
        if hasattr(self._wrapped, name):
            return getattr(self._wrapped, name)

        # Fall back to attribute access (this will raise AttributeError if not found)
        return getattr(self._wrapped, name)

    def __setattr__(self, name: str, value: Any) -> None:
        """Set attributes on wrapper or delegate to wrapped entity."""
        if name.startswith("_") or name in self.__dict__ or hasattr(type(self), name):
            super().__setattr__(name, value)
        else:
            setattr(self._wrapped, name, value)

    def __getitem__(self, key, /):
        """Delegate item access to the wrapped entity."""
        if key in self._key_mapping:
            key = self._key_mapping[key]
        return self._wrapped[key]  # type: ignore[attr-defined]

    def __setitem__(self, key: Any, value: Any) -> None:
        """Set items on the wrapped entity."""
        if key in self._key_mapping:
            key = self._key_mapping[key]
        self._wrapped[key] = value  # type: ignore[attr-defined]

    def __contains__(self, key: Any) -> bool:
        """Check if the wrapped entity contains a key."""
        if key in self._key_mapping:
            key = self._key_mapping[key]
        return key in self._wrapped

    def get_stack(self) -> list[Any]:
        """
        Return the list of objects from the outermost wrapper
        all the way down to the final un­wrapped entity.
        """
        chain = []
        obj = self
        while isinstance(obj, Wrapper):
            chain.append(obj)
            obj = obj.unwrap()
        chain.append(obj)  # the final, non‐wrapper object
        return chain

    def wrapper_types(self) -> list[str]:
        """
        Return the list of class‐names for each layer in the wrapper chain,
        from outermost to innermost.
        """
        return [type(o).__name__ for o in self.get_stack()]

    def wrapper_depth(self) -> int:
        """
        How many Wrapper‐layers are there?  (Does *not* count the final unwrapped.)
        """
        # subtract 1 if you don’t want to count the raw object
        return len(self.get_stack()) - 1


class Spatial(Wrapper):
    """
    Wrapper class providing spatial operations for entities.

    Defines spatial functionality for entities that have xyz coordinates
    and can perform geometric transformations.
    """

    def __post_init__(self):
        """Post-initialization hook for spatial operations."""
        # Can be no-op for now, or add any spatial-specific initialization
        pass

    @property
    def xyz(self) -> np.ndarray:
        """Get the xyz coordinates as a numpy array."""
        atoms = self.atoms
        return np.array([atom["xyz"] for atom in atoms], dtype=float)

    @xyz.setter
    def xyz(self, value: ArrayLike) -> None:
        """Set the xyz coordinates from an array-like object."""
        new_xyz = np.array(value, dtype=float)
        if new_xyz.shape[0] != len(self.atoms):
            raise ValueError("New xyz coordinates must match the number of atoms.")
        for i, atom in enumerate(self.atoms):
            atom["xyz"] = new_xyz[i]

    @property
    def positions(self) -> np.ndarray:
        """Get the positions as a numpy array (alias for xyz)."""
        return self.xyz

    @positions.setter
    def positions(self, value: ArrayLike) -> None:
        """Set the positions from an array-like object (alias for xyz)."""
        self.xyz = value

    @property
    def symbols(self) -> list[str]:
        """Get the symbols of all atoms in the structure."""
        return [atom.get("symbol", atom.get("name", "")) for atom in self.atoms]

    def distance_to(self, other) -> float:
        """
        Calculate the distance to another spatial entity.

        Args:
            other: Another spatial entity (can be wrapped or unwrapped)

        Returns:
            The Euclidean distance between the two entities
        """
        other_xyz = other.xyz if hasattr(other, "xyz") else other
        return float(np.linalg.norm(self.xyz - other_xyz))

    def move(self, vector: ArrayLike) -> None:
        """
        Translate the entity by a given vector.

        Args:
            vector: Translation vector
        """
        self.xyz = self.xyz + np.array(vector)

    def move_to(
        self, target_position: ArrayLike, reference_atom_index: int = 0
    ) -> None:
        """
        Move the entity so that a reference atom is at the target position.
        All other atoms maintain their relative positions.

        Args:
            target_position: Target position for the reference atom
            reference_atom_index: Index of the atom to use as reference (default: 0)
        """
        target_position = np.array(target_position, dtype=float)

        # Get current position of reference atom
        current_reference_pos = self.atoms[reference_atom_index]["xyz"]

        # Calculate translation vector
        translation_vector = target_position - current_reference_pos

        # Apply translation
        self.move(translation_vector)

        # Verify the reference atom is now at target position
        assert np.allclose(
            self.atoms[reference_atom_index]["xyz"], target_position
        ), f"Reference atom not moved to target position. Expected {target_position}, got {self.atoms[reference_atom_index]['xyz']}"

    def scale(self, factor: float | ArrayLike, origin: ArrayLike | None = None) -> None:
        """
        Scale the entity by a given factor around an origin.

        Args:
            factor: Scaling factor (scalar or per-axis)
            origin: Origin point for scaling (defaults to current position)
        """
        if origin is None:
            origin = self.xyz
        else:
            origin = np.array(origin)

        factor = np.array(factor)
        self.xyz = origin + factor * (self.xyz - origin)

    def rotate(
        self, axis: ArrayLike, angle: float, origin: ArrayLike | None = None
    ) -> None:
        """
        Rotate the entity around an axis by a given angle.

        Args:
            axis: Rotation axis vector
            angle: Rotation angle in radians
            origin: Origin point for rotation (defaults to [0,0,0])
        """
        if origin is None:
            origin = np.zeros(3)
        else:
            origin = np.array(origin)
        axis = np.array(axis, dtype=float)
        translated = self.xyz - origin
        rotated = rotate_by_rodrigues(translated, axis, angle)
        if isinstance(rotated, np.ndarray) and rotated.shape == (1, 3):
            rotated = rotated[0]
        rotated = np.where(np.abs(rotated) < 1e-12, 0.0, rotated)
        self.xyz = origin + rotated

    def rotate_with_matrix(self, R: np.ndarray):
        """
        Rotate the entity using a rotation matrix.

        Args:
            R: 3x3 rotation matrix
        """
        self.xyz = np.dot(self.xyz, R)

    def __call__(self, **kwargs):
        """
        Create a new instance of the wrapped entity with optional modifications.

        This method enables wrapped entities to be used as factory functions,
        creating copies of themselves with potentially modified properties.

        Args:
            **kwargs: Properties to pass to the wrapped entity's constructor

        Returns:
            A new instance of the wrapped entity with copied data
        """
        # If the wrapped entity has a __call__ method, use it
        if hasattr(self._wrapped, "__call__"):
            new_wrapped = self._wrapped(**kwargs)
            # Return a new Spatial wrapping the new instance
            return Spatial(new_wrapped)

        # Otherwise, create a new instance using the wrapped entity's class
        wrapped_class = type(self._wrapped)

        # Try to create a new instance with kwargs
        new_instance = wrapped_class(**kwargs)

        # Copy relevant data from the current wrapped entity
        if hasattr(self._wrapped, "items"):
            for key, value in self._wrapped.items():
                if key not in kwargs:  # Don't override kwargs
                    if hasattr(value, "copy"):
                        new_instance[key] = value.copy()
                    else:
                        new_instance[key] = value

        return new_instance


class HierarchyWrapper(Wrapper):
    """
    Wrapper class providing hierarchical structure functionality.

    Allows entities to have parent-child relationships and provides
    methods for navigating and manipulating the hierarchy.
    """

    def __init__(self, wrapped):
        """
        Initialize hierarchy wrapper.

        Args:
            wrapped: The entity to wrap
        """
        super().__init__(wrapped)
        if not hasattr(self._wrapped, "_parent"):
            self._wrapped._parent = None
        if not hasattr(self._wrapped, "_children"):
            self._wrapped._children = []

    @property
    def parent(self):
        """Get the parent entity."""
        return self._wrapped._parent

    @property
    def children(self):
        """Get the list of child entities."""
        return (
            self._wrapped._children.copy()
        )  # Return a copy to prevent external modification

    @property
    def is_root(self) -> bool:
        """Check if this entity is a root (has no parent)."""
        return self._wrapped._parent is None

    @property
    def is_leaf(self) -> bool:
        """Check if this entity is a leaf (has no children)."""
        return len(self._wrapped._children) == 0

    @property
    def depth(self) -> int:
        """Get the depth of this entity in the hierarchy (root = 0)."""
        if self.is_root or self.parent is None:
            return 0
        # Handle case where parent might also be wrapped
        parent_depth = self.parent.depth if hasattr(self.parent, "depth") else 0
        return parent_depth + 1

    def add_child(self, child) -> None:
        """
        Add a child entity to this entity.

        Args:
            child: The child entity to add (can be wrapped or unwrapped)
        """
        # Get the underlying entity if child is wrapped
        child_entity = child.unwrap() if hasattr(child, "unwrap") else child

        if child_entity not in self._wrapped._children:
            self._wrapped._children.append(child_entity)
            child_entity._parent = self._wrapped

    def remove_child(self, child) -> None:
        """
        Remove a child entity from this entity.

        Args:
            child: The child entity to remove (can be wrapped or unwrapped)
        """
        # Get the underlying entity if child is wrapped
        child_entity = child.unwrap() if hasattr(child, "unwrap") else child

        if child_entity in self._wrapped._children:
            self._wrapped._children.remove(child_entity)
            child_entity._parent = None

    def get_root(self):
        """
        Get the root entity of the hierarchy.

        Returns:
            The root entity
        """
        if self.is_root or self.parent is None:
            return self._wrapped
        # Handle case where parent might also be wrapped
        parent_root = (
            self.parent.get_root() if hasattr(self.parent, "get_root") else self.parent
        )
        return parent_root

    def get_ancestors(self) -> list:
        """
        Get all ancestors of this entity (from parent to root).

        Returns:
            List of ancestor entities
        """
        ancestors = []
        current = self.parent
        while current is not None:
            ancestors.append(current)
            current = current.parent if hasattr(current, "parent") else None
        return ancestors

    def get_descendants(self) -> list:
        """
        Get all descendants of this entity.

        Returns:
            List of descendant entities
        """
        descendants = []
        for child in self._wrapped._children:
            descendants.append(child)
            # Handle case where child might be wrapped
            if hasattr(child, "get_descendants"):
                descendants.extend(child.get_descendants())
            elif hasattr(child, "_children"):
                # Recursively get descendants for unwrapped children
                child_wrapper = HierarchyWrapper(child)
                descendants.extend(child_wrapper.get_descendants())
        return descendants

    def find_by_condition(self, condition: Callable[[Any], bool]):
        """
        Find the first descendant (including self) that satisfies a condition.

        Args:
            condition: A function that takes an entity and returns True if it matches

        Returns:
            The first matching entity, or None
        """
        if condition(self._wrapped):
            return self._wrapped
        for child in self._wrapped._children:
            # Check the child directly
            if condition(child):
                return child
            # Recursively search child's descendants
            if hasattr(child, "_children"):
                child_wrapper = HierarchyWrapper(child)
                found = child_wrapper.find_by_condition(condition)
                if found:
                    return found
        return None
