"""
Wrapper classes for the molpy framework.

This module provides wrapper classes that can be composed to add functionality
to molecular entities. Replaces the old mixin-based approach with a more
flexible composition pattern.
"""

import copy
from typing import Any, Callable, Generic, Self, TypeVar

import numpy as np
from numpy.typing import ArrayLike

from ..op import rotate_by_rodrigues
from .entity import Struct

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

        # Find other Wrapper bases from MRO (excluding the current class and base Wrapper)
        wrapper_bases = []
        for base in cls.__mro__[1:]:  # Skip cls itself
            if (
                issubclass(base, Wrapper)
                and base != Wrapper
                and base not in wrapper_bases  # Avoid duplicates
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
        self._wrapped = wrapped
        self._key_mapping = key_mapping if key_mapping else {}

        # Run __post_init__ hooks for auto-composition
        if not self._wrappers_initialized:
            remaining_props = self._run_post_init_hooks(**props)

            # If no wrapped entity, create a Struct with remaining props
            if self._wrapped is None:
                self._wrapped = Struct(**remaining_props)
            # Otherwise, assign remaining props to the innermost Struct
            elif remaining_props:
                innermost = self._get_innermost()
                for k, v in remaining_props.items():
                    innermost[k] = v

            self._wrappers_initialized = True

    def _run_post_init_hooks(self, **props):
        """
        Run __post_init__ hooks for auto-composition.

        Each wrapper's __post_init__ can consume (pop) kwargs it needs.
        Remaining kwargs are passed down and eventually assigned to innermost Struct.

        The __post_init__ methods are called following MRO order (most derived first).
        Each can pop kwargs and return the remaining dict.

        Returns:
            dict: Remaining props after all __post_init__ calls
        """
        remaining_props = dict(props)

        # Get all __post_init__ methods in MRO order
        # We want to call them from most derived to least derived
        post_inits = []
        for klass in type(self).__mro__:
            if klass is object:
                continue
            # Only get __post_init__ defined directly on this class
            if "__post_init__" in klass.__dict__:
                post_inits.append(klass.__post_init__)

        # Call each __post_init__ in order
        for post_init_method in post_inits:
            result = post_init_method(self, **remaining_props)
            if result is not None:
                remaining_props = result

        return remaining_props

    def __post_init__(self, **props):
        """
        Post-initialization hook for wrapper setup.

        Override this method in subclasses to perform one-time initialization
        that depends on the wrapped entity being fully constructed.

        Subclasses should:
        1. Pop any kwargs they need from props
        2. Process those kwargs to set up wrapper-specific state
        3. Return the remaining props dict (or None to keep all remaining)

        Args:
            **props: Properties passed during initialization

        Returns:
            dict or None: Remaining props to pass down, or None to keep all
        """
        return props

    def unwrap(self) -> T:
        """Get the wrapped entity."""
        return self._wrapped

    def _get_innermost(self) -> Any:
        """Get the innermost (unwrapped) entity in the wrapper chain."""
        try:
            obj = object.__getattribute__(self, "_wrapped")
        except AttributeError:
            return None

        if obj is None:
            return None

        while isinstance(obj, Wrapper):
            obj = obj._wrapped
        return obj

    def __getattr__(self, name: str):
        """
        Get attribute by searching through wrapper chain.

        Searches from outer to inner layers, raising AttributeError if not found.
        """
        sentinel = object()

        # Access _wrapped using object.__getattribute__ to avoid recursion
        try:
            obj = object.__getattribute__(self, "_wrapped")
        except AttributeError:
            raise AttributeError(f"{type(self).__name__} has no attribute {name!r}")

        if obj is None:
            raise AttributeError(f"{type(self).__name__} has no attribute {name!r}")

        seen = set()

        while True:
            oid = id(obj)
            if oid in seen:
                raise RuntimeError("cycle detected in _wrapped chain")
            seen.add(oid)

            try:
                return obj[name]  # dict / dict-like
            except (KeyError, TypeError, AttributeError):
                pass

            val = getattr(obj, name, sentinel)
            if val is not sentinel:
                return val

            next_obj = getattr(obj, "_wrapped", sentinel)
            if next_obj is sentinel:
                break
            obj = next_obj

        raise AttributeError(
            f"{type(self).__name__} (and its wrapped chain) has no attribute {name!r}"
        )

    def __setattr__(self, name: str, value: Any) -> None:
        """
        Set attribute by searching from outer to inner layers.

        Private attributes (_*) are set on the wrapper itself.
        For other attributes:
        - If exists in any layer, update it there
        - If doesn't exist, assign to innermost Struct
        """
        # Private attributes and class properties go to wrapper
        if name.startswith("_") or hasattr(type(self), name):
            super().__setattr__(name, value)
            return

        # If attribute already exists in wrapper's __dict__, update it
        try:
            wrapper_dict = object.__getattribute__(self, "__dict__")
            if name in wrapper_dict:
                super().__setattr__(name, value)
                return
        except AttributeError:
            pass

        # Get _wrapped using object.__getattribute__ to avoid recursion
        try:
            obj = object.__getattribute__(self, "_wrapped")
        except AttributeError:
            # During __init__, _wrapped might not exist yet
            super().__setattr__(name, value)
            return

        # If _wrapped is None (during initialization), set on wrapper
        if obj is None:
            super().__setattr__(name, value)
            return

        # Search through wrapper chain for existing attribute
        seen = set()

        while obj is not None:
            oid = id(obj)
            if oid in seen:
                raise RuntimeError("cycle detected in _wrapped chain")
            seen.add(oid)

            # Check if attribute exists in current layer
            if hasattr(obj, name):
                setattr(obj, name, value)
                return

            # Check if it's dict-like and contains the key
            if hasattr(obj, "__contains__") and name in obj:
                obj[name] = value
                return

            # Move to next layer
            obj = getattr(obj, "_wrapped", None)

        # Attribute not found anywhere, assign to innermost
        innermost = self._get_innermost()
        if innermost is not None:
            if hasattr(innermost, "__setitem__"):
                innermost[name] = value
            else:
                setattr(innermost, name, value)
        else:
            # If no innermost exists (during initialization), set on wrapper
            super().__setattr__(name, value)

    def __getitem__(self, key, /):
        """
        Get item by searching from outer to inner layers.

        Searches through the wrapper chain from outermost to innermost,
        raising KeyError if not found anywhere.
        """
        if key in self._key_mapping:
            key = self._key_mapping[key]

        # Search through wrapper chain from outer to inner
        obj = self
        seen = set()

        while True:
            oid = id(obj)
            if oid in seen:
                raise RuntimeError("cycle detected in _wrapped chain")
            seen.add(oid)

            # Try direct __dict__ access first (for wrapper's own attributes)
            if hasattr(obj, "__dict__") and key in obj.__dict__:
                return obj.__dict__[key]

            # Try dict-like access
            if hasattr(obj, "_wrapped"):
                wrapped = obj._wrapped
                if hasattr(wrapped, "__getitem__"):
                    try:
                        return wrapped[key]
                    except (KeyError, TypeError):
                        pass

                # Move to next layer
                if isinstance(wrapped, Wrapper):
                    obj = wrapped
                    continue
                else:
                    # Reached innermost non-wrapper, try one more time
                    if hasattr(wrapped, "__getitem__"):
                        return wrapped[key]  # Let it raise KeyError
                    break
            else:
                break

        raise KeyError(f"Key {key!r} not found in wrapper chain")

    def __setitem__(self, key: Any, value: Any) -> None:
        """
        Set item by searching from outer to inner layers.

        If key exists in any layer, update it there.
        If key doesn't exist anywhere, assign to innermost Struct.
        """
        if key in self._key_mapping:
            key = self._key_mapping[key]

        # Search through wrapper chain to find existing key
        obj = self
        seen = set()

        while True:
            oid = id(obj)
            if oid in seen:
                raise RuntimeError("cycle detected in _wrapped chain")
            seen.add(oid)

            # Check if key exists in current layer's __dict__
            if hasattr(obj, "__dict__") and key in obj.__dict__:
                obj.__dict__[key] = value
                return

            # Check if key exists in wrapped entity
            if hasattr(obj, "_wrapped"):
                wrapped = obj._wrapped
                if hasattr(wrapped, "__contains__") and key in wrapped:
                    wrapped[key] = value
                    return

                # Move to next layer
                if isinstance(wrapped, Wrapper):
                    obj = wrapped
                    continue
                else:
                    # Reached innermost, check one more time
                    if hasattr(wrapped, "__contains__") and key in wrapped:
                        wrapped[key] = value
                        return
                    break
            else:
                break

        # Key not found anywhere, assign to innermost
        innermost = self._get_innermost()
        innermost[key] = value

    def __contains__(self, key: Any) -> bool:
        """
        Check if key exists anywhere in the wrapper chain.

        Searches from outer to inner layers.
        """
        if key in self._key_mapping:
            key = self._key_mapping[key]

        # Search through wrapper chain
        obj = self
        seen = set()

        while True:
            oid = id(obj)
            if oid in seen:
                return False
            seen.add(oid)

            # Check wrapper's __dict__
            if hasattr(obj, "__dict__") and key in obj.__dict__:
                return True

            # Check wrapped entity
            if hasattr(obj, "_wrapped"):
                wrapped = obj._wrapped
                if hasattr(wrapped, "__contains__") and key in wrapped:
                    return True

                # Move to next layer
                if isinstance(wrapped, Wrapper):
                    obj = wrapped
                    continue
                else:
                    # Reached innermost
                    return hasattr(wrapped, "__contains__") and key in wrapped
            else:
                return False

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

    def copy(self: Self) -> Self:
        """
        Create a deep copy of the wrapper and its wrapped entity.

        Always returns a complete deep copy, including all wrapper layers
        and the innermost wrapped entity.

        Returns:
            A new instance with the same wrapper hierarchy and deep-copied data
        """
        return copy.deepcopy(self)

    def __call__(self):
        return self.copy()


class Spatial(Wrapper):
    """
    Wrapper class providing spatial operations for entities.

    Defines spatial functionality for entities that have xyz coordinates
    and can perform geometric transformations.
    """

    def __post_init__(self, **props):
        """
        Initialize Spatial wrapper.

        Can consume 'anchor' prop to specify default anchor point for move_to.
        """
        # Store anchor if provided (can be int index or coordinate)
        self._anchor = props.pop("anchor", None)
        return props

    def distance_to(self, other) -> float:
        """
        Calculate the distance to another spatial entity.

        Args:
            other: Another spatial entity (can be wrapped or unwrapped)

        Returns:
            The Euclidean distance between the two entities
        """
        return float(np.linalg.norm(self.xyz - other.xyz))

    def move(self, vector: ArrayLike) -> Self:
        """
        Translate the entity by a given vector.

        Args:
            vector: Translation vector
        """
        self.xyz = self.xyz + np.array(vector)

        return self

    def move_to(self, target: ArrayLike, anchor: ArrayLike | int | None = None) -> Self:
        """
        Move the entity so that a reference atom is at the target position.
        All other atoms maintain their relative positions.

        Args:
            target: Target position for the reference atom
            anchor: Either:
                   - None: use self._anchor if set, otherwise use first atom (xyz[0])
                   - int: use xyz[anchor] as reference point
                   - ArrayLike: use this coordinate as reference point
        """
        target = np.array(target, dtype=float)

        # Determine which anchor to use
        if anchor is None and hasattr(self, "_anchor"):
            anchor = self._anchor

        # Get current position of reference atom
        xyz = np.atleast_2d(self.xyz)  # Ensure 2D array for consistent indexing
        if anchor is None:
            anchor_pos = xyz[0]
        elif isinstance(anchor, int):
            anchor_pos = xyz[anchor]
        else:
            anchor_pos = np.array(anchor, dtype=float)

        # Calculate translation vector
        translation_vector = target - anchor_pos

        # Apply translation
        self.move(translation_vector)

        # Verify the reference atom is now at target position
        # (check the position that was used as anchor)
        xyz_after = np.atleast_2d(self.xyz)
        if isinstance(anchor, int):
            new_pos = xyz_after[anchor]
        elif anchor is None:
            new_pos = xyz_after[0]
        else:
            new_pos = anchor_pos + translation_vector

        assert np.allclose(
            new_pos, target
        ), f"Reference atom not moved to target position. Expected {target}, got {new_pos}"

        return self

    def scale(self, factor: float | ArrayLike, origin: ArrayLike | None = None) -> Self:
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
        return self

    def rotate(
        self, axis: ArrayLike, angle: float, origin: ArrayLike | None = None
    ) -> Self:
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
        return self

    def rotate_to(
        self, target_direction: ArrayLike, current_direction: ArrayLike
    ) -> Self:
        """
        Rotate the entity so that a current direction aligns with a target direction.
        Args:
            target_direction: The desired direction vector after rotation.
            current_direction: The current direction vector before rotation.
        """
        target_direction = np.array(target_direction, dtype=float)
        current_direction = np.array(current_direction, dtype=float)
        target_direction /= np.linalg.norm(target_direction)
        current_direction /= np.linalg.norm(current_direction)
        rotation_axis = np.cross(current_direction, target_direction)
        if np.linalg.norm(rotation_axis) < 1e-12:
            # Directions are parallel, no rotation needed
            return self
        rotation_axis /= np.linalg.norm(rotation_axis)
        cos_angle = np.clip(np.dot(current_direction, target_direction), -1.0, 1.0)
        angle = np.arccos(cos_angle)
        self.rotate(rotation_axis, angle)
        return self

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


class Hierarchy(Wrapper):
    """
    Wrapper class providing hierarchical structure functionality.

    Allows entities to have parent-child relationships and provides
    methods for navigating and manipulating the hierarchy.
    """

    _parent: Wrapper | None
    _children: list[Wrapper]

    def __init__(self, wrapped=None, **props):
        """
        Initialize hierarchy wrapper.

        Args:
            wrapped: The entity to wrap
            **props: Properties to pass to parent __init__
        """
        super().__init__(wrapped, **props)

    def __post_init__(self, **props):
        """
        Initialize Hierarchy wrapper.

        Hierarchy doesn't consume any specific props, so pass them all down.
        Also initializes parent/children if needed.
        """
        if not hasattr(self._wrapped, "_parent"):
            self._wrapped._parent = None
        if not hasattr(self._wrapped, "_children"):
            self._wrapped._children = []
        return props

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

    def add_child(self, child) -> Self:
        """
        Add a child entity to this entity.

        Args:
            child: The child entity to add (can be wrapped or unwrapped)
        """
        self._children.append(child)
        self.merge(child)
        if isinstance(child, Hierarchy):
            child._parent = self
        return child

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
                child_wrapper = Hierarchy(child)
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
                child_wrapper = Hierarchy(child)
                found = child_wrapper.find_by_condition(condition)
                if found:
                    return found
        return None
