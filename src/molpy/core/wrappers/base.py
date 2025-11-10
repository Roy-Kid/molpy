from typing import Any, Generic, TypeVar, Self
from ..entity import Assembly
import copy

T = TypeVar("T")

class Wrapper(Generic[T]):
    """
    Base class for all wrappers.

    Wrappers provide a composable way to add functionality to entities
    without using inheritance mixins.

    Supports both explicit wrapping and multiple inheritance auto-composition.
    """
    _wrapped: T

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
        wrapped: T,
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
    
    def __repr__(self) -> str:
        """
        Simple representation showing wrapper type and wrapped object type.
        
        For detailed tree view, use tree_repr().
        """
        wrapped = self.unwrap()
        wrapped_type = type(wrapped).__name__
        return f"{type(self).__name__}({wrapped_type})"
    
    def tree_repr(self, _indent: str = "", _is_last: bool = True) -> str:
        """
        Generate a tree-style representation of the wrapper stack.
        
        Returns:
            Multi-line string showing the wrapper hierarchy
            
        Example:
            Monomer(Atomistic)
            └── Atomistic: 2 atoms, 1 bond
        """
        # Get current wrapper info
        lines = []
        
        # Prefix for current level
        if _indent == "":
            # Root level - no prefix
            prefix = ""
        else:
            # Child level
            prefix = _indent + ("└── " if _is_last else "├── ")
        
        # Current wrapper line
        wrapper_info = self._repr_info()
        lines.append(f"{prefix}{type(self).__name__}: {wrapper_info}")
        
        # Get wrapped object
        wrapped = self.unwrap()
        
        # Prepare indent for next level
        if _indent == "":
            next_indent = ""
        else:
            next_indent = _indent + ("    " if _is_last else "│   ")
        
        # Recursively print wrapped object
        if isinstance(wrapped, Wrapper):
            # Another wrapper - recurse
            lines.append(wrapped.tree_repr(_indent=next_indent, _is_last=True))
        else:
            # Final wrapped object
            final_prefix = next_indent + "└── "
            if hasattr(wrapped, '__repr__'):
                final_repr = repr(wrapped)
            else:
                final_repr = f"{type(wrapped).__name__}"
            lines.append(f"{final_prefix}{final_repr}")
        
        return "\n".join(lines)
    
    def _repr_info(self) -> str:
        """
        Override this to provide wrapper-specific info in tree_repr.
        
        Returns:
            String describing this wrapper's state
        """
        return "wrapper"
