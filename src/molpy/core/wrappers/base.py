"""
Base wrapper class for Assembly objects.

Provides recursive unwrapping and transparent delegation to wrapped assemblies.
"""

from typing import Any, Generic, TypeVar, Self
from ..entity import Assembly

T = TypeVar("T", bound=Assembly)


class Wrapper(Generic[T]):
    """
    Base wrapper class for Assembly objects.
    
    Provides:
    - Recursive unwrapping: automatically unwraps nested wrappers to innermost Assembly
    - Transparent delegation: forwards attribute/method access to wrapped object
    - Type preservation: wrappers maintain identity semantics
    
    Type parameter:
        T: The type of Assembly being wrapped (bound to Assembly)
    """
    
    inner: Assembly  # The innermost wrapped Assembly (runtime type)
    
    def __init__(self, wrapped: T | "Wrapper[T]", **props):
        """Initialize wrapper by recursively unwrapping to innermost Assembly.
        
        Args:
            wrapped: Either a concrete Assembly instance or another Wrapper
            **props: Additional properties (passed to __post_init__)
        
        The wrapper automatically unwraps nested wrappers to find the innermost
        Assembly, storing it in `self.inner`.
        """
        # Recursively unwrap to innermost Assembly
        current = wrapped
        while isinstance(current, Wrapper):
            current = current.inner
        
        # Store innermost instance
        self.inner = current
        
        # Call post-init hook for subclass initialization
        self.__post_init__(**props)
    
    def __post_init__(self, **props):
        """Post-initialization hook for subclass setup.
        
        Override this in subclasses to:
        - Initialize wrapper-specific state
        - Consume relevant kwargs from props
        - Return remaining props dict (or None)
        
        Args:
            **props: Properties passed during initialization
            
        Returns:
            dict or None: Remaining props to pass down, or None
        """
        return props
    
    def unwrap(self) -> T:
        """Return the innermost wrapped Assembly.
        
        Returns:
            The innermost Assembly instance (after all unwrapping)
        """
        return self.inner  # type: ignore[return-value]
    
    def __getattr__(self, name: str):
        """Delegate attribute/method access to inner object.
        
        This enables transparent delegation of all undefined methods/attributes
        to the wrapped Assembly. Only methods explicitly defined in wrapper
        subclasses will not be delegated.
        
        Args:
            name: Attribute or method name to access
            
        Returns:
            The attribute/method from inner object
            
        Raises:
            AttributeError: If inner object doesn't have the attribute
        """
        # First try dict-like access (for Assembly props)
        if hasattr(self.inner, '__getitem__') and hasattr(self.inner, '__contains__'):
            if name in self.inner:
                return self.inner[name]
        
        # Then try regular attribute access
        try:
            return getattr(self.inner, name)
        except AttributeError:
            raise AttributeError(
                f"{type(self).__name__} has no attribute {name!r}"
            )
    
    def __setattr__(self, name: str, value: Any) -> None:
        """Set attribute on wrapper or delegate to inner.
        
        Private attributes (starting with _) and wrapper-specific attributes
        are set on the wrapper itself. Other attributes are delegated to inner.
        
        Args:
            name: Attribute name
            value: Value to set
        """
        # Private attributes and 'inner' always go to wrapper
        if name.startswith("_") or name == "inner":
            super().__setattr__(name, value)
            return
        
        # If attribute exists in wrapper's __dict__, update it
        try:
            wrapper_dict = object.__getattribute__(self, "__dict__")
            if name in wrapper_dict:
                super().__setattr__(name, value)
                return
        except AttributeError:
            pass
        
        # Try to set in inner (dict-like first, then attribute)
        if hasattr(self.inner, '__setitem__'):
            self.inner[name] = value
        else:
            setattr(self.inner, name, value)
    
    def __getitem__(self, key: str):
        """Delegate dict-style access to inner object.
        
        Args:
            key: Key to access
            
        Returns:
            Value from inner object
        """
        return self.inner[key]
    
    def __setitem__(self, key: str, value: Any) -> None:
        """Delegate dict-style write to inner object.
        
        Args:
            key: Key to set
            value: Value to assign
        """
        self.inner[key] = value
    
    def __contains__(self, key: str) -> bool:
        """Check if key exists in inner object.
        
        Args:
            key: Key to check
            
        Returns:
            True if key exists in inner
        """
        return key in self.inner
    
    def get(self, key: str, default=None):
        """Get value from inner with default.
        
        Args:
            key: Key to get
            default: Default value if key not found
            
        Returns:
            Value from inner or default
        """
        return self.inner.get(key, default)
    
    def copy(self: Self) -> Self:
        """Create a deep copy of the wrapper and its wrapped entity.
        
        Returns:
            A new wrapper instance with deep-copied inner Assembly
        """
        import copy as copy_module
        
        # Deep copy the inner
        new_inner = copy_module.deepcopy(self.inner)
        
        # Create new wrapper with copied inner
        new_wrapper = type(self)(new_inner)  # type: ignore[arg-type]
        
        # Copy wrapper-specific attributes
        for key, value in self.__dict__.items():
            if key != 'inner':
                new_wrapper.__dict__[key] = copy_module.deepcopy(value)
        
        return new_wrapper
    
    def __call__(self) -> Self:
        """Calling a wrapper creates a deep copy.
        
        Returns:
            New wrapper with copied structure
        """
        return self.copy()
    
    def __repr__(self) -> str:
        """Simple representation showing wrapper type and wrapped object."""
        return f"<{type(self).__name__} wrapping {self.inner!r}>"
