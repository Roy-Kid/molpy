"""Registry for managing parameter formatters with MRO-based lookup.

This module provides a thread-safe registry for parameter formatters that
supports inheritance-based lookup through the Method Resolution Order (MRO).
"""

from threading import RLock
from typing import Any, Callable, Optional


class ParameterFormatterRegistry:
    """Thread-safe registry for parameter formatters.

    The registry maps Style classes to formatter functions. When looking up
    a formatter, it walks the style class's MRO to find the most specific
    formatter available.

    Examples:
        >>> registry = ParameterFormatterRegistry()
        >>> registry.register(BondHarmonicStyle, format_bond_harmonic)
        >>> formatter = registry.get_formatter(BondHarmonicStyle)
        >>> params = formatter(bond_type)

        Using the registry to format a type:
        >>> params = registry.format(bond_type, bond_style)
    """

    def __init__(self):
        """Initialize an empty formatter registry."""
        self._formatters: dict[type, Callable[[Any], list[float]]] = {}
        self._lock = RLock()

    def register(
        self, style_class: type, formatter: Callable[[Any], list[float]]
    ) -> None:
        """Register a formatter for a style class.

        Args:
            style_class: The Style class to register the formatter for
                        (e.g., BondHarmonicStyle, AngleHarmonicStyle)
            formatter: Callable that takes a Type object and returns list[float]

        Examples:
            >>> registry.register(BondHarmonicStyle, format_bond_harmonic)
            >>> registry.register(MyCustomStyle, lambda typ: [typ.params.kwargs['k']])
        """
        with self._lock:
            self._formatters[style_class] = formatter

    def get_formatter(
        self, style_class: type
    ) -> Optional[Callable[[Any], list[float]]]:
        """Get formatter for a style class with MRO-based lookup.

        Looks up the formatter for the given style class. If no exact match
        is found, walks the MRO to find a formatter registered for a parent class.

        Args:
            style_class: The Style class to get the formatter for

        Returns:
            The formatter function if found, None otherwise

        Examples:
            >>> formatter = registry.get_formatter(BondHarmonicStyle)
            >>> if formatter:
            ...     params = formatter(bond_type)
        """
        with self._lock:
            # Direct lookup
            if style_class in self._formatters:
                return self._formatters[style_class]

            # MRO lookup - find formatter for parent class
            for base in style_class.__mro__[1:]:
                if base in self._formatters:
                    return self._formatters[base]

            return None

    def format(self, typ: Any, style: Any) -> list[float]:
        """Format a type using the registered formatter for its style.

        Args:
            typ: Type object (BondType, AngleType, etc.)
            style: Style object that contains this type

        Returns:
            List of numeric parameters formatted for output

        Raises:
            ValueError: If no formatter is registered for the style class

        Examples:
            >>> params = registry.format(bond_type, bond_style)
        """
        formatter = self.get_formatter(type(style))
        if formatter is None:
            raise ValueError(
                f"No formatter registered for style class {type(style).__name__}. "
                f"Available formatters: {list(self._formatters.keys())}"
            )
        return formatter(typ)

    def copy(self) -> "ParameterFormatterRegistry":
        """Create a shallow copy of this registry.

        The copy contains the same formatter mappings but is independent,
        so registering new formatters in the copy won't affect the original.

        Returns:
            A new ParameterFormatterRegistry with the same formatters

        Examples:
            >>> custom_registry = default_registry.copy()
            >>> custom_registry.register(MyStyle, my_formatter)
            >>> # default_registry is unchanged
        """
        new_registry = ParameterFormatterRegistry()
        with self._lock:
            new_registry._formatters = self._formatters.copy()
        return new_registry

    def list_formatters(self) -> list[type]:
        """List all registered style classes.

        Returns:
            List of style classes that have formatters registered

        Examples:
            >>> registered = registry.list_formatters()
            >>> print(f"Formatters registered for: {registered}")
        """
        with self._lock:
            return list(self._formatters.keys())
