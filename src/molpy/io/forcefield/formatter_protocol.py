"""Protocol definition for parameter formatters.

This module defines the interface that parameter formatters must implement.
Formatters convert Type objects (BondType, AngleType, etc.) into lists of
numeric parameters suitable for force field file output.
"""

from typing import Any, Protocol, runtime_checkable


@runtime_checkable
class ParameterFormatter(Protocol):
    """Protocol for parameter formatters.

    Formatters are callables that take a Type object and return a list of
    numeric parameters. They can be implemented as functions or callable classes.

    Examples:
        Function-based formatter:
        >>> def format_bond_harmonic(typ) -> list[float]:
        ...     return [typ.params.kwargs['k'], typ.params.kwargs['r0']]

        Class-based formatter:
        >>> class BondFormatter:
        ...     def __init__(self, unit_conversion: float = 1.0):
        ...         self.unit_conversion = unit_conversion
        ...
        ...     def __call__(self, typ) -> list[float]:
        ...         k = typ.params.kwargs['k'] * self.unit_conversion
        ...         r0 = typ.params.kwargs['r0']
        ...         return [k, r0]
    """

    def __call__(self, typ: Any) -> list[float]:
        """Format a Type object into numeric parameters.

        Args:
            typ: Type object (BondType, AngleType, DihedralType, etc.)
                 containing parameters in typ.params.kwargs

        Returns:
            List of numeric parameters in the format expected by the
            force field file writer

        Raises:
            KeyError: If required parameters are missing from typ.params.kwargs
            TypeError: If parameters cannot be converted to float
            ValueError: If parameters are invalid for this formatter
        """
        ...
