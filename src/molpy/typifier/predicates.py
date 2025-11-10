"""Predicate factories for SMARTS pattern matching.

This module provides predicate functions for vertex and edge matching
in SMARTS patterns. Predicates are pure callable objects that check
specific properties against atom/bond attributes.
"""

from typing import Any, Callable, Protocol
from dataclasses import dataclass


# ===================================================================
#                     Predicate Protocol
# ===================================================================

class VertexPredicate(Protocol):
    """Protocol for vertex (atom) predicates."""
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        """Check if vertex attributes satisfy the predicate.
        
        Args:
            attrs: Vertex attributes dict with keys like:
                - element: str (e.g., "C", "N", "O")
                - is_aromatic: bool
                - charge: int
                - degree: int (number of bonds)
                - hyb: int | None (1=sp, 2=sp2, 3=sp3)
                - in_ring: bool
                - Any other custom attributes
        
        Returns:
            True if predicate is satisfied
        """
        ...


class EdgePredicate(Protocol):
    """Protocol for edge (bond) predicates."""
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        """Check if edge attributes satisfy the predicate.
        
        Args:
            attrs: Edge attributes dict with keys like:
                - order: int | str (1, 2, 3, ":")
                - is_aromatic: bool
                - is_in_ring: bool
                - Any other custom attributes
        
        Returns:
            True if predicate is satisfied
        """
        ...


# ===================================================================
#                     Predicate Metadata
# ===================================================================

@dataclass
class PredicateMeta:
    """Metadata for predicate scoring and debugging."""
    name: str
    weight: int = 1  # Weight for specificity scoring


# ===================================================================
#                  Vertex Predicate Factories
# ===================================================================

class ElementPredicate:
    """Predicate for element matching."""
    
    def __init__(self, element: str):
        """
        Args:
            element: Element symbol (e.g., "C", "N", "O", "*")
        """
        self.element = element.upper() if element != "*" else "*"
        self.meta = PredicateMeta(name=f"element={element}", weight=0)
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        if self.element == "*":
            return True
        attr_element = attrs.get("element", "")
        if isinstance(attr_element, str):
            return attr_element.upper() == self.element
        # Handle atomic number
        attr_num = attrs.get("number", None)
        if attr_num is not None:
            from molpy.core.element import Element
            try:
                return Element(self.element).number == attr_num
            except:
                return False
        return False
    
    def __repr__(self) -> str:
        return f"is_element({self.element!r})"


class AromaticPredicate:
    """Predicate for aromaticity."""
    
    def __init__(self, is_aromatic: bool = True):
        self.is_aromatic = is_aromatic
        self.meta = PredicateMeta(
            name=f"is_aromatic={is_aromatic}",
            weight=2
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("is_aromatic", False) == self.is_aromatic
    
    def __repr__(self) -> str:
        return f"is_aromatic({self.is_aromatic})"


class ChargePredicate:
    """Predicate for formal charge."""
    
    def __init__(self, charge: int):
        self.charge = charge
        self.meta = PredicateMeta(
            name=f"charge={charge}",
            weight=1
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("charge", 0) == self.charge
    
    def __repr__(self) -> str:
        return f"charge({self.charge})"


class DegreePredicate:
    """Predicate for connectivity degree."""
    
    def __init__(self, degree: int):
        self.degree = degree
        self.meta = PredicateMeta(
            name=f"degree={degree}",
            weight=1
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("degree", 0) == self.degree
    
    def __repr__(self) -> str:
        return f"degree({self.degree})"


class HybridizationPredicate:
    """Predicate for hybridization state."""
    
    def __init__(self, hyb: int):
        """
        Args:
            hyb: Hybridization (1=sp, 2=sp2, 3=sp3)
        """
        self.hyb = hyb
        self.meta = PredicateMeta(
            name=f"hyb={hyb}",
            weight=1
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("hyb", None) == self.hyb
    
    def __repr__(self) -> str:
        return f"hyb({self.hyb})"


class InRingPredicate:
    """Predicate for ring membership."""
    
    def __init__(self, in_ring: bool = True):
        self.in_ring = in_ring
        self.meta = PredicateMeta(
            name=f"in_ring={in_ring}",
            weight=2
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("in_ring", False) == self.in_ring
    
    def __repr__(self) -> str:
        return f"in_ring({self.in_ring})"


class RingSizePredicate:
    """Predicate for ring size."""
    
    def __init__(self, size: int):
        self.size = size
        self.meta = PredicateMeta(
            name=f"ring_size={size}",
            weight=3
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        cycles = attrs.get("cycles", set())
        for cycle in cycles:
            if len(cycle) == self.size:
                return True
        return False
    
    def __repr__(self) -> str:
        return f"ring_size({self.size})"


class RingCountPredicate:
    """Predicate for number of rings."""
    
    def __init__(self, count: int):
        self.count = count
        self.meta = PredicateMeta(
            name=f"ring_count={count}",
            weight=2
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        cycles = attrs.get("cycles", set())
        return len(cycles) == self.count
    
    def __repr__(self) -> str:
        return f"ring_count({self.count})"


class CustomPredicate:
    """Custom predicate with user-provided function."""
    
    def __init__(self, func: Callable[[dict[str, Any]], bool], 
                 name: str = "custom", weight: int = 4):
        self.func = func
        self.meta = PredicateMeta(name=name, weight=weight)
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return self.func(attrs)
    
    def __repr__(self) -> str:
        return f"custom({self.meta.name})"


# ===================================================================
#                  Edge Predicate Factories
# ===================================================================

class BondOrderPredicate:
    """Predicate for bond order."""
    
    def __init__(self, order: int | str):
        """
        Args:
            order: Bond order (1, 2, 3, or ":" for aromatic)
        """
        self.order = order
        self.meta = PredicateMeta(
            name=f"bond_order={order}",
            weight=3
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        attr_order = attrs.get("order", 1)
        # Handle both int and string representations
        if isinstance(attr_order, str) and isinstance(self.order, str):
            return attr_order == self.order
        if isinstance(attr_order, (int, float)) and isinstance(self.order, (int, float)):
            return int(attr_order) == int(self.order)
        # Mixed types - convert to string for comparison
        return str(attr_order) == str(self.order)
    
    def __repr__(self) -> str:
        return f"bond_order({self.order})"


class BondAromaticPredicate:
    """Predicate for aromatic bonds."""
    
    def __init__(self, is_aromatic: bool = True):
        self.is_aromatic = is_aromatic
        self.meta = PredicateMeta(
            name=f"bond_aromatic={is_aromatic}",
            weight=2
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("is_aromatic", False) == self.is_aromatic
    
    def __repr__(self) -> str:
        return f"bond_aromatic({self.is_aromatic})"


class BondInRingPredicate:
    """Predicate for bonds in rings."""
    
    def __init__(self, is_in_ring: bool = True):
        self.is_in_ring = is_in_ring
        self.meta = PredicateMeta(
            name=f"bond_in_ring={is_in_ring}",
            weight=2
        )
    
    def __call__(self, attrs: dict[str, Any]) -> bool:
        return attrs.get("is_in_ring", False) == self.is_in_ring
    
    def __repr__(self) -> str:
        return f"bond_in_ring({self.is_in_ring})"


# ===================================================================
#                  Convenience Functions
# ===================================================================

def is_element(element: str) -> ElementPredicate:
    """Create element predicate.
    
    Example:
        >>> pred = is_element("C")
        >>> pred({"element": "C"})
        True
    """
    return ElementPredicate(element)


def is_aromatic(value: bool = True) -> AromaticPredicate:
    """Create aromaticity predicate."""
    return AromaticPredicate(value)


def charge(value: int) -> ChargePredicate:
    """Create charge predicate."""
    return ChargePredicate(value)


def degree(value: int) -> DegreePredicate:
    """Create degree predicate."""
    return DegreePredicate(value)


def hyb(value: int) -> HybridizationPredicate:
    """Create hybridization predicate."""
    return HybridizationPredicate(value)


def in_ring(value: bool = True) -> InRingPredicate:
    """Create ring membership predicate."""
    return InRingPredicate(value)


def ring_size(size: int) -> RingSizePredicate:
    """Create ring size predicate."""
    return RingSizePredicate(size)


def ring_count(count: int) -> RingCountPredicate:
    """Create ring count predicate."""
    return RingCountPredicate(count)


def bond_order(order: int | str) -> BondOrderPredicate:
    """Create bond order predicate."""
    return BondOrderPredicate(order)


def bond_aromatic(value: bool = True) -> BondAromaticPredicate:
    """Create bond aromaticity predicate."""
    return BondAromaticPredicate(value)


def bond_in_ring(value: bool = True) -> BondInRingPredicate:
    """Create bond in ring predicate."""
    return BondInRingPredicate(value)


def custom_vertex(func: Callable[[dict[str, Any]], bool], 
                  name: str = "custom", weight: int = 4) -> CustomPredicate:
    """Create custom vertex predicate."""
    return CustomPredicate(func, name, weight)


def wildcard() -> ElementPredicate:
    """Create wildcard predicate that matches any atom."""
    return ElementPredicate("*")
