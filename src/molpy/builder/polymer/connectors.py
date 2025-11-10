"""
Connector abstraction for topology-only polymer assembly.

Connectors decide which ports to connect between adjacent monomers
in a linear sequence. This is purely topological - no geometry involved.
"""

from typing import Literal, Any, Mapping, Tuple, Callable, Iterable
from ...core.wrappers.monomer import Monomer, Port
from ..errors import (
    AmbiguousPortsError,
    MissingConnectorRule,
    NoCompatiblePortsError,
    BondKindConflictError,
)


BondKind = Literal["-", "=", "#", ":"]


class ConnectorContext(dict[str, Any]):
    """
    Shared context passed to connectors during linear build.
    
    Contains information like:
    - step: int (current connection step index)
    - sequence: str (full sequence being built)
    - left_label: str (label of left monomer)
    - right_label: str (label of right monomer)
    - audit: list (accumulated connection records)
    """
    pass


class Connector:
    """
    Abstract base for port selection between two adjacent monomers.
    
    This is topology-only: connectors decide WHICH ports to connect,
    not HOW to position them geometrically.
    """
    
    def select_ports(
        self,
        left: Monomer,
        right: Monomer,
        left_ports: Mapping[str, Port],   # unconsumed ports
        right_ports: Mapping[str, Port],  # unconsumed ports
        ctx: ConnectorContext,
    ) -> Tuple[str, str, BondKind | None]:
        """
        Select which ports to connect between left and right monomers.
        
        Args:
            left: Left monomer in the sequence
            right: Right monomer in the sequence
            left_ports: Available (unconsumed) ports on left monomer
            right_ports: Available (unconsumed) ports on right monomer
            ctx: Shared context with step info, sequence, etc.
            
        Returns:
            Tuple of (left_port_name, right_port_name, optional_bond_kind_override)
            
        Raises:
            AmbiguousPortsError: Cannot uniquely determine ports
            NoCompatiblePortsError: No valid port pair found
            MissingConnectorRule: Required rule not found (TableConnector)
        """
        raise NotImplementedError("Subclasses must implement select_ports()")


class AutoConnector(Connector):
    """
    BigSMILES-guided automatic port selection.
    
    Strategy:
    1. If left has port with role='right' and right has role='left' → use those
    2. Else if each side has exactly one unconsumed port → use that pair
    3. Else raise AmbiguousPortsError
    
    This implements the common case where:
    - BigSMILES uses [<] for "left" role and [>] for "right" role
    - We connect left's "right" port to right's "left" port
    """
    
    def select_ports(
        self,
        left: Monomer,
        right: Monomer,
        left_ports: Mapping[str, Port],
        right_ports: Mapping[str, Port],
        ctx: ConnectorContext,
    ) -> Tuple[str, str, BondKind | None]:
        """Select ports using BigSMILES role heuristics."""
        
        # Strategy 1: Try role-based selection (BigSMILES < and >)
        # Case 1a: left has role='right', right has role='left' (normal chain extension)
        left_right_role = [name for name, p in left_ports.items() 
                          if p.data.get('role') == 'right']
        right_left_role = [name for name, p in right_ports.items() 
                          if p.data.get('role') == 'left']
        
        if len(left_right_role) == 1 and len(right_left_role) == 1:
            return (left_right_role[0], right_left_role[0], None)
        
        # Case 1b: left has role='left', right has role='right' (terminus connection)
        # E.g., HO[*:t] (left terminus) connects to CC[>] (right-directed monomer)
        left_left_role = [name for name, p in left_ports.items()
                         if p.data.get('role') == 'left']
        right_right_role = [name for name, p in right_ports.items()
                           if p.data.get('role') == 'right']
        
        if len(left_left_role) == 1 and len(right_right_role) == 1:
            return (left_left_role[0], right_right_role[0], None)
        
        # Strategy 2: If exactly one port on each side, use those
        if len(left_ports) == 1 and len(right_ports) == 1:
            left_name = next(iter(left_ports.keys()))
            right_name = next(iter(right_ports.keys()))
            return (left_name, right_name, None)
        
        # Strategy 3: Ambiguous - cannot decide
        raise AmbiguousPortsError(
            f"Cannot auto-select ports between {ctx.get('left_label')} and {ctx.get('right_label')}: "
            f"left has {len(left_ports)} available ports {list(left_ports.keys())}, "
            f"right has {len(right_ports)} available ports {list(right_ports.keys())}. "
            "Use TableConnector or CallbackConnector to specify explicit rules."
        )


class TableConnector(Connector):
    """
    Rule-based port selection using a lookup table.
    
    Maps (left_label, right_label) → (left_port, right_port [, bond_kind])
    
    Example:
        rules = {
            ("A", "B"): ("1", "2"),
            ("B", "A"): ("3", "1", "="),  # with bond kind override
            ("T", "A"): ("t", "1"),
        }
        connector = TableConnector(rules, fallback=AutoConnector())
    """
    
    def __init__(
        self,
        rules: Mapping[tuple[str, str], tuple[str, str] | tuple[str, str, BondKind]],
        fallback: Connector | None = None,
    ):
        """
        Initialize table connector.
        
        Args:
            rules: Mapping from (left_label, right_label) to port specifications
            fallback: Optional connector to try if pair not in rules
        """
        self.rules = dict(rules)  # Convert to dict for internal use
        self.fallback = fallback
    
    def select_ports(
        self,
        left: Monomer,
        right: Monomer,
        left_ports: Mapping[str, Port],
        right_ports: Mapping[str, Port],
        ctx: ConnectorContext,
    ) -> Tuple[str, str, BondKind | None]:
        """Select ports using table lookup."""
        
        left_label = ctx.get('left_label', '')
        right_label = ctx.get('right_label', '')
        key = (left_label, right_label)
        
        if key in self.rules:
            rule = self.rules[key]
            if len(rule) == 2:
                return (rule[0], rule[1], None)
            else:
                return (rule[0], rule[1], rule[2])
        
        # Try fallback if available
        if self.fallback is not None:
            return self.fallback.select_ports(left, right, left_ports, right_ports, ctx)
        
        # No rule and no fallback
        raise MissingConnectorRule(
            f"No rule found for ({left_label}, {right_label}) and no fallback connector"
        )


class CallbackConnector(Connector):
    """
    User-defined callback for port selection.
    
    Example:
        def my_selector(left, right, left_ports, right_ports, ctx):
            # Custom logic here
            return ("port_out", "port_in", "-")
        
        connector = CallbackConnector(my_selector)
    """
    
    def __init__(
        self, 
        fn: Callable[
            [Monomer, Monomer, Mapping[str, Port], Mapping[str, Port], ConnectorContext],
            Tuple[str, str] | Tuple[str, str, BondKind]
        ]
    ):
        """
        Initialize callback connector.
        
        Args:
            fn: Callable that takes (left, right, left_ports, right_ports, ctx)
                and returns (left_port, right_port [, bond_kind])
        """
        self.fn = fn
    
    def select_ports(
        self,
        left: Monomer,
        right: Monomer,
        left_ports: Mapping[str, Port],
        right_ports: Mapping[str, Port],
        ctx: ConnectorContext,
    ) -> Tuple[str, str, BondKind | None]:
        """Select ports using user callback."""
        
        result = self.fn(left, right, left_ports, right_ports, ctx)
        
        if len(result) == 2:
            return (result[0], result[1], None)
        else:
            return (result[0], result[1], result[2])


class ChainConnector(Connector):
    """
    Try a list of connectors in order; first one that succeeds wins.
    
    Example:
        connector = ChainConnector([
            TableConnector(specific_rules),
            AutoConnector(),
        ])
    """
    
    def __init__(self, connectors: Iterable[Connector]):
        """
        Initialize chain connector.
        
        Args:
            connectors: List of connectors to try in order
        """
        self.connectors = list(connectors)
    
    def select_ports(
        self,
        left: Monomer,
        right: Monomer,
        left_ports: Mapping[str, Port],
        right_ports: Mapping[str, Port],
        ctx: ConnectorContext,
    ) -> Tuple[str, str, BondKind | None]:
        """Try connectors in order until one succeeds."""
        
        errors = []
        for connector in self.connectors:
            try:
                return connector.select_ports(left, right, left_ports, right_ports, ctx)
            except (AmbiguousPortsError, MissingConnectorRule, NoCompatiblePortsError) as e:
                errors.append(f"{type(connector).__name__}: {e}")
                continue
        
        # All connectors failed
        raise AmbiguousPortsError(
            f"All connectors failed for ({ctx.get('left_label')}, {ctx.get('right_label')}): "
            f"{'; '.join(errors)}"
        )
