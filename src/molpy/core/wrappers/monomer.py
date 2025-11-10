from ..entity import Entity, Assembly
from .base import Wrapper
from typing import TypeVar, Self, Generic

T = TypeVar("T", bound=Assembly)


class Port(Entity):
    """
    Port represents a connection point on a monomer.
    
    A Port is an Entity that wraps a target entity (e.g., an atom)
    and provides a named interface for connecting monomers.
    
    Attributes:
        name: Port identifier (e.g., 'in', 'out', 'port_1')
        target: The underlying entity this port connects to
        role: Optional role ('left'/'right' from BigSMILES </>)
        bond_kind: Optional bond type ('-', '=', '#', ':')
        compat: Compatibility spec (set of port names, or '*' for any)
        multiplicity: How many times this port can be used (default 1)
        priority: Selection priority when multiple ports available (default 0)
    """
    
    def __init__(
        self, 
        name: str, 
        target: Entity,
        *,
        role: str | None = None,
        bond_kind: str | None = None,
        compat: set[str] | str | None = None,
        multiplicity: int = 1,
        priority: int = 0,
    ):
        super().__init__()
        self.data['name'] = name
        self.data['target'] = target
        self.data['role'] = role
        self.data['bond_kind'] = bond_kind
        self.data['compat'] = compat if compat is not None else '*'
        self.data['multiplicity'] = multiplicity
        self.data['priority'] = priority
    
    @property
    def name(self) -> str:
        return self.data['name']
    
    @property
    def target(self) -> Entity:
        return self.data['target']
    
    @property
    def role(self) -> str | None:
        return self.data.get('role')
    
    @property
    def bond_kind(self) -> str | None:
        return self.data.get('bond_kind')
    
    @property
    def compat(self) -> set[str] | str:
        return self.data.get('compat', '*')
    
    @property
    def multiplicity(self) -> int:
        return self.data.get('multiplicity', 1)
    
    @multiplicity.setter
    def multiplicity(self, value: int) -> None:
        self.data['multiplicity'] = value
    
    @property
    def priority(self) -> int:
        return self.data.get('priority', 0)
    
    def orientation(self) -> list[float] | None:
        """
        Get explicit orientation vector for this port.
        
        Returns orientation if stored in port or target entity data,
        otherwise None (caller should infer from connectivity).
        
        Returns:
            Unit orientation vector [x, y, z] or None
        """
        # Check port data first
        if 'orientation' in self.data:
            return self.data['orientation']
        
        # Check target entity
        if hasattr(self.target, 'data') and 'orientation' in self.target.data:
            return self.target.data['orientation']
        
        return None
    
    def separation(self) -> float | None:
        """
        Get explicit separation distance for this port.
        
        Returns separation if stored in port or target entity data,
        otherwise None (caller should use VDW-based distance).
        
        Returns:
            Separation distance in Angstroms or None
        """
        # Check port data first
        if 'separation' in self.data:
            return self.data['separation']
        
        # Check target entity
        if hasattr(self.target, 'data') and 'separation' in self.target.data:
            return self.target.data['separation']
        
        return None
    
    def __repr__(self) -> str:
        parts = [f"name={self.name!r}"]
        if self.role:
            parts.append(f"role={self.role!r}")
        if self.bond_kind:
            parts.append(f"bond={self.bond_kind!r}")
        if self.multiplicity != 1:
            parts.append(f"mult={self.multiplicity}")
        return f"Port({', '.join(parts)})"


class Monomer(Wrapper[T], Generic[T]):
    """
    Monomer wraps an Assembly and adds port management.
    
    Ports define connection points for polymer assembly.
    The wrapped Assembly contains the molecular structure,
    and ports reference specific entities within that structure.
    """
    
    def __init__(self, wrapped: T):
        super().__init__(wrapped)
        # Store ports as a private dict on the wrapper layer
        self._ports: dict[str, Port] = {}
    
    def set_port(
        self, 
        name: str, 
        target: Entity,
        *,
        role: str | None = None,
        bond_kind: str | None = None,
        compat: set[str] | str | None = None,
        multiplicity: int = 1,
        priority: int = 0,
    ) -> None:
        """
        Create or update a port pointing to a target entity.
        
        Args:
            name: Port identifier
            target: Entity to connect (typically from wrapped Assembly)
            role: Optional BigSMILES role ('left'/'right')
            bond_kind: Optional bond type specification
            compat: Compatibility spec (set of names or '*')
            multiplicity: How many times port can be used
            priority: Selection priority
        """
        port = Port(
            name, 
            target, 
            role=role,
            bond_kind=bond_kind,
            compat=compat,
            multiplicity=multiplicity,
            priority=priority,
        )
        self._ports[name] = port
    
    def get_port(self, name: str) -> Port | None:
        """Get port by name, or None if not found."""
        return self._ports.get(name)
    
    def port_names(self) -> list[str]:
        """Return list of all port names."""
        return list(self._ports.keys())
    
    @property
    def ports(self) -> dict[str, Port]:
        """Access ports dictionary (read-only view recommended)."""
        return self._ports
    
    def copy(self) -> Self:
        """
        Deep copy the monomer, including wrapped Assembly and ports.
        
        Returns:
            New Monomer with copied structure and remapped ports
        """
        # Copy the wrapped assembly
        wrapped = self.unwrap()
        new_wrapped = wrapped.copy()
        
        # Create new monomer
        new_monomer = type(self)(new_wrapped)
        
        # Get entity mapping from copy operation
        # We need to find the mapping between old and new entities
        emap = self._build_entity_map(wrapped, new_wrapped)
        
        # Resolve and set ports on new monomer
        for name, port in self._ports.items():
            if port.target in emap:
                new_monomer.set_port(
                    name, 
                    emap[port.target],
                    role=port.role,
                    bond_kind=port.bond_kind,
                    compat=port.compat,
                    multiplicity=port.multiplicity,
                    priority=port.priority,
                )
        
        return new_monomer
    
    def _build_entity_map(self, old_assembly: Assembly, new_assembly: Assembly) -> dict[Entity, Entity]:
        """
        Build mapping between entities in old and new assemblies.
        
        This is a heuristic based on iteration order consistency.
        For more robust mapping, consider adding entity IDs.
        """
        emap: dict[Entity, Entity] = {}
        
        # Iterate through all entity types
        for entity_type in old_assembly.entities.classes():
            old_entities = old_assembly.entities.bucket(entity_type)
            new_entities = new_assembly.entities.bucket(entity_type)
            
            # Assume same iteration order after copy
            if len(old_entities) == len(new_entities):
                for old_ent, new_ent in zip(old_entities, new_entities):
                    emap[old_ent] = new_ent
        
        return emap
    
    def resolve_ports(self, emap: dict[Entity, Entity]) -> dict[str, Entity]:
        """
        Map ports through an entity mapping.
        
        Useful when merging monomers into polymers.
        
        Args:
            emap: Mapping from old entities to new entities
            
        Returns:
            Dictionary mapping port names to resolved entities
        """
        resolved: dict[str, Entity] = {}
        for name, port in self._ports.items():
            if port.target in emap:
                resolved[name] = emap[port.target]
        return resolved
    
    def _repr_info(self) -> str:
        """Provide monomer-specific info for tree representation."""
        port_names = list(self._ports.keys())
        if port_names:
            ports_str = ", ".join(port_names)
            return f"{len(self._ports)} port(s) [{ports_str}]"
        return "no ports"
