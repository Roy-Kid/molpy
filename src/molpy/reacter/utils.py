"""
Utility functions for assembly manipulation in reactions.

This module provides helper functions for finding neighbors
and other common operations needed for reactions.
"""

from copy import deepcopy
from molpy.core.entity import Entity
from molpy.core.atomistic import Atomistic, Atom, Bond


def find_neighbors(
    assembly: Atomistic,
    atom: Entity,
    *,
    element: str | None = None,
) -> list[Entity]:
    """
    Find neighboring atoms of a given atom.
    
    Args:
        assembly: Atomistic assembly containing the atom
        atom: Atom entity to find neighbors of
        element: Optional element symbol to filter by (e.g., 'H', 'C')
    
    Returns:
        List of neighboring atom entities
    
    Example:
        >>> h_neighbors = find_neighbors(asm, carbon_atom, element='H')
        >>> all_neighbors = find_neighbors(asm, carbon_atom)
    """
    neighbors: list[Entity] = []
    
    # Look through all bonds
    for bond in assembly.links.bucket(Bond):
        # Use identity check (is) not equality check (==)
        if any(ep is atom for ep in bond.endpoints):
            # Found a bond involving this atom
            for endpoint in bond.endpoints:
                if endpoint is not atom:
                    # Filter by element if specified
                    if element is None or endpoint.get('symbol') == element:
                        neighbors.append(endpoint)
    
    return neighbors


def get_bond_between(
    assembly: Atomistic,
    i: Entity,
    j: Entity,
) -> Bond | None:
    """
    Find existing bond between two atoms.
    
    Args:
        assembly: Atomistic assembly containing the atoms
        i: First atom entity
        j: Second atom entity
    
    Returns:
        Bond entity if found, None otherwise
    
    Example:
        >>> bond = get_bond_between(asm, atom1, atom2)
        >>> if bond:
        ...     print(f"Bond order: {bond.get('order', 1)}")
    """
    for bond in assembly.links.bucket(Bond):
        endpoints = bond.endpoints
        # Use identity check (is) not equality check (==)
        if any(ep is i for ep in endpoints) and any(ep is j for ep in endpoints):
            return bond
    return None


def count_bonds(assembly: Atomistic, atom: Entity) -> int:
    """
    Count the number of bonds connected to an atom.
    
    Args:
        assembly: Atomistic assembly containing the atom
        atom: Atom entity to count bonds for
    
    Returns:
        Number of bonds
    
    Example:
        >>> valence = count_bonds(asm, carbon_atom)
        >>> print(f"Carbon has {valence} bonds")
    """
    count = 0
    for bond in assembly.links.bucket(Bond):
        # Use identity check (is) not equality check (==)
        if any(ep is atom for ep in bond.endpoints):
            count += 1
    return count


def remove_dummy_atoms(assembly: Atomistic) -> list[Entity]:
    """
    Remove all dummy atoms (element '*' or symbol '*') from assembly.
    
    Args:
        assembly: Atomistic assembly to clean
    
    Returns:
        List of removed dummy atoms
    
    Example:
        >>> removed = remove_dummy_atoms(merged_asm)
        >>> print(f"Removed {len(removed)} dummy atoms")
    """
    dummy_atoms: list[Entity] = []
    
    for atom in assembly.entities.bucket(Atom):
        symbol = atom.get('symbol', '')
        element = atom.get('element', '')
        if symbol == '*' or element == '*':
            dummy_atoms.append(atom)
    
    if dummy_atoms:
        assembly.remove_entity(*dummy_atoms, drop_incident_links=True)
    
    return dummy_atoms
