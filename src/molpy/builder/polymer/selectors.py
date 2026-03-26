"""Leaving group selector functions for polymer building.

Most selectors have been consolidated into ``molpy.reacter.selectors``.
This module retains builder-specific utilities and re-exports for
backward compatibility.
"""

from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Callable

if TYPE_CHECKING:
    from molpy.core.atomistic import Atom, Atomistic

# Type alias — also available from molpy.reacter.selectors
LeavingGroupSelector = Callable[["Atomistic", "Atom"], list["Atom"]]


# ============================================================================
# Port marker processing utilities  (builder-specific, stays here)
# ============================================================================


def process_port_markers(monomer: "Atomistic") -> "Atomistic":
    """Process port markers in monomer structure.

    Converts port marker notation [>] and [<] to port attributes on real atoms:
    1. Find wildcard atoms (*) with port attribute
    2. Find atoms bonded to these wildcards
    3. Transfer port attribute to bonded atoms
    4. Remove wildcard atoms

    Args:
        monomer: Atomistic structure with port markers from parse_smiles("[>]...[<]")

    Returns:
        Atomistic structure with port attributes on real atoms, wildcards removed

    Example:
        >>> from molpy.parser import parse_molecule
        >>> monomer = parse_molecule("[>]CCOCC[<]")
        >>> # parse_molecule handles IR conversion internally
        >>> monomer = process_port_markers(monomer)
        >>> # After processing: C - C - O - C - C (port on first and last C)
    """
    from molpy.core.atomistic import Atomistic as _Atomistic

    # Find wildcard atoms with port attribute
    port_markers: list[tuple["Atom", str]] = []
    for atom in monomer.atoms:
        if atom.get("symbol") == "*" and atom.get("port"):
            port_markers.append((atom, atom.get("port")))

    if not port_markers:
        return monomer

    # Find atoms connected to each port marker
    port_transfers: list[tuple["Atom", "Atom", str]] = []
    for wildcard, port_value in port_markers:
        connected_atom = None
        for bond in monomer.bonds:
            if bond.itom == wildcard:
                connected_atom = bond.jtom
                break
            elif bond.jtom == wildcard:
                connected_atom = bond.itom
                break

        if connected_atom:
            port_transfers.append((wildcard, connected_atom, port_value))

    # Create new Atomistic without wildcard atoms
    new_struct = _Atomistic()

    atom_mapping = {}
    for atom in monomer.atoms:
        if atom.get("symbol") != "*" or not atom.get("port"):
            atom_data = dict(atom.items())
            new_atom = new_struct.def_atom(**atom_data)
            atom_mapping[atom] = new_atom

    # Transfer port attributes to connected atoms
    for wildcard, real_atom, port_value in port_transfers:
        if real_atom in atom_mapping:
            new_atom = atom_mapping[real_atom]
            new_atom["port"] = port_value

    # Add bonds (excluding bonds to wildcards)
    for bond in monomer.bonds:
        if bond.itom in atom_mapping and bond.jtom in atom_mapping:
            new_itom = atom_mapping[bond.itom]
            new_jtom = atom_mapping[bond.jtom]
            new_struct.def_bond(
                new_itom,
                new_jtom,
                order=bond.get("order"),
                kind=bond.get("kind"),
            )

    return new_struct


# ============================================================================
# Deprecated re-exports — import from molpy.reacter.selectors instead
# ============================================================================


def _deprecated_reexport(name: str):
    warnings.warn(
        f"Importing {name} from molpy.builder.polymer.selectors is deprecated. "
        f"Use molpy.reacter.selectors instead.",
        DeprecationWarning,
        stacklevel=2,
    )


def select_bonded_hydrogens(
    monomer: "Atomistic", connection_atom: "Atom"
) -> list["Atom"]:
    """Deprecated: use ``select_hydrogens(n=None)`` from ``molpy.reacter.selectors``."""
    _deprecated_reexport("select_bonded_hydrogens")
    from molpy.reacter.selectors import select_hydrogens

    return select_hydrogens(n=None)(monomer, connection_atom)


def select_by_smarts(
    smarts_pattern: str,
    max_distance: int | None = None,
) -> LeavingGroupSelector:
    """Create a selector that matches atoms by SMARTS pattern (stub)."""
    _deprecated_reexport("select_by_smarts")

    def selector(monomer: "Atomistic", connection_atom: "Atom") -> list["Atom"]:
        return []

    return selector


def get_port_atoms(monomer: "Atomistic") -> dict[str, "Atom"]:
    """Deprecated: use ``get_ports`` from ``molpy.reacter.selectors``."""
    _deprecated_reexport("get_port_atoms")
    from molpy.reacter.selectors import get_ports

    return get_ports(monomer)
