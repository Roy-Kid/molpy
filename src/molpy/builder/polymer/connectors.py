"""
Connector for polymer assembly.

The Connector decides which ports to connect between adjacent monomers
and executes the chemical reaction via a Reacter.
"""

from collections.abc import Mapping
from typing import Any

from molpy.core.atomistic import Atomistic
from molpy.core.entity import Entity
from molpy.reacter.base import Reacter
from molpy.typifier.atomistic import TypifierBase

from .errors import AmbiguousPortsError
from .port_utils import PortInfo
from .types import ConnectionMetadata, ConnectionResult


class ConnectorContext(dict[str, Any]):
    """Shared context passed to the connector during linear build.

    Keys:
    - step: int (current connection step index)
    - left_label: str (label of left monomer)
    - right_label: str (label of right monomer)
    - sequence: list[str] (full sequence being built)
    """

    pass


class Connector:
    """Select ports and execute reactions between adjacent monomers.

    Port selection strategy (applied in order):
    1. Explicit port_map lookup for (left_label, right_label)
    2. Role-based: left's role='right' port ↔ right's role='left' port
    3. Single-port: each side has exactly one unconsumed port
    4. Common name: both sides share a port name
    5. Raise AmbiguousPortsError

    The Reacter handles the actual bond formation and leaving group removal.
    """

    def __init__(
        self,
        reacter: Reacter,
        *,
        port_map: dict[tuple[str, str], tuple[str, str]] | None = None,
        overrides: dict[tuple[str, str], Reacter] | None = None,
    ):
        """Initialize the connector.

        Args:
            reacter: Default Reacter for bond formation.
            port_map: Optional explicit mapping from (left_label, right_label)
                to (left_port_name, right_port_name).
            overrides: Optional mapping from (left_label, right_label) to
                a specialized Reacter for that pair.
        """
        self.default = reacter
        self.port_map = port_map or {}
        self.overrides = overrides or {}
        self._history: list = []

    def get_reacter(self, left_type: str, right_type: str) -> Reacter:
        """Get the appropriate Reacter for a structure pair."""
        return self.overrides.get((left_type, right_type), self.default)

    def select_ports(
        self,
        left: Atomistic,
        right: Atomistic,
        left_ports: Mapping[str, list[PortInfo]],
        right_ports: Mapping[str, list[PortInfo]],
        ctx: ConnectorContext,
    ) -> tuple[str, int, str, int, None]:
        """Select which ports to connect.

        Args:
            left: Left Atomistic structure.
            right: Right Atomistic structure.
            left_ports: Available ports on left (name → list[PortInfo]).
            right_ports: Available ports on right (name → list[PortInfo]).
            ctx: Context with step info and labels.

        Returns:
            (left_port_name, left_idx, right_port_name, right_idx, None)

        Raises:
            AmbiguousPortsError: Cannot determine which ports to connect.
        """
        left_label = ctx.get("left_label", "")
        right_label = ctx.get("right_label", "")

        # 1. Explicit port_map
        key = (left_label, right_label)
        if key in self.port_map:
            port_L, port_R = self.port_map[key]
            if port_L not in left_ports:
                raise AmbiguousPortsError(
                    f"Port '{port_L}' not found on left ({left_label})"
                )
            if port_R not in right_ports:
                raise AmbiguousPortsError(
                    f"Port '{port_R}' not found on right ({right_label})"
                )
            return (port_L, 0, port_R, 0, None)

        # 2. Role-based: left[role=right] ↔ right[role=left]
        left_right = [
            (name, idx)
            for name, plist in left_ports.items()
            for idx, p in enumerate(plist)
            if p.role == "right"
        ]
        right_left = [
            (name, idx)
            for name, plist in right_ports.items()
            for idx, p in enumerate(plist)
            if p.role == "left"
        ]
        if len(left_right) == 1 and len(right_left) == 1:
            return (
                left_right[0][0],
                left_right[0][1],
                right_left[0][0],
                right_left[0][1],
                None,
            )

        # Also try left[role=left] ↔ right[role=right] (terminus)
        left_left = [
            (name, idx)
            for name, plist in left_ports.items()
            for idx, p in enumerate(plist)
            if p.role == "left"
        ]
        right_right = [
            (name, idx)
            for name, plist in right_ports.items()
            for idx, p in enumerate(plist)
            if p.role == "right"
        ]
        if len(left_left) == 1 and len(right_right) == 1:
            return (
                left_left[0][0],
                left_left[0][1],
                right_right[0][0],
                right_right[0][1],
                None,
            )

        # 3. Single-port on each side
        if len(left_ports) == 1 and len(right_ports) == 1:
            ln = next(iter(left_ports))
            rn = next(iter(right_ports))
            return (ln, 0, rn, 0, None)

        # 4. Common port name
        common = set(left_ports) & set(right_ports)
        if common:
            name = next(iter(common))
            return (name, 0, name, 0, None)

        # 5. Give up
        raise AmbiguousPortsError(
            f"Cannot auto-select ports between {left_label} and {right_label}: "
            f"left ports={list(left_ports.keys())}, right ports={list(right_ports.keys())}."
        )

    def connect(
        self,
        left: Atomistic,
        right: Atomistic,
        left_type: str,
        right_type: str,
        port_atom_L: Entity,
        port_atom_R: Entity,
        typifier: TypifierBase | None = None,
    ) -> ConnectionResult:
        """Execute the chemical reaction between two structures.

        Args:
            left: Left Atomistic structure.
            right: Right Atomistic structure.
            left_type: Label of left monomer.
            right_type: Label of right monomer.
            port_atom_L: Port atom on left.
            port_atom_R: Port atom on right.
            typifier: Optional typifier for incremental re-typing.

        Returns:
            ConnectionResult with product and metadata.
        """
        from molpy.reacter.base import ReactionResult

        reacter = self.get_reacter(left_type, right_type)

        product_set: ReactionResult = reacter.run(
            left,
            right,
            port_atom_L=port_atom_L,
            port_atom_R=port_atom_R,
            compute_topology=True,
            typifier=typifier,
        )

        self._history.append(product_set)

        metadata = ConnectionMetadata(
            port_L=port_atom_L.get("port", "unknown"),
            port_R=port_atom_R.get("port", "unknown"),
            reaction_name=reacter.name,
            formed_bonds=product_set.topology_changes.new_bonds,
            new_angles=product_set.topology_changes.new_angles,
            new_dihedrals=product_set.topology_changes.new_dihedrals,
            modified_atoms=product_set.topology_changes.modified_atoms,
            requires_retype=product_set.metadata.requires_retype,
            entity_maps=product_set.metadata.entity_maps,
        )

        return ConnectionResult(
            product=product_set.product_info.product, metadata=metadata
        )
