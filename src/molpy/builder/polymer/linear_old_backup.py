"""
Polymer assembly with support for chemical reactions.

This module provides the PolymerBuilder class with a linear() method that:
- Assembles monomers in sequence based on connector rules
- Supports both topology-only and chemical reaction connectors
- Deep-copies monomer topology (no geometry)
- Manages port consumption (used ports are removed)
- Provides detailed audit trail of connections
- Optionally applies geometric transforms via Placer
- Supports automatic retypification after reactions
"""

from typing import Literal, Any, Mapping, cast, TYPE_CHECKING
from ...core.wrappers.monomer import Monomer, Port
from ...core.wrappers.polymer import Polymer
from ...core.atomistic import Bond, Atom, Atomistic
from ...core.entity import Entity, Link
from ...core.ops.geometry import _vec_add, _vec_scale
from .connectors import Connector, ConnectorContext, BondKind, ReacterConnector
from ..errors import (
    SequenceError,
    AmbiguousPortsError,
    NoCompatiblePortsError,
    BondKindConflictError,
    PortReuseError,
)
from .geom_utils import Placer
from molpy.typifier.atomistic import OplsAtomisticTypifier


class PolymerBuilder:
    """
    Polymer constructed from linear sequence of monomers.
    
    Phase 1: Topology only (connection decisions), no geometry.
    """
    
    @classmethod
    def linear(
        cls,
        *,
        sequence: str,
        library: Mapping[str, Monomer],
        connector: Connector,
        consume_ports: bool = True,
        bond_kind: BondKind | None = None,
        stereo: str | None = None,
        details: dict | None = None,
        # Geometry layer
        placer: "Placer | None" = None,
        geom_ctx: dict[str, Any] | None = None,
        # NEW: Typification layer
        typifier: "OplsAtomisticTypifier | None" = None,
        auto_retypify: bool = True,
    ) -> Polymer:
        """
        Assemble linear polymer from sequence (topology + optional geometry).
        
        Args:
            sequence: String of monomer labels (e.g., "TABBAAT")
                     Termini must be included explicitly
            library: Mapping from label to Monomer template
            connector: Strategy for selecting which ports to connect
            consume_ports: If True, remove used ports after connection (default True)
            bond_kind: Global bond kind override (None = use port defaults)
            stereo: Global stereochemistry override (Phase 1: accepted but not used)
            details: Optional dict to fill with audit info
            placer: Optional geometry strategy for rigid transforms (Phase 1.5)
            geom_ctx: Optional dict for logging placement details
            
        Returns:
            Polymer as a Monomer (wrapped assembly with remaining ports)
            
        Raises:
            SequenceError: Invalid sequence (too short, unknown labels)
            AmbiguousPortsError: Cannot determine port selection
            NoCompatiblePortsError: No compatible ports found
            BondKindConflictError: Conflicting bond kind specifications
            PortReuseError: Attempt to use consumed port
            
        Example:
            >>> A = bigsmiles_to_monomer("[<]CCO[>]")
            >>> B = bigsmiles_to_monomer("[<]CC(C)O[>]")
            >>> T = bigsmiles_to_monomer("HO[*:t]")
            >>> # Topology only
            >>> poly = PolymerBuilder.linear(
            ...     sequence="TABBAT",
            ...     library={"A": A, "B": B, "T": T},
            ...     connector=AutoConnector(),
            ... )
            >>> # With geometry
            >>> from molpy.builder.polymer.geom_utils import DockPlacer
            >>> placer = DockPlacer(vdw_scale=0.80, vdw_scales_by_bond={"-": 0.80, "=": 0.72})
            >>> poly = PolymerBuilder.linear(
            ...     sequence="TABBAT",
            ...     library={"A": A, "B": B, "T": T},
            ...     connector=AutoConnector(),
            ...     placer=placer,
            ...     geom_ctx={"placement_log": []},
            ... )
        """
        # Validate sequence
        if len(sequence) < 2:
            raise SequenceError(f"Sequence must have at least 2 labels, got: {sequence!r}")
        
        # Validate all labels exist in library
        for label in sequence:
            if label not in library:
                raise SequenceError(f"Label {label!r} not found in library")
        
        # Initialize audit trail
        audit = []
        
        # Start with deep copy of first monomer
        first_copy = library[sequence[0]].copy()
        # Create polymer wrapper - Polymer() already handles port transfer
        current_polymer = Polymer(first_copy)
        
        # Transfer ports from the first monomer copy into the polymer
        # Important: must map port targets through the copy operation
        for pname, p in first_copy.ports.items():
            # Port target should already be in the wrapped assembly
            current_polymer.set_port(
                pname, p.target, 
                role=p.role, 
                bond_kind=p.bond_kind, 
                compat=p.compat, 
                priority=p.priority
            )
        
        # Build chain by iteratively adding monomers
        for step in range(len(sequence) - 1):
            left_label = sequence[step]
            right_label = sequence[step + 1]
            
            # Get fresh copy of right monomer
            right_monomer = library[right_label].copy()
            
            # Get unconsumed ports
            left_ports = current_polymer.ports  # ref not copy
            right_ports = right_monomer.ports  # ref not copy
            
            if not left_ports:
                raise NoCompatiblePortsError(
                    f"Step {step}: Left monomer ({left_label}) has no available ports"
                )
            if not right_ports:
                raise NoCompatiblePortsError(
                    f"Step {step}: Right monomer ({right_label}) has no available ports"
                )
            
            # Build connector context
            ctx = ConnectorContext(
                step=step,
                sequence=sequence,
                left_label=left_label,
                right_label=right_label,
                audit=audit,
            )
            
            # Check if connector is a ReacterConnector (chemical reaction mode)
            # TODO: too much redundancy args for select_ports
            #       remove connectorContext!
            use_reacter = isinstance(connector, ReacterConnector)
            
            if use_reacter:
                # === CHEMICAL REACTION MODE ===
                # Use ReacterConnector.connect() for full reaction
                
                # 1. Select ports using connector (reacter-guided)
                left_port_name, right_port_name, connector_bond_kind = connector.select_ports(
                    cast(Monomer, current_polymer),
                    right_monomer,
                    left_ports,
                    right_ports,
                    ctx,
                )
                
                # 2. Apply geometric transform if placer provided
                # TODO: make geom_ctx as args of placer!
                if placer is not None:
                    if geom_ctx is None:
                        geom_ctx = {}
                    
                    left_port = left_ports[left_port_name]
                    right_port = right_ports[right_port_name]
                    bk: BondKind = cast(BondKind, connector_bond_kind or "-")
                    
                    transform = placer.place_pair(
                        cast(Monomer, current_polymer),
                        right_monomer,
                        left_port,
                        right_port,
                        bond_kind=bk,
                        ctx=geom_ctx,
                    )
                    
                    R = transform["R_right"]
                    t = transform["t_right"]
                    _apply_transform(right_monomer, R, t)
                
                # 3. Execute chemical reaction
                product_assembly, metadata = connector.connect(
                    cast(Monomer, current_polymer),
                    right_monomer,
                    left_label,
                    right_label,
                    left_port_name,
                    right_port_name,
                )
                
                # 4. Auto-retypify if requested
                if auto_retypify and typifier and metadata.get('needs_retypification'):
                    _retypify_from_metadata(product_assembly, metadata, typifier)
                
                # 5. Update current_polymer to product
                current_polymer = Polymer(product_assembly)
                
                # 6. Transfer unconsumed ports from both monomers
                # CRITICAL: Since Reacter uses direct merge (no copy), port targets
                # are ALREADY pointing to the correct atoms in the product!
                # NO entity_maps needed!
                
                # 6a. Transfer LEFT monomer's unused ports (from current_polymer)
                for port_name, port in left_ports.items():
                    if port_name == left_port_name:
                        continue  # This port was consumed, skip it
                    # Direct transfer - port.target is already in product
                    current_polymer.set_port(
                        port_name,
                        port.target,
                        role=port.role,
                        bond_kind=port.bond_kind,
                        compat=port.compat,
                        priority=port.priority,
                    )
                
                # 6b. Transfer RIGHT monomer's unused ports
                for port_name, port in right_ports.items():
                    if port_name == right_port_name:
                        continue  # This port was consumed, skip it
                    # Direct transfer - port.target is already in product
                    current_polymer.set_port(
                        port_name,
                        port.target,
                        role=port.role,
                        bond_kind=port.bond_kind,
                        compat=port.compat,
                        priority=port.priority,
                    )
                
                # Record in audit
                audit.append({
                    "step": step,
                    "left_label": left_label,
                    "right_label": right_label,
                    "left_port": left_port_name,
                    "right_port": right_port_name,
                    "bond_kind": metadata.get('bond_kind', '-'),
                    "reaction": metadata.get('reaction_name', 'unknown'),
                })
                
            else:
                # === TOPOLOGY-ONLY MODE (Original logic) ===
                
                # Select ports using connector
                left_port_name, right_port_name, connector_bond_kind = connector.select_ports(
                    cast(Monomer, current_polymer),
                    right_monomer,
                    left_ports,
                    right_ports,
                    ctx,
                )
            
            # Validate port selection
            if left_port_name not in left_ports:
                raise NoCompatiblePortsError(
                    f"Step {step}: Connector selected unavailable left port {left_port_name!r}"
                )
            if right_port_name not in right_ports:
                raise NoCompatiblePortsError(
                    f"Step {step}: Connector selected unavailable right port {right_port_name!r}"
                )
            
            left_port = left_ports[left_port_name]
            right_port = right_ports[right_port_name]
            
            # Check port compatibility (if compat is defined)
            if not _ports_compatible(left_port, right_port):
                raise NoCompatiblePortsError(
                    f"Step {step}: Ports {left_port_name!r} and {right_port_name!r} are not compatible"
                )
            
            # Resolve bond kind with priority: global > connector > port defaults
            resolved_bond_kind = _resolve_bond_kind(
                bond_kind, 
                connector_bond_kind,
                left_port.bond_kind,
                right_port.bond_kind,
            )
            
            # Apply geometric transform if placer is provided
            if placer is not None:
                # Create geometry context if not provided
                if geom_ctx is None:
                    geom_ctx = {}
                
                # Convert resolved_bond_kind to BondKind literal
                bk: BondKind = cast(BondKind, resolved_bond_kind or "-")
                
                # Compute rigid transform
                transform = placer.place_pair(
                    cast(Monomer, current_polymer),
                    right_monomer,
                    left_port,
                    right_port,
                    bond_kind=bk,
                    ctx=geom_ctx,
                )
                
                # Apply transform to right monomer atoms BEFORE merge
                R = transform["R_right"]
                t = transform["t_right"]
                _apply_transform(right_monomer, R, t)
            
            # Merge right monomer into polymer (deep copy into polymer assembly)
            wrapped_polymer = current_polymer
            wrapped_right = right_monomer
            emap = wrapped_polymer.merge(wrapped_right)
            
            # Create bond between the two port targets
            left_target = left_port.target
            right_target_old = right_port.target
            # Map right target through entity map produced by merge
            right_target = emap[right_target_old]
            
            # Add bond to polymer
            bond_attrs = {}
            if resolved_bond_kind:
                bond_attrs['kind'] = resolved_bond_kind
            if stereo:
                bond_attrs['stereo'] = stereo
            
            # Use Bond if targets are Atoms, otherwise use generic Link
            if isinstance(left_target, Atom) and isinstance(right_target, Atom):
                bond = Bond(left_target, right_target, **bond_attrs)
            else:
                bond = Link([left_target, right_target], **bond_attrs)
            wrapped_polymer.links.add(bond)
            
            # Consume ports (remove used ports from polymer)
            if consume_ports:
                # Remove the left port that was just used
                if left_port_name in current_polymer.ports:
                    del current_polymer.ports[left_port_name]

                # Transfer right_monomer's unused ports into polymer with remapped targets
                # Skip the right port that was just used
                for r_name, r_port in right_monomer.ports.items():
                    if r_name == right_port_name:
                        continue  # This port was consumed, skip it
                    if r_port.target in emap:
                        new_target = emap[r_port.target]
                        # Add unused port to polymer
                        current_polymer.set_port(
                            r_name,
                            new_target,
                            role=r_port.role,
                            bond_kind=r_port.bond_kind,
                            compat=r_port.compat,
                            priority=r_port.priority,
                        )
            
            # Record in audit trail
            audit.append({
                "step": step,
                "left_label": left_label,
                "right_label": right_label,
                "left_port": left_port_name,
                "right_port": right_port_name,
                "bond_kind": resolved_bond_kind or "-",
            })
        
        # Fill details if provided
        if details is not None:
            details["connections"] = audit
            details["leftover_ports"] = list(current_polymer.ports.keys())

        return current_polymer


def _ports_compatible(p1, p2) -> bool:
    """
    Check if two ports are compatible based on compat specs.
    
    Rules:
    - If either has compat='*', they're compatible
    - Otherwise, p2.name must be in p1.compat AND p1.name must be in p2.compat
    """
    c1 = p1.compat
    c2 = p2.compat
    
    # If either compat is None, assume compatible
    if c1 is None or c2 is None:
        return True
    
    if c1 == '*' or c2 == '*':
        return True
    
    # Both have specific compat sets
    if isinstance(c1, set) and isinstance(c2, set):
        return p2.name in c1 and p1.name in c2
    
    # Mixed types - treat strings as single-element sets
    c1_set = {c1} if isinstance(c1, str) else c1
    c2_set = {c2} if isinstance(c2, str) else c2
    
    return p2.name in c1_set and p1.name in c2_set


def _resolve_bond_kind(
    global_kind: BondKind | None,
    connector_kind: BondKind | None,
    left_kind: str | None,
    right_kind: str | None,
) -> str:
    """
    Resolve bond kind with priority: global > connector > port defaults.
    
    Raises:
        BondKindConflictError: If port defaults conflict
    """
    # Priority 1: Global override
    if global_kind is not None:
        return global_kind
    
    # Priority 2: Connector override
    if connector_kind is not None:
        return connector_kind
    
    # Priority 3: Port defaults
    if left_kind and right_kind:
        if left_kind != right_kind:
            raise BondKindConflictError(
                f"Conflicting bond kinds: left={left_kind!r}, right={right_kind!r}"
            )
        return left_kind
    elif left_kind:
        return left_kind
    elif right_kind:
        return right_kind
    else:
        # Default to single bond
        return "-"


def _apply_transform(assembly, R: list[list[float]], t: list[float]) -> None:
    """
    Apply rigid transform (rotation + translation) to all atoms in assembly.
    
    Args:
        assembly: Assembly or Atomistic object with atoms
        R: 3x3 rotation matrix
        t: Translation vector [x, y, z]
    """
    # mat_vec_mult is local helper, vec_add from core/ops/geometry
    def _mat_vec_mult(M: list[list[float]], v: list[float]) -> list[float]:
        return [
            M[0][0] * v[0] + M[0][1] * v[1] + M[0][2] * v[2],
            M[1][0] * v[0] + M[1][1] * v[1] + M[1][2] * v[2],
            M[2][0] * v[0] + M[2][1] * v[1] + M[2][2] * v[2],
        ]
    
    # Get all atoms
    atoms = assembly.atoms if hasattr(assembly, 'atoms') else assembly.entities.bucket(Atom)
    
    for atom in atoms:
        # Get current position
        pos = atom.data.get('pos') or atom.data.get('xyz')
        if pos is None:
            continue  # Skip atoms without positions
        
        # Apply rotation then translation: p' = R*p + t
        pos_rotated = _mat_vec_mult(R, pos)
        pos_new = _vec_add(pos_rotated, t)
        
        # Update position
        atom.data['pos'] = pos_new
        atom.data['xyz'] = pos_new  # Also update xyz for compatibility


def _retypify_from_metadata(
    assembly: Atomistic,
    metadata: dict[str, Any],
    typifier: "OplsAtomTypifier",
) -> None:
    """
    Retypify topology after chemical reaction based on metadata.
    
    Args:
        assembly: The atomistic assembly to retypify
        metadata: Metadata from ReacterConnector containing:
            - modified_atoms: set of atoms whose types may have changed
            - new_bonds: list of newly created bonds
            - new_angles: list of newly created angles
            - new_dihedrals: list of newly created dihedrals
        typifier: OPLS atom typifier
    """
    # 1. Retypify modified atoms
    modified_atoms = metadata.get('modified_atoms', set())
    for atom in modified_atoms:
        try:
            typifier.typify(atom)
        except Exception:
            pass  # Skip atoms that can't be typed
    
    # 2. Typify new bonds
    new_bonds = metadata.get('new_bonds', [])
    for bond in new_bonds:
        try:
            typifier.typify(bond)
        except Exception:
            pass
    
    # 3. Typify new angles
    new_angles = metadata.get('new_angles', [])
    for angle in new_angles:
        try:
            typifier.typify(angle)
        except Exception:
            pass
    
    # 4. Typify new dihedrals
    new_dihedrals = metadata.get('new_dihedrals', [])
    for dihedral in new_dihedrals:
        try:
            typifier.typify(dihedral)
        except Exception:
            pass
