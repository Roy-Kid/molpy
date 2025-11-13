"""
Geometry utilities for polymer assembly.

Provides Rodrigues rotation, VDW-based docking, and placer strategies.
Directly uses core.ops.geometry functions - no wrapper layer.
"""

import math
from typing import Protocol, Literal, Any, TYPE_CHECKING

# Direct import from core ops - no wrapper
from ...core.ops.geometry import (
    _dot,
    _cross,
    _norm,
    _unit,
    _vec_sub,
    _vec_add,
    _vec_scale,
)

from ..errors import OrientationUnavailableError, PositionMissingError

if TYPE_CHECKING:
    from ...core.entity import Entity
    from ...core.wrappers.monomer import Monomer, Port
    from ...core.atomistic import Atom


BondKind = Literal["-", "=", "#", ":"]


# ============================================================================
# Core Math Functions
# ============================================================================

def _mat_vec_mult(R: list[list[float]], v: list[float]) -> list[float]:
    """Multiply 3x3 matrix by 3D vector: R @ v."""
    return [
        R[0][0] * v[0] + R[0][1] * v[1] + R[0][2] * v[2],
        R[1][0] * v[0] + R[1][1] * v[1] + R[1][2] * v[2],
        R[2][0] * v[0] + R[2][1] * v[1] + R[2][2] * v[2],
    ]


def rodrigues(a: list[float], b: list[float]) -> list[list[float]]:
    """
    Compute rotation matrix that maps unit vector 'a' to unit vector 'b'.
    
    Uses Rodrigues' rotation formula. Handles the antiparallel case robustly.
    
    Args:
        a: Source unit vector (should be normalized)
        b: Target unit vector (should be normalized)
        
    Returns:
        3x3 rotation matrix as list[list[float]]
    """
    # Normalize inputs
    n_a = _norm(a)
    n_b = _norm(b)
    if n_a < 1e-10 or n_b < 1e-10:
        raise ValueError(f"Cannot compute rotation with zero vector")
    
    a_norm = _unit(a)
    b_norm = _unit(b)
    
    # Compute dot product
    c = _dot(a_norm, b_norm)
    
    # Check if already aligned
    if c > 1.0 - 1e-10:
        return [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    
    # Check if antiparallel (180° rotation)
    if c < -1.0 + 1e-10:
        # Find an orthogonal axis
        if abs(a_norm[0]) < abs(a_norm[1]) and abs(a_norm[0]) < abs(a_norm[2]):
            orth = [1.0, 0.0, 0.0]
        elif abs(a_norm[1]) < abs(a_norm[2]):
            orth = [0.0, 1.0, 0.0]
        else:
            orth = [0.0, 0.0, 1.0]
        
        k = _unit(_cross(a_norm, orth))
        
        # 180° rotation around k
        return [
            [2 * k[0] * k[0] - 1, 2 * k[0] * k[1], 2 * k[0] * k[2]],
            [2 * k[1] * k[0], 2 * k[1] * k[1] - 1, 2 * k[1] * k[2]],
            [2 * k[2] * k[0], 2 * k[2] * k[1], 2 * k[2] * k[2] - 1],
        ]
    
    # General case: Rodrigues formula
    k = _unit(_cross(a_norm, b_norm))
    s = math.sqrt(1 - c * c)
    
    K = [[0.0, -k[2], k[1]], [k[2], 0.0, -k[0]], [-k[1], k[0], 0.0]]
    K2 = [
        [k[0] * k[0] - 1, k[0] * k[1], k[0] * k[2]],
        [k[1] * k[0], k[1] * k[1] - 1, k[1] * k[2]],
        [k[2] * k[0], k[2] * k[1], k[2] * k[2] - 1],
    ]
    
    R = [
        [1.0 + s * K[0][0] + (1 - c) * K2[0][0],
         s * K[0][1] + (1 - c) * K2[0][1],
         s * K[0][2] + (1 - c) * K2[0][2]],
        [s * K[1][0] + (1 - c) * K2[1][0],
         1.0 + s * K[1][1] + (1 - c) * K2[1][1],
         s * K[1][2] + (1 - c) * K2[1][2]],
        [s * K[2][0] + (1 - c) * K2[2][0],
         s * K[2][1] + (1 - c) * K2[2][1],
         1.0 + s * K[2][2] + (1 - c) * K2[2][2]],
    ]
    
    return R


def infer_orientation(port: "Port") -> list[float]:
    """
    Infer orientation vector for a port.
    
    Rules:
    1. If port has explicit orientation, use it
    2. If target has one heavy neighbor, use direction target→neighbor
    3. Else average normalized vectors from target to all heavy neighbors
    4. Fallback: [1.0, 0.0, 0.0]
    
    Args:
        port: Port object with target entity
        
    Returns:
        Unit orientation vector [x, y, z]
    """
    # Check for explicit orientation in port data
    if hasattr(port, 'data') and 'orientation' in port.data:
        orient = port.data['orientation']
        if orient is not None and len(orient) == 3:
            return _unit(orient)
    
    # Fallback: arbitrary direction (proper connectivity inference needs Struct API)
    return [1.0, 0.0, 0.0]


def get_vdw_radius(entity: "Entity") -> float:
    """
    Get van der Waals radius for an entity.
    
    Priority:
    1. Explicit 'vdw' in entity.data
    2. Element class lookup by symbol
    3. Default to carbon (1.70 Å)
    
    Args:
        entity: Entity with element information
        
    Returns:
        VDW radius in Angstroms
    """
    if 'vdw' in entity.data:
        return float(entity.data['vdw'])
    
    try:
        from ...core.element import Element
        symbol = entity.data.get('symbol') or entity.data.get('type', 'C')
        elem = Element(symbol)
        return elem.vdw
    except (ImportError, KeyError, AttributeError):
        return 1.70


# ============================================================================
# Placer Protocol and Implementations
# ============================================================================

class Placer(Protocol):
    """
    Strategy interface for computing rigid transforms to dock monomers.
    Topology is already decided by Connector.
    """
    
    def place_pair(
        self,
        left: "Monomer",
        right: "Monomer",
        left_port: "Port",
        right_port: "Port",
        *,
        bond_kind: BondKind,
        ctx: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Compute rigid transform to dock right monomer onto left.
        
        Returns:
            Dict with "R_right" (3x3 matrix), "t_right" (3D vector), "distance" (float)
        """
        ...


class NoOpPlacer:
    """No-op placer returning identity transform."""
    
    def place_pair(
        self,
        left: "Monomer",
        right: "Monomer",
        left_port: "Port",
        right_port: "Port",
        *,
        bond_kind: BondKind,
        ctx: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        return {
            "R_right": [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
            "t_right": [0.0, 0.0, 0.0],
            "distance": 0.0,
        }


class DockPlacer:
    """
    VDW-based docking placer using Rodrigues rotation.
    
    Automatically infers port orientations using two strategies:
    1. Strategy 1 (from_coords_nearest_heavy): Orient from port anchor toward nearest heavy atom
    2. Strategy 2 (from_role_main_chain): Use main-chain direction from left/right role ports
    
    Args:
        vdw_scale: Global VDW scale factor (default: 0.80)
        vdw_scales_by_bond: Per-bond-kind scale overrides (e.g., {"-": 0.80, "=": 0.72})
        allow_role_direction: Enable strategy 2 (main-chain from roles) (default: True)
    
    Example:
        >>> placer = DockPlacer(
        ...     vdw_scale=0.80,
        ...     vdw_scales_by_bond={"-": 0.80, "=": 0.72, "#": 0.66},
        ...     allow_role_direction=True
        ... )
    """
    
    def __init__(
        self,
        vdw_scale: float = 0.80,
        vdw_scales_by_bond: dict[BondKind, float] | None = None,
        allow_role_direction: bool = True,
    ):
        self.vdw_scale = vdw_scale
        self.vdw_scales_by_bond = vdw_scales_by_bond or {}
        self.allow_role_direction = allow_role_direction
    
    def place_pair(
        self,
        left: "Monomer",
        right: "Monomer",
        left_port: "Port",
        right_port: "Port",
        *,
        bond_kind: BondKind,
        ctx: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Compute rigid transform to dock right monomer onto left.
        
        Args:
            left: Left monomer (already positioned)
            right: Right monomer (to be transformed)
            left_port: Port on left monomer
            right_port: Port on right monomer
            bond_kind: Bond type ("-", "=", "#", ":")
            ctx: Optional context dict for logging
            
        Returns:
            Dict with "R_right" (3x3 rotation matrix), "t_right" (translation), "distance"
            
        Raises:
            PositionMissingError: If port anchor has no coordinates
            OrientationUnavailableError: If cannot infer orientation
        """
        if ctx is None:
            ctx = {}
        
        # Get port positions
        pL = self._get_position(left_port.target)
        pR = self._get_position(right_port.target)
        
        # Infer orientations using strategies
        uL = self._infer_orientation(left_port, left)
        uR = self._infer_orientation(right_port, right)
        
        # Compute docking distance from VDW radii
        rvdw_L = get_vdw_radius(left_port.target)
        rvdw_R = get_vdw_radius(right_port.target)
        
        # Get scale factor (per-bond override or global)
        scale = self.vdw_scales_by_bond.get(bond_kind, self.vdw_scale)
        
        d = scale * (rvdw_L + rvdw_R)
        
        # Rotation: align uR to -uL
        neg_uL = _vec_scale(uL, -1.0)
        R = rodrigues(uR, neg_uL)
        
        # Target position for right port: pL - d*uL
        target_pos = _vec_sub(pL, _vec_scale(uL, d))
        
        # Current position after rotation
        pR_rotated = _mat_vec_mult(R, pR)
        
        # Translation
        t = _vec_sub(target_pos, pR_rotated)
        
        # Log if requested
        if "placement_log" in ctx:
            ctx["placement_log"].append({
                "left_port_pos": pL,
                "right_port_pos": pR,
                "left_orientation": uL,
                "right_orientation": uR,
                "rvdw_L": rvdw_L,
                "rvdw_R": rvdw_R,
                "scale": scale,
                "distance": d,
                "rotation": R,
                "translation": t,
            })
        
        return {"R_right": R, "t_right": t, "distance": d}
    
    def _get_position(self, entity: "Entity") -> list[float]:
        """
        Get 3D position from entity.
        
        Raises:
            PositionMissingError: If no position data found
        """
        pos = entity.data.get('pos') or entity.data.get('xyz')
        if pos is None:
            raise PositionMissingError(f"Entity has no position data: {entity}")
        if len(pos) != 3:
            raise PositionMissingError(f"Position must be 3D, got: {pos}")
        return list(pos)
    
    def _infer_orientation(self, port: "Port", monomer: "Monomer") -> list[float]:
        """
        Infer orientation for a port using fixed strategy order:
        1. Strategy 1: Nearest heavy neighbor
        2. Strategy 2: Role-based main-chain direction (if allow_role_direction=True)
        
        Args:
            port: Port to infer orientation for
            monomer: Monomer containing the port
            
        Returns:
            Unit orientation vector [x, y, z]
            
        Raises:
            OrientationUnavailableError: If neither strategy produces valid orientation
        """
        # Check for explicit orientation in port or target
        explicit_orient = self._get_explicit_orientation(port)
        if explicit_orient is not None:
            return explicit_orient
        
        # Strategy 1: Nearest heavy neighbor
        try:
            orient = self._orient_from_nearest_heavy(port, monomer)
            if orient is not None:
                return orient
        except Exception:
            pass
        
        # Strategy 2: Role-based main-chain direction
        if self.allow_role_direction:
            try:
                orient = self._orient_from_role_main_chain(port, monomer)
                if orient is not None:
                    return orient
            except Exception:
                pass
        
        # Both strategies failed
        raise OrientationUnavailableError(
            f"Cannot infer orientation for port {port.name!r}: "
            f"no heavy neighbors and no role-based main chain available"
        )
    
    def _get_explicit_orientation(self, port: "Port") -> list[float] | None:
        """Check for explicit orientation in port or target data."""
        # Check port data
        if hasattr(port, 'data') and 'orientation' in port.data:
            orient = port.data['orientation']
            if orient is not None and len(orient) == 3:
                n = _norm(orient)
                if n > 1e-10:
                    return _unit(orient)
        
        # Check target entity data
        if hasattr(port.target, 'data') and 'orientation' in port.target.data:
            orient = port.target.data['orientation']
            if orient is not None and len(orient) == 3:
                n = _norm(orient)
                if n > 1e-10:
                    return _unit(orient)
        
        return None
    
    def _orient_from_nearest_heavy(self, port: "Port", monomer: "Monomer") -> list[float] | None:
        """
        Strategy 1: Orient from port anchor toward nearest heavy atom.
        
        Args:
            port: Port with target anchor atom
            monomer: Monomer to search for neighbors
            
        Returns:
            Unit orientation vector or None if no heavy neighbors
        """
        from ...core.atomistic import Atom, Bond
        
        anchor = port.target
        if not isinstance(anchor, Atom):
            return None
        
        # Get anchor position
        try:
            anchor_pos = self._get_position(anchor)
        except PositionMissingError:
            return None
        
        # Find all bonded neighbors
        assembly = monomer.unwrap()
        neighbors = []
        
        for bond in assembly.bonds:
            if bond.itom == anchor:
                neighbors.append(bond.jtom)
            elif bond.jtom == anchor:
                neighbors.append(bond.itom)
        
        if not neighbors:
            return None
        
        # Filter to heavy atoms (non-hydrogen)
        heavy_neighbors = [
            n for n in neighbors
            if n.data.get('symbol', 'C') != 'H'
        ]
        
        if not heavy_neighbors:
            return None
        
        # Find nearest heavy neighbor
        min_dist = float('inf')
        nearest = None
        
        for neighbor in heavy_neighbors:
            try:
                n_pos = self._get_position(neighbor)
                dist = _norm(_vec_sub(n_pos, anchor_pos))
                if dist < min_dist:
                    min_dist = dist
                    nearest = neighbor
            except PositionMissingError:
                continue
        
        if nearest is None:
            return None
        
        # Direction: anchor → nearest heavy
        nearest_pos = self._get_position(nearest)
        direction = _vec_sub(nearest_pos, anchor_pos)
        n = _norm(direction)
        
        if n < 1e-10:
            return None
        
        return _unit(direction)
    
    def _orient_from_role_main_chain(self, port: "Port", monomer: "Monomer") -> list[float] | None:
        """
        Strategy 2: Use main-chain direction from left/right role ports.
        
        Main chain vector g = pos(right) - pos(left)
        - Ports with role='right' use +g
        - Ports with role='left' use -g
        
        Args:
            port: Port to orient
            monomer: Monomer containing ports with roles
            
        Returns:
            Unit orientation vector or None if roles not available
        """
        # Check if port has a role
        port_role = port.role
        if port_role not in ('left', 'right'):
            return None
        
        # Find left and right role ports
        left_port = None
        right_port = None
        
        for name, p in monomer.ports.items():
            if p.role == 'left' and left_port is None:
                left_port = p
            elif p.role == 'right' and right_port is None:
                right_port = p
        
        if left_port is None or right_port is None:
            return None
        
        # Get positions
        try:
            left_pos = self._get_position(left_port.target)
            right_pos = self._get_position(right_port.target)
        except PositionMissingError:
            return None
        
        # Compute main-chain direction: right - left
        chain_dir = _vec_sub(right_pos, left_pos)
        n = _norm(chain_dir)
        
        if n < 1e-10:
            return None
        
        chain_unit = _unit(chain_dir)
        
        # Apply based on port role
        if port_role == 'right':
            return chain_unit  # +g
        else:  # port_role == 'left'
            return _vec_scale(chain_unit, -1.0)  # -g
