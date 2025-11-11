"""
Tests for Placer strategies (NoOpPlacer, DockPlacer).
"""

import pytest
import math
from molpy.builder.polymer import NoOpPlacer, DockPlacer
from molpy.core.ops.geometry import _dot, _norm, _vec_sub, _unit
from molpy.core.atomistic import Atomistic, Atom
from molpy.core.wrappers.monomer import Monomer


def create_simple_monomer(name: str, port_pos: list[float], port_dir: list[float] | None = None):
    """
    Create a simple monomer with one atom at origin and one port atom.
    
    Args:
        name: Monomer label
        port_pos: Position of port atom
        port_dir: Optional explicit orientation
    """
    atomistic = Atomistic()
    
    # Central atom at origin
    center = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    
    # Port atom
    port_atom = atomistic.add_atom(symbol="C", pos=port_pos)
    
    # Add bond
    atomistic.add_bond(center, port_atom, order=1)
    
    # Create monomer
    monomer = Monomer(atomistic)
    
    # Add port
    monomer.set_port(
        f"port_{name}",
        port_atom,
        role="left" if name.startswith("L") else "right",
    )
    
    # Set explicit orientation if provided
    if port_dir is not None:
        port_atom.data['orientation'] = port_dir
    
    return monomer


class TestNoOpPlacer:
    """Test NoOp placer that returns identity transform."""
    
    def test_returns_identity(self):
        """Test that NoOpPlacer returns identity transform."""
        placer = NoOpPlacer()
        
        # Create dummy monomers
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0])
        right = create_simple_monomer("R1", [1.0, 0.0, 0.0])
        
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        assert left_port is not None
        assert right_port is not None
        
        ctx = GeometryContext()
        result = placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        
        # Should be identity
        assert result["R_right"] == [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
        assert result["t_right"] == [0.0, 0.0, 0.0]
        assert result["distance"] == 0.0


class TestDockPlacer:
    """Test VDW-based docking placer."""
    
    def test_simple_alignment_x_axis(self):
        """Test alignment of two monomers along x-axis."""
        placer = DockPlacer()
        
        # Left monomer: port at +x pointing in +x direction
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0], port_dir=[1.0, 0.0, 0.0])
        
        # Right monomer: port at -x pointing in -x direction (will flip to point in +x)
        right = create_simple_monomer("R1", [-1.0, 0.0, 0.0], port_dir=[-1.0, 0.0, 0.0])
        
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        
        ctx = GeometryContext({"vdw_scale": 0.80})
        result = placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        
        R = result["R_right"]
        t = result["t_right"]
        d = result["distance"]
        
        # Check that rotation aligns right orientation (-x) to -left orientation (-x)
        # i.e., align [-1, 0, 0] to [-1, 0, 0], which should be identity
        right_dir = [-1.0, 0.0, 0.0]
        left_dir = [1.0, 0.0, 0.0]
        neg_left_dir = [-1.0, 0.0, 0.0]
        
        # Apply rotation to right direction
        rotated_dir = [
            R[0][0] * right_dir[0] + R[0][1] * right_dir[1] + R[0][2] * right_dir[2],
            R[1][0] * right_dir[0] + R[1][1] * right_dir[1] + R[1][2] * right_dir[2],
            R[2][0] * right_dir[0] + R[2][1] * right_dir[1] + R[2][2] * right_dir[2],
        ]
        
        # Should align with -left_dir
        assert _dot(rotated_dir, neg_left_dir) == pytest.approx(1.0, abs=1e-6)
        
        # Distance should be VDW-based (2 * 1.70 * 0.80 = 2.72)
        assert d == pytest.approx(0.80 * 2 * 1.70, abs=0.01)
    
    def test_vdw_distance_carbon_carbon(self):
        """Test that distance is correctly computed from VDW radii."""
        placer = DockPlacer()
        
        # Create monomers with explicit VDW
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0])
        right = create_simple_monomer("R1", [-1.0, 0.0, 0.0])
        
        # Set explicit VDW radii (carbon)
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        left_port.target.data['vdw'] = 1.70
        right_port.target.data['vdw'] = 1.70
        
        ctx = GeometryContext({"vdw_scale": 0.80})
        result = placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        
        # Distance should be 0.80 * (1.70 + 1.70) = 2.72
        assert result["distance"] == pytest.approx(2.72, abs=0.01)
    
    def test_bond_kind_scale_override(self):
        """Test per-bond-kind scale factor override."""
        placer = DockPlacer()
        
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0])
        right = create_simple_monomer("R1", [-1.0, 0.0, 0.0])
        
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        
        # Set per-bond scales
        ctx = GeometryContext({
            "vdw_scale": 0.80,  # Default
            "vdw_scales_by_bond": {
                "-": 0.80,
                "=": 0.72,
                "#": 0.66,
            }
        })
        
        # Single bond
        result_single = placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        assert result_single["distance"] == pytest.approx(0.80 * 2 * 1.70, abs=0.01)
        
        # Double bond
        result_double = placer.place_pair(left, right, left_port, right_port, bond_kind="=", ctx=ctx)
        assert result_double["distance"] == pytest.approx(0.72 * 2 * 1.70, abs=0.01)
        
        # Triple bond
        result_triple = placer.place_pair(left, right, left_port, right_port, bond_kind="#", ctx=ctx)
        assert result_triple["distance"] == pytest.approx(0.66 * 2 * 1.70, abs=0.01)
    
    def test_placement_log(self):
        """Test that placement details are logged to context."""
        placer = DockPlacer()
        
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0])
        right = create_simple_monomer("R1", [-1.0, 0.0, 0.0])
        
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        
        ctx = GeometryContext({"placement_log": []})
        placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        
        # Check that log was populated
        assert len(ctx["placement_log"]) == 1
        log_entry = ctx["placement_log"][0]
        
        assert "left_port_pos" in log_entry
        assert "right_port_pos" in log_entry
        assert "distance" in log_entry
        assert "rotation" in log_entry
        assert "translation" in log_entry
    
    def test_90_degree_rotation(self):
        """Test placement requiring 90° rotation."""
        placer = DockPlacer()
        
        # Left monomer: port at +x pointing in +x
        left = create_simple_monomer("L1", [1.0, 0.0, 0.0], port_dir=[1.0, 0.0, 0.0])
        
        # Right monomer: port at +y pointing in +y (needs to rotate to -x)
        right = create_simple_monomer("R1", [0.0, 1.0, 0.0], port_dir=[0.0, 1.0, 0.0])
        
        left_port = left.get_port("port_L1")
        assert left_port is not None
        right_port = right.get_port("port_R1")
        assert right_port is not None
        
        ctx = GeometryContext({"vdw_scale": 0.80})
        result = placer.place_pair(left, right, left_port, right_port, bond_kind="-", ctx=ctx)
        
        R = result["R_right"]
        
        # Check rotation: should map +y to -x
        right_dir = [0.0, 1.0, 0.0]
        neg_left_dir = [-1.0, 0.0, 0.0]
        
        rotated_dir = [
            R[0][0] * right_dir[0] + R[0][1] * right_dir[1] + R[0][2] * right_dir[2],
            R[1][0] * right_dir[0] + R[1][1] * right_dir[1] + R[1][2] * right_dir[2],
            R[2][0] * right_dir[0] + R[2][1] * right_dir[1] + R[2][2] * right_dir[2],
        ]
        
        assert _dot(rotated_dir, neg_left_dir) == pytest.approx(1.0, abs=1e-6)
