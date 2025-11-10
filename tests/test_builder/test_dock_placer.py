"""
Tests for DockPlacer with automatic orientation strategies.

Tests cover:
1. Strategy 1: Nearest heavy neighbor orientation
2. Strategy 2: Role-based main-chain orientation
3. Fallback handling and error cases
4. Per-bond VDW scaling
5. Integration with PolymerBuilder.linear
"""

import pytest
import math
from molpy.builder.polymer.geom_utils import DockPlacer, rodrigues, get_vdw_radius
from molpy.builder.errors import OrientationUnavailableError, PositionMissingError
from molpy.core.atomistic import Atomistic, Atom
from molpy.core.wrappers.monomer import Monomer


def test_rodrigues_aligned():
    """Test Rodrigues rotation for already aligned vectors."""
    a = [1.0, 0.0, 0.0]
    b = [1.0, 0.0, 0.0]
    R = rodrigues(a, b)
    
    # Should be identity
    expected = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]]
    for i in range(3):
        for j in range(3):
            assert abs(R[i][j] - expected[i][j]) < 1e-10


def test_rodrigues_antiparallel():
    """Test Rodrigues rotation for 180° rotation."""
    a = [1.0, 0.0, 0.0]
    b = [-1.0, 0.0, 0.0]
    R = rodrigues(a, b)
    
    # Apply rotation
    from molpy.builder.polymer.geom_utils import _mat_vec_mult
    result = _mat_vec_mult(R, a)
    
    # Should align with b
    for i in range(3):
        assert abs(result[i] - b[i]) < 1e-6


def test_rodrigues_90_degrees():
    """Test Rodrigues rotation for 90° rotation."""
    a = [1.0, 0.0, 0.0]
    b = [0.0, 1.0, 0.0]
    R = rodrigues(a, b)
    
    from molpy.builder.polymer.geom_utils import _mat_vec_mult
    result = _mat_vec_mult(R, a)
    
    # Should align with b
    for i in range(3):
        assert abs(result[i] - b[i]) < 1e-6


def test_vdw_radius_from_element():
    """Test VDW radius lookup from Element class."""
    atom = Atom(name="C1", symbol="C", xyz=[0.0, 0.0, 0.0])
    r = get_vdw_radius(atom)
    
    # Carbon VDW radius
    assert abs(r - 1.70) < 0.01
    
    # Hydrogen
    atom_h = Atom(name="H1", symbol="H", xyz=[0.0, 0.0, 0.0])
    r_h = get_vdw_radius(atom_h)
    assert abs(r_h - 1.20) < 0.01


def test_vdw_radius_explicit():
    """Test VDW radius from explicit data."""
    atom = Atom(name="X", symbol="X", xyz=[0.0, 0.0, 0.0])
    atom.data['vdw'] = 2.5
    
    r = get_vdw_radius(atom)
    assert abs(r - 2.5) < 1e-10


def test_nearest_heavy_orientation():
    """Test Strategy 1: Nearest heavy neighbor orientation."""
    # Create monomer with C-C-H structure
    # C1 at origin, C2 at (1.5, 0, 0), H at (0, 1.0, 0)
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c2 = atomistic.add_atom(symbol="C", pos=[1.5, 0.0, 0.0])  # Nearest heavy
    h1 = atomistic.add_atom(symbol="H", pos=[0.0, 1.0, 0.0])
    
    atomistic.add_bond(c1, c2, order=1)
    atomistic.add_bond(c1, h1, order=1)
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    # Test orientation inference
    placer = DockPlacer()
    port = monomer.ports["p1"]
    
    orient = placer._orient_from_nearest_heavy(port, monomer)
    
    # Should point from C1 toward C2: [1, 0, 0]
    assert orient is not None
    assert abs(orient[0] - 1.0) < 1e-6
    assert abs(orient[1]) < 1e-6
    assert abs(orient[2]) < 1e-6


def test_nearest_heavy_multiple_neighbors():
    """Test Strategy 1 picks nearest heavy when multiple exist."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c2_near = atomistic.add_atom(symbol="C", pos=[1.0, 0.0, 0.0])  # Nearest
    c3_far = atomistic.add_atom(symbol="C", pos=[0.0, 3.0, 0.0])   # Farther
    
    atomistic.add_bond(c1, c2_near, order=1)
    atomistic.add_bond(c1, c3_far, order=1)
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer()
    port = monomer.ports["p1"]
    
    orient = placer._orient_from_nearest_heavy(port, monomer)
    
    # Should point toward c2_near: [1, 0, 0]
    assert orient is not None
    assert abs(orient[0] - 1.0) < 1e-6
    assert abs(orient[1]) < 1e-6


def test_role_main_chain_orientation():
    """Test Strategy 2: Role-based main-chain orientation."""
    atomistic = Atomistic()
    c_left = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_mid = atomistic.add_atom(symbol="C", pos=[1.5, 0.0, 0.0])
    c_right = atomistic.add_atom(symbol="C", pos=[3.0, 0.0, 0.0])
    
    atomistic.add_bond(c_left, c_mid, order=1)
    atomistic.add_bond(c_mid, c_right, order=1)
    
    monomer = Monomer(atomistic)
    monomer.set_port("left", c_left, role="left", bond_kind="-", compat="*", multiplicity=1, priority=0)
    monomer.set_port("right", c_right, role="right", bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer(allow_role_direction=True)
    
    # Test left port (should use -g)
    left_port = monomer.ports["left"]
    orient_left = placer._orient_from_role_main_chain(left_port, monomer)
    
    # Main chain g = right - left = [3, 0, 0] - [0, 0, 0] = [3, 0, 0] -> unit [1, 0, 0]
    # Left port uses -g = [-1, 0, 0]
    assert orient_left is not None
    assert abs(orient_left[0] - (-1.0)) < 1e-6
    assert abs(orient_left[1]) < 1e-6
    assert abs(orient_left[2]) < 1e-6
    
    # Test right port (should use +g)
    right_port = monomer.ports["right"]
    orient_right = placer._orient_from_role_main_chain(right_port, monomer)
    
    assert orient_right is not None
    assert abs(orient_right[0] - 1.0) < 1e-6
    assert abs(orient_right[1]) < 1e-6
    assert abs(orient_right[2]) < 1e-6


def test_orientation_unavailable_error():
    """Test OrientationUnavailableError when no strategy succeeds."""
    # Create monomer with isolated atom (no neighbors, no roles)
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer(allow_role_direction=False)
    port = monomer.ports["p1"]
    
    with pytest.raises(OrientationUnavailableError):
        placer._infer_orientation(port, monomer)


def test_explicit_orientation_override():
    """Test that explicit orientation takes precedence."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c1.data['orientation'] = [0.0, 1.0, 0.0]  # Explicit Y direction
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer()
    port = monomer.ports["p1"]
    
    orient = placer._infer_orientation(port, monomer)
    
    # Should use explicit orientation
    assert abs(orient[0]) < 1e-6
    assert abs(orient[1] - 1.0) < 1e-6
    assert abs(orient[2]) < 1e-6


def test_per_bond_vdw_scaling():
    """Test per-bond VDW scaling."""
    # Create simple monomers with ports
    left_atomistic = Atomistic()
    c_left = left_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_left.data['orientation'] = [1.0, 0.0, 0.0]
    
    right_atomistic = Atomistic()
    c_right = right_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_right.data['orientation'] = [1.0, 0.0, 0.0]
    
    left_mono = Monomer(left_atomistic)
    left_mono.set_port("out", c_left, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    right_mono = Monomer(right_atomistic)
    right_mono.set_port("in", c_right, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    # Test single bond (scale 0.80)
    placer_single = DockPlacer(
        vdw_scale=0.90,  # Default
        vdw_scales_by_bond={"-": 0.80, "=": 0.72}
    )
    
    result_single = placer_single.place_pair(
        left_mono, right_mono,
        left_mono.ports["out"], right_mono.ports["in"],
        bond_kind="-",
        ctx={}
    )
    
    # VDW radius for carbon: ~1.70
    # Distance should be 0.80 * (1.70 + 1.70) = 2.72
    expected_dist_single = 0.80 * 2 * 1.70
    assert abs(result_single["distance"] - expected_dist_single) < 0.01
    
    # Test double bond (scale 0.72)
    result_double = placer_single.place_pair(
        left_mono, right_mono,
        left_mono.ports["out"], right_mono.ports["in"],
        bond_kind="=",
        ctx={}
    )
    
    expected_dist_double = 0.72 * 2 * 1.70
    assert abs(result_double["distance"] - expected_dist_double) < 0.01
    
    # Verify different distances
    assert abs(result_single["distance"] - result_double["distance"]) > 0.1


def test_docking_alignment():
    """Test that docking aligns orientations correctly (uR → -uL)."""
    # Left monomer at origin, orientation [1, 0, 0]
    left_atomistic = Atomistic()
    c_left = left_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_left.data['orientation'] = [1.0, 0.0, 0.0]
    
    # Right monomer at [5, 0, 0], orientation [1, 0, 0]
    right_atomistic = Atomistic()
    c_right = right_atomistic.add_atom(symbol="C", pos=[5.0, 0.0, 0.0])
    c_right.data['orientation'] = [1.0, 0.0, 0.0]
    
    left_mono = Monomer(left_atomistic)
    left_mono.set_port("out", c_left, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    right_mono = Monomer(right_atomistic)
    right_mono.set_port("in", c_right, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer(vdw_scale=0.80)
    
    result = placer.place_pair(
        left_mono, right_mono,
        left_mono.ports["out"], right_mono.ports["in"],
        bond_kind="-",
        ctx={}
    )
    
    R = result["R_right"]
    t = result["t_right"]
    
    # uL = [1, 0, 0], uR = [1, 0, 0]
    # Should rotate uR to -uL = [-1, 0, 0] (180° rotation)
    from molpy.builder.polymer.geom_utils import _mat_vec_mult
    uR_rotated = _mat_vec_mult(R, [1.0, 0.0, 0.0])
    
    # Should be [-1, 0, 0]
    assert abs(uR_rotated[0] - (-1.0)) < 1e-6
    assert abs(uR_rotated[1]) < 1e-6
    assert abs(uR_rotated[2]) < 1e-6


def test_position_missing_error():
    """Test PositionMissingError when atom has no coordinates."""
    atomistic = Atomistic()
    # Use add_atom to create atom without position
    # (don't provide pos argument)
    # Actually, we need to test when atom exists but has no pos
    c1 = atomistic.add_atom(symbol="C")
    # Remove pos if it was set
    if 'pos' in c1.data:
        del c1.data['pos']
    if 'xyz' in c1.data:
        del c1.data['xyz']
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer()
    
    with pytest.raises(PositionMissingError):
        placer._get_position(c1)


def test_placement_logging():
    """Test that placement details are logged when ctx has placement_log."""
    left_atomistic = Atomistic()
    c_left = left_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_left.data['orientation'] = [1.0, 0.0, 0.0]
    
    right_atomistic = Atomistic()
    c_right = right_atomistic.add_atom(symbol="C", pos=[5.0, 0.0, 0.0])
    c_right.data['orientation'] = [1.0, 0.0, 0.0]
    
    left_mono = Monomer(left_atomistic)
    left_mono.set_port("out", c_left, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    right_mono = Monomer(right_atomistic)
    right_mono.set_port("in", c_right, role=None, bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    placer = DockPlacer(vdw_scale=0.80)
    ctx = {"placement_log": []}
    
    result = placer.place_pair(
        left_mono, right_mono,
        left_mono.ports["out"], right_mono.ports["in"],
        bond_kind="-",
        ctx=ctx
    )
    
    # Check log was updated
    assert len(ctx["placement_log"]) == 1
    log_entry = ctx["placement_log"][0]
    
    assert "distance" in log_entry
    assert "rvdw_L" in log_entry
    assert "rvdw_R" in log_entry
    assert "scale" in log_entry
    assert abs(log_entry["scale"] - 0.80) < 1e-10


def test_integration_with_linear_builder():
    """Test DockPlacer integration with PolymerBuilder.linear."""
    from molpy.builder.polymer.linear import PolymerBuilder
    from molpy.builder.polymer.connectors import AutoConnector
    
    # Create simple monomers with explicit orientations
    # Monomer A: C at origin, orientation +X
    a_atomistic = Atomistic()
    c_a = a_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_a.data['orientation'] = [1.0, 0.0, 0.0]
    
    mono_a = Monomer(a_atomistic)
    mono_a.set_port("out", c_a, role="right", bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    # Monomer B: C at origin, orientation +X
    b_atomistic = Atomistic()
    c_b1 = b_atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c_b1.data['orientation'] = [-1.0, 0.0, 0.0]
    c_b2 = b_atomistic.add_atom(symbol="C", pos=[1.5, 0.0, 0.0])
    c_b2.data['orientation'] = [1.0, 0.0, 0.0]
    b_atomistic.add_bond(c_b1, c_b2, order=1)
    
    mono_b = Monomer(b_atomistic)
    mono_b.set_port("in", c_b1, role="left", bond_kind="-", compat="*", multiplicity=1, priority=0)
    mono_b.set_port("out", c_b2, role="right", bond_kind="-", compat="*", multiplicity=1, priority=0)
    
    # Build polymer with geometry
    placer = DockPlacer(vdw_scale=0.80)
    ctx = {"placement_log": []}
    
    poly = PolymerBuilder.linear(
        sequence="AB",
        library={"A": mono_a, "B": mono_b},
        connector=AutoConnector(),
        placer=placer,
        geom_ctx=ctx,
    )
    
    # Check that placement occurred
    assert len(ctx["placement_log"]) == 1
    
    # Check that polymer has correct atom count
    atoms = list(poly.unwrap().atoms)
    assert len(atoms) == 3  # 1 from A + 2 from B
    
    # Check that bond was created
    bonds = list(poly.unwrap().bonds)
    assert len(bonds) == 2  # 1 original in B + 1 new connection


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
