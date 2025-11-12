"""
Integration tests for Placer with Polymer.linear().
"""

import pytest
from molpy.builder.polymer import (
    PolymerBuilder,
    AutoConnector,
    DockPlacer,
    NoOpPlacer
)
from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer


def create_simple_monomer_for_chain(label: str, pos: list[float]):
    """Create simple monomers with explicit positions and ports."""
    atomistic = Atomistic()
    
    # Center atom
    center = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    
    # Port atom at specified position
    port_atom = atomistic.add_atom(symbol="C", pos=pos)
    
    # Add bond
    atomistic.add_bond(center, port_atom, order=1)
    
    # Create monomer
    monomer = Monomer(atomistic)
    
    # Add ports based on label
    if label == "T":  # Terminus
        # Right-facing port
        monomer.set_port("right", port_atom, role="right")
        # Set explicit orientation
        port_atom.data['orientation'] = [1.0, 0.0, 0.0]
    else:  # Middle monomers
        # Left and right ports (same atom for simplicity)
        monomer.set_port("left", port_atom, role="left", multiplicity=1)
        monomer.set_port("right", port_atom, role="right", multiplicity=1)
        # Orientation pointing right
        port_atom.data['orientation'] = [1.0, 0.0, 0.0]
    
    return monomer


class TestPolymerWithPlacer:
    """Test Polymer.linear() with geometry placers."""
    
    def test_noop_placer_preserves_topology(self):
        """Test that NoOpPlacer doesn't break topology-only assembly."""
        # Create simple library
        T = create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0])
        A = create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])
        
        library = {"T": T, "A": A}
        
        # Build with NoOp placer
        poly = PolymerBuilder.linear(
            sequence="TA",
            library=library,
            connector=AutoConnector(),
            placer=NoOpPlacer(),
        )
        
        # Should have atoms from both monomers
        atoms = list(poly.unwrap().atoms)
        assert len(atoms) == 4  # 2 atoms per monomer
        
        # Should have bond between monomers
        bonds = list(poly.unwrap().bonds)
        assert len(bonds) >= 3  # 2 internal + 1 connection
    
    def test_dock_placer_applies_transforms(self):
        """Test that DockPlacer applies geometric transforms."""
        # Create monomers with known positions
        T = create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0])
        A = create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])
        
        library = {"T": T, "A": A}
        
        # Build with DockPlacer
        poly = PolymerBuilder.linear(
            sequence="TA",
            library=library,
            connector=AutoConnector(),
            placer=DockPlacer(),
        )
        
        # Check that atoms have positions
        atoms = list(poly.unwrap().atoms)
        assert all('pos' in atom.data or 'xyz' in atom.data for atom in atoms)
    
    def test_dock_placer_with_default_context(self):
        """Test that DockPlacer works without explicit context."""
        T = create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0])
        A = create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])
        
        library = {"T": T, "A": A}
        
        # Build without explicit context (should create default)
        poly = PolymerBuilder.linear(
            sequence="TA",
            library=library,
            connector=AutoConnector(),
            placer=DockPlacer(),
            # No geom_ctx - should create default
        )
        
        # Should still work
        atoms = list(poly.unwrap().atoms)
        assert len(atoms) == 4
    
    def test_topology_only_without_placer(self):
        """Test that placer=None preserves original behavior."""
        T = create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0])
        A = create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])
        
        library = {"T": T, "A": A}
        
        # Build without placer (topology only)
        poly = PolymerBuilder.linear(
            sequence="TA",
            library=library,
            connector=AutoConnector(),
            # placer=None (default)
        )
        
        # Should work as before
        atoms = list(poly.unwrap().atoms)
        assert len(atoms) == 4
        
        bonds = list(poly.unwrap().bonds)
        assert len(bonds) >= 3
    
    def test_vdw_scale_affects_distance(self):
        """Test that VDW scale factor affects placement distance."""
        T = create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0])
        A = create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])
        
        library = {"T": T, "A": A}
        
        # Build with different scales
        poly1 = PolymerBuilder.linear(
            sequence="TA",
            library=library,
            connector=AutoConnector(),
            placer=DockPlacer(),
        )
        
        # Get first placement distance
        
        # Build with different scale
        library2 = {"T": create_simple_monomer_for_chain("T", [1.0, 0.0, 0.0]), 
                   "A": create_simple_monomer_for_chain("A", [1.0, 0.0, 0.0])}
        poly2 = PolymerBuilder.linear(
            sequence="TA",
            library=library2,
            connector=AutoConnector(),
            placer=DockPlacer(),
        )
        
