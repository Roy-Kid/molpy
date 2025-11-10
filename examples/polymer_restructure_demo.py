"""
Example usage of the reorganized polymer module.

After restructuring:
- builder/polymer/linear.py: PolymerBuilder
- builder/polymer/connectors.py: Connector strategies
- builder/polymer/geom_utils.py: Geometry utilities and Placer strategies
"""

from molpy.builder.polymer import (
    PolymerBuilder,
    AutoConnector,
    DockPlacer,
    GeometryContext,
)
from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer


def create_example_monomer(name: str) -> Monomer:
    """Create a simple monomer for demonstration."""
    atomistic = Atomistic()
    
    # Add atoms
    c1 = atomistic.add_atom(symbol="C", pos=[0.0, 0.0, 0.0])
    c2 = atomistic.add_atom(symbol="C", pos=[1.5, 0.0, 0.0])
    
    # Add bond
    atomistic.add_bond(c1, c2, order=1)
    
    # Create monomer
    mono = Monomer(atomistic)
    
    # Add ports with orientations
    mono.set_port("left", c1, role="left")
    c1.data['orientation'] = [-1.0, 0.0, 0.0]  # Points left
    
    mono.set_port("right", c2, role="right")
    c2.data['orientation'] = [1.0, 0.0, 0.0]  # Points right
    
    return mono


# Example 1: Topology-only assembly (original behavior)
print("=== Example 1: Topology-only ===")
A = create_example_monomer("A")
B = create_example_monomer("B")

poly1 = PolymerBuilder.linear(
    sequence="ABBA",
    library={"A": A, "B": B},
    connector=AutoConnector(),
    # No placer - topology only
)
print(f"Atoms: {len(list(poly1.unwrap().atoms))}")
print(f"Bonds: {len(list(poly1.unwrap().bonds))}")


# Example 2: With geometry (VDW-based docking)
print("\n=== Example 2: With geometry ===")
A2 = create_example_monomer("A")
B2 = create_example_monomer("B")

ctx = GeometryContext({
    "vdw_scale": 0.80,  # 80% of VDW sum
    "placement_log": [],
})

poly2 = PolymerBuilder.linear(
    sequence="ABBA",
    library={"A": A2, "B": B2},
    connector=AutoConnector(),  # Topology
    placer=DockPlacer(),        # Geometry
    geom_ctx=ctx,
)

print(f"Atoms: {len(list(poly2.unwrap().atoms))}")
print(f"Bonds: {len(list(poly2.unwrap().bonds))}")
print(f"Placements: {len(ctx['placement_log'])}")

# Inspect placements
for i, entry in enumerate(ctx["placement_log"]):
    print(f"  Step {i+1}: distance={entry['distance']:.2f} Å, "
          f"VDW={entry['rvdw_L']:.2f}+{entry['rvdw_R']:.2f}")


print("\n=== New structure ===")
print("✅ builder/polymer/")
print("   ├── __init__.py      # Public API")
print("   ├── linear.py        # PolymerBuilder")
print("   ├── connectors.py    # Connector strategies")
print("   └── geom_utils.py    # Geometry + Placer (uses core/ops/geometry)")
print("\n✅ Direct use of core/ops/geometry - no wrapper layer!")
