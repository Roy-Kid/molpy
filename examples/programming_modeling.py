"""
Demonstration of MolPy's programming-based modeling capabilities.

This example shows how to use core primitives to build complex systems.
"""

from molpy import Atomistic
import math
import random


# ========== Define Molecular Templates ==========

class Water(Atomistic):
    """Simple water molecule"""
    def __post_init__(self, **props):
        o = self.add_atom(symbol="O", xyz=[0, 0, 0])
        h1 = self.add_atom(symbol="H", xyz=[0.757, 0.586, 0])
        h2 = self.add_atom(symbol="H", xyz=[-0.757, 0.586, 0])
        self.add_bond(o, h1)
        self.add_bond(o, h2)


class Methane(Atomistic):
    """Simple methane molecule"""
    def __post_init__(self, **props):
        c = self.add_atom(symbol="C", xyz=[0, 0, 0])
        h1 = self.add_atom(symbol="H", xyz=[0.629, 0.629, 0.629])
        h2 = self.add_atom(symbol="H", xyz=[-0.629, -0.629, 0.629])
        h3 = self.add_atom(symbol="H", xyz=[-0.629, 0.629, -0.629])
        h4 = self.add_atom(symbol="H", xyz=[0.629, -0.629, -0.629])
        for h in [h1, h2, h3, h4]:
            self.add_bond(c, h)


# ========== Example 1: Basic System Building ==========

def example_basic_building():
    """Build system by adding molecules with transformations"""
    print("\n=== Example 1: Basic System Building ===")
    
    system = Atomistic()
    
    # Add molecules at different positions
    system += Water()
    system += Water().move([5, 0, 0])
    system += Water().move([10, 0, 0])
    
    # Add with chained transformations
    system += Methane().move([2.5, 5, 0]).rotate([0, 0, 1], 45)
    
    print(f"System: {system}")
    print(f"Total atoms: {len(system)}")
    

# ========== Example 2: Linear Array ==========

def example_linear_array():
    """Create molecules in a line using replicate()"""
    print("\n=== Example 2: Linear Array ===")
    
    # Method 1: Using loop
    system1 = Atomistic()
    for i in range(5):
        system1 += Water().move([i * 5, 0, 0])
    
    print(f"Method 1 (loop): {system1}")
    
    # Method 2: Using replicate
    system2 = Water().replicate(5, lambda mol, i: mol.move([i * 5, 0, 0]))
    
    print(f"Method 2 (replicate): {system2}")


# ========== Example 3: 2D Grid ==========

def example_grid():
    """Create a 2D grid of molecules"""
    print("\n=== Example 3: 2D Grid ===")
    
    nx, ny = 4, 4
    spacing = 5.0
    
    n = nx * ny
    grid = Water().replicate(
        n,
        lambda mol, idx: mol.move([
            (idx % nx) * spacing,
            (idx // nx) * spacing,
            0
        ])
    )
    
    print(f"Grid: {grid}")


# ========== Example 4: Ring Structure ==========

def example_ring():
    """Create molecules arranged in a ring"""
    print("\n=== Example 4: Ring Structure ===")
    
    n_molecules = 8
    radius = 10.0
    
    ring = Atomistic()
    for i in range(n_molecules):
        angle = (i / n_molecules) * 2 * math.pi
        
        ring += Water().move([
            radius * math.cos(angle),
            radius * math.sin(angle),
            0
        ]).rotate([0, 0, 1], math.degrees(angle))
    
    print(f"Ring: {ring}")


# ========== Example 5: Random Assembly ==========

def example_random():
    """Create system with random placement and rotation"""
    print("\n=== Example 5: Random Assembly ===")
    
    system = Atomistic()
    
    molecules = [Water, Methane]
    weights = [0.7, 0.3]
    
    for i in range(20):
        # Random selection
        MolClass = random.choices(molecules, weights=weights)[0]
        
        # Random position
        pos = [random.uniform(0, 30) for _ in range(3)]
        
        # Random rotation
        axis = [random.random(), random.random(), random.random()]
        angle = random.uniform(0, 360)
        
        system += MolClass().move(pos).rotate(axis, angle)
    
    print(f"Random system: {system}")


# ========== Example 6: Cumulative Array ==========

def example_cumulative():
    """MolTemplate-style cumulative transformations"""
    print("\n=== Example 6: Cumulative Array ===")
    
    system = Atomistic()
    current = Water()
    system += current
    
    for i in range(4):
        # Each copy is transformed relative to previous
        current = current.copy().move([3, 0, 0]).rotate([0, 0, 1], 30)
        system += current
    
    print(f"Cumulative array: {system}")


# ========== Example 7: Hierarchical Assembly ==========

def example_hierarchical():
    """Build complex structures from sub-assemblies"""
    print("\n=== Example 7: Hierarchical Assembly ===")
    
    def water_cluster(n=4):
        """Create a small cluster of water molecules"""
        cluster = Atomistic()
        for i in range(n):
            angle = (i / n) * 360
            cluster += Water().move([3, 0, 0]).rotate([0, 0, 1], angle)
        return cluster
    
    # Build system from clusters
    system = Atomistic()
    for i in range(3):
        system += water_cluster(4).move([i * 10, 0, 0])
    
    print(f"Hierarchical system: {system}")


# ========== Example 8: Combining Systems ==========

def example_combination():
    """Combine multiple systems using + operator"""
    print("\n=== Example 8: System Combination ===")
    
    # Create separate layers
    water_layer = Water().replicate(9, lambda mol, i: mol.move([
        (i % 3) * 5,
        (i // 3) * 5,
        0
    ]))
    
    methane_layer = Methane().replicate(4, lambda mol, i: mol.move([
        i * 5,
        0,
        10
    ]))
    
    # Combine using + operator
    system = water_layer + methane_layer
    
    print(f"Combined system: {system}")
    print(f"Water layer: {water_layer}")
    print(f"Methane layer: {methane_layer}")


# ========== Main ==========

if __name__ == "__main__":
    print("=" * 60)
    print("MolPy Programming-Based Modeling Examples")
    print("=" * 60)
    
    example_basic_building()
    example_linear_array()
    example_grid()
    example_ring()
    example_random()
    example_cumulative()
    example_hierarchical()
    example_combination()
    
    print("\n" + "=" * 60)
    print("All examples completed successfully!")
    print("=" * 60)
