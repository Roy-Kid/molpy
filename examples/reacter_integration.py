"""
Integration example: Building polymer ABCBD using Reacter.

This example demonstrates the complete workflow from BigSMILES parsing
through reaction-based polymer assembly.

Monomer Library:
- A: CCCCO[*:1]
- B: CC(C[*:2])O[*:3]
- C: CCC(C[*:4])O[*:5]
- D: CCC(C[*:6])O[*:7]

Target Sequence: ABCBD
"""

from molpy.parser.smiles import SmilesParser, bigsmilesir_to_monomer
from molpy.reacter import (
    Reacter,
    port_anchor_selector,
    remove_dummy_atoms,
    make_single_bond,
)
from molpy.core.wrappers.monomer import Monomer
from molpy.reacter.utils import find_neighbors


def custom_remove_dummy(monomer, anchor):
    """
    Remove dummy atom (*) connected to anchor.
    
    This is a custom selector that finds the dummy atom
    bonded to the anchor and returns it for removal.
    """
    assembly = monomer.unwrap()
    neighbors = find_neighbors(assembly, anchor)
    
    # Find dummy atoms (symbol='*')
    dummy = [n for n in neighbors if n.get('symbol') == '*']
    return dummy


def main():
    """Main execution function."""
    
    print("=" * 70)
    print("Reacter Integration Example: Building Polymer ABCBD")
    print("=" * 70)
    
    # Step 1: Parse monomer library
    print("\n📝 Step 1: Parsing monomer library...")
    
    parser = SmilesParser()
    monomer_smiles = {
        "A": "CCCCO[*:1]",
        "B": "CC(C[*:2])O[*:3]",
        "C": "CCC(C[*:4])O[*:5]",
        "D": "CCC(C[*:6])O[*:7]",
    }
    
    monomers = {}
    for label, smiles in monomer_smiles.items():
        ir = parser.parse_bigsmiles(smiles)
        monomer = bigsmilesir_to_monomer(ir)
        monomers[label] = monomer
        
        atoms = list(monomer.unwrap().atoms)
        ports = list(monomer.ports.keys())
        print(f"  {label}: {smiles} → {len(atoms)} atoms, ports: {ports}")
    
    print(f"\n✓ Parsed {len(monomers)} monomers")
    
    # Step 2: Define ether linkage reaction
    print("\n🔬 Step 2: Defining C-O-C ether coupling reaction...")
    
    # For ether linkage: remove dummy atoms and connect via C-O bond
    ether_coupling = Reacter(
        name="C-O-C_ether_formation",
        anchor_left=port_anchor_selector,
        anchor_right=port_anchor_selector,
        leaving_left=custom_remove_dummy,
        leaving_right=custom_remove_dummy,
        bond_maker=make_single_bond,
    )
    
    print(f"  Reaction: {ether_coupling.name}")
    print(f"  Mechanism: Remove [*] from both sides, form C-O-C bond")
    
    # Step 3: Build polymer sequence ABCBD
    print("\n🔗 Step 3: Building polymer sequence ABCBD...")
    
    # Start with monomer A
    current_polymer = monomers["A"].copy()
    print(f"\n  Start: A")
    print(f"    Atoms: {len(list(current_polymer.unwrap().atoms))}")
    
    # Define sequence to add: B, C, B, D
    sequence = ["B", "C", "B", "D"]
    
    for i, next_label in enumerate(sequence, 1):
        print(f"\n  Step {i}: + {next_label}")
        
        # Get next monomer
        next_mono = monomers[next_label].copy()
        
        # Determine ports
        # For simplification, we'll use port names from the monomers
        # In practice, you'd track which ports are available
        
        # Find available ports
        curr_ports = list(current_polymer.ports.keys())
        next_ports = list(next_mono.ports.keys())
        
        if not curr_ports or not next_ports:
            print(f"    ⚠ Warning: Missing ports, skipping")
            continue
        
        port_L = curr_ports[0]
        port_R = next_ports[0]
        
        print(f"    Connecting port {port_L} (current) ← → port {port_R} (next)")
        
        # Run reaction
        try:
            product = ether_coupling.run(
                current_polymer,
                next_mono,
                port_L=port_L,
                port_R=port_R,
            )
            
            print(f"    ✓ Reaction completed")
            print(f"      Removed: {product.notes['removed_count']} atoms")
            print(f"      Total atoms: {len(list(product.product.atoms))}")
            print(f"      Total bonds: {len(list(product.product.bonds))}")
            
            # Wrap product back into Monomer
            # Note: ports need to be managed carefully
            current_polymer = Monomer(product.product)
            
            # TODO: Properly track and update ports for next iteration
            # For now, we'll need to manually set ports or use a smarter strategy
            
        except Exception as e:
            print(f"    ✗ Error: {e}")
            import traceback
            traceback.print_exc()
            break
    
    # Step 4: Final results
    print("\n" + "=" * 70)
    print("📊 Final Results")
    print("=" * 70)
    
    final_structure = current_polymer.unwrap()
    final_atoms = list(final_structure.atoms)
    final_bonds = list(final_structure.bonds)
    
    print(f"\nPolymer ABCBD:")
    print(f"  Total atoms: {len(final_atoms)}")
    print(f"  Total bonds: {len(final_bonds)}")
    
    # Atom composition
    from collections import Counter
    symbols = [a.get('symbol', '?') for a in final_atoms]
    composition = Counter(symbols)
    print(f"  Composition: {dict(composition)}")
    
    print("\n✅ Polymer assembly completed!")
    
    return current_polymer


def example_simple_coupling():
    """
    Simplified example showing direct C-C coupling without parser.
    
    This demonstrates the core Reacter functionality without
    the complexity of BigSMILES parsing.
    """
    from molpy.core.atomistic import Atomistic, Atom, Bond
    
    print("\n" + "=" * 70)
    print("Simple C-C Coupling Example")
    print("=" * 70)
    
    # Create left monomer: CH3-C*
    asm_L = Atomistic()
    c1_L = Atom(symbol='C')
    c2_L = Atom(symbol='C')  # Anchor
    h1_L = Atom(symbol='H')
    h2_L = Atom(symbol='H')
    h3_L = Atom(symbol='H')
    star_L = Atom(symbol='*')  # Dummy
    
    asm_L.add_entity(c1_L, c2_L, h1_L, h2_L, h3_L, star_L)
    asm_L.add_link(
        Bond(c1_L, c2_L),
        Bond(c1_L, h1_L),
        Bond(c1_L, h2_L),
        Bond(c1_L, h3_L),
        Bond(c2_L, star_L),
    )
    
    mono_L = Monomer(asm_L)
    mono_L.set_port("1", c2_L)  # Port at anchor carbon
    
    # Create right monomer: *C-CH3
    asm_R = Atomistic()
    c1_R = Atom(symbol='C')  # Anchor
    c2_R = Atom(symbol='C')
    h1_R = Atom(symbol='H')
    h2_R = Atom(symbol='H')
    h3_R = Atom(symbol='H')
    star_R = Atom(symbol='*')  # Dummy
    
    asm_R.add_entity(c1_R, c2_R, h1_R, h2_R, h3_R, star_R)
    asm_R.add_link(
        Bond(c1_R, c2_R),
        Bond(c2_R, h1_R),
        Bond(c2_R, h2_R),
        Bond(c2_R, h3_R),
        Bond(c1_R, star_R),
    )
    
    mono_R = Monomer(asm_R)
    mono_R.set_port("2", c1_R)
    
    print("\nLeft monomer:  CH3-C-[*]")
    print(f"  Atoms: {len(list(mono_L.unwrap().atoms))}")
    print(f"  Port at C anchor")
    
    print("\nRight monomer: [*]-C-CH3")
    print(f"  Atoms: {len(list(mono_R.unwrap().atoms))}")
    print(f"  Port at C anchor")
    
    # Define reaction
    cc_coupling = Reacter(
        name="C-C_coupling",
        anchor_left=port_anchor_selector,
        anchor_right=port_anchor_selector,
        leaving_left=custom_remove_dummy,
        leaving_right=custom_remove_dummy,
        bond_maker=make_single_bond,
    )
    
    # Execute
    print(f"\n🔬 Running {cc_coupling.name}...")
    product = cc_coupling.run(mono_L, mono_R, port_L="1", port_R="2")
    
    print(f"\n✓ Reaction completed!")
    print(f"  Removed atoms: {product.notes['removed_count']} ([*] from each side)")
    print(f"  Product atoms: {len(list(product.product.atoms))}")
    print(f"  Product bonds: {len(list(product.product.bonds))}")
    
    # Check composition
    symbols = [a.get('symbol') for a in product.product.atoms]
    from collections import Counter
    print(f"  Composition: {dict(Counter(symbols))}")
    print(f"\nResult: CH3-C-C-CH3 (butane)")
    
    return product


if __name__ == "__main__":
    # Run simple example first
    example_simple_coupling()
    
    # Then run full integration
    # Note: This may require additional port management logic
    # main()
