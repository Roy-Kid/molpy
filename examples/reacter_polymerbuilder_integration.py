"""
Complete example: Building ABCBD polymer with Reacter and retypification.

This example demonstrates:
1. Using ReacterConnector with default and override reacters
2. Building sequential polymer with explicit port selection
3. Tracking topology changes (angles, dihedrals)
4. Batch retypification of modified atoms

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
    make_single_bond,
)
from molpy.reacter.selectors import remove_dummy_atoms
from molpy.reacter.utils import find_neighbors
from molpy.builder.polymer import PolymerBuilder, ReacterConnector


def custom_remove_dummy(monomer, anchor):
    """Remove dummy atom (*) connected to anchor."""
    assembly = monomer.unwrap()
    neighbors = find_neighbors(assembly, anchor)
    return [n for n in neighbors if n.get('symbol') == '*']


# Port selection strategy: use first available port
def first_available_ports(left, right, left_ports, right_ports, ctx):
    """Select the first available port from each monomer."""
    if not left_ports or not right_ports:
        raise ValueError("No available ports")
    return list(left_ports.keys())[0], list(right_ports.keys())[0]



def main():
    """Main workflow demonstration."""
    
    print("=" * 70)
    print("Reacter + PolymerBuilder Integration Example")
    print("Building Polymer: ABCBD")
    print("=" * 70)
    
    # Step 1: Parse monomer library
    print("\n📝 Step 1: Parse monomer library")
    print("-" * 70)
    
    parser = SmilesParser()
    monomer_smiles = {
        "A": "CCCCO[*:1]",
        "B": "CC(C[*:2])O[*:3]",
        "C": "CCC(C[*:4])O[*:5]",
        "D": "CCC(C[*:6])O[*:7]",
    }
    
    library = {}
    for label, smiles in monomer_smiles.items():
        ir = parser.parse_bigsmiles(smiles)
        monomer = bigsmilesir_to_monomer(ir)
        library[label] = monomer
        
        atoms = list(monomer.unwrap().atoms)
        ports = list(monomer.ports.keys())
        print(f"  {label}: {smiles}")
        print(f"     → {len(atoms)} atoms, ports: {ports}")
    
    # Step 2: Define reactions
    print("\n🔬 Step 2: Define reactions")
    print("-" * 70)
    
    # Default: C-O-C ether linkage
    default_ether = Reacter(
        name="C-O-C_ether_formation",
        anchor_left=port_anchor_selector,
        anchor_right=port_anchor_selector,
        leaving_left=custom_remove_dummy,
        leaving_right=custom_remove_dummy,
        bond_maker=make_single_bond,
    )
    print(f"  Default: {default_ether.name}")
    
    # Special reaction for B-C connection (example)
    # In this case, we'll use the same reaction, but you could customize it
    special_BC = Reacter(
        name="B-C_special_ether",
        anchor_left=port_anchor_selector,
        anchor_right=port_anchor_selector,
        leaving_left=custom_remove_dummy,
        leaving_right=custom_remove_dummy,
        bond_maker=make_single_bond,
    )
    print(f"  Override (B-C): {special_BC.name}")
    
    # Step 3: Create connector
    print("\n🔗 Step 3: Create connector with overrides")
    print("-" * 70)
    
    connector = ReacterConnector(
        default=default_ether,
        port_strategy=first_available_ports,  # Use explicit port selection
        overrides={
            ('B', 'C'): special_BC,  # Use special reacter for B-C
        },
    )
    print("  ✓ ReacterConnector initialized")
    print("    - Default reacter for most connections")
    print("    - Override: B-C uses special reacter")
    print("    - Port strategy: first_available_ports")
    
    # Step 4: Build polymer using PolymerBuilder
    print("\n🏗️  Step 4: Build polymer sequence")
    print("-" * 70)
    
    try:
        polymer = PolymerBuilder.linear(
            sequence='ABCBD',
            library=library,
            connector=connector,
        )
        
        print("  ✓ Polymer built successfully!")
        print(f"    Total atoms: {len(list(polymer.unwrap().atoms))}")
        print(f"    Total bonds: {len(list(polymer.unwrap().bonds))}")
        
    except Exception as e:
        print(f"  ⚠ Error during build: {e}")
        import traceback
        traceback.print_exc()
        polymer = None
    
    # Step 5: Analyze reaction metadata
    print("\n📊 Step 5: Analyze reaction metadata")
    print("-" * 70)
    
    history = connector.get_history()
    print(f"  Total reactions performed: {len(history)}")
    
    total_removed = 0
    total_new_bonds = 0
    total_new_angles = 0
    total_new_dihedrals = 0
    
    for i, product in enumerate(history, 1):
        notes = product.notes
        total_removed += notes['removed_count']
        total_new_bonds += len(notes.get('new_bonds', []))
        total_new_angles += len(notes.get('new_angles', []))
        total_new_dihedrals += len(notes.get('new_dihedrals', []))
        
        print(f"\n  Reaction {i}: {notes['reaction_name']}")
        print(f"    Removed: {notes['removed_count']} atoms")
        print(f"    New bonds: {len(notes.get('new_bonds', []))}")
        print(f"    New angles: {len(notes.get('new_angles', []))}")
        print(f"    New dihedrals: {len(notes.get('new_dihedrals', []))}")
        print(f"    Modified atoms: {len(notes.get('modified_atoms', []))}")
    
    print(f"\n  Summary:")
    print(f"    Total removed: {total_removed} atoms")
    print(f"    Total new bonds: {total_new_bonds}")
    print(f"    Total new angles: {total_new_angles}")
    print(f"    Total new dihedrals: {total_new_dihedrals}")
    
    # Step 6: Retypification
    print("\n🔄 Step 6: Prepare for retypification")
    print("-" * 70)
    
    if connector.needs_retypification():
        modified_atoms = connector.get_all_modified_atoms()
        print(f"  ⚠ Retypification needed!")
        print(f"  Modified atoms: {len(modified_atoms)}")
        
        # Simulate retypification
        print("\n  Retypification workflow:")
        print("    1. Collect all modified atoms")
        print("    2. Clear old atom types (if any)")
        print("    3. Run OplsAtomTypifier on modified atoms")
        print("    4. Re-match bond/angle/dihedral parameters")
        
        print("\n  Example code:")
        print("    ```python")
        print("    from molpy.typifier.atomistic import OplsAtomTypifier")
        print("    from molpy.io import read_xml_forcefield")
        print("")
        print("    ff = read_xml_forcefield('oplsaa.xml')")
        print("    typifier = OplsAtomTypifier(ff)")
        print("")
        print("    # Retypify only modified atoms")
        print("    for atom in modified_atoms:")
        print("        # Clear old type")
        print("        atom.data.pop('type', None)")
        print("")
        print("    # Run typifier on whole assembly")
        print("    typifier.typify(polymer.unwrap())")
        print("    ```")
    else:
        print("  ✓ No retypification needed")
    
    # Step 7: Final summary
    print("\n" + "=" * 70)
    print("✅ Workflow Complete!")
    print("=" * 70)
    
    if polymer:
        print(f"\nFinal Polymer Statistics:")
        print(f"  Sequence: ABCBD")
        print(f"  Atoms: {len(list(polymer.unwrap().atoms))}")
        print(f"  Bonds: {len(list(polymer.unwrap().bonds))}")
        
        # Composition
        from collections import Counter
        symbols = [a.get('symbol', '?') for a in polymer.unwrap().atoms]
        composition = Counter(symbols)
        print(f"  Composition: {dict(composition)}")
    
    print("\n📚 Key Features Demonstrated:")
    print("  ✓ Default reacter for standard connections")
    print("  ✓ Override reacters for special monomer pairs")
    print("  ✓ Explicit port selection strategy")
    print("  ✓ Topology tracking (angles, dihedrals)")
    print("  ✓ Atom type change tracking")
    print("  ✓ Retypification metadata")
    print("  ✓ Reaction history for debugging")


if __name__ == "__main__":
    main()
