"""
Example: Generate 3D coordinates for monomers using RDKit.

This example shows how to:
1. Parse a SMILES string to get a molecular structure
2. Convert it to a Monomer wrapper
3. Generate 3D coordinates using RDKit's ETKDG method
4. Access the generated coordinates
"""

from molpy.adapter.converter import convert
from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer
from molpy.parser.smiles import SmilesParser


def create_monomer_with_3d(smiles: str, optimize: bool = True) -> Monomer:
    """
    Create a Monomer from SMILES and generate 3D coordinates.
    
    Args:
        smiles: SMILES string representation
        optimize: Whether to optimize geometry with MMFF94
        
    Returns:
        Monomer with 3D coordinates
    """
    # Parse SMILES
    parser = SmilesParser()
    smiles_ir = parser.parse_smiles(smiles)
    
    # Convert to Atomistic via registered converters
    atomistic = convert(smiles_ir, Atomistic)
    
    # Wrap in Monomer
    monomer = Monomer(atomistic)
    
    # Generate 3D coordinates
    view = monomer.as_rdkit()
    view.embed_3d(optimize=optimize)
    view.flush_coords_back()
    
    return monomer


if __name__ == "__main__":
    # Example 1: Ethanol
    print("=" * 60)
    print("Example 1: Ethanol (CCO)")
    print("=" * 60)
    
    ethanol = create_monomer_with_3d("CCO")
    print(f"Monomer: {ethanol}")
    
    atoms = list(ethanol.unwrap().atoms)
    print(f"\nAtoms and coordinates:")
    for i, atom in enumerate(atoms):
        symbol = atom.data.get("symbol", "?")
        pos = atom.data.get("pos", [0, 0, 0])
        print(f"  Atom {i} ({symbol}): [{pos[0]:8.4f}, {pos[1]:8.4f}, {pos[2]:8.4f}]")
    
    # Example 2: Benzene
    print("\n" + "=" * 60)
    print("Example 2: Benzene (c1ccccc1)")
    print("=" * 60)
    
    benzene = create_monomer_with_3d("c1ccccc1")
    print(f"Monomer: {benzene}")
    
    atoms = list(benzene.unwrap().atoms)
    print(f"\nAtoms and coordinates:")
    for i, atom in enumerate(atoms):
        symbol = atom.data.get("symbol", "?")
        pos = atom.data.get("pos", [0, 0, 0])
        print(f"  Atom {i} ({symbol}): [{pos[0]:8.4f}, {pos[1]:8.4f}, {pos[2]:8.4f}]")
    
    # Example 3: Monomer with port (for polymer assembly)
    print("\n" + "=" * 60)
    print("Example 3: Monomer with connection port")
    print("=" * 60)
    
    styrene = create_monomer_with_3d("C=Cc1ccccc1")  # Styrene monomer
    
    # Add a port at the vinyl carbon for polymerization
    atoms = list(styrene.unwrap().atoms)
    vinyl_carbon = atoms[0]  # First carbon in C=C
    
    styrene.set_port(
        "vinyl",
        vinyl_carbon,
        role="left",
        bond_kind="-",
        multiplicity=1
    )
    
    print(f"Monomer: {styrene}")
    print(f"Ports: {styrene.ports}")
    print(f"\nVinyl carbon position: {vinyl_carbon.data.get('pos')}")
    
    print("\n✓ All examples completed successfully!")
