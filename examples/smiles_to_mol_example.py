"""
Example demonstrating SmilesIR → RDKit Mol conversion and visualization.

This example shows how to:
1. Parse SMILES strings into IR
2. Convert IR to RDKit Mol objects
3. Visualize molecules using draw_molecule
"""

from molpy.parser.smiles import SmilesParser, smilesir_to_mol
from molpy.adapter.rdkit import draw_molecule

# Create parser
parser = SmilesParser()

# Example 1: Simple molecules
print("=== Example 1: Simple Molecules ===")
for smiles in ["CCO", "C=C", "C#N"]:
    ir = parser.parse_smiles(smiles)
    mol = smilesir_to_mol(ir)
    print(f"{smiles}: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")

# Example 2: Aromatic molecules
print("\n=== Example 2: Aromatic Molecules ===")
benzene_ir = parser.parse_smiles("c1ccccc1")
benzene_mol = smilesir_to_mol(benzene_ir)
print(f"Benzene: {benzene_mol.GetNumAtoms()} atoms, all aromatic: {all(a.GetIsAromatic() for a in benzene_mol.GetAtoms())}")

# Example 3: Charged molecules
print("\n=== Example 3: Charged Molecules ===")
ammonium_ir = parser.parse_smiles("[NH4+]")
ammonium_mol = smilesir_to_mol(ammonium_ir)
n_atom = ammonium_mol.GetAtomWithIdx(0)
print(f"Ammonium: charge={n_atom.GetFormalCharge()}, H count={n_atom.GetTotalNumHs()}")

# Example 4: Isotopes
print("\n=== Example 4: Isotopes ===")
c13_ir = parser.parse_smiles("[13C]")
c13_mol = smilesir_to_mol(c13_ir)
c_atom = c13_mol.GetAtomWithIdx(0)
print(f"Carbon-13: isotope={c_atom.GetIsotope()}")

# Example 5: Visualization
print("\n=== Example 5: Visualization ===")
print("Drawing molecules...")

# Draw from IR
aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
aspirin_ir = parser.parse_smiles(aspirin_smiles)
img1 = draw_molecule(aspirin_ir, size=(400, 300))
print(f"Aspirin from IR: {img1.size}")

# Draw from SMILES string directly
img2 = draw_molecule("c1ccccc1", size=(300, 300), show_atom_idx=True)
print(f"Benzene with atom indices: {img2.size}")

# Draw with highlights
img3 = draw_molecule("CCO", highlight_atoms=[1, 2])
print(f"Ethanol with highlights: {img3.size}")

print("\n✓ All examples completed successfully!")
