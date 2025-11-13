"""
Example demonstrating SmilesIR → RDKit Mol conversion and visualization.

This example shows how to:
1. Parse SMILES strings into IR
2. Convert IR to RDKit Mol objects
3. Visualize molecules using draw_molecule
"""

from rdkit.Chem import Draw

from molpy.adapter.converter import convert
from molpy.adapter.rdkit_adapter import RDKitWrapper
from molpy.parser.smiles import SmilesParser, smilesir_to_mol

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

# Draw from IR using RDKitWrapper
aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
aspirin_ir = parser.parse_smiles(aspirin_smiles)
aspirin_view = convert(aspirin_ir, RDKitWrapper)
aspirin_img = Draw.MolToImage(aspirin_view.mol, size=(400, 300))
print(f"Aspirin from IR: {aspirin_img.size}")

# Draw from SMILES string directly
benzene_ir = parser.parse_smiles("c1ccccc1")
benzene_view = convert(benzene_ir, RDKitWrapper)
benzene_img = Draw.MolToImage(benzene_view.mol, size=(300, 300))
print(f"Benzene (2D depiction): {benzene_img.size}")

# Draw with highlights
ethanol_ir = parser.parse_smiles("CCO")
ethanol_view = convert(ethanol_ir, RDKitWrapper)
drawer = Draw.MolDraw2DCairo(300, 300)
Draw.rdMolDraw2D.PrepareMolForDrawing(ethanol_view.mol)
drawer.DrawMolecule(ethanol_view.mol, highlightAtoms=[1, 2])
drawer.FinishDrawing()
print(f"Ethanol highlighted image bytes: {len(drawer.GetDrawingText())}")

print("\n✓ All examples completed successfully!")
