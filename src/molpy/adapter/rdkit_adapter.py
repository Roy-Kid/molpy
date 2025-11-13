"""RDKit adapter - explicit conversion functions."""

from typing import Optional
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from molpy.core.atomistic import Atomistic, Atom, Bond
from molpy.core.wrappers.monomer import Monomer


def mol_to_atomistic(mol: Chem.Mol) -> Atomistic:
    """Convert RDKit Mol to Atomistic."""
    atomistic = Atomistic()
    atom_map = {}
    
    # Add atoms
    for rdkit_atom in mol.GetAtoms():
        attrs = {
            "symbol": rdkit_atom.GetSymbol(),
            "charge": rdkit_atom.GetFormalCharge(),
        }
        
        # Add coordinates if available
        if mol.GetNumConformers() > 0:
            conf = mol.GetConformer()
            pos = conf.GetAtomPosition(rdkit_atom.GetIdx())
            attrs["xyz"] = [pos.x, pos.y, pos.z]
        
        atom = atomistic.add_atom(**attrs)
        atom_map[rdkit_atom.GetIdx()] = atom
    
    # Add bonds
    for rdkit_bond in mol.GetBonds():
        a = atom_map[rdkit_bond.GetBeginAtomIdx()]
        b = atom_map[rdkit_bond.GetEndAtomIdx()]
        
        bond_order_map = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 1.5,
        }
        order = bond_order_map.get(rdkit_bond.GetBondType(), 1)
        atomistic.add_bond(a, b, order=order)
    
    return atomistic


def atomistic_to_mol(atomistic: Atomistic) -> Chem.Mol:
    """Convert Atomistic to RDKit Mol."""
    mol = Chem.RWMol()
    atom_map = {}
    
    # Add atoms
    for atom in atomistic.atoms:
        symbol = atom.get("symbol", "C")
        charge = atom.get("charge", 0)
        if charge is None:
            charge = 0
        
        rdkit_atom = Chem.Atom(symbol)
        rdkit_atom.SetFormalCharge(int(charge))
        idx = mol.AddAtom(rdkit_atom)
        atom_map[atom] = idx
    
    # Add bonds
    for bond in atomistic.bonds:
        a, b = bond.endpoints
        begin_idx = atom_map[a]
        end_idx = atom_map[b]
        
        order = bond.get("order", 1)
        if order == 1:
            bond_type = Chem.BondType.SINGLE
        elif order == 2:
            bond_type = Chem.BondType.DOUBLE
        elif order == 3:
            bond_type = Chem.BondType.TRIPLE
        elif order == 1.5:
            bond_type = Chem.BondType.AROMATIC
        else:
            bond_type = Chem.BondType.SINGLE
        
        mol.AddBond(begin_idx, end_idx, bond_type)
    
    # Add coordinates if available
    has_coords = all("xyz" in atom or "pos" in atom for atom in atomistic.atoms)
    if has_coords:
        conf = Chem.Conformer(len(atomistic.atoms))
        for atom, idx in atom_map.items():
            xyz = atom.get("xyz", atom.get("pos", [0, 0, 0]))
            conf.SetAtomPosition(idx, xyz)
        mol.AddConformer(conf)
    
    return mol.GetMol()


def monomer_to_mol(monomer: Monomer) -> Chem.Mol:
    """Convert Monomer to RDKit Mol."""
    # Monomer wraps Atomistic, so just convert the wrapped object
    return atomistic_to_mol(monomer)


def draw_molecule(
    mol: Chem.Mol,
    width: int = 400,
    height: int = 300,
    highlight_atoms: Optional[list] = None,
) -> str:
    """Draw molecule and return SVG string."""
    drawer = Draw.rdMolDraw2D.MolDraw2DSVG(width, height)
    
    if highlight_atoms:
        drawer.DrawMolecule(mol, highlightAtoms=highlight_atoms)
    else:
        drawer.DrawMolecule(mol)
    
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def generate_3d_coords(
    monomer: Monomer,
    add_hydrogens: bool = True,
    optimize: bool = True,
) -> Monomer:
    """Generate 3D coordinates for monomer."""
    # Convert to RDKit mol
    mol = monomer_to_mol(monomer)
    
    # Add hydrogens if requested
    if add_hydrogens:
        mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(mol, randomSeed=42)
    
    # Optimize if requested
    if optimize:
        AllChem.MMFFOptimizeMolecule(mol)
    
    # Convert back to atomistic, then wrap in monomer
    atomistic = mol_to_atomistic(mol)
    return Monomer(atomistic)
