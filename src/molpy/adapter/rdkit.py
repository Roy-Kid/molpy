from .registry import REG, convert
from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor
from IPython.display import display, SVG
from typing import Callable
from typing_extensions import Concatenate, ParamSpec, TypeVar
from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer
from rdkit.Chem.Draw import rdMolDraw2D

P = ParamSpec("P")
R = TypeVar("R")


def mol_to_atomistic(mol: Chem.Mol) -> "Atomistic":
    """
    Convert RDKit Mol to Atomistic representation.
    
    Args:
        mol: RDKit Mol object
        
    Returns:
        Atomistic object with atoms and bonds
        
    Example:
        >>> from rdkit import Chem
        >>> mol = Chem.MolFromSmiles("CCO")
        >>> atomistic = mol_to_atomistic(mol)
    """
    from molpy.core.atomistic import Atomistic
    
    atomistic = Atomistic()
    
    # Store atom mapping: RDKit atom idx -> Atomistic Atom object
    atom_map = {}
    
    # Add atoms
    for rdkit_atom in mol.GetAtoms():
        idx = rdkit_atom.GetIdx()
        symbol = rdkit_atom.GetSymbol()
        
        # Get 3D coordinates if available
        conformer = mol.GetConformer() if mol.GetNumConformers() > 0 else None
        if conformer is not None:
            pos = conformer.GetAtomPosition(idx)
            atom = atomistic.add_atom(
                symbol=symbol,
                pos=[pos.x, pos.y, pos.z],
                atomic_num=rdkit_atom.GetAtomicNum(),
                formal_charge=rdkit_atom.GetFormalCharge(),
            )
        else:
            atom = atomistic.add_atom(
                symbol=symbol,
                atomic_num=rdkit_atom.GetAtomicNum(),
                formal_charge=rdkit_atom.GetFormalCharge(),
            )
        
        atom_map[idx] = atom
    
    # Add bonds
    for rdkit_bond in mol.GetBonds():
        begin_idx = rdkit_bond.GetBeginAtomIdx()
        end_idx = rdkit_bond.GetEndAtomIdx()
        bond_type = rdkit_bond.GetBondType()
        
        # Map RDKit bond type to bond order
        bond_order_map = {
            Chem.BondType.SINGLE: 1,
            Chem.BondType.DOUBLE: 2,
            Chem.BondType.TRIPLE: 3,
            Chem.BondType.AROMATIC: 1.5,
        }
        
        atomistic.add_bond(
            atom_map[begin_idx],
            atom_map[end_idx],
            order=bond_order_map.get(bond_type, 1),
            type=str(bond_type),
        )
    
    return atomistic

REG.register(Chem.Mol, Atomistic, mol_to_atomistic)


def atomistic_to_mol(atomistic: Atomistic) -> Chem.Mol:
    """
    Convert Atomistic to RDKit Mol.
    
    Args:
        atomistic: Atomistic object with atoms and bonds
        
    Returns:
        RDKit Mol object
        
    Example:
        >>> from molpy.core.atomistic import Atomistic
        >>> atomistic = Atomistic()
        >>> c1 = atomistic.add_atom(symbol="C")
        >>> c2 = atomistic.add_atom(symbol="C")
        >>> atomistic.add_bond(c1, c2, order=1)
        >>> mol = atomistic_to_mol(atomistic)
    """
    from rdkit import Chem
    
    # Create editable mol
    mol = Chem.RWMol()
    
    # Store atom mapping: Atomistic Atom -> RDKit atom idx
    atom_map = {}
    
    # Add atoms
    for i, atom in enumerate(atomistic.atoms):
        rdkit_atom = Chem.Atom(atom.data.get("symbol", "C"))
        
        # Set properties if available
        if "atomic_num" in atom.data:
            rdkit_atom.SetAtomicNum(atom.data["atomic_num"])
        if "formal_charge" in atom.data:
            rdkit_atom.SetFormalCharge(atom.data["formal_charge"])
        
        idx = mol.AddAtom(rdkit_atom)
        mol.GetAtomWithIdx(idx).SetIntProp("molpy_idx", i)
        atom_map[atom] = idx
    
    # Add bonds
    for bond in atomistic.bonds:
        atom1 = bond.itom
        atom2 = bond.jtom
        
        order = bond.data.get("order", 1)
        
        # Map bond order to RDKit bond type
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
        
        mol.AddBond(atom_map[atom1], atom_map[atom2], bond_type)
    
    # Add coordinates if available
    if any("pos" in atom.data for atom in atomistic.atoms):
        from rdkit.Chem import AllChem
        conformer = Chem.Conformer(len(atom_map))
        
        for atom, idx in atom_map.items():
            pos = atom.data.get("pos", [0.0, 0.0, 0.0])
            conformer.SetAtomPosition(idx, pos)
        
        mol.AddConformer(conformer)
    
    return mol.GetMol()


REG.register(Atomistic, Chem.Mol, atomistic_to_mol)


def monomer_to_mol(monomer: "Monomer") -> Chem.Mol:
    """
    Convert Monomer to RDKit Mol by unwrapping to Atomistic first.
    
    Args:
        monomer: Monomer wrapper object
        
    Returns:
        RDKit Mol object
    """
    from molpy.core.wrappers.monomer import Monomer as MonomerCls
    
    if not isinstance(monomer, MonomerCls):
        raise TypeError(f"Expected Monomer, got {type(monomer)}")
    
    atomistic = monomer.unwrap()
    return atomistic_to_mol(atomistic)


# Register via predicate to avoid circular import at module level
from molpy.core.wrappers.monomer import Monomer as MonomerCls
REG.register(MonomerCls, Chem.Mol, monomer_to_mol)


def _accepts_rdkit_src(
    func: Callable[Concatenate[Chem.Mol, P], R]
) -> Callable[Concatenate[object, P], R]:
    def wrapper(src: object, /, *args: P.args, **kwargs: P.kwargs) -> R:
        if isinstance(src, Chem.Mol):
            mol = src
        elif isinstance(src, str):
            mol = Chem.MolFromSmiles(src)
            if mol is None:
                raise ValueError(f"Invalid SMILES string: {src}")
        else:
            mol = convert(src, Chem.Mol)
        return func(mol, *args, **kwargs)
    return wrapper


@_accepts_rdkit_src
def draw_molecule(
    mol: Chem.Mol,
    *,
    size: tuple[int, int] = (320, 260),
    show_indices: bool = True,
    show_explicit_H: bool = False,
    kekulize: bool = True,
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
    title: str | None = None,
    show: bool = True,
):
    """
    RDKit SVG drawing with nicer layout:
    - rdCoordGen 2D coords
    - extra padding so labels/indices don't clip
    - fixed and slightly larger font
    - optional: convert implicit H to explicit for display
    """
    # 0) coords
    if not mol.GetNumConformers():
        rdDepictor.SetPreferCoordGen(True)  # use rdCoordGen backend
        rdDepictor.Compute2DCoords(mol)

    # 1) make a copy to draw
    dm = Chem.Mol(mol)
    if show_explicit_H:
        dm = Chem.AddHs(dm)
        # 重新布局一次，避免加H后重叠
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(dm)

    # 2) prepare options
    w, h = size
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = drawer.drawOptions()

    # —— 关键：给足边距，避免裁剪（默认≈0.05，太小）
    opts.padding = 0.10  # 0.10–0.18 之间可调
    # 原子标签周围再给点额外留白
    opts.additionalAtomLabelPadding = 0.06

    # 字体/线宽（固定值更稳定，不会因为缩放挤掉 H/C）
    opts.fixedFontSize = 13
    opts.minFontSize = 9
    opts.bondLineWidth = 2

    # 索引：不覆盖元素符号（只打开开关）
    opts.addAtomIndices = bool(show_indices)

    # 可读性
    # opts.useBWAtomPalette(False)         # 彩色元素
    opts.addStereoAnnotation = False
    opts.explicitMethyl = True           # -CH3 画成 C-H3 更清楚

    # 3) draw
    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        dm,
        highlightAtoms=highlight_atoms or [],
        highlightBonds=highlight_bonds or [],
        legend=(title or "")
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    if show:
        display(SVG(svg))
    return svg

  

to_rdkit = convert[Chem.Mol]


def generate_3d_coords(
    monomer: "Monomer",
    *,
    optimize: bool = True,
    random_seed: int = 42,
    force_tol: float = 0.001,
    max_iters: int = 200,
) -> "Monomer":
    """
    Generate 3D coordinates for a Monomer using RDKit's ETKDG method.
    
    This function:
    1. Converts Monomer -> RDKit Mol
    2. Generates 3D coordinates via ETKDG
    3. Optionally optimizes geometry with MMFF94
    4. Binds coordinates back to the original Monomer atoms
    5. Returns the Monomer with updated positions
    
    Args:
        monomer: Monomer object to generate 3D coords for
        optimize: Whether to optimize geometry with MMFF94 (default: True)
        random_seed: Random seed for reproducibility (default: 42)
        force_tol: Force field convergence tolerance (default: 0.001)
        max_iters: Max iterations for optimization (default: 200)
        
    Returns:
        The same Monomer object with updated 3D coordinates
        
    Raises:
        RuntimeError: If 3D coordinate generation fails
        
    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> from molpy.core.wrappers.monomer import Monomer
        >>> parser = SmilesParser()
        >>> smiles_ir = parser.parse_smiles("CCO")
        >>> atomistic = smiles_ir.to_atomistic()
        >>> monomer = Monomer(atomistic)
        >>> monomer_3d = generate_3d_coords(monomer)
        >>> # Now monomer has 3D coordinates
    """
    from rdkit.Chem import AllChem
    from molpy.core.wrappers.monomer import Monomer as MonomerCls
    
    if not isinstance(monomer, MonomerCls):
        raise TypeError(f"Expected Monomer, got {type(monomer)}")
    
    # Convert Monomer -> RDKit Mol
    mol = convert(monomer, Chem.Mol)
    
    # Sanitize the molecule (calculate implicit valence, aromaticity, etc.)
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        raise RuntimeError(f"Failed to sanitize molecule: {e}")
    
    # Store atom mapping: Atomistic Atom -> RDKit atom idx
    # We need to track which atom in the monomer corresponds to which RDKit atom
    atomistic = monomer.unwrap()
    atom_list = list(atomistic.atoms)
    
    # Add hydrogens if needed
    mol = Chem.AddHs(mol)
    
    # Generate 3D coordinates using ETKDG
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    
    result = AllChem.EmbedMolecule(mol, params)
    if result == -1:
        raise RuntimeError("Failed to generate 3D coordinates for molecule")
    
    # Optimize geometry with MMFF94 if requested
    if optimize:
        try:
            AllChem.MMFFOptimizeMolecule(
                mol, 
                maxIters=max_iters,
                nonBondedThresh=force_tol,
            )
        except Exception as e:
            # MMFF might fail for some molecules, continue with unoptimized coords
            import warnings
            warnings.warn(f"MMFF optimization failed: {e}. Using unoptimized coordinates.")
    
    # Get the conformer with 3D coordinates
    conformer = mol.GetConformer()
    
    # Bind coordinates back to original Monomer atoms
    # Note: RDKit might have added explicit hydrogens, so we only update
    # the heavy atoms that were in the original structure
    for i, atom in enumerate(atom_list):
        if i >= mol.GetNumAtoms():
            break
        pos = conformer.GetAtomPosition(i)
        atom.data["pos"] = [pos.x, pos.y, pos.z]
        atom.data["xyz"] = [pos.x, pos.y, pos.z]  # Also set xyz for compatibility
    
    return monomer