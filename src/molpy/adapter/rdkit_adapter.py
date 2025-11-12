from __future__ import annotations

from typing import Any, Callable, TYPE_CHECKING
from typing_extensions import Concatenate, ParamSpec, TypeVar

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG

from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer

# Avoid circular import: delay SmilesIR import
from molpy.parser.smiles import SmilesIR

P = ParamSpec("P")
R = TypeVar("R")

# ------------------------- helpers & constants -------------------------

# 用于“回传映射”的稳定标签：优先使用实体内置 id，否则用添加顺序
MP_ID = "mp_id"

_BOND_ORDER_TO_RDKIT: dict[float, Chem.BondType] = {
    1.0: Chem.BondType.SINGLE,
    2.0: Chem.BondType.DOUBLE,
    3.0: Chem.BondType.TRIPLE,
    1.5: Chem.BondType.AROMATIC,  # 仅用于可视化/近似，真正芳香度靠 RDKit 判定
}
_RDKIT_TO_BOND_ORDER: dict[Chem.BondType, float] = {
    Chem.BondType.SINGLE: 1.0,
    Chem.BondType.DOUBLE: 2.0,
    Chem.BondType.TRIPLE: 3.0,
    Chem.BondType.AROMATIC: 1.5,
}

def _ensure_2d(mol: Chem.Mol) -> None:
    """Compute 2D coordinates in-place if not present."""
    if not mol.GetNumConformers():
        rdDepictor.SetPreferCoordGen(True)
        rdDepictor.Compute2DCoords(mol)

def _label_rdkit_heavy_atoms(rwmol: Chem.RWMol, mon_atoms: list[Any]) -> None:
    """
    给 RDKit 重原子打回传标签，用于把坐标/氢回填回原结构。
    假设 rwmol 中的前 len(mon_atoms) 个是重原子，顺序与构建时一致。
    """
    n_heavy = len(mon_atoms)
    for i in range(n_heavy):
        # 优先用实体自带 id（如果有）
        ent = mon_atoms[i]
        ent_id = ent.get("id", i)
        rwmol.GetAtomWithIdx(i).SetIntProp(MP_ID, int(ent_id))

def _build_mon_heavy_index(mon_atoms: list[Any]) -> dict[int, Any]:
    """根据实体 id 构建：heavy_id -> monomer 的重原子实体。"""
    table: dict[int, Any] = {}
    for i, ent in enumerate(mon_atoms):
        hid = ent.get("id", i)
        table[int(hid)] = ent
    return table

def _rdkit_bond_type(order: float) -> Chem.BondType:
    return _BOND_ORDER_TO_RDKIT.get(float(order), Chem.BondType.SINGLE)

def _order_from_rdkit(bt: Chem.BondType) -> float:
    return _RDKIT_TO_BOND_ORDER.get(bt, 1.0)

# ------------------------- converters -------------------------
def mol_to_atomistic(mol: Chem.Mol) -> Atomistic:
    """
    RDKit Mol -> Atomistic
    - 保留（若有）坐标
    - 复制 formal_charge、atomic_num
    - bond order 做基础映射（1/2/3/1.5）
    """
    atomistic = Atomistic()
    atom_map: dict[int, Any] = {}

    conf = mol.GetConformer(0) if mol.GetNumConformers() else None

    # atoms
    for a in mol.GetAtoms():
        idx = a.GetIdx()
        data: dict[str, Any] = {
            "symbol": a.GetSymbol(),
            "atomic_num": a.GetAtomicNum(),
            "formal_charge": a.GetFormalCharge(),
        }
        if conf is not None:
            p = conf.GetAtomPosition(idx)
            data["pos"] = [float(p.x), float(p.y), float(p.z)]
        atom = atomistic.add_atom(**data)
        atom_map[idx] = atom

    # bonds
    for b in mol.GetBonds():
        i = b.GetBeginAtomIdx()
        j = b.GetEndAtomIdx()
        order = _order_from_rdkit(b.GetBondType())
        atomistic.add_bond(atom_map[i], atom_map[j], order=order, type=str(b.GetBondType()))

    return atomistic


def atomistic_to_mol(atomistic: Atomistic) -> Chem.Mol:
    """
    Atomistic -> RDKit Mol
    - 给每个 RDKit 重原子设置 MP_ID，便于往返映射
    - 如存在 'pos'，写入 conformer
    """
    rwmol = Chem.RWMol()
    atom_map: dict[Any, int] = {}
    mon_atoms = list(atomistic.atoms)

    # atoms
    for i, ent in enumerate(mon_atoms):
        sym = ent.get("symbol", "C")
        rd_atom = Chem.Atom(sym)
        if "atomic_num" in ent:
            rd_atom.SetAtomicNum(int(ent["atomic_num"]))
        if "formal_charge" in ent:
            rd_atom.SetFormalCharge(int(ent["formal_charge"]))

        ridx = rwmol.AddAtom(rd_atom)
        rwmol.GetAtomWithIdx(ridx).SetIntProp(MP_ID, int(ent.get("id", i)))
        atom_map[ent] = ridx

    # bonds
    for b in atomistic.bonds:
        i = atom_map[b.itom]
        j = atom_map[b.jtom]
        bt = _rdkit_bond_type(b.get("order", 1))
        rwmol.AddBond(i, j, bt)

    mol = rwmol.GetMol()

    # coordinates
    if any("pos" in ent for ent in mon_atoms):
        conf = Chem.Conformer(len(mon_atoms))
        for ent, ridx in atom_map.items():
            x, y, z = ent.get("pos", [0.0, 0.0, 0.0])
            conf.SetAtomPosition(ridx, (float(x), float(y), float(z)))
        mol.AddConformer(conf, assignId=True)

    return mol


def monomer_to_mol(monomer: Monomer) -> Chem.Mol:
    """Monomer -> RDKit Mol（unwrap 后复用 atomistic_to_mol）"""
    return atomistic_to_mol(monomer.unwrap())


def smilesir_to_mol(ir: "SmilesIR") -> Chem.Mol:
    """
    Convert SmilesIR to RDKit Mol by directly constructing the molecule graph.

    This approach preserves IR-specific information and supports extended syntax
    (BigSMILES, G-BigSMILES) where explicit topology is essential.

    Args:
        ir: SmilesIR instance with atoms and bonds

    Returns:
        RDKit Mol object

    Raises:
        ValueError: if IR contains invalid molecular data

    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> parser = SmilesParser()
        >>> ir = parser.parse_smiles("CCO")
        >>> mol = smilesir_to_mol(ir)
        >>> mol.GetNumAtoms()
        3
    """
    # Import here to avoid circular dependency
    from molpy.parser.smiles import SmilesIR, AtomIR
    
    assert isinstance(ir, SmilesIR), "Input must be a SmilesIR instance"

    if not ir.atoms:
        # Empty molecule
        return Chem.Mol()

    # Bond type mapping
    bond_type_map = {
        "-": Chem.BondType.SINGLE,
        "=": Chem.BondType.DOUBLE,
        "#": Chem.BondType.TRIPLE,
        ":": Chem.BondType.AROMATIC,
        "/": Chem.BondType.SINGLE,  # Stereochemistry, treat as single for now
        "\\": Chem.BondType.SINGLE,  # Stereochemistry, treat as single for now
    }

    # Create editable molecule
    mol = Chem.RWMol()

    # Map AtomIR -> RDKit atom index (using object identity)
    atom_to_idx: dict[int, int] = {}

    # Add atoms
    for atom_ir in ir.atoms:
        # Skip non-AtomIR entities (e.g., BondDescriptorIR in BigSMILES)
        if not isinstance(atom_ir, AtomIR):
            continue
            
        # Handle aromatic symbols (lowercase in SMILES → uppercase + aromatic flag)
        symbol = atom_ir.symbol.upper() if atom_ir.symbol.islower() else atom_ir.symbol
        is_aromatic = atom_ir.symbol.islower()

        # Create RDKit atom
        rdkit_atom = Chem.Atom(symbol)

        # Set properties
        if atom_ir.charge is not None:
            rdkit_atom.SetFormalCharge(atom_ir.charge)

        if atom_ir.isotope is not None:
            rdkit_atom.SetIsotope(atom_ir.isotope)

        if atom_ir.h_count is not None:
            rdkit_atom.SetNumExplicitHs(atom_ir.h_count)

        # Handle chirality
        if atom_ir.chiral is not None:
            if atom_ir.chiral == "@":
                rdkit_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CCW)
            elif atom_ir.chiral == "@@":
                rdkit_atom.SetChiralTag(Chem.ChiralType.CHI_TETRAHEDRAL_CW)
            # Other chiral tags can be added as needed

        # Set aromaticity
        if is_aromatic:
            rdkit_atom.SetIsAromatic(True)

        # Add atom and store mapping (use id() for object identity)
        atom_idx = mol.AddAtom(rdkit_atom)
        atom_to_idx[id(atom_ir)] = atom_idx

    # Add bonds
    for bond_ir in ir.bonds:
        # Skip bonds involving non-AtomIR entities
        from molpy.parser.smiles import AtomIR
        if not (isinstance(bond_ir.start, AtomIR) and isinstance(bond_ir.end, AtomIR)):
            continue
            
        start_idx = atom_to_idx.get(id(bond_ir.start))
        end_idx = atom_to_idx.get(id(bond_ir.end))

        if start_idx is None or end_idx is None:
            continue  # Skip bonds to filtered atoms

        # Determine bond type (upgrade single bonds between aromatic atoms to aromatic)
        bond_type_str = bond_ir.bond_type
        if (
            bond_type_str == "-"
            and bond_ir.start.symbol.islower()
            and bond_ir.end.symbol.islower()
        ):
            # Single bond between aromatic atoms → aromatic bond
            bond_type = Chem.BondType.AROMATIC
        else:
            bond_type = bond_type_map.get(bond_type_str)
            if bond_type is None:
                raise ValueError(f"Unknown bond type: {bond_type_str}")

        mol.AddBond(start_idx, end_idx, bond_type)

    # Convert to immutable Mol
    final_mol = mol.GetMol()

    # Sanitize molecule (compute aromaticity, implicit Hs, etc.)
    try:
        Chem.SanitizeMol(final_mol)
    except Exception as e:
        # If sanitization fails, return unsanitized molecule with warning
        import warnings

        warnings.warn(
            f"Molecule sanitization failed: {e}. Returning unsanitized molecule."
        )

    return final_mol


def bigsmilesir_to_mol(ir: "BigSmilesIR") -> Chem.Mol:
    """
    Convert BigSmilesIR to RDKit Mol (chemical structure only, no BigSMILES markers).

    Uses the degenerate() method to strip bond descriptors and get clean chemistry.

    Args:
        ir: BigSmilesIR instance

    Returns:
        RDKit Mol object with only the chemical structure

    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> parser = SmilesParser()
        >>> ir = parser.parse_bigsmiles("{[<]CC[>]}")
        >>> mol = bigsmilesir_to_mol(ir)
        >>> mol.GetNumAtoms()  # Only 2 carbons, no descriptors
        2
    """
    from molpy.parser.smiles import BigSmilesIR
    
    # Get clean chemical structure (no bond descriptors)
    clean_ir = ir.degenerate()
    return smilesir_to_mol(clean_ir)


def bigsmilesir_to_monomer(ir: "BigSmilesIR") -> Monomer[Atomistic]:
    """
    Convert BigSmilesIR to Monomer with 3D coordinates and hydrogens.

    This is a convenience function that:
    1. Extracts monomer topology from BigSmilesIR
    2. Converts to RDKit Mol
    3. Generates 3D coordinates with ETKDG
    4. Adds explicit hydrogens
    5. Transfers coordinates back to Monomer

    Args:
        ir: BigSmilesIR from parser (must contain exactly ONE repeat unit)

    Returns:
        Monomer[Atomistic] with:
        - Ports set from bond descriptors
        - 3D coordinates on all atoms
        - Explicit hydrogens added

    Raises:
        ValueError: If IR contains multiple repeat units

    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> parser = SmilesParser()
        >>> ir = parser.parse_bigsmiles("{[<]CC[>]}")
        >>> monomer = bigsmilesir_to_monomer(ir)
        >>> len(list(monomer.unwrap().atoms))  # Has C + H atoms
        8
        >>> monomer.port_names()
        ['in', 'out']
    """
    from molpy.parser.smiles import bigsmilesir_to_monomer as extract_monomer
    
    # Extract topology-only monomer
    monomer = extract_monomer(ir)
    
    # Generate 3D coordinates and add hydrogens
    return generate_3d_coords(monomer)


def bigsmilesir_to_polymerspec(ir: "BigSmilesIR") -> "PolymerSpec":
    """
    Convert BigSmilesIR to PolymerSpec with 3D monomers.

    This is a convenience function that:
    1. Extracts polymer specification from BigSmilesIR
    2. For each monomer, generates 3D coordinates and adds hydrogens

    Args:
        ir: BigSmilesIR from parser

    Returns:
        PolymerSpec with all monomers having 3D coordinates and explicit H

    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> parser = SmilesParser()
        >>> ir = parser.parse_bigsmiles("{[<]CC[>]}{[<]OCC[>]}")
        >>> spec = bigsmilesir_to_polymerspec(ir)
        >>> spec.topology
        'block_copolymer'
        >>> len(spec.all_monomers)
        2
    """
    from molpy.parser.smiles import bigsmilesir_to_polymerspec as extract_spec
    
    # Extract topology-only spec
    spec = extract_spec(ir)
    
    # Generate 3D coords for all monomers
    for segment in spec.segments:
        for i, monomer in enumerate(segment.monomers):
            segment.monomers[i] = generate_3d_coords(monomer)
        for i, eg_monomer in enumerate(segment.end_groups):
            segment.end_groups[i] = generate_3d_coords(eg_monomer)
    
    # Regenerate all_monomers from updated segments
    spec.all_monomers = [
        monomer for segment in spec.segments for monomer in segment.monomers
    ]
    
    return spec


# Note: SmilesIR/BigSmilesIR converter functions are defined here
# to avoid circular import issues

# ------------------------- drawing -------------------------

@_accepts_rdkit_src
def draw_molecule(
    mol: Chem.Mol,
    *,
    size: tuple[int, int] = (320, 260),
    show_indices: bool = True,
    show_explicit_H: bool = False,
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
    title: str | None = None,
    show: bool = True,
) -> str:
    """
    以 SVG 方式绘制，不清空前输出；保留元素符号并在旁标注索引。
    """
    _ensure_2d(mol)
    dm = Chem.Mol(mol)
    if show_explicit_H:
        dm = Chem.AddHs(dm)
        _ensure_2d(dm)

    w, h = size
    drawer = rdMolDraw2D.MolDraw2DSVG(w, h)
    opts = drawer.drawOptions()
    opts.padding = 0.12
    opts.additionalAtomLabelPadding = 0.06
    opts.fixedFontSize = 13
    opts.minFontSize = 9
    opts.bondLineWidth = 2
    opts.addAtomIndices = bool(show_indices)
    opts.addStereoAnnotation = False
    opts.explicitMethyl = True

    rdMolDraw2D.PrepareAndDrawMolecule(
        drawer,
        dm,
        highlightAtoms=highlight_atoms or [],
        highlightBonds=highlight_bonds or [],
        legend=(title or ""),
    )
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    if show:
        display(SVG(svg))
    return svg

# ------------------------- 3D generation & back-transfer -------------------------

def transfer_coords_and_h_to_monomer(monomer: Monomer, mol_with_h: Chem.Mol) -> Monomer:
    """
    将 RDKit 全氢分子的坐标 + 新增氢回填到 Monomer（unwrap 的 Atomistic）：
    - 重原子：通过 MP_ID 做一一映射（不依赖索引顺序）
    - 氢：通过其唯一邻居的 MP_ID 定位要挂接的重原子
    """
    atomistic = monomer.unwrap()
    mon_atoms = list(atomistic.atoms)
    heavy_index = _build_mon_heavy_index(mon_atoms)

    if mol_with_h.GetNumConformers() == 0:
        raise ValueError("RDKit mol has no conformer")
    conf = mol_with_h.GetConformer()

    # 回填重原子坐标
    for a in mol_with_h.GetAtoms():
        if a.GetAtomicNum() == 1:
            continue
        if not a.HasProp(MP_ID):
            # 若对方缺少标签，尝试按序 fallback（不建议发生）
            hid = a.GetIdx()
        else:
            hid = int(a.GetProp(MP_ID))
        ent = heavy_index.get(hid)
        if ent is None:
            raise ValueError(f"No monomer heavy atom for {MP_ID}={hid}")
        p = conf.GetAtomPosition(a.GetIdx())
        ent["pos"] = [float(p.x), float(p.y), float(p.z)]

    # 添加氢并连到其邻居重原子
    for a in mol_with_h.GetAtoms():
        if a.GetAtomicNum() != 1:
            continue
        bonds = list(a.GetBonds())
        if not bonds:
            continue
        other = bonds[0].GetOtherAtom(a)
        if not other.HasProp(MP_ID):
            continue
        hid = int(other.GetProp(MP_ID))
        heavy_ent = heavy_index.get(hid)
        if heavy_ent is None:
            continue

        p = conf.GetAtomPosition(a.GetIdx())
        h_ent = atomistic.add_atom(symbol="H", atomic_num=1, pos=[float(p.x), float(p.y), float(p.z)])
        atomistic.add_bond(heavy_ent, h_ent, order=1.0)

    return monomer

def generate_3d_coords(
    monomer: Monomer,
    *,
    optimize: bool = True,
    random_seed: int = 42,
    max_iters: int = 200,
) -> Monomer:
    """
    使用 RDKit ETKDG + (可选)MMFF 对 Monomer 生成 3D，并把坐标/氢回填。
    - Monomer -> RDKit (重原子) 并打 MP_ID
    - AddHs/Embed/MMFF
    - 回填（坐标 + 新氢）
    """
    if not isinstance(monomer, Monomer):
        raise TypeError(f"Expected Monomer, got {type(monomer)}")

    # 1) Monomer -> RDKit（重原子）
    atomistic = monomer.unwrap()
    mon_atoms = list(atomistic.atoms)
    rwmol = Chem.RWMol()

    ridx_map: dict[Any, int] = {}
    for i, ent in enumerate(mon_atoms):
        rd_atom = Chem.Atom(ent.get("symbol", "C"))
        if "atomic_num" in ent:
            rd_atom.SetAtomicNum(int(ent["atomic_num"]))
        if "formal_charge" in ent:
            rd_atom.SetFormalCharge(int(ent["formal_charge"]))
        ridx = rwmol.AddAtom(rd_atom)
        rwmol.GetAtomWithIdx(ridx).SetIntProp(MP_ID, int(ent.get("id", i)))
        ridx_map[ent] = ridx

    for b in atomistic.bonds:
        i = ridx_map[b.itom]
        j = ridx_map[b.jtom]
        rwmol.AddBond(i, j, _rdkit_bond_type(b.get("order", 1.0)))

    mol = rwmol.GetMol()

    # 2) AddHs + ETKDG + (opt)MMFF
    molH = Chem.AddHs(mol, addCoords=True)
    params = AllChem.ETKDGv3()
    params.randomSeed = int(random_seed)
    # 更稳：允许随机初值
    params.useRandomCoords = True
    if AllChem.EmbedMolecule(molH, params) == -1:
        # 再试一次
        params.useRandomCoords = True
        if AllChem.EmbedMolecule(molH, params) == -1:
            raise RuntimeError("ETKDG embedding failed")

    if optimize:
        try:
            AllChem.MMFFOptimizeMolecule(molH, maxIters=int(max_iters))
        except Exception as e:
            # 某些分子失败也无妨，保留嵌入结构
            import warnings
            warnings.warn(f"MMFF optimization failed: {e}")

    # 3) 回填坐标 + H
    return transfer_coords_and_h_to_monomer(monomer, molH)



# ===================================================================
#   4. Converter: SmartsIR -> RDKit Mol
# ===================================================================

def smartsir_to_mol(ir: SmartsIR) -> "Chem.Mol":
    """
    Convert SmartsIR to RDKit Mol query object.
    
    This creates a molecule object that can be used for substructure searching.
    The resulting Mol has query atoms that match the SMARTS pattern.
    
    Args:
        ir: SmartsIR instance with atoms and bonds
        
    Returns:
        RDKit Mol object configured as a query
        
    Raises:
        ImportError: if RDKit is not available
        ValueError: if IR contains invalid pattern data
        
    Example:
        >>> parser = SmartsParser()
        >>> ir = parser.parse_smarts("[#6]")
        >>> mol = smartsir_to_mol(ir)
        >>> # Use mol for substructure search
    """
    
    if not ir.atoms:
        return Chem.Mol()
    
    # For now, convert to SMARTS string and use RDKit's parser
    # TODO: Direct construction from IR for better control
    smarts_str = _ir_to_smarts_string(ir)
    mol = Chem.MolFromSmarts(smarts_str)
    
    if mol is None:
        raise ValueError(f"Failed to create RDKit Mol from SMARTS: {smarts_str}")
    
    return mol
