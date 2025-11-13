from __future__ import annotations

from typing import Any, Callable
from typing_extensions import Concatenate, ParamSpec, TypeVar

from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import display, SVG

from molpy.core.atomistic import Atomistic
from molpy.core.wrappers.monomer import Monomer

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

# ------------------------- drawing -------------------------

def _accepts_rdkit_src(func: Callable[Concatenate[Chem.Mol, P], R]) -> Callable[Concatenate[object, P], R]:
    def wrapper(src: object, /, *args: P.args, **kwargs: P.kwargs) -> R:
        if isinstance(src, Chem.Mol):
            mol = src
        elif isinstance(src, str):
            mol = Chem.MolFromSmiles(src)
            if mol is None:
                raise ValueError(f"Invalid SMILES: {src}")
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
