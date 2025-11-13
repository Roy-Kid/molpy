"""通用分子表示转换器。"""

from __future__ import annotations

from typing import Any, Callable, TypeVar, overload

from rdkit import Chem
from rdkit.Chem.rdchem import ChiralType

from molpy.core.atomistic import Atomistic
from molpy.parser.smiles import SmilesIR

T = TypeVar("T")
U = TypeVar("U")


def smilesir_to_mol(ir: SmilesIR) -> Chem.Mol:
    """将 `SmilesIR` 转换为 RDKit `Chem.Mol`（内部工具）。"""
    mol = Chem.RWMol()
    atom_map: dict[int, int] = {}
    aromatic_flags: dict[int, bool] = {}

    for atom_ir in ir.atoms:
        symbol = atom_ir.symbol or "C"
        # Aromatic atoms are encoded as lower-case in SMILES
        is_aromatic = symbol.islower()
        if len(symbol) == 1:
            rd_symbol = symbol.upper()
        else:
            rd_symbol = symbol[0].upper() + symbol[1:].lower()

        rdkit_atom = Chem.Atom(rd_symbol)
        if is_aromatic:
            rdkit_atom.SetIsAromatic(True)

        if atom_ir.charge is not None:
            rdkit_atom.SetFormalCharge(int(atom_ir.charge))
        if getattr(atom_ir, "isotope", None) is not None:
            rdkit_atom.SetIsotope(int(atom_ir.isotope))
        if getattr(atom_ir, "h_count", None) is not None:
            rdkit_atom.SetNumExplicitHs(int(atom_ir.h_count))
        chiral = getattr(atom_ir, "chiral", None)
        if isinstance(chiral, str):
            if chiral == "@":
                rdkit_atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CCW)
            elif chiral == "@@":
                rdkit_atom.SetChiralTag(ChiralType.CHI_TETRAHEDRAL_CW)

        idx = mol.AddAtom(rdkit_atom)
        key = id(atom_ir)
        atom_map[key] = idx
        aromatic_flags[key] = is_aromatic

    bond_type_map: dict[str, Chem.BondType] = {
        "-": Chem.BondType.SINGLE,
        "=": Chem.BondType.DOUBLE,
        "#": Chem.BondType.TRIPLE,
        ":": Chem.BondType.AROMATIC,
    }

    for bond_ir in ir.bonds:
        start_key = id(bond_ir.start)
        end_key = id(bond_ir.end)
        begin = atom_map.get(start_key)
        end = atom_map.get(end_key)
        if begin is None or end is None or begin == end:
            continue
        if mol.GetBondBetweenAtoms(begin, end) is not None:
            continue

        bond_symbol = bond_ir.bond_type or "-"
        # If both atoms are aromatic and no explicit aromatic bond is set, use aromatic bond
        start_is_aromatic = aromatic_flags.get(start_key, False)
        end_is_aromatic = aromatic_flags.get(end_key, False)
        if bond_symbol == "-" and start_is_aromatic and end_is_aromatic:
            bond_symbol = ":"

        bond_type = bond_type_map.get(bond_symbol, Chem.BondType.SINGLE)
        mol.AddBond(begin, end, bond_type)

    try:
        Chem.SanitizeMol(mol)
    except Exception:
        pass

    return mol.GetMol()


class ConverterRegistry:
    """类型到类型的转换注册表。"""

    def __init__(self) -> None:
        self._converters: dict[tuple[type[Any], type[Any]], Callable[[Any], Any]] = {}

    def register(self, source: type[Any], target: type[Any], converter: Callable[[Any], Any]) -> None:
        self._converters[(source, target)] = converter

    def resolve(self, source_obj: object, target: type[Any]) -> Callable[[Any], Any] | None:
        for source_type in type(source_obj).__mro__:
            converter = self._converters.get((source_type, target))
            if converter is not None:
                return converter
        return None


REG = ConverterRegistry()


def register(source: type[Any], target: type[Any]) -> Callable[[Callable[[Any], Any]], Callable[[Any], Any]]:
    """将函数注册到转换表的装饰器。"""

    def decorator(func: Callable[[Any], Any]) -> Callable[[Any], Any]:
        REG.register(source, target, func)
        return func

    return decorator


@overload
def convert(obj: T, target: type[T]) -> T: ...


@overload
def convert(obj: T, target: type[U]) -> U: ...


def convert(obj: T, target: type[U]) -> U:
    """根据注册表执行类型转换。"""
    if isinstance(obj, target):  # type: ignore[arg-type]
        return obj  # type: ignore[return-value]

    converter = REG.resolve(obj, target)
    if converter is not None:
        return converter(obj)  # type: ignore[return-value]

    raise TypeError(f"无法将 {type(obj)} 转换为 {target}")


@register(SmilesIR, Atomistic)
def _smilesir_to_atomistic(ir: SmilesIR) -> Atomistic:
    from molpy.adapter.rdkit_adapter import RDKitWrapper

    wrapper = RDKitWrapper.from_mol(smilesir_to_mol(ir))
    return wrapper.to_atomistic(sync_coords=True)

