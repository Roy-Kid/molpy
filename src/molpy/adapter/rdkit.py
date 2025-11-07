from .registry import REG
from rdkit import Chem


class to_rdkit:

    def __call__(self, src: object) -> Chem.Mol:
        converter = REG.resolve(src, Chem.Mol)
        if converter is None:
            raise TypeError(f"No converter registered for {type(src)} -> Chem.Mol")
        return converter(src)
    
    def __class_getitem__(cls, src: type) -> type:
        return REG.get_converter(src, Chem.Mol)
    

class from_rdkit:

    def __call__(self, src: Chem.Mol) -> object:
        converter = REG.resolve(src, object)
        if converter is None:
            raise TypeError(f"No converter registered for Chem.Mol -> {object}")
        return converter(src)
    
    def __class_getitem__(cls, dst: type) -> type:
        return REG.get_converter(Chem.Mol, dst)


def draw_molecule(
    src,
    *,
    size: tuple[int, int] = (200, 200),
    show_atom_idx: bool = False,
    highlight_atoms: list[int] | None = None,
    highlight_bonds: list[int] | None = None,
    title: str | None = None,
):
    """
    Draw molecular structure using RDKit.
    
    Args:
        src: SmilesIR, RDKit Mol, or SMILES string
        size: Image size (width, height)
        show_atom_idx: Show atom indices on atoms
        highlight_atoms: List of atom indices to highlight
        highlight_bonds: List of bond indices to highlight
        
    Returns:
        PIL Image object (displays inline in Jupyter notebooks)
        
    Example:
        >>> from molpy.parser.smiles import SmilesParser
        >>> parser = SmilesParser()
        >>> ir = parser.parser_smiles("c1ccccc1")
        >>> img = draw_molecule(ir, show_atom_idx=True)
    """
    from rdkit.Chem import Draw
    
    # Convert to Mol if needed
    if isinstance(src, str):
        mol = Chem.MolFromSmiles(src)
        if mol is None:
            raise ValueError(f"Invalid SMILES string: {src}")
    elif isinstance(src, Chem.Mol):
        mol = src
    else:
        # Try to convert using registered converter
        mol = to_rdkit()(src)
    
    # Add atom labels if requested
    if show_atom_idx:
        for atom in mol.GetAtoms():
            atom.SetProp('atomLabel', str(atom.GetIdx()))
    
    # Draw using RDKit's built-in method
    return Draw.MolToImage(
        mol,
        size=size,
        highlightAtoms=highlight_atoms or [],
        highlightBonds=highlight_bonds or [],
        legend=title or "",
    )

