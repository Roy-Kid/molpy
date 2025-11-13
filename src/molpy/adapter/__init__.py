from .converter import REG, convert, smilesir_to_mol
from .rdkit_adapter import RDKitWrapper

__all__ = [
    "RDKitWrapper",
    "convert",
    "smilesir_to_mol",
    "REG",
]

# Register parser converters (after imports to avoid circular dependency)
def _register_parser_converters():
    """Register converters from parser module."""
    try:
        from molpy.parser.smiles import (
            BigSmilesIR,
            SmilesIR,
            bigsmilesir_to_monomer,
            bigsmilesir_to_polymerspec,
            smilesir_to_mol as parser_smilesir_to_mol,
        )
        from molpy.core.wrappers.monomer import Monomer
        from molpy.core.atomistic import Atomistic
        from rdkit import Chem
        
        # Register SmilesIR -> Mol conversion
        REG.register(SmilesIR, Chem.Mol, parser_smilesir_to_mol)
        
        # Register BigSmilesIR conversions
        REG.register(BigSmilesIR, Monomer[Atomistic], bigsmilesir_to_monomer)
        
        # Note: PolymerSpec registration would go here if needed
        # REG.register(BigSmilesIR, PolymerSpec, bigsmilesir_to_polymerspec)
        
    except ImportError:
        pass  # Parser or dependencies not available

# Auto-register on import
_register_parser_converters()