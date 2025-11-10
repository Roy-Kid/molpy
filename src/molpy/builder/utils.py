"""
Helper function to convert BigSMILES/SMILES to Monomer for polymer assembly.
"""

from ..core.wrappers.monomer import Monomer
from ..core.atomistic import Atomistic
from ..parser.smiles import SmilesParser


def bigsmiles_to_monomer(smiles: str) -> Monomer[Atomistic]:
    """
    Convert BigSMILES or plain SMILES with port markers to Monomer.
    
    Supports:
    1. BigSMILES with stochastic objects: {[<]CC[>]} 
    2. Plain SMILES with atom class ports: CCCCO[*:1]
    
    Port roles from BigSMILES:
    - [<] creates port with role='left'
    - [>] creates port with role='right'
    - [*:n] creates port with name='port_n' or just 'n'
    
    Args:
        smiles: BigSMILES or SMILES string
        
    Returns:
        Monomer with ports set (topology only, no coordinates)
        
    Example:
        >>> mono = bigsmiles_to_monomer("[<]CCO[>]")
        >>> mono.port_names()
        ['left', 'right']  # or ['in', 'out'] depending on implementation
    """
    from ..parser.smiles import bigsmilesir_to_monomer
    
    parser = SmilesParser()
    ir = parser.parse_bigsmiles(smiles)
    monomer = bigsmilesir_to_monomer(ir)
    
    return monomer
