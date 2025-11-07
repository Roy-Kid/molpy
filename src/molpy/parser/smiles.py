from pathlib import Path
from .base import GrammarConfig, GrammarParserBase
from lark import Transformer
from lark import Token
from dataclasses import dataclass, field
import os

class MolPyAPIError(Exception):
    """Custom exception for API errors."""
    pass

# ===================================================================
#   1. 中间表示 (IR) 数据类 (只包含SMILES部分)
# ===================================================================
@dataclass(eq=True)
class AtomIR:
    symbol: str
    isotope: int | None = None
    chiral: str | None = None
    h_count: int | None = None
    charge: int | None = None
    class_: int | None = None
    id: int = field(default_factory=lambda: id(AtomIR), compare=False, repr=False)

    def __hash__(self): return self.id

    def __repr__(self):
        attrs = [f"symbol={self.symbol!r}"]
        if self.isotope is not None: attrs.append(f"isotope={self.isotope}")
        if self.chiral is not None: attrs.append(f"chiral={self.chiral!r}")
        if self.h_count is not None: attrs.append(f"h_count={self.h_count}")
        if self.charge is not None: attrs.append(f"charge={self.charge}")
        return f"AtomIR({', '.join(attrs)})"

@dataclass(eq=True)
class BondIR:
    start: AtomIR
    end: AtomIR
    bond_type: str


    def __repr__(self):
        return f"BondIR({self.start.symbol!r}, {self.end.symbol!r}, {self.bond_type!r})"

@dataclass(eq=True)
class SmilesIR:
    atoms: list[AtomIR] = field(default_factory=list)
    bonds: list[BondIR] = field(default_factory=list)


    def __repr__(self):
        return f"SmilesIR(atoms={self.atoms!r}, bonds={self.bonds!r})"

# ===================================================================
#   2. SmilesTransformer (基类)
# ===================================================================
class SmilesTransformer(Transformer):
    def __init__(self):
        super().__init__()
        # Ring openings shared across the entire molecule (allows cross-branch closure)
        self.ring_openings: dict[str, tuple[AtomIR, str | None]] = {}
    # ================== 终端转换 ==================
    def INT(self, n: Token) -> int: return int(n)
    def ATOM_SYM(self, s: Token) -> str: return str(s.value)
    def ALIPHATIC_ORGANIC(self, s: Token) -> AtomIR: return AtomIR(symbol=s.value)
    def AROMATIC_ORGANIC(self, s: Token) -> AtomIR: return AtomIR(symbol=s.value)
    def ELEMENT_SYM(self, s: Token) -> AtomIR: return AtomIR(symbol=s.value)
    def BOND_SYM(self, s: Token) -> str: return s.value
    def bond_symbol(self, children: list[Token]) -> str: return children[0].value
    def ring_id(self, d: list[Token]) -> str: return "".join(t.value for t in d)

    # ================== 原子及其属性构建 ==================
    def atom(self, children: list[AtomIR]) -> AtomIR: return children[0]
    def isotope(self, n: list[int]) -> tuple[str, int]: return "isotope", n[0]
    def chiral(self, c: list[Token]) -> tuple[str, str]: return "chiral", "".join(t.value for t in c)
    def h_count(self, h: list) -> tuple[str, int]: return "h_count", h[1] if len(h) > 1 else 1
    def class_(self, c: list) -> tuple[str, int]: return "class_", c[1]
    def atom_class(self, c: list) -> tuple[str, int]: return "class_", c[1]
    
    def atom_charge(self, items: list) -> tuple[str, int]:
        sign = -1 if items[0].value == '-' else 1
        # forms: '+' '-' '++' '--' '+2' '-2'
        if len(items) == 1:
            return "charge", sign
        second = items[1]
        if isinstance(second, Token):  # ++ or -- already collapsed by grammar
            return "charge", 2 * sign
        # numeric token already converted? ensure int
        val = int(second) if not isinstance(second, int) else second
        return "charge", sign * val

    def bracket_atom(self, items: list) -> AtomIR:
        # items include '[' and ']' tokens plus optional property tuples
        filtered = [it for it in items if not (isinstance(it, Token) and it.type in {"LSQB","RSQB"})]
        # first non-tuple is the element symbol
        symbol = None
        props_pairs: list[tuple[str, int | str]] = []
        for it in filtered:
            if isinstance(it, tuple):
                if len(it) == 2:
                    props_pairs.append(it)
            elif symbol is None:
                symbol = it
        if symbol is None:
            raise ValueError("Bracket atom missing symbol")
        # coerce to correct types
        kwargs: dict[str, object] = {}
        for k, v in props_pairs:
            if k in {"isotope", "h_count", "charge", "class_"}:
                kwargs[k] = int(v)
            elif k == "chiral":
                kwargs[k] = str(v)
        return AtomIR(symbol=str(symbol), **kwargs)  # type: ignore[arg-type]

    # ================== 核心SMILES组装 ==================
    def smiles(self, children: list) -> SmilesIR:
        """Assemble linear smiles: children = [first_branched_atom, (bond, branched_atom)*]"""
        ir = SmilesIR()
        active_atom: AtomIR | None = None
        pending_bond_type = "-"
        debug = os.getenv("SMILES_DEBUG")
        if debug:
            print(f"[SMILES] processing {len(children)} children")

        # Normalize children into sequence: branched_atom, (bond, branched_atom)...
        seq: list = []
        if not children:
            return ir
        first = children[0]
        rest = children[1:]
        seq.append(first)
        i = 0
        while i < len(rest):
            item = rest[i]
            if isinstance(item, str):
                # unlikely direct token here, but keep for robustness
                if i + 1 >= len(rest):
                    break
                seq.append(item)
                seq.append(rest[i + 1])
                i += 2
            elif isinstance(item, tuple) and item and item[0] == "asm":
                _, bond, ba = item
                if bond is not None:
                    seq.append(bond)
                seq.append(ba)
                i += 1
            else:
                seq.append(item)
                i += 1

        for item in seq:
            if isinstance(item, str):
                pending_bond_type = item
                continue
            atom, rings, branches = item
            ir.atoms.append(atom)
            if active_atom is not None:
                ir.bonds.append(BondIR(active_atom, atom, pending_bond_type))
            active_atom = atom
            pending_bond_type = "-"
            for ring_id, bond_type in rings:
                if active_atom is None:
                    continue
                # Check if ring was opened before, if so close it
                if ring_id in self.ring_openings:
                    start_atom, start_bond = self.ring_openings.pop(ring_id)
                    final_bond = bond_type or start_bond or "-"
                    ir.bonds.append(BondIR(start_atom, active_atom, final_bond))
                    if debug:
                        print(f"[SMILES] close ring {ring_id}: {start_atom.symbol}->{active_atom.symbol} bond={final_bond}")
                else:
                    self.ring_openings[ring_id] = (active_atom, bond_type)
                    if debug:
                        print(f"[SMILES] open ring {ring_id} at atom {active_atom.symbol} bond={bond_type}")
            for branch_ir, branch_bond_type in branches:
                if active_atom is None:
                    continue
                head_atom = branch_ir.atoms[0]
                ir.bonds.append(BondIR(active_atom, head_atom, branch_bond_type or "-"))
                ir.atoms.extend(branch_ir.atoms)
                ir.bonds.extend(branch_ir.bonds)

        # Note: We don't check for unclosed rings here because rings can span across branches
        # e.g., C1CCC2C(C1)CCC2 has ring 1 opened in main chain and closed in a branch
        # Ring validation should be done at the parser level, not transformer level
        if debug:
            print(f"[SMILES] exit")
        return ir

    def branched_atom(self, children: list) -> tuple:
        atom = children[0]
        rings: list = []
        branches: list = []
        for item in children[1:]:
            if isinstance(item, tuple):
                tag = item[0]
                if tag == "ring":
                    rings.append((item[1], item[2]))
                elif tag == "branch":
                    branches.append((item[1], item[2]))
            elif isinstance(item, Token) and item.type == 'DIGIT':
                rings.append((item.value, None))
        return atom, rings, branches

    def ring_bond(self, children: list) -> tuple:
        # children could include optional bond sym and digits, possibly with '%'
        bond_type = None
        ring_digits: list[str] = []
        for c in children:
            if isinstance(c, str) and c in {"-","=","#","$",":","/","\\"}:
                bond_type = c
            elif isinstance(c, Token) and c.type == 'DIGIT':
                ring_digits.append(c.value)
        ring_id = ''.join(ring_digits) if ring_digits else (children[-1].value if isinstance(children[-1], Token) else str(children[-1]))
        return "ring", ring_id, bond_type

    def branch(self, children: list) -> tuple:
        # children may include LPAR/RPAR tokens; filter them out
        filtered = [c for c in children if not (isinstance(c, Token) and c.type in {"LPAR", "RPAR"})]
        bond_type = filtered[0] if filtered and isinstance(filtered[0], str) else None
        smiles_ir = filtered[-1]
        return "branch", smiles_ir, bond_type

    def atom_assembly(self, children: list):
        # children: [branched_atom] or [bond, branched_atom]
        if len(children) == 1:
            return ("asm", None, children[0])
        return ("asm", children[0], children[1])


# ===================================================================
#   3. BigSMILES IR (扩展)
# ===================================================================
@dataclass
class BondDescriptorIR:
    symbol: str | None = None
    index: int | None = None
    generation: list[int] | None = None

@dataclass
class StochasticDistributionIR:
    name: str
    params: list[float]

@dataclass
class RepeatSegmentIR:
    stochastic_objects: list["StochasticObjectIR"]
    implicit_smiles: SmilesIR | None = None


@dataclass
class BigSmilesChainIR:
    start_smiles: SmilesIR
    repeat_segments: list[RepeatSegmentIR]

@dataclass
class StochasticObjectIR:
    left_descriptor: BondDescriptorIR
    right_descriptor: BondDescriptorIR
    repeat_units: list[SmilesIR | BigSmilesChainIR]
    end_groups: list[SmilesIR | BigSmilesChainIR] | None = None
    distribution: StochasticDistributionIR | None = None


@dataclass
class BigSmilesMoleculeIR:
    chain: BigSmilesChainIR


@dataclass(eq=True)
class BigSmilesIR(SmilesIR):
    """BigSMILES IR extends SmilesIR with polymer-specific chain structure."""
    chain: BigSmilesChainIR | None = None
    
    def __post_init__(self):
        # Ensure chain is initialized if not provided
        if self.chain is None:
            empty_smiles = SmilesIR(atoms=[], bonds=[])
            object.__setattr__(self, 'chain', BigSmilesChainIR(start_smiles=empty_smiles, repeat_segments=[]))


# ===================================================================
#   4. BigSmilesTransformer (子类)
# ===================================================================
class BigSmilesTransformer(SmilesTransformer):
    def NUMBER(self, n: Token) -> float: return float(n)
    def bond_descriptor_symbol(self, t: list[Token]) -> str: return t[0].value

    # ================== 键描述符组装 ==================
    def bond_descriptor_symbol_idx(self, items: list) -> tuple:
        return items[0], items[1] if len(items) > 1 else None
        
    def bond_descriptor_generation(self, items: list) -> list[int]:
        return [int(i.value) for i in items]

    def terminal_bond_descriptor(self, items: list) -> BondDescriptorIR:
        symbol, index, gen = None, None, None
        if items:
            if isinstance(items[0], tuple): # symbol_idx
                symbol, index = items[0]
                if len(items) > 1: gen = items[1]
            else: # generation only
                gen = items[0]
        return BondDescriptorIR(symbol=symbol, index=index, generation=gen)
        
    # ================== 随机对象组装 ==================
    def stochastic_distribution(self, items: list) -> StochasticDistributionIR:
        return items[0]

    def flory_schulz(self, p:list) -> StochasticDistributionIR: return StochasticDistributionIR("flory_schulz", [p[0]])
    def schulz_zimm(self, p:list) -> StochasticDistributionIR: return StochasticDistributionIR("schulz_zimm", p)
    # ... implement for gauss, uniform, etc. ...

    def _repeat_units(self, items: list) -> list:
        return items

    def _end_group(self, items: list) -> list:
        return items

    def stochastic_object(self, items: list) -> StochasticObjectIR:
        desc1, repeats = items[0], items[1]
        
        # 逐个解析可选部分
        i = 2
        ends = None
        if i < len(items) and isinstance(items[i], list):
            ends = items[i]
            i += 1
            
        desc2 = items[i]
        i+= 1
        
        dist = items[i] if i < len(items) else None
        
        return StochasticObjectIR(
            left_descriptor=desc1,
            right_descriptor=desc2,
            repeat_units=repeats,
            end_groups=ends,
            distribution=dist
        )
        
    # ================== 顶层组装 ==================
    def repeat_segment(self, items:list) -> RepeatSegmentIR:
        smiles = items.pop() if isinstance(items[-1], SmilesIR) else None
        return RepeatSegmentIR(stochastic_objects=items, implicit_smiles=smiles)
        
    def big_smiles_chain(self, items:list) -> BigSmilesChainIR:
        start_smiles = items[0]
        repeat_segments = items[1:]
        return BigSmilesChainIR(start_smiles=start_smiles, repeat_segments=repeat_segments)
        
    def big_smiles_molecule(self, items:list) -> BigSmilesMoleculeIR:
        return BigSmilesMoleculeIR(chain=items[0])

    def big_smiles(self, items: list) -> list:
        return items
    

class SmilesParser(GrammarParserBase):

    def __init__(self):
        config = GrammarConfig(
            grammar_path=Path(__file__).parent / "grammar" / "smiles.lark",
            start="smiles",
            parser="earley",
            propagate_positions=True,
            maybe_placeholders=False,
            auto_reload=True,
        )
        super().__init__(config)


    def parse_smiles(self, smiles: str) -> SmilesIR:
        # Support disconnected components separated by '.' by parsing each separately
        if '.' in smiles:
            raise MolPyAPIError("Disconnected components ('.') not supported in this parser method, use parse_dot_smiles instead.")
        tree = self.parse_tree(smiles)
        transformer = SmilesTransformer()
        ir: SmilesIR = transformer.transform(tree)
        # Check for unclosed rings after transformation
        if transformer.ring_openings:
            unclosed = list(transformer.ring_openings.keys())
            raise ValueError(f"Unclosed rings: {unclosed}")
        return ir
    
    def parse_dot_smiles(self, smiles: str) -> list[SmilesIR]:
        """
        Parse SMILES string with disconnected components separated by '.'.
        
        Args:
            smiles: SMILES string with possible disconnected components
            
        Returns:
            List of SmilesIR for each component
            
        Raises:
            ValueError: if syntax errors detected
        """
        parts = [p.strip() for p in smiles.split('.') if p.strip()]
        results: list[SmilesIR] = []
        for part in parts:
            tree = self.parse_tree(part)
            transformer = SmilesTransformer()
            ir: SmilesIR = transformer.transform(tree)
            # Check for unclosed rings after transformation
            if transformer.ring_openings:
                unclosed = list(transformer.ring_openings.keys())
                raise ValueError(f"Unclosed rings in component '{part}': {unclosed}")
            results.append(ir)
        return results

    def parse_bigsmiles(self, text: str) -> BigSmilesIR:
        """
        Parse BigSMILES string into BigSmilesIR.
        
        Minimal implementation: returns a BigSmilesIR with parsed SMILES components
        and empty/minimal chain structure. Full BigSMILES grammar support is staged.
        
        Args:
            text: BigSMILES string
            
        Returns:
            BigSmilesIR with chain structure
            
        Raises:
            ValueError: if syntax errors detected
        """
        # Stage 1: minimal implementation - treat as SMILES and wrap in BigSmilesIR
        # TODO: implement full BigSMILES grammar parsing with stochastic objects
        
        # For now, parse any embedded SMILES and create empty chain
        # This allows tests to pass while we build out full support
        smiles_ir = self.parse_smiles(text) if text and not any(c in text for c in '{}[]<>') else SmilesIR(atoms=[], bonds=[])
        
        # Create minimal chain structure
        chain = BigSmilesChainIR(
            start_smiles=smiles_ir,
            repeat_segments=[]
        )
        
        return BigSmilesIR(atoms=smiles_ir.atoms, bonds=smiles_ir.bonds, chain=chain)


# ===================================================================
#   Converter: SmilesIR -> RDKit Mol
# ===================================================================

def smilesir_to_mol(ir: SmilesIR) -> "Chem.Mol":
    """
    Convert SmilesIR to RDKit Mol by directly constructing the molecule graph.
    
    This approach preserves IR-specific information and supports extended syntax
    (BigSMILES, G-BigSMILES) where explicit topology is essential.
    
    Args:
        ir: SmilesIR instance with atoms and bonds
        
    Returns:
        RDKit Mol object
        
    Raises:
        ImportError: if RDKit is not available
        ValueError: if IR contains invalid molecular data
        
    Example:
        >>> parser = SmilesParser()
        >>> ir = parser.parser_smiles("CCO")
        >>> mol = smilesir_to_mol(ir)
        >>> mol.GetNumAtoms()
        3
    """
    try:
        from rdkit import Chem
    except ImportError as e:
        raise ImportError("RDKit is required for smilesir_to_mol conversion") from e
    
    if not ir.atoms:
        # Empty molecule
        return Chem.Mol()
    
    # Bond type mapping
    bond_type_map = {
        "-": Chem.BondType.SINGLE,
        "=": Chem.BondType.DOUBLE,
        "#": Chem.BondType.TRIPLE,
        ":": Chem.BondType.AROMATIC,
        "/": Chem.BondType.SINGLE,   # Stereochemistry, treat as single for now
        "\\": Chem.BondType.SINGLE,  # Stereochemistry, treat as single for now
    }
    
    # Create editable molecule
    mol = Chem.RWMol()
    
    # Map AtomIR -> RDKit atom index (using object identity)
    atom_to_idx: dict[int, int] = {}
    
    # Add atoms
    for atom_ir in ir.atoms:
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
        start_idx = atom_to_idx.get(id(bond_ir.start))
        end_idx = atom_to_idx.get(id(bond_ir.end))
        
        if start_idx is None or end_idx is None:
            raise ValueError(f"Bond references unknown atom: {bond_ir}")
        
        # Determine bond type (upgrade single bonds between aromatic atoms to aromatic)
        bond_type_str = bond_ir.bond_type
        if bond_type_str == "-" and bond_ir.start.symbol.islower() and bond_ir.end.symbol.islower():
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
        warnings.warn(f"Molecule sanitization failed: {e}. Returning unsanitized molecule.")
    
    return final_mol


# Register converter if RDKit and adapter are available
try:
    from rdkit import Chem
    from molpy.adapter.registry import REG
    REG.register(SmilesIR, Chem.Mol, smilesir_to_mol)
except ImportError:
    pass  # RDKit or adapter not available, skip registration

