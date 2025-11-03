# smiles_parser.py
"""
SMILES -> IR single-file parser using Lark (Python 3.11+)

Features:
- Atoms (including bracket form with isotope/charge/H count/chirality/class)
- Bonds: implicit, -, =, #, $, :, /, \
- Branches (...) and ring closures (1..9 and %NN)
- Aromatic small letters treated as aromatic by writing convention

Limitations:
- No aromaticity perception beyond "as-written"
- Stereochemistry kept as up/down for / and \\ (no full E/Z inference)
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum, auto
from pathlib import Path
from typing import Iterable

from lark import Tree, Token

from .base import GrammarConfig, GrammarParserBase


# =========================
# IR (enums + data classes)
# =========================

class BondOrder(Enum):
    """Standardized bond order labels."""
    SINGLE = auto()
    DOUBLE = auto()
    TRIPLE = auto()
    QUAD = auto()
    AROMATIC = auto()
    UNKNOWN = auto()


class BondStereo(Enum):
    """Stereochemical decoration for bonds."""
    NONE = auto()
    UP = auto()      # "/" direction as written
    DOWN = auto()    # "\" direction as written
    E = auto()       # placeholder (not inferred here)
    Z = auto()       # placeholder (not inferred here)
    UNKNOWN = auto()


AROMATIC_ELEMENTS: set[str] = {"b", "c", "n", "o", "s", "p", "se", "as"}


@dataclass
class AtomIR:
    """
    Atom node in the molecular graph IR.
    - element: kept as written ("c" means aromatic carbon)
    - h_count_explicit: explicit H count in bracket form (None if absent)
    - atom_map: value from [X:12] if present
    """
    id: int
    element: str
    aromatic: bool = False
    isotope: int | None = None
    charge: int = 0
    h_count_explicit: int | None = None
    chiral: str | None = None
    atom_map: int | None = None
    # Future extension (e.g., BigSMILES ports)
    port_tag: str | None = None
    meta: dict[str, str] = field(default_factory=dict)


@dataclass
class BondIR:
    """Undirected edge between atoms u and v."""
    u: int
    v: int
    order: BondOrder = BondOrder.SINGLE
    stereo: BondStereo = BondStereo.NONE
    aromatic: bool = False


@dataclass
class MolGraphIR:
    """
    Molecule graph IR without coordinates.
    - atoms: list of AtomIR
    - bonds: list of BondIR
    - total_charge: integer sum of atomic charges
    - has_disconnected: should remain False for plain SMILES (no ".")
    """
    atoms: list[AtomIR] = field(default_factory=list)
    bonds: list[BondIR] = field(default_factory=list)
    total_charge: int = 0
    has_disconnected: bool = False
    source: str | None = None


# =========================
# SmilesParser (from file)
# =========================

# Map a bond symbol to (BondOrder, BondStereo)
BOND_MAP: dict[str, tuple[BondOrder, BondStereo]] = {
    "":   (BondOrder.SINGLE,   BondStereo.NONE),   # implicit
    "-":  (BondOrder.SINGLE,   BondStereo.NONE),
    "=":  (BondOrder.DOUBLE,   BondStereo.NONE),
    "#":  (BondOrder.TRIPLE,   BondStereo.NONE),
    "$":  (BondOrder.QUAD,     BondStereo.NONE),
    ":":  (BondOrder.AROMATIC, BondStereo.NONE),
    "/":  (BondOrder.SINGLE,   BondStereo.UP),
    "\\": (BondOrder.SINGLE,   BondStereo.DOWN),
}


class SmilesParser(GrammarParserBase[MolGraphIR]):
    """
    SMILES parser that loads a Lark grammar from a .lark file.

    Usage:
        parser = SmilesParser()  # looks for ./smiles.lark by default
        ir = parser.parse("C(=O)O")
    """

    def __init__(
        self,
        grammar_path: Path | None = None,
        *,
        start: str = "smiles",
        auto_reload: bool = True,
    ):
        if grammar_path is None:
            grammar_dir = Path(__file__).parent / "grammar"
            new_path = grammar_dir / "smiles.lark"
            legacy_path = grammar_dir / "base_smiles.lark"
            grammar_path = new_path if new_path.exists() else legacy_path
        cfg = GrammarConfig(
            grammar_path=grammar_path,
            start=start,
            parser="lalr",
            propagate_positions=False,
            maybe_placeholders=False,
            auto_reload=auto_reload,
        )
        # 在 GrammarParserBase 里给 Lark(...) 传 lexer="contextual"
        # 如果没暴露参数，就在 _compile_grammar 里加：lexer="contextual"

        super().__init__(cfg)
        # runtime state per-parse
        self._ring_open: dict[int, tuple[int, str | None]] = {}
        self._next_atom_id: int = 0

    # ---------- public API ----------

    def build(self, tree: Tree) -> MolGraphIR:
        """Convert the Lark parse tree to MolGraphIR."""
        self._ring_open.clear()
        self._next_atom_id = 0
        self._ir = MolGraphIR(source=None)
        self._parse_smiles(tree)
        return self._ir

    # ---------- primitives ----------

    def _new_atom(
        self,
        element: str,
        aromatic: bool = False,
        isotope: int | None = None,
        charge: int = 0,
        h_explicit: int | None = None,
        chiral: str | None = None,
        atom_map: int | None = None,
    ) -> int:
        aid = self._next_atom_id
        self._next_atom_id += 1
        self._ir.atoms.append(
            AtomIR(
                id=aid,
                element=element,
                aromatic=aromatic,
                isotope=isotope,
                charge=charge,
                h_count_explicit=h_explicit,
                chiral=chiral,
                atom_map=atom_map,
            )
        )
        self._ir.total_charge += charge
        return aid

    def _add_bond(self, u: int, v: int, symbol: str | None):
        sym = symbol or ""
        order, stereo = BOND_MAP.get(sym, (BondOrder.UNKNOWN, BondStereo.NONE))

        # Implicit bond between aromatic atoms: treat as aromatic
        if sym == "" and self._ir.atoms[u].aromatic and self._ir.atoms[v].aromatic:
            order = BondOrder.AROMATIC

        self._ir.bonds.append(
            BondIR(u=u, v=v, order=order, stereo=stereo, aromatic=(order == BondOrder.AROMATIC))
        )

    # ---------- walkers ----------

    def _parse_smiles(self, node: Tree) -> tuple[int, int]:
        """
        smiles: chain (DOT chain)*
        Returns (first_atom_id, last_atom_id).
        """
        assert node.data == "smiles"
        assert node.children, "SMILES grammar should yield at least one chain"
        first_id, last_id = self._parse_chain(node.children[0])

        i = 1
        while i < len(node.children):
            child = node.children[i]
            if isinstance(child, Token) and child.type == "DOT":
                self._ir.has_disconnected = True
                if i + 1 >= len(node.children):
                    raise ValueError("Dangling disconnected component without a chain")
                next_chain = node.children[i + 1]
                if not isinstance(next_chain, Tree) or next_chain.data != "chain":
                    raise ValueError("Expected chain after dot in SMILES")
                _, last_id = self._parse_chain(next_chain)
                i += 2
                continue
            if isinstance(child, Tree) and child.data == "chain":
                _, last_id = self._parse_chain(child)
                i += 1
                continue
            raise ValueError(f"Unexpected node in smiles: {child!r}")

        return first_id, last_id

    def _parse_chain(self, node: Tree) -> tuple[int, int]:
        """
        chain: branched_atom (ring_bond | atom_assembly)*
        """
        assert node.data == "chain"
        first_id, last_id = self._parse_branched_atom(node.children[0])

        for child in node.children[1:]:
            if isinstance(child, Tree):
                if child.data == "atom_assembly":
                    # atom_assembly: BOND_SYM? branched_atom
                    bsym: str | None = None
                    ba: Tree | None = None
                    for c in child.children:
                        if isinstance(c, Token) and c.type in ("BOND_SYM", "CHAIN_BOND_SYM"):
                            bsym = c.value
                        elif isinstance(c, Tree) and c.data == "branched_atom":
                            ba = c
                    assert ba is not None
                    bf, bl = self._parse_branched_atom(ba)
                    self._add_bond(last_id, bf, bsym)
                    last_id = bl
                elif child.data == "ring_bond":
                    self._handle_ring_bond(last_id, child)
                else:
                    raise ValueError(f"Unexpected node in chain: {child.data}")
            else:
                raise ValueError(f"Unexpected token in chain: {child!r}")

        return first_id, last_id

    def _parse_branched_atom(self, node: Tree) -> tuple[int, int]:
        """
        branched_atom: atom branch* ring_bond*
        Returns (first_atom_id, last_atom_id).
        """
        assert node.data == "branched_atom"
        atom_node = node.children[0]
        cur = self._parse_atom(atom_node)
        first = cur

        i = 1
        # branch*
        while i < len(node.children) and isinstance(node.children[i], Tree) and node.children[i].data == "branch":
            br = node.children[i]
            # branch: "(" BOND_SYM? smiles ")"
            bsym: str | None = None
            subtree: Tree | None = None
            for bc in br.children:
                if isinstance(bc, Token) and bc.type in ("BOND_SYM", "CHAIN_BOND_SYM"):
                    bsym = bc.value
                elif isinstance(bc, Tree) and bc.data == "smiles":
                    subtree = bc
            assert subtree is not None
            bf, _ = self._parse_smiles(subtree)
            self._add_bond(cur, bf, bsym)
            i += 1

        # ring_bond*
        while i < len(node.children) and isinstance(node.children[i], Tree) and node.children[i].data == "ring_bond":
            self._handle_ring_bond(cur, node.children[i])
            i += 1

        return first, cur

    def _handle_ring_bond(self, cur: int, node: Tree):
        """
        ring_bond: RING_BOND_SYM? RING_NUM
        """
        bsym: str | None = None
        idx_str = ""

        for c in node.children:
            if isinstance(c, Token) and c.type == "RING_BOND_SYM":
                bsym = c.value
            elif isinstance(c, Token) and c.type == "RING_NUM":
                val = c.value
                if val.startswith("%"):
                    idx_str = val[1:]
                else:
                    idx_str = val

        if not idx_str:
            return
        idx = int(idx_str)

        if idx not in self._ring_open:
            self._ring_open[idx] = (cur, bsym)
        else:
            other, bsym0 = self._ring_open.pop(idx)
            use_sym = bsym if bsym is not None else bsym0
            self._add_bond(other, cur, use_sym)

    def _parse_atom(self, node: Tree | Token) -> int:
        if isinstance(node, Token):
            # ATOM_SIMPLE or STAR as token
            el = node.value
            if el == "*":
                return self._new_atom("*", False)
            aromatic = (el in AROMATIC_ELEMENTS)
            return self._new_atom(el, aromatic)

        if node.data == "atom":
            return self._parse_atom(node.children[0])

        if node.data == "atom_simple":
            t = node.children[0]
            if isinstance(t, Token):
                if t.type == "STAR":
                    return self._new_atom("*", False)
                elif t.type in ("ALIPHATIC_ORGANIC", "ELEMENT_SYM"):
                    el = t.value
                    aromatic = (el in AROMATIC_ELEMENTS)
                    return self._new_atom(el, aromatic)
                elif t.type in ("AROMATIC_ORGANIC", "AROMATIC_SYM"):
                    el = t.value
                    aromatic = True
                    return self._new_atom(el, aromatic)
            # Fallback
            raise ValueError(f"Unhandled atom_simple node: {node}")

        if node.data == "bracket_atom":
            # "[" isotope? atom_symbol chiral? h_count? atom_charge? atom_class? "]"
            isotope: int | None = None
            symbol: str | None = None
            chiral: str | None = None
            h_explicit: int | None = None
            charge: int = 0
            atom_map: int | None = None

            for c in node.children:
                if isinstance(c, Tree):
                    if c.data == "isotope":
                        isotope = int(c.children[0].value)
                    elif c.data == "atom_symbol":
                        symbol = self._atom_symbol(c)
                    elif c.data == "chiral":
                        chiral = c.children[0].value if c.children else "@"
                    elif c.data == "h_count":
                        h_explicit = 1 if len(c.children) == 0 else int(c.children[-1].value)
                    elif c.data == "atom_charge":
                        charge = self._parse_charge(c)
                    elif c.data == "atom_class":
                        atom_map = int(c.children[-1].value)

            el = symbol or "*"
            aromatic = (el.lower() in AROMATIC_ELEMENTS and el.islower())
            return self._new_atom(
                element=el,
                aromatic=aromatic,
                isotope=isotope,
                charge=charge,
                h_explicit=h_explicit,
                chiral=chiral,
                atom_map=atom_map,
            )
        
        elif node.data == "ALIPHATIC_ORGANIC":
            # Unused with current grammar (wrapped by atom_simple)
            el = node.children[0].value
            aromatic = (el in AROMATIC_ELEMENTS)
            return self._new_atom(el, aromatic)

        elif node.data == "AROMATIC_ORGANIC":
            # Unused with current grammar (wrapped by atom_simple)
            el = node.children[0].value
            aromatic = (el in AROMATIC_ELEMENTS)
            return self._new_atom(el, aromatic)

        raise ValueError(f"Unhandled atom node: {node}")

    # ---------- helpers ----------

    @staticmethod
    def _atom_symbol(node: Tree) -> str:
        """Return symbol from atom_symbol rule."""
        c = node.children[0]
        return c.value if isinstance(c, Token) else c.children[0].value  # aromatic_symbol tree

    @staticmethod
    def _parse_charge(node: Tree) -> int:
        """
        Parse atom_charge:
          "-", "+", "--", "++", "-2", "+3"
        """
        if not node.children:
            return 0
        t0 = node.children[0]
        if isinstance(t0, Token):
            ttype = t0.type
            if ttype == "CHG_DMINUS":
                return -2
            if ttype == "CHG_DPLUS":
                return +2
            if ttype in ("CHG_MINUS", "CHG_PLUS"):
                n = 1
                if len(node.children) > 1 and isinstance(node.children[1], Token) and node.children[1].type == "INT":
                    n = int(node.children[1].value)
                return -n if ttype == "CHG_MINUS" else +n
        return 0

    @staticmethod
    def _extract_bond_symbol(node: Tree | Token | None) -> str | None:
        """
        Extracts a bond symbol from a `bond_symbol` rule node or returns None.
        """
        if node is None:
            return None
        if isinstance(node, Token):
            # Not expected with current grammar; kept for robustness.
            return node.value
        if isinstance(node, Tree) and node.data == "bond_symbol":
            # bond_symbol is a rule consisting of a single literal token.
            # e.g., children: [Token('-')] or [Token('/')] ...
            for ch in node.children:
                if isinstance(ch, Token):
                    return ch.value
        return None


def parse_smiles(text: str, grammar_path: Path | None = None) -> MolGraphIR:
    """
    Convenience function if you don't want to instantiate the class yourself.
    """
    parser = SmilesParser(grammar_path=grammar_path)
    return parser.parse(text)
