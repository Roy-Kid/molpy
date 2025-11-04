from __future__ import annotations

from pathlib import Path
from typing import Sequence, cast, Literal

from lark import Lark, Token, Tree, UnexpectedInput

from .smiles_ir import (
    Span,
    IRAtom,
    IRBond,
    IRRingBond,
    IRBranch,
    IRDiagnostic,
    # SMILES
    IRSmilesDot,
    IRSmilesMolecule,
    IRSmilesDocument,
    IRSmilesPart,
    # BigSMILES
    IRBigGeneration,
    IRBigIndexContext,
    IRBigIndexExpr,
    IRBigNonCovalent,
    IRBigConnector,
    IRBigDistribution,
    IRBigSystemDot,
    IRBigStochasticBlock,
    IRBigFragmentRef,
    IRBigFragmentDef,
    IRBigPart,
    IRBigSmilesMolecule,
    IRBigSmilesDocument,
    # G-BigSMILES
    IRGPart,
    IRGBigSmilesMolecule,
    IRGBigSmilesDocument,
)

__all__ = ["parse_smiles", "parse_bigsmiles", "parse_gbigsmiles"]

_GRAMMAR_PATH = Path(__file__).with_name("grammar").joinpath("smiles.lark")
_GRAMMAR = _GRAMMAR_PATH.read_text(encoding="utf-8")
_PARSER = Lark(
    _GRAMMAR,
    start="document",
    parser="earley",
    lexer="dynamic",
    propagate_positions=True,
    maybe_placeholders=False,
    import_paths=[str(_GRAMMAR_PATH.parent)],
)

def _span_from_exception(exc: UnexpectedInput) -> Span | None:
    if exc.pos_in_stream is None:
        return None
    column = getattr(exc, "column", None) or 1
    line = getattr(exc, "line", None) or 1
    return Span(
        start=exc.pos_in_stream,
        end=exc.pos_in_stream,
        line=line,
        column=column,
    )

def parse_smiles(smiles: str) -> IRSmilesDocument:
    parser = SmilesParser()
    return cast(IRSmilesDocument, parser.parse(smiles, flavor="smiles"))

def parse_bigsmiles(smiles: str) -> IRBigSmilesDocument:
    parser = SmilesParser()
    return cast(IRBigSmilesDocument, parser.parse(smiles, flavor="bigsmiles"))

def parse_gbigsmiles(smiles: str) -> IRGBigSmilesDocument:
    parser = SmilesParser()
    return cast(IRGBigSmilesDocument, parser.parse(smiles, flavor="gbigsmiles"))


class SmilesParser:

    def __init__(self):
        self._parser = _PARSER
        self.features: set[str] = set()
        self.errors: list[IRDiagnostic] = []
        # Holds current input text during a parse; used for node text slicing
        self._input_text: str | None = None

    def parse(self, smiles: str, flavor) -> IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument:
        """
        Parse input string and convert to IR via build().
        """
        try:
            # Stash the current input for downstream helpers that need original text
            self._input_text = smiles
            tree = _PARSER.parse(smiles)
            return self.build(tree, flavor)
        except UnexpectedInput as exc:
            self.errors.append(IRDiagnostic(level="error", message=str(exc), span=_span_from_exception(exc)))
            match flavor:
                case "smiles":
                    return IRSmilesDocument(
                        molecules=[],
                        features=set(),
                        errors=self.errors,
                    )
                case "bigsmiles":
                    return IRBigSmilesDocument(
                        molecules=[],
                        fragments=[],
                        features=set(),
                        errors=self.errors,
                    )
                case "gbigsmiles":
                    return IRGBigSmilesDocument(
                        molecules=[],
                        fragments=[],
                        features=set(),
                        errors=self.errors,
                    )
                case _:
                    return IRSmilesDocument(
                        molecules=[],
                        features=set(),
                        errors=self.errors,
                    )
        finally:
            # Ensure no lingering references between parses
            self._input_text = None

    def build(self, tree: Tree, flavor: str) -> IRSmilesDocument | IRBigSmilesDocument | IRGBigSmilesDocument:
        big_molecules: list[IRBigSmilesMolecule] = []
        fragments: list[IRBigFragmentDef] = []

        for child in tree.children:
            if isinstance(child, Tree):
                if child.data == "molecule":
                    big_molecules.append(self._build_molecule_big(child))
                elif child.data == "fragment_definition":
                    fragments.append(self._build_fragment_definition(child))

        if big_molecules:
            self.features.add("smiles")

        if flavor == "smiles":
            smiles_molecules: list[IRSmilesMolecule] = []
            for m in big_molecules:
                parts_smiles: list[IRSmilesPart] = [p for p in m.parts if isinstance(p, IRAtom)]
                smiles_molecules.append(IRSmilesMolecule(parts=parts_smiles, span=m.span))
            return IRSmilesDocument(
                flavor="smiles",
                molecules=smiles_molecules,
                features=self.features | {"smiles"},
                errors=self.errors,
            )

        if flavor == "gbigsmiles":
            return IRGBigSmilesDocument(
                flavor="gbigsmiles",
                molecules=[IRGBigSmilesMolecule(parts=cast(list[IRGPart], list(m.parts)), generation=m.generation, span=m.span) for m in big_molecules],
                fragments=fragments,
                features=self.features,
                errors=self.errors,
            )

        return IRBigSmilesDocument(
            flavor="bigsmiles",
            molecules=big_molecules,
            fragments=fragments,
            features=self.features,
            errors=self.errors,
        )

    # ---------- Molecules ----------

    def _build_molecule_big(self, node: Tree) -> IRBigSmilesMolecule:
        parts: list[IRBigPart] = []
        generation: IRBigSystemDot | None = None

        for child in node.children:
            if isinstance(child, Tree):
                if child.data == "molecule_body":
                    parts.extend(self._build_molecule_body(child))
                elif child.data == "dot_generation":
                    generation = self._build_dot_generation(child)

        return IRBigSmilesMolecule(parts=parts, generation=generation, span=self._span(node))

    def _build_molecule_body(self, node: Tree) -> list[IRBigPart]:
        parts: list[IRBigPart] = []

        for child in node.children:
            if not isinstance(child, Tree):
                continue
            if child.data == "smiles":
                parts.extend(self._build_smiles(child))
            elif child.data == "repeat_block":
                parts.extend(self._build_repeat_block(child))
            elif child.data == "stochastic_object":
                parts.append(self._build_stochastic_object(child))
            elif child.data == "bond_descriptor":
                parts.append(self._build_bond_descriptor(child))
            elif child.data == "fragment_declaration":
                parts.append(self._build_fragment_ref(child))

        return parts

    def _build_smiles(self, node: Tree) -> list[IRBigPart]:
        parts: list[IRBigPart] = []
        children = [child for child in node.children if isinstance(child, Tree)]
        if not children:
            return parts

        first = children[0]
        parts.append(cast(IRBigPart, self._build_branched_atom(first, leading_bond=None)))

        for atom_assembly in children[1:]:
            if atom_assembly.data != "atom_assembly":
                continue
            bond: IRBond | None = None
            target: Tree | None = None
            for child in atom_assembly.children:
                if isinstance(child, Tree) and child.data == "bond_symbol":
                    bond = self._build_bond(child)
                elif isinstance(child, Tree):
                    target = child
            if target is None:
                continue
            parts.append(cast(IRBigPart, self._build_branched_atom(target, leading_bond=bond)))

        return parts

    def _build_repeat_block(self, node: Tree) -> list[IRBigPart]:
        parts: list[IRBigPart] = []
        for child in node.children:
            if isinstance(child, Tree):
                if child.data == "stochastic_object":
                    parts.append(self._build_stochastic_object(child))
                elif child.data == "smiles":
                    parts.extend(self._build_smiles(child))
                elif child.data == "molecule_body":
                    parts.extend(self._build_molecule_body(child))
        return parts

    # ---------- Branched atoms ----------

    def _build_branched_atom(self, node: Tree, *, leading_bond: IRBond | None) -> object:
        atom_entry: Tree | None = None
        ring_nodes: list[Tree] = []
        branch_nodes: list[Tree] = []

        for child in node.children:
            if not isinstance(child, Tree):
                continue
            if child.data == "atom_entry":
                atom_entry = child
            elif child.data == "ring_bond":
                ring_nodes.append(child)
            elif child.data == "branch":
                branch_nodes.append(child)

        part = self._build_atom_entry(atom_entry, node)

        if isinstance(part, IRAtom):
            part.bond_to_prev = leading_bond
            part.rings.extend(self._build_ring_bonds(ring_nodes))
            part.branches.extend(self._build_branches(branch_nodes))
        return part

    def _build_atom_entry(self, node: Tree | None, parent: Tree) -> object:
        if node is None:
            return self._make_placeholder_atom(parent)

        child = next((c for c in node.children if isinstance(c, Tree)), None)
        if child is None:
            return self._make_placeholder_atom(parent)

        if child.data == "atom":
            return self._build_atom(child, parent)
        if child.data == "stochastic_object":
            return self._build_stochastic_object(child)
        if child.data == "bond_descriptor":
            descriptor = self._build_bond_descriptor(child)
            self.features.add("bigsmiles")
            if descriptor.connector == ":" or descriptor.mode == "non_covalent" or descriptor.noncovalent:
                self.features.add("gbigsmiles")
            return descriptor
        if child.data == "fragment_declaration":
            self.features.add("bigsmiles")
            return self._build_fragment_ref(child)

        return self._make_placeholder_atom(parent)

    def _build_atom(self, node: Tree, parent: Tree) -> IRAtom:
        atom_kind: Literal['bracket','aliphatic','aromatic','star'] = "aliphatic"
        symbol: str | None = None
        isotope: int | None = None
        chiral: str | None = None
        hcount: int | None = None
        charge: int | None = None
        clazz: int | None = None

        for child in node.children:
            if not isinstance(child, Tree):
                continue
            if child.data == "bracket_atom":
                (
                    atom_kind,
                    symbol,
                    isotope,
                    chiral,
                    hcount,
                    charge,
                    clazz,
                ) = self._build_bracket_atom(child)
            elif child.data == "aliphatic_organic":
                atom_kind = "aliphatic"
                symbol = self._node_text(child)
            elif child.data == "aromatic_organic":
                atom_kind = "aromatic"
                symbol = self._node_text(child)
            elif child.data == "STAR":
                atom_kind = "star"
                symbol = "*"

        self.features.add("smiles")

        return IRAtom(
            kind=atom_kind,
            symbol=symbol,
            isotope=isotope,
            chiral=chiral,
            hcount=hcount,
            charge=charge,
            clazz=clazz,
            rings=[],
            branches=[],
            span=self._span(parent),
        )

    def _build_bracket_atom(
        self, node: Tree
    ) -> tuple[Literal['bracket'], str | None, int | None, str | None, int | None, int | None, int | None]:
        symbol: str | None = None
        isotope: int | None = None
        chiral: str | None = None
        hcount: int | None = None
        charge: int | None = None
        clazz: int | None = None

        for child in node.children:
            if isinstance(child, Tree):
                if child.data == "isotope":
                    value = self._find_token(child, "INT")
                    isotope = int(value) if value else None
                elif child.data == "atom_symbol":
                    text = self._node_text(child)
                    symbol = text if text else None
                elif child.data == "chiral":
                    chiral = self._node_text(child)
                elif child.data == "h_count":
                    hcount = self._parse_h_count(child)
                elif child.data == "atom_charge":
                    charge = self._parse_charge(child)
                elif child.data == "atom_class":
                    value = self._find_token(child, "INT")
                    clazz = int(value) if value else None

        return "bracket", symbol, isotope, chiral, hcount, charge, clazz

    def _parse_h_count(self, node: Tree) -> int:
        token = self._find_token(node, "INT")
        return int(token) if token else 1

    def _parse_charge(self, node: Tree) -> int | None:
        value = self._node_text(node)
        if not value:
            sign = "+"
            token = self._find_token(node, "INT")
            if token is None:
                return None
            amount = int(token)
            return amount if sign == "+" else -amount

        if value in ("+", "-"):
            return 1 if value == "+" else -1
        if value in ("++", "--"):
            return 2 if value == "++" else -2

        sign = 1 if value[0] == "+" else -1
        magnitude = int(value[1:]) if len(value) > 1 else 1
        return sign * magnitude

    def _build_ring_bonds(self, nodes: Sequence[Tree]) -> list[IRRingBond]:
        return [self._build_ring_bond(node) for node in nodes]

    def _build_ring_bond(self, node: Tree) -> IRRingBond:
        bond: IRBond | None = None
        digits: list[str] = []

        for child in node.children:
            if isinstance(child, Tree) and child.data == "bond_symbol":
                bond = self._build_bond(child)
            elif isinstance(child, Token) and child.type in {"DIGIT", "INT"}:
                digits.append(child.value)

        if not digits:
            text = self._node_text(node)
            digits = [ch for ch in text if ch.isdigit()]

        index = int("".join(digits)) if digits else 0

        return IRRingBond(
            bond=bond,
            index=index,
            span=self._span(node),
        )

    def _build_branches(self, nodes: Sequence[Tree]) -> list[IRBranch]:
        return [self._build_branch(node) for node in nodes]

    def _build_branch(self, node: Tree) -> IRBranch:
        bond: IRBond | None = None
        content: list[object] = []

        for child in node.children:
            if isinstance(child, Tree) and child.data == "bond_symbol":
                bond = self._build_bond(child)
            elif isinstance(child, Tree) and child.data == "molecule_body":
                content.extend(self._build_molecule_body(child))

        branch = IRBranch(bond=bond, content=content, span=self._span(node))
        first_atom = next((part for part in content if isinstance(part, IRAtom)), None)
        if isinstance(first_atom, IRAtom) and first_atom.bond_to_prev is None:
            first_atom.bond_to_prev = bond
        return branch

    # ---------- Stochastic objects ----------

    def _build_stochastic_object(self, node: Tree) -> IRBigStochasticBlock:
        left_descriptor: IRBigConnector | None = None
        right_descriptor: IRBigConnector | None = None
        repeat_units: list[IRBigSmilesMolecule] = []
        end_group: list[IRBigSmilesMolecule] | None = None
        distribution: IRBigDistribution | None = None

        for child in node.children:
            if not isinstance(child, Tree):
                continue
            if child.data == "terminal_bond_descriptor":
                descriptor = self._build_terminal_descriptor(child)
                if left_descriptor is None:
                    left_descriptor = descriptor
                else:
                    right_descriptor = descriptor
            elif child.data == "repeat_units":
                repeat_units = self._build_repeat_units(child)
            elif child.data == "end_group":
                end_group = self._build_end_group(child)
            elif child.data == "stochastic_generation":
                distribution = self._build_stochastic_distribution(child)

        if left_descriptor is None:
            left_descriptor = IRBigConnector(
                mode="simple",
                connector=None,
                index=None,
                span=self._span(node),
                generation=None,
                inner=None,
                noncovalent=None,
            )
        if right_descriptor is None:
            right_descriptor = IRBigConnector(
                mode="simple",
                connector=None,
                index=None,
                span=self._span(node),
                generation=None,
                inner=None,
                noncovalent=None,
            )

        self.features.add("bigsmiles")
        if distribution is not None:
            self.features.add("gbigsmiles")

        return IRBigStochasticBlock(
            left_terminal=left_descriptor,
            right_terminal=right_descriptor,
            units=repeat_units,
            end_group=end_group,
            distribution=distribution,
            span=self._span(node),
        )

    def _build_repeat_units(self, node: Tree) -> list[IRBigSmilesMolecule]:
        molecules: list[IRBigSmilesMolecule] = []
        for child in node.children:
            if isinstance(child, Tree):
                if child.data == "smiles":
                    molecules.append(self._inline_molecule_from_smiles(child))
                elif child.data == "molecule_body":
                    parts = self._build_molecule_body(child)
                    molecules.append(IRBigSmilesMolecule(parts=parts, span=self._span(child)))
                elif child.data == "monomer_list":
                    molecules.extend(self._build_repeat_units(child))
        return molecules

    def _build_end_group(self, node: Tree) -> list[IRBigSmilesMolecule]:
        molecules: list[IRBigSmilesMolecule] = []
        for child in node.children:
            if isinstance(child, Tree):
                if child.data == "smiles":
                    molecules.append(self._inline_molecule_from_smiles(child))
                elif child.data == "molecule_body":
                    parts = self._build_molecule_body(child)
                    molecules.append(IRBigSmilesMolecule(parts=parts, span=self._span(child)))
                elif child.data == "monomer_list":
                    molecules.extend(self._build_repeat_units(child))
        return molecules

    def _inline_molecule_from_smiles(self, node: Tree) -> IRBigSmilesMolecule:
        parts = self._build_smiles(node)
        return IRBigSmilesMolecule(parts=parts, span=self._span(node))

    # ---------- Bond descriptors ----------

    def _build_terminal_descriptor(self, node: Tree) -> IRBigConnector:
        connector: str | None = None
        index: int | None = None
        generation: IRBigGeneration | None = None

        for child in node.children:
            if isinstance(child, Tree) and child.data == "bond_descriptor_symbol_idx":
                symbol = self._find_token(child, "INT")
                raw = self._node_text(child)
                connector = raw[:-len(symbol)] if symbol else raw
                if symbol:
                    index = int(symbol)
                else:
                    index = None
                if connector:
                    connector = connector.strip()
                if connector == "":
                    connector = None
            elif isinstance(child, Tree) and child.data == "bond_descriptor_generation":
                generation = self._build_generation(child)

        descriptor = IRBigConnector(
            mode="simple",
            connector=self._as_connector(connector),
            index=index,
            generation=generation,
            inner=None,
            noncovalent=None,
            span=self._span(node),
        )
        if connector in {"$", "<", ">"}:
            self.features.add("bigsmiles")
        return descriptor

    def _build_bond_descriptor(self, node: Tree) -> IRBigConnector:
        if node.data == "bond_descriptor":
            child = next((c for c in node.children if isinstance(c, Tree)), None)
            if child is None:
                return IRBigConnector(
                    mode="simple",
                    connector=None,
                    index=None,
                    generation=None,
                    inner=None,
                    noncovalent=None,
                    span=self._span(node),
                )
            return self._build_bond_descriptor(child)

        if node.data == "simple_bond_descriptor":
            inner = next((c for c in node.children if isinstance(c, Tree)), None)
            return self._build_simple_descriptor(inner, node)

        if node.data == "non_covalent_bond_descriptor":
            inner = next((c for c in node.children if isinstance(c, Tree)), None)
            descriptor = self._build_non_covalent_descriptor(inner, node)
            self.features.add("gbigsmiles")
            return descriptor

        if node.data == "ladder_bond_descriptor":
            parts: list[IRBigConnector] = []
            index: int | None = None
            for child in node.children:
                if isinstance(child, Tree):
                    if child.data in {
                        "simple_bond_descriptor",
                        "non_covalent_bond_descriptor",
                        "bond_descriptor",
                        "inner_bond_descriptor",
                        "inner_non_covalent_descriptor",
                        "inner_ambi_covalent_descriptor",
                    }:
                        parts.append(self._build_bond_descriptor(child))
                    elif child.data == "bond_descriptor_symbol_idx":
                        token = self._find_token(child, "INT")
                        index = int(token) if token else None
                elif isinstance(child, Token) and child.type == "INT":
                    index = int(child.value)

            descriptor = IRBigConnector(
                mode="ladder",
                connector=None,
                index=index,
                generation=None,
                inner=parts,
                noncovalent=None,
                span=self._span(node),
            )
            self.features.add("bigsmiles")
            if any(part.noncovalent for part in parts):
                self.features.add("gbigsmiles")
            return descriptor

        if node.data in {"inner_bond_descriptor", "bond_descriptor_symbol_idx"}:
            return self._build_simple_descriptor(node, node)

        if node.data == "inner_non_covalent_descriptor":
            descriptor = self._build_non_covalent_descriptor(node, node)
            self.features.add("gbigsmiles")
            return descriptor

        return IRBigConnector(
            mode="simple",
            connector=None,
            index=None,
            generation=None,
            inner=None,
            noncovalent=None,
            span=self._span(node),
        )

    def _build_simple_descriptor(self, node: Tree | None, span_node: Tree) -> IRBigConnector:
        connector: str | None = None
        index: int | None = None
        generation: IRBigGeneration | None = None

        if node is not None:
            for child in node.children:
                if isinstance(child, Tree) and child.data == "bond_descriptor_symbol_idx":
                    text = self._node_text(child)
                    token = self._find_token(child, "INT")
                    connector = text[:-len(token)] if token else text
                    if connector:
                        connector = connector.strip()
                    index = int(token) if token else None
                elif isinstance(child, Tree) and child.data == "bond_descriptor_generation":
                    generation = self._build_generation(child)

        descriptor = IRBigConnector(
            mode="simple",
            connector=self._as_connector(connector),
            index=index,
            generation=generation,
            inner=None,
            noncovalent=None,
            span=self._span(span_node),
        )

        if descriptor.connector in {"$", "<", ">"}:
            self.features.add("bigsmiles")
        return descriptor

    def _build_non_covalent_descriptor(self, node: Tree | None, span_node: Tree) -> IRBigConnector:
        connector = None
        index: int | None = None
        label: int | None = None
        context: IRBigIndexContext | None = None

        if node is not None:
            for child in node.children:
                if isinstance(child, Tree):
                    if child.data == "bond_descriptor_symbol":
                        connector = self._node_text(child) or ":"
                    elif child.data == "INT":
                        label = int(self._node_text(child))
                    elif child.data == "non_covalent_context":
                        context = self._build_non_covalent_context(child)
                elif isinstance(child, Token) and child.type == "INT":
                    label = int(child.value)

        index = label if label is not None else index

        noncovalent = IRBigNonCovalent(
            label=label,
            context=context,
            span=self._span(node) if node is not None else self._span(span_node),
        )

        descriptor = IRBigConnector(
            mode="non_covalent",
            connector=":",
            index=index,
            generation=None,
            inner=None,
            noncovalent=noncovalent,
            span=self._span(span_node),
        )
        self.features.add("gbigsmiles")
        return descriptor

    def _build_non_covalent_context(self, node: Tree) -> IRBigIndexContext:
        expr: IRBigIndexExpr | None = None
        kv: dict[str, str] = {}

        for child in node.children:
            if isinstance(child, Tree) and child.data == "index_expression" and expr is None:
                expr = self._build_index_expression(child)
            elif isinstance(child, Tree) and child.data == "non_covalent_key_value_pair":
                key = self._read_text_block(child, "non_covalent_key")
                value = self._read_text_block(child, "non_covalent_value")
                if key:
                    kv[key] = value

        if expr is None:
            expr = IRBigIndexExpr(left=0, op=None, right=None, unary=None, span=self._span(node))

        if kv:
            self.features.add("gbigsmiles")

        return IRBigIndexContext(expr=expr, kv=kv, span=self._span(node))

    def _read_text_block(self, node: Tree, name: str) -> str:
        target = next((c for c in node.children if isinstance(c, Tree) and c.data == name), None)
        if target is None:
            return ""
        chars: list[str] = []
        for child in target.children:
            if isinstance(child, Token):
                chars.append(child.value)
        if not chars:
            text = self._node_text(target)
            if text:
                return text
        return "".join(chars)

    # ---------- Index expressions ----------

    def _build_index_expression(self, node: Tree) -> IRBigIndexExpr:
        # Robust text-based parser for index expressions: handles parentheses, !, ~, &
        text = self._node_text(node)
        return self._parse_index_expr_text(text, self._span(node))

    def _build_index_statement(self, node: Tree | Token) -> IRBigIndexExpr | int:
        if isinstance(node, Token):
            if node.type == "INT":
                return int(node.value)
            text = node.value if hasattr(node, "value") else ""
            return int(text) if text.isdigit() else 0

        children = [child for child in node.children if isinstance(child, (Tree, Token))]
        if not children:
            return 0

        first = children[0]
        if isinstance(first, Tree) and first.data == "unary_index_operator":
            unary = self._node_text(first)
            rest = self._build_index_statement(children[1])
            return IRBigIndexExpr(left=rest, op=None, right=None, unary=("!" if unary == "!" else None), span=self._span(node))

        # handle binary op when operator token is directly present
        if len(children) >= 3 and isinstance(children[1], Token) and children[1].value in {"&", "~"}:
            return self._build_index_binary(node)

        if len(children) == 1 and isinstance(children[0], Token):
            return int(children[0].value)

        if len(children) == 1 and isinstance(children[0], Tree):
            child = children[0]
            if child.data in {"branched_index_expression", "unbranched_index_expression"}:
                return self._build_index_binary(child)
            if child.data == "index_statement":
                return self._build_index_statement(child)

        if isinstance(first, Tree) and first.data in {"branched_index_expression", "unbranched_index_expression"}:
            return self._build_index_binary(first)

        if len(children) >= 3 and isinstance(children[1], Tree) and children[1].data == "binary_index_operator":
            return self._build_index_binary(node)

        token = self._find_token(node, "INT")
        return int(token) if token else 0

    def _build_index_binary(self, node: Tree) -> IRBigIndexExpr:
        items = [child for child in node.children if isinstance(child, (Tree, Token))]
        if len(items) < 3:
            return IRBigIndexExpr(left=0, op=None, right=None, unary=None, span=self._span(node))

        left = self._build_index_statement(items[0])
        op_token = items[1]
        op_text = self._node_text(op_token)
        op: Literal['~','&'] | None = "~" if op_text == "~" else ("&" if op_text == "&" else None)
        right = self._build_index_statement(items[2])

        return IRBigIndexExpr(
            left=left,
            op=op,
            right=right,
            unary=None,
            span=self._span(node),
        )

    # ---------- Text parsing for index expressions ----------

    def _parse_index_expr_text(self, s: str, span: Span | None) -> IRBigIndexExpr:
        i = 0

        def peek() -> str:
            return s[i] if i < len(s) else ""

        def consume() -> str:
            nonlocal i
            ch = peek()
            i += 1
            return ch

        def skip_ws() -> None:
            nonlocal i
            while i < len(s) and s[i].isspace():
                i += 1

        def parse_int() -> int:
            nonlocal i
            skip_ws()
            start = i
            while i < len(s) and s[i].isdigit():
                i += 1
            if start == i:
                return 0
            return int(s[start:i])

        def parse_term():
            skip_ws()
            if peek() == "!":
                consume()
                inner = parse_term()
                return IRBigIndexExpr(left=inner, op=None, right=None, unary="!", span=span)
            if peek() == "(":
                consume()
                expr = parse_expr()
                if peek() == ")":
                    consume()
                return expr
            # number
            return parse_int()

        def parse_expr():
            left = parse_term()
            skip_ws()
            while True:
                opch = peek()
                if opch not in {"~", "&"}:
                    break
                consume()
                right = parse_term()
                op_lit: Literal['~','&'] | None = "~" if opch == "~" else "&"
                left = IRBigIndexExpr(left=left, op=op_lit, right=right, unary=None, span=span)
                skip_ws()
            return left

        # strip surrounding pipes if present (defensive)
        if s and s[0] == "|" and s[-1:] == "|":
            s = s[1:-1]
        result = parse_expr()
        if isinstance(result, int):
            return IRBigIndexExpr(left=result, op=None, right=None, unary=None, span=span)
        return result

    # ---------- Generations & distributions ----------

    def _build_generation(self, node: Tree) -> IRBigGeneration:
        values: list[float] = []
        for child in node.children:
            if isinstance(child, Token) and child.type == "NUMBER":
                values.append(float(child.value))
        return IRBigGeneration(values=values, span=self._span(node))

    def _build_stochastic_distribution(self, node: Tree) -> IRBigDistribution | None:
        distribution = next((c for c in node.children if isinstance(c, Tree) and c.data == "stochastic_distribution"), None)
        if distribution is None:
            return None
        variant = next((c for c in distribution.children if isinstance(c, Tree)), None)
        target = variant or distribution
        kind = str(variant.data) if variant is not None else str(distribution.data)
        params = [
            float(tok.value)
            for tok in target.children
            if isinstance(tok, Token) and tok.type == "NUMBER"
        ]
        p1 = params[0] if params else None
        p2 = params[1] if len(params) > 1 else None

        self.features.add("gbigsmiles")

        return IRBigDistribution(
            kind=self._as_distribution_kind(kind),
            p1=p1,
            p2=p2,
            span=self._span(target),
        )

    # ---------- Fragments & dots ----------

    def _build_fragment_definition(self, node: Tree) -> IRBigFragmentDef:
        name_node = next((c for c in node.children if isinstance(c, Tree) and c.data == "fragment_name"), None)
        name = self._node_text(name_node) if name_node is not None else ""
        molecules: list[IRBigSmilesMolecule] = []
        for child in node.children:
            if isinstance(child, Tree) and child.data == "molecule":
                molecules.append(self._build_molecule_big(child))

        self.features.add("bigsmiles")

        return IRBigFragmentDef(name=name, molecules=molecules, span=self._span(node))

    def _build_fragment_ref(self, node: Tree) -> IRBigFragmentRef:
        name_node = next((c for c in node.children if isinstance(c, Tree) and c.data == "fragment_name"), None)
        name = self._node_text(name_node) if name_node is not None else ""
        return IRBigFragmentRef(name=name, span=self._span(node))

    def _build_dot_generation(self, node: Tree) -> IRBigSystemDot:
        system_size: float | None = None
        size_node = next((c for c in node.children if isinstance(c, Tree) and c.data == "dot_system_size"), None)
        if size_node is not None:
            token = self._find_token(size_node, "NUMBER")
            system_size = float(token) if token else None
        return IRBigSystemDot(system_size=system_size, span=self._span(node))

    # ---------- Helpers ----------

    def _build_bond(self, node: Tree | Token | None) -> IRBond | None:
        if node is None:
            return None
        text = self._node_text(node)
        text = text.strip()
        if not text:
            return None
        return IRBond(kind=self._as_bond_kind(text))

    def _make_placeholder_atom(self, node: Tree) -> IRAtom:
        return IRAtom(
            kind="aliphatic",
            symbol=None,
            isotope=None,
            chiral=None,
            hcount=None,
            charge=None,
            clazz=None,
            rings=[],
            branches=[],
            span=self._span(node),
        )

    def _find_token(self, node: Tree, token_type: str) -> str | None:
        for child in node.children:
            if isinstance(child, Token) and child.type == token_type:
                return child.value
        return None

    def _node_text(self, node: Tree | Token | None) -> str:
        if node is None:
            return ""
        if isinstance(node, Token):
            value = getattr(node, "value", None)
            return str(value) if value is not None else ""
        meta = getattr(node, "meta", None)
        if meta and meta.start_pos is not None and meta.end_pos is not None:
            # Prefer slicing from the current input text if available
            if isinstance(self._input_text, str) and meta.end_pos <= len(self._input_text):
                return self._input_text[meta.start_pos : meta.end_pos]
            # Fallback: return empty string if input text isn't available
            return ""
        return ""

    def _span(self, node: Tree | Token | None) -> Span | None:
        if node is None:
            return None
        if isinstance(node, Token):
            start = getattr(node, "pos_in_stream", None)
            end = getattr(node, "end_pos", None)
            line = getattr(node, "line", None)
            column = getattr(node, "column", None)
            if None in (start, end, line, column):
                return None
            assert isinstance(start, int) and isinstance(end, int) and isinstance(line, int) and isinstance(column, int)
            return Span(start=start, end=end, line=line, column=column)
        meta = getattr(node, "meta", None)
        if meta is None:
            return None
        start = getattr(meta, "start_pos", None)
        end = getattr(meta, "end_pos", None)
        line = getattr(meta, "line", None)
        column = getattr(meta, "column", None)
        if None in (start, end, line, column):
            return None
        assert isinstance(start, int) and isinstance(end, int) and isinstance(line, int) and isinstance(column, int)
        return Span(start=start, end=end, line=line, column=column)

    # ---------- Type adaptation helpers ----------

    def _as_connector(self, value: str | None) -> Literal['$', '<', '>', ':'] | None:
        if value in {"$", "<", ">", ":"}:
            return cast(Literal['$', '<', '>', ':'], value)
        return None

    def _as_bond_kind(self, text: str) -> Literal['-', '=', '#', ';', ':', '/', '\\']:
        if text in {"-", "=", "#", ";", ":", "/", "\\"}:
            return cast(Literal['-', '=', '#', ';', ':', '/', '\\'], text)
        # Fallback to single bond
        return "-"

    def _as_distribution_kind(self, name: str) -> "Literal['flory_schulz', 'schulz_zimm', 'gauss', 'uniform', 'log_normal', 'poisson']":
        n = name.lower()
        if "flory" in n:
            return "flory_schulz"
        if "schulz" in n:
            return "schulz_zimm"
        if "gauss" in n:
            return "gauss"
        if "uniform" in n:
            return "uniform"
        if "log" in n:
            return "log_normal"
        if "poisson" in n:
            return "poisson"
        return "uniform"
