"""
Utilities for parsing (G-)BigSMILES strings using the validated Lark grammar
in ``grammar/big_smiles.lark``.  Rather than constructing a full molecular IR,
we expose a light-weight summary that is convenient for tests and downstream
tools: repeat units extracted from stochastic objects and the metadata carried
by non-covalent descriptors.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

from lark import Tree, Token

from .base import GrammarConfig, GrammarParserBase


# --------------------------------------------------------------------------- #
# Dataclasses
# --------------------------------------------------------------------------- #


@dataclass(slots=True)
class NonCovalentDescriptor:
    """Non-covalent descriptor such as ``[>:|1,Type=HB|]``."""

    symbol: str
    index: int | None
    expression: str
    attributes: dict[str, str] = field(default_factory=dict)


@dataclass(slots=True)
class StochasticObjectSummary:
    """Summary of a stochastic object ``{ ... }`` block."""

    left_symbol: str | None
    right_symbol: str | None
    repeat_units: list[str] = field(default_factory=list)
    has_end_group: bool = False
    distribution: str | None = None


@dataclass(slots=True)
class GBigSmilesIR:
    """Aggregate information extracted from a G-BigSMILES string."""

    source: str
    stochastic_objects: list[StochasticObjectSummary] = field(default_factory=list)
    non_covalent_descriptors: list[NonCovalentDescriptor] = field(default_factory=list)


# --------------------------------------------------------------------------- #
# Parser implementation
# --------------------------------------------------------------------------- #


class GBigSmilesParser(GrammarParserBase[GBigSmilesIR]):
    """
    Parser wrapper around the ``big_smiles.lark`` grammar.

    The grammar itself cannot currently accept non-covalent descriptors that
    appear in terminal positions (e.g. ``[>:|1,Type=HB|]``).  To work around
    this without modifying the grammar, we sanitise the input before parsing
    by removing the ``:...`` part and recording the metadata ourselves.
    """

    def __init__(
        self,
        grammar_path: Path | None = None,
        *,
        start: str = "big_smiles",
        auto_reload: bool = True,
    ):
        if grammar_path is None:
            grammar_path = Path(__file__).parent / "grammar/big_smiles.lark"
        config = GrammarConfig(
            grammar_path=grammar_path,
            start=start,
            parser="earley",
            propagate_positions=True,
            maybe_placeholders=False,
            auto_reload=auto_reload,
        )
        super().__init__(config)
        self._source: str = ""
        self._sanitised_text: str = ""
        self._pending_contexts: list[NonCovalentDescriptor] = []

    # Public API ------------------------------------------------------------ #

    def parse(self, text: str) -> GBigSmilesIR:
        self._source = text
        sanitised, contexts = self._sanitise_text(text)
        self._sanitised_text = sanitised
        self._pending_contexts = contexts
        tree = self.parse_tree(sanitised)
        return self.build(tree)

    def build(self, tree: Tree) -> GBigSmilesIR:
        contexts = list(self._pending_contexts)
        stochastic_objects: list[StochasticObjectSummary] = []

        for node in tree.iter_subtrees_topdown():
            if node.data == "stochastic_object":
                stochastic_objects.append(self._summarise_stochastic_object(node))

        return GBigSmilesIR(
            source=self._source,
            stochastic_objects=stochastic_objects,
            non_covalent_descriptors=contexts,
        )

    # ------------------------------------------------------------------ #
    # Sanitisation (pre-parsing)
    # ------------------------------------------------------------------ #

    def _sanitise_text(self, text: str) -> tuple[str, list[NonCovalentDescriptor]]:
        """
        Remove the ``:`` + context segments from non-covalent descriptors so the
        existing grammar can parse the remainder.  The removed metadata is
        returned in a list (left-to-right order).
        """
        contexts: list[NonCovalentDescriptor] = []
        out: list[str] = []
        i = 0
        length = len(text)

        while i < length:
            ch = text[i]
            if ch == "[" and i + 1 < length:
                symbol_char = text[i + 1]
                if symbol_char in (">", "<", "$"):
                    j = i + 2
                    while j < length and text[j].isspace():
                        j += 1
                    if j < length and text[j] == ":":
                        j += 1
                        idx_start = j
                        while j < length and text[j].isdigit():
                            j += 1
                        index = int(text[idx_start:j]) if j > idx_start else None

                        expression = ""
                        attributes: dict[str, str] = {}
                        if j < length and text[j] == "|":
                            ctx_end = self._find_matching_bar(text, j)
                            context_body = text[j + 1 : ctx_end]
                            expression, attributes = self._parse_context_body(context_body)
                            j = ctx_end + 1

                        contexts.append(
                            NonCovalentDescriptor(
                                symbol=symbol_char,
                                index=index,
                                expression=expression,
                                attributes=attributes,
                            )
                        )

                        while j < length and text[j] != "]":
                            j += 1
                        out.append("[")
                        out.append(symbol_char)
                        i = j
                        continue

            out.append(ch)
            i += 1

        return "".join(out), contexts

    @staticmethod
    def _find_matching_bar(text: str, start: int) -> int:
        """Return the index of the closing ``|`` matching the one at ``start``."""
        assert text[start] == "|"
        pos = start + 1
        while pos < len(text):
            if text[pos] == "|" and text[pos - 1] != "\\":
                return pos
            pos += 1
        raise ValueError("Unterminated non-covalent context")

    @staticmethod
    def _parse_context_body(body: str) -> tuple[str, dict[str, str]]:
        parts = [segment.strip() for segment in body.split(",") if segment.strip()]
        expression = parts[0] if parts else ""
        attributes: dict[str, str] = {}
        for part in parts[1:]:
            if "=" in part:
                key, value = part.split("=", 1)
                attributes[key.strip()] = value.strip()
        return expression, attributes

    # ------------------------------------------------------------------ #
    # Tree walking helpers
    # ------------------------------------------------------------------ #

    def _summarise_stochastic_object(self, node: Tree) -> StochasticObjectSummary:
        terminals = [
            child for child in node.children
            if isinstance(child, Tree) and child.data == "terminal_bond_descriptor"
        ]

        left_symbol = self._extract_terminal_symbol(terminals[0]) if terminals else None
        right_symbol = (
            self._extract_terminal_symbol(terminals[-1]) if len(terminals) >= 2 else left_symbol
        )

        repeat_units: list[str] = []
        for child in node.children:
            if isinstance(child, Tree):
                if child.data in {"smiles", "big_smiles_molecule"}:
                    repeat_units.append(self._slice_from_sanitised(child).strip())
                elif child.data == "_monomer_list":
                    repeat_units.extend(self._collect_repeat_units(child))

        has_end_group = any(
            isinstance(child, Tree) and child.data == "_end_group" for child in node.children
        )

        distribution = None
        for child in node.children:
            if isinstance(child, Tree) and child.data == "stochastic_generation":
                for grand in child.children:
                    if isinstance(grand, Tree):
                        distribution = grand.data
                        break

        return StochasticObjectSummary(
            left_symbol=left_symbol,
            right_symbol=right_symbol,
            repeat_units=repeat_units,
            has_end_group=has_end_group,
            distribution=distribution,
        )

    def _extract_terminal_symbol(self, node: Tree) -> str | None:
        # Prefer direct symbol index if present
        target = next(
            (child for child in node.children if isinstance(child, Tree) and child.data == "bond_descriptor_symbol_idx"),
            None,
        )
        if target is not None:
            snippet = self._slice_from_sanitised(target).strip()
            return snippet[0] if snippet else None

        # Handle ladder/non-covalent/simple variants nested under terminal_bond_descriptor
        # Search depth-first for either bond_descriptor_symbol_idx or bond_descriptor_symbol
        def find_symbol_subtree(t: Tree) -> Tree | None:
            for ch in t.children:
                if isinstance(ch, Tree):
                    if ch.data in {"bond_descriptor_symbol_idx", "bond_descriptor_symbol"}:
                        return ch
                    found = find_symbol_subtree(ch)
                    if found is not None:
                        return found
            return None

        sym_tree = find_symbol_subtree(node)
        if sym_tree is not None:
            snippet = self._slice_from_sanitised(sym_tree).strip()
            return snippet[0] if snippet else None
        return None

    def _collect_repeat_units(self, node: Tree) -> list[str]:
        units: list[str] = []
        for child in node.children:
            if isinstance(child, Tree):
                if child.data in {"smiles", "big_smiles_molecule"}:
                    units.append(self._slice_from_sanitised(child).strip())
                elif child.data == "_monomer_list":
                    units.extend(self._collect_repeat_units(child))
        return units

    def _stringify(self, node: Tree | Token | Iterable[Tree | Token]) -> str:
        if isinstance(node, Token):
            return node.value
        if isinstance(node, Tree):
            return "".join(self._stringify(child) for child in node.children)
        return "".join(self._stringify(child) for child in node)

    def _slice_from_sanitised(self, node: Tree | Token | None) -> str:
        if node is None:
            return ""
        meta = getattr(node, "meta", None)
        if meta and meta.start_pos is not None and meta.end_pos is not None:
            return self._sanitised_text[meta.start_pos:meta.end_pos]
        return self._stringify(node)


def parse_g_big_smiles(text: str, grammar_path: Path | None = None) -> GBigSmilesIR:
    """Convenience wrapper mirroring the behaviour of ``GBigSmilesParser.parse``."""
    parser = GBigSmilesParser(grammar_path=grammar_path)
    return parser.parse(text)
