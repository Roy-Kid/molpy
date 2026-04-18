"""Hand-written tokenizer + recursive-descent parser for MolTemplate (.lt).

Design note:
    Real moltemplate files interleave structured directives (``ClassName {``,
    ``write_once("Section") {``, ``new Foo.move(...)``) with free-form LAMMPS
    coeff lines inside the block bodies. A pure Lark grammar is awkward
    because the body grammar differs per section. Instead we parse structure
    via a tokenizer that tracks brace depth, capture block bodies as raw text
    lines, and let ``builder.py`` regex-parse them per section.
"""

from __future__ import annotations

import re
from pathlib import Path

from .ir import (
    ArrayDim,
    ClassDef,
    Document,
    ImportStmt,
    NewStmt,
    Statement,
    Transform,
    WriteBlock,
    WriteOnceBlock,
)


# ---------------------------------------------------------------------------
# Tokenizer
# ---------------------------------------------------------------------------

_TOKEN_RE = re.compile(
    r"""
      (?P<WS>\s+)
    | (?P<COMMENT>\#[^\n]*)
    | (?P<STRING>"[^"]*")
    | (?P<LBRACE>\{)
    | (?P<RBRACE>\})
    | (?P<LPAREN>\()
    | (?P<RPAREN>\))
    | (?P<LBRACK>\[)
    | (?P<RBRACK>\])
    | (?P<COMMA>,)
    | (?P<DOT>\.)
    | (?P<EQ>=)
    | (?P<NUMBER>-?\d+\.\d+(?:[eE][+-]?\d+)?|-?\d+(?:[eE][+-]?\d+)?)
    | (?P<IDENT>[A-Za-z_][A-Za-z_0-9/\-:@$]*)
    | (?P<OTHER>.)
    """,
    re.VERBOSE,
)

_STRUCTURED_KEYWORDS = {"write", "write_once", "import", "new", "inherits"}


class Token:
    __slots__ = ("kind", "value", "line")

    def __init__(self, kind: str, value: str, line: int):
        self.kind = kind
        self.value = value
        self.line = line

    def __repr__(self) -> str:
        return f"Token({self.kind!r}, {self.value!r}, line={self.line})"


def tokenize(source: str) -> list[Token]:
    """Lex ``source`` into a list of tokens (whitespace/comments dropped)."""
    tokens: list[Token] = []
    line = 1
    for m in _TOKEN_RE.finditer(source):
        kind = m.lastgroup or "OTHER"
        value = m.group()
        if kind == "WS":
            line += value.count("\n")
            continue
        if kind == "COMMENT":
            continue
        if kind == "STRING":
            value = value[1:-1]
        tokens.append(Token(kind, value, line))
    tokens.append(Token("EOF", "", line))
    return tokens


# ---------------------------------------------------------------------------
# Brace-aware body extraction
# ---------------------------------------------------------------------------

def _extract_block_body(source: str, start_idx: int) -> tuple[list[str], int]:
    """Extract raw text lines inside a ``{ ... }`` block starting at ``start_idx``.

    ``source[start_idx]`` must be '{'. Returns (lines, end_idx_exclusive) where
    ``source[end_idx_exclusive - 1] == '}'`` and ``lines`` is the stripped-
    line body (no leading/trailing whitespace on each line, empty lines kept
    if user-written).
    """
    assert source[start_idx] == "{", f"expected '{{' at {start_idx}"
    depth = 1
    i = start_idx + 1
    body_start = i
    while i < len(source) and depth > 0:
        c = source[i]
        if c == "{":
            depth += 1
        elif c == "}":
            depth -= 1
            if depth == 0:
                break
        elif c == "#":
            # skip to end of line
            nl = source.find("\n", i)
            if nl == -1:
                i = len(source)
                break
            i = nl
            continue
        elif c == '"':
            # skip over string literal
            close = source.find('"', i + 1)
            if close == -1:
                raise SyntaxError(
                    f"unterminated string starting at offset {i}"
                )
            i = close + 1
            continue
        i += 1
    if depth != 0:
        raise SyntaxError(f"unbalanced '{{' starting at {start_idx}")
    body = source[body_start:i]
    lines = [
        ln.strip() for ln in body.splitlines() if ln.strip() and not ln.strip().startswith("#")
    ]
    return lines, i + 1


# ---------------------------------------------------------------------------
# Parser
# ---------------------------------------------------------------------------

class MolTemplateParser:
    """Parse MolTemplate source into a ``Document`` IR.

    The parser uses a two-pass strategy:
      1. Scan the source character-by-character to find top-level statements,
         extracting block bodies via :func:`_extract_block_body`.
      2. Tokenize the structural prefix of each statement to identify
         keywords and names.

    This keeps brace-depth handling simple and avoids grammar explosion
    for the free-form content inside ``write``/``write_once`` blocks.
    """

    def parse(self, source: str) -> Document:
        doc = Document()
        doc.statements = self._parse_block(source, 0, len(source))
        return doc

    def _parse_block(self, source: str, start: int, end: int) -> list[Statement]:
        """Parse statements in ``source[start:end]`` (top-level or inside a class)."""
        statements: list[Statement] = []
        i = start
        while i < end:
            # skip whitespace
            while i < end and source[i].isspace():
                i += 1
            if i >= end:
                break
            # comments
            if source[i] == "#":
                nl = source.find("\n", i)
                i = end if nl == -1 else nl
                continue
            # import "file.lt"
            if source.startswith("import", i) and _is_word_boundary(source, i, i + 6):
                stmt, i = self._parse_import(source, i + 6, end)
                statements.append(stmt)
                continue
            # write(...) { ... }
            if source.startswith("write_once", i) and _is_word_boundary(source, i, i + 10):
                stmt, i = self._parse_write(source, i + 10, end, once=True)
                statements.append(stmt)
                continue
            if source.startswith("write", i) and _is_word_boundary(source, i, i + 5):
                stmt, i = self._parse_write(source, i + 5, end, once=False)
                statements.append(stmt)
                continue
            # instance = new ClassName
            # Look ahead: IDENT '=' 'new' ...
            m = re.match(r"([A-Za-z_][\w:$@]*)\s*=\s*new\b", source[i:end])
            if m:
                stmt, i = self._parse_new(source, i, end)
                statements.append(stmt)
                continue
            # class definition: IDENT (inherits ...)? { ... }
            m = re.match(r"([A-Za-z_][\w]*)\s*(?:inherits\b([^{]*))?\s*\{", source[i:end])
            if m and source[i + m.end() - 1] == "{":
                stmt, i = self._parse_class(source, i, end)
                statements.append(stmt)
                continue
            # Unknown construct — skip to next line to stay robust
            nl = source.find("\n", i)
            i = end if nl == -1 else nl + 1
        return statements

    def _parse_import(self, source: str, i: int, end: int) -> tuple[ImportStmt, int]:
        # Expect: whitespace then "path"
        while i < end and source[i].isspace():
            i += 1
        if i >= end or source[i] != '"':
            raise SyntaxError(f"import: expected quoted path near offset {i}")
        close = source.find('"', i + 1)
        if close == -1:
            raise SyntaxError("import: unterminated string")
        path = source[i + 1:close]
        return ImportStmt(path=path), close + 1

    def _parse_write(
        self, source: str, i: int, end: int, *, once: bool
    ) -> tuple[WriteBlock | WriteOnceBlock, int]:
        # Expect: '(' "section name (may contain parens)" ')' '{' ... '}'
        while i < end and source[i].isspace():
            i += 1
        if source[i] != "(":
            raise SyntaxError(f"write: expected '(' at offset {i}")
        # Skip the opening '('
        k = i + 1
        # Skip leading whitespace
        while k < end and source[k].isspace():
            k += 1
        if k < end and source[k] in ('"', "'"):
            # Quoted section: find closing quote, then the matching ')'
            quote = source[k]
            close_q = source.find(quote, k + 1)
            if close_q == -1:
                raise SyntaxError("write: unterminated quoted section name")
            section = source[k + 1:close_q]
            # Now find the next ')' after the closing quote
            close_paren = source.find(")", close_q + 1)
            if close_paren == -1:
                raise SyntaxError("write: unterminated '('")
        else:
            # Unquoted: find next ')' at paren-depth 0
            close_paren = source.find(")", k)
            if close_paren == -1:
                raise SyntaxError("write: unterminated '('")
            section = source[k:close_paren].strip()
        j = close_paren + 1
        while j < end and source[j].isspace():
            j += 1
        if j >= end or source[j] != "{":
            raise SyntaxError(f"write: expected '{{' at offset {j}")
        body, end_idx = _extract_block_body(source, j)
        if once:
            return WriteOnceBlock(section=section, body_lines=body), end_idx
        return WriteBlock(section=section, body_lines=body), end_idx

    def _parse_new(self, source: str, i: int, end: int) -> tuple[NewStmt, int]:
        m = re.match(
            r"([A-Za-z_][\w:$@]*)\s*=\s*new\s*(?:\[\s*(\d+)\s*\])?\s*([A-Za-z_][\w]*)",
            source[i:end],
        )
        if m is None:
            raise SyntaxError(f"new: parse failure near offset {i}")
        instance = m.group(1)
        count = int(m.group(2)) if m.group(2) else 1
        cls_name = m.group(3)
        cursor = i + m.end()
        transforms: list[Transform] = []
        arrays: list[ArrayDim] = []
        while cursor < end:
            # skip whitespace (incl. newlines — chains may wrap)
            k = cursor
            while k < end and source[k].isspace():
                k += 1
            if k >= end:
                break
            c = source[k]
            # Post-class array dimension: [N] or [N].move(...)
            if c == "[":
                am = re.match(
                    r"\[\s*(\d+)\s*\]\s*(?:\.\s*([A-Za-z_][\w]*)\s*\(\s*([^)]*)\)\s*)?",
                    source[k:end],
                )
                if am is None:
                    break
                n = int(am.group(1))
                tr: Transform | None = None
                if am.group(2) is not None:
                    op = am.group(2)
                    args_str = (am.group(3) or "").strip()
                    args: list[float] = []
                    if args_str:
                        for part in args_str.split(","):
                            part = part.strip()
                            if part:
                                args.append(float(part))
                    tr = Transform(op=op, args=args)
                arrays.append(ArrayDim(count=n, transform=tr))
                cursor = k + am.end()
                continue
            # Per-instance transform chain: .move(...) / .rot(...) / etc.
            if c == ".":
                tm = re.match(
                    r"\.\s*([A-Za-z_][\w]*)\s*\(\s*([^)]*)\)",
                    source[k:end],
                )
                if tm is None:
                    break
                op = tm.group(1)
                args_str = tm.group(2).strip()
                args = []
                if args_str:
                    for part in args_str.split(","):
                        part = part.strip()
                        if part:
                            args.append(float(part))
                transforms.append(Transform(op=op, args=args))
                cursor = k + tm.end()
                continue
            break
        return NewStmt(
            instance_name=instance,
            class_name=cls_name,
            count=count,
            transforms=transforms,
            arrays=arrays,
        ), cursor

    def _parse_class(self, source: str, i: int, end: int) -> tuple[ClassDef, int]:
        m = re.match(
            r"([A-Za-z_][\w]*)\s*(?:inherits\s+([^{]+))?\s*\{",
            source[i:end],
        )
        if m is None:
            raise SyntaxError(f"class: parse failure near offset {i}")
        name = m.group(1)
        bases_raw = m.group(2) or ""
        bases = [
            b.strip() for b in re.split(r"[,\s]+", bases_raw.strip()) if b.strip()
        ]
        brace_idx = i + m.end() - 1
        # Find matching '}' while respecting nested braces
        depth = 1
        j = brace_idx + 1
        while j < end and depth > 0:
            c = source[j]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
                if depth == 0:
                    break
            elif c == "#":
                nl = source.find("\n", j)
                j = end if nl == -1 else nl
                continue
            elif c == '"':
                close = source.find('"', j + 1)
                if close == -1:
                    raise SyntaxError("class: unterminated string")
                j = close + 1
                continue
            j += 1
        if depth != 0:
            raise SyntaxError(f"class '{name}': unbalanced braces")
        body_start = brace_idx + 1
        body_end = j
        inner_stmts = self._parse_block(source, body_start, body_end)
        return ClassDef(name=name, bases=bases, statements=inner_stmts), j + 1


def _is_word_boundary(source: str, start: int, end: int) -> bool:
    """Return True if ``source[end]`` is not an identifier continuation char."""
    if end >= len(source):
        return True
    c = source[end]
    return not (c.isalnum() or c == "_")


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def parse_string(source: str) -> Document:
    """Parse MolTemplate source text into a ``Document`` IR tree."""
    return MolTemplateParser().parse(source)


def parse_file(path: str | Path) -> Document:
    """Parse a ``.lt`` file from disk."""
    text = Path(path).read_text()
    return parse_string(text)
