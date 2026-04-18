"""Intermediate representation (IR) dataclasses for MolTemplate parsing.

The parser produces a tree of these nodes; the builder walks the tree to
materialise ``ForceField`` and ``Atomistic`` objects.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Union


@dataclass
class Transform:
    """A chained transform call on a ``new`` instance (``.move(...)``, etc.)."""

    op: str  # "move" | "rot" | "rotvv" | "scale" | "quat"
    args: list[float] = field(default_factory=list)


@dataclass
class NewStmt:
    """``inst = new ClassName[.move(...)...]`` — molecule instantiation."""

    instance_name: str
    class_name: str
    count: int = 1  # for array form: inst = new [N] ClassName
    transforms: list[Transform] = field(default_factory=list)


@dataclass
class WriteBlock:
    """``write("Section Name") { ... body lines ... }``."""

    section: str
    body_lines: list[str] = field(default_factory=list)


@dataclass
class WriteOnceBlock:
    """``write_once("Section Name") { ... }`` — identical to WriteBlock semantically."""

    section: str
    body_lines: list[str] = field(default_factory=list)


@dataclass
class ImportStmt:
    """``import "file.lt"``."""

    path: str


# Forward declaration for recursive ClassDef.statements
Statement = Union[
    "ClassDef", NewStmt, WriteBlock, WriteOnceBlock, ImportStmt
]


@dataclass
class ClassDef:
    """``ClassName [inherits Base1, Base2] { ... statements ... }``."""

    name: str
    bases: list[str] = field(default_factory=list)
    statements: list[Statement] = field(default_factory=list)


@dataclass
class Document:
    """Top-level parsed document — a sequence of statements."""

    statements: list[Statement] = field(default_factory=list)
