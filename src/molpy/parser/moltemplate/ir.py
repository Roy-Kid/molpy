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
class ArrayDim:
    """One dimension of the ``new Cls [N].move(dx, dy, dz)`` array form.

    ``count`` is N (number of copies along this dimension); ``transform`` is
    applied repeatedly — copy ``k`` is shifted by ``k * args`` of the
    transform. If ``transform`` is ``None`` the copies are placed at the
    same position (rare but legal per moltemplate).
    """

    count: int
    transform: Transform | None = None


@dataclass
class NewStmt:
    """``inst = new ClassName[.move(...)...]`` — molecule instantiation.

    Supports post-class multi-dimensional arrays, e.g.::

        m = new Butane [12].move(0, 0, 5.2) [12].move(0, 5.2, 0) [6].move(10.4, 0, 0)

    which produces ``12 * 12 * 6 = 864`` copies on a regular grid.
    ``transforms`` is the per-instance transform chain applied *before*
    array expansion. ``arrays`` is the list of dimensions (empty when no
    ``[N]`` form was written).
    """

    instance_name: str
    class_name: str
    count: int = 1  # legacy single-count form: `new [N] Cls`
    transforms: list[Transform] = field(default_factory=list)
    arrays: list[ArrayDim] = field(default_factory=list)


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
