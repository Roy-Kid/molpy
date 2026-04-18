"""Build ``ForceField`` and ``Atomistic`` objects from a MolTemplate IR.

Section directives we understand:

* ``Data Masses``          -- @atom:X  mass  [# comment]
* ``Data Atoms``           -- $atom:id  $mol:m  @atom:T  charge  x  y  z
* ``Data Bonds``           -- $bond:id  @bond:T  $atom:a  $atom:b
* ``Data Angles``          -- $angle:id  @angle:T  $atom:a  $atom:b  $atom:c
* ``Data Dihedrals``       -- $dihedral:id  @dihedral:T  $atom:...x4
* ``Data Impropers``       -- $improper:id  @improper:T  $atom:...x4
* ``Data Charges``         -- set type @atom:X charge q
* ``In Charges``           -- alias of Data Charges
* ``In Settings``          -- pair_coeff / bond_coeff / angle_coeff / dihedral_coeff / improper_coeff
* ``In Init``              -- plain lammps init lines; stored in ForceField.metadata.

All other sections are preserved as raw text in ``ForceField.metadata[section]``
(list of lines) so that downstream emitters can re-export them verbatim.
"""

from __future__ import annotations

import re
from copy import deepcopy
from pathlib import Path
from typing import Any

from molpy.core.atomistic import Atom, Atomistic, Bond, Angle, Dihedral
from molpy.core.forcefield import (
    AtomStyle,
    AtomType,
    ForceField,
    Style,
    Type,
)

from .ir import (
    ArrayDim,
    ClassDef,
    Document,
    ImportStmt,
    NewStmt,
    Transform,
    WriteBlock,
    WriteOnceBlock,
)
from .parser import parse_file


# ---------------------------------------------------------------------------
# Helpers: strip `@` / `$` prefixes
# ---------------------------------------------------------------------------

_AT_RE = re.compile(r"^@(?:atom|bond|angle|dihedral|improper|pair|mol):(.+)$")
_DOLLAR_RE = re.compile(r"^\$(?:atom|bond|angle|dihedral|improper|mol):(.+)$")


def _strip_prefix(token: str) -> str:
    """Strip ``@atom:``/``$atom:`` (and siblings) from a token."""
    m = _AT_RE.match(token) or _DOLLAR_RE.match(token)
    return m.group(1) if m else token


# ---------------------------------------------------------------------------
# Flatten class tree into a flat stream of statements, respecting `inherits`
# ---------------------------------------------------------------------------

def _flatten(
    doc: Document, *, include_path: Path | None = None
) -> list[WriteBlock | WriteOnceBlock | NewStmt]:
    """Collapse the IR into a flat list of blocks and new-statements.

    ``inherits`` relations are resolved by merging parent class statements
    into derived classes (left-to-right, later parents override earlier ones).
    """
    class_map: dict[str, ClassDef] = {}

    def collect_classes(stmts):
        for s in stmts:
            if isinstance(s, ClassDef):
                class_map[s.name] = s
                collect_classes(s.statements)

    collect_classes(doc.statements)

    # Resolve inheritance: expand each class's statements by prepending its
    # bases' statements (in order).
    resolved_bodies: dict[str, list] = {}

    def resolve(name: str, stack=()) -> list:
        if name in resolved_bodies:
            return resolved_bodies[name]
        if name in stack:
            raise ValueError(f"Circular inheritance detected: {stack + (name,)}")
        cls = class_map.get(name)
        if cls is None:
            return []
        body: list = []
        for base in cls.bases:
            body.extend(resolve(base, stack + (name,)))
        # Classes inside classes still contribute their own flatten-visible
        # statements (write/write_once/new/imports); nested class defs stay
        # scoped but their bodies are available when that class is inherited.
        for s in cls.statements:
            if isinstance(s, ClassDef):
                continue
            body.append(s)
        resolved_bodies[name] = body
        return body

    for n in list(class_map):
        resolve(n)

    # Walk top-level statements, keeping order; when we see a ClassDef we
    # attach its resolved body *immediately* so that any later reference to
    # its name (via `new`) can be expanded by the builder.
    flat: list = []
    for s in doc.statements:
        if isinstance(s, ClassDef):
            flat.extend(resolved_bodies.get(s.name, []))
        else:
            flat.append(s)

    return flat


def _moltemplate_ff_search_paths() -> list[Path]:
    """Return extra directories to search for ``import`` statements.

    Checks (in order): the installed ``moltemplate`` package's
    ``force_fields/`` directory and any colon-separated paths in the
    ``MOLTEMPLATE_PATH`` env var.
    """
    import os

    dirs: list[Path] = []
    try:
        import moltemplate  # type: ignore

        root = Path(moltemplate.__file__).parent
        ff_dir = root / "force_fields"
        if ff_dir.is_dir():
            dirs.append(ff_dir)
    except ImportError:
        pass
    for raw in os.environ.get("MOLTEMPLATE_PATH", "").split(":"):
        raw = raw.strip()
        if raw:
            dirs.append(Path(raw))
    return dirs


# Tracks files already imported (keyed by resolved absolute path) so that
# recursive cross-imports don't loop or duplicate the same content.
_ImportCache = set


def _resolve_imports(
    doc: Document,
    base_dir: Path,
    *,
    seen: set[Path] | None = None,
    search_paths: list[Path] | None = None,
) -> Document:
    """Inline ``import "path"`` statements recursively.

    Missing imports **warn** rather than raising — moltemplate itself also
    continues parsing if an optional include is absent, and we want this
    implementation to be robust on partial FF installations.
    """
    import warnings

    if seen is None:
        seen = set()
    if search_paths is None:
        search_paths = _moltemplate_ff_search_paths()

    resolved = Document()
    for s in doc.statements:
        if isinstance(s, ImportStmt):
            imp_path = _locate_import(s.path, base_dir, search_paths)
            if imp_path is None:
                warnings.warn(
                    f"moltemplate import not found: {s.path!r} (searched "
                    f"{base_dir} and {search_paths})",
                    stacklevel=2,
                )
                continue
            if imp_path in seen:
                continue
            seen.add(imp_path)
            sub_doc = parse_file(imp_path)
            sub_doc = _resolve_imports(
                sub_doc,
                imp_path.parent,
                seen=seen,
                search_paths=search_paths,
            )
            resolved.statements.extend(sub_doc.statements)
        else:
            resolved.statements.append(s)
    return resolved


def _locate_import(
    raw_path: str, base_dir: Path, search_paths: list[Path]
) -> Path | None:
    """Find an imported .lt file in base_dir or any search path."""
    candidate = (base_dir / raw_path).resolve()
    if candidate.exists():
        return candidate
    for sp in search_paths:
        candidate = (sp / raw_path).resolve()
        if candidate.exists():
            return candidate
        # basename-only fallback (some files use relative paths)
        candidate = sp / Path(raw_path).name
        if candidate.exists():
            return candidate.resolve()
    return None


# ---------------------------------------------------------------------------
# Force-field construction
# ---------------------------------------------------------------------------

_COEFF_RE = re.compile(
    r"""
    ^\s*(?P<kind>pair|bond|angle|dihedral|improper)_coeff\s+
    (?P<body>.+?)\s*$
    """,
    re.VERBOSE,
)

_SET_CHARGE_RE = re.compile(
    r"^\s*set\s+type\s+@atom:(\S+)\s+charge\s+(-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s*$"
)

_AT_IDENT_RE = re.compile(r"@atom:(\S+)")
_BOND_IDENT_RE = re.compile(r"@bond:(\S+)")
_ANGLE_IDENT_RE = re.compile(r"@angle:(\S+)")
_DIHE_IDENT_RE = re.compile(r"@dihedral:(\S+)")
_IMPR_IDENT_RE = re.compile(r"@improper:(\S+)")


def _ensure_atomtype(ff: ForceField, name: str) -> AtomType:
    style = ff.get_style_by_name("full", AtomStyle)
    if style is None:
        style = ff.def_style(AtomStyle("full"))
    # Fast path: per-style dict cache to avoid O(n) scans. Lives as a
    # private attribute on the AtomStyle instance; refreshed lazily.
    cache: dict[str, AtomType] = getattr(style, "_mt_type_cache", None)
    if cache is None:
        cache = {t.name: t for t in style.types.bucket(AtomType)}
        style._mt_type_cache = cache  # type: ignore[attr-defined]
    existing = cache.get(name)
    if existing is not None:
        return existing
    at = style.def_type(name=name, mass=0.0, charge=0.0, type_=name)
    cache[name] = at
    return at


def _split_line(line: str) -> list[str]:
    # Strip inline comments
    if "#" in line:
        line = line.split("#", 1)[0]
    return line.split()


def _parse_data_masses(ff: ForceField, lines: list[str]) -> None:
    for ln in lines:
        toks = _split_line(ln)
        if len(toks) < 2 or not toks[0].startswith("@atom:"):
            continue
        name = _strip_prefix(toks[0])
        try:
            mass = float(toks[1])
        except ValueError:
            continue
        at = _ensure_atomtype(ff, name)
        at["mass"] = mass


def _parse_data_charges(ff: ForceField, lines: list[str]) -> None:
    for ln in lines:
        m = _SET_CHARGE_RE.match(ln)
        if m:
            name, charge = m.group(1), float(m.group(2))
            at = _ensure_atomtype(ff, name)
            at["charge"] = charge


def _parse_in_settings(ff: ForceField, lines: list[str]) -> None:
    from molpy.potential import angle as _a  # registers styles
    from molpy.potential import bond as _b
    from molpy.potential import dihedral as _d
    from molpy.potential import improper as _i
    from molpy.potential import pair as _p
    # Access _kernel_registry via ForceField

    for raw in lines:
        m = _COEFF_RE.match(raw)
        if not m:
            continue
        kind = m.group("kind")
        body = m.group("body")
        toks = _split_line(body)
        if len(toks) < 2:
            continue

        if kind == "pair":
            # pair_coeff @atom:I @atom:J [style_name] p1 p2 ...
            if not (toks[0].startswith("@atom:") and toks[1].startswith("@atom:")):
                continue
            name_i = _strip_prefix(toks[0])
            name_j = _strip_prefix(toks[1])
            # Skip wildcard-only lines -- they're applied via pattern match
            # by moltemplate but we don't yet model that resolver.
            if "*" in name_i or "*" in name_j:
                continue
            rest = toks[2:]
            style_name, params = _split_style_params("pair", rest)
            style = _get_or_def_style(ff, "pair", style_name)
            if style is None:
                continue
            at1 = _ensure_atomtype(ff, name_i)
            at2 = _ensure_atomtype(ff, name_j)
            _call_def_type(style, "pair", at1, at2, params)

        elif kind == "bond":
            if not toks[0].startswith("@bond:"):
                continue
            tname = _strip_prefix(toks[0])
            rest = toks[1:]
            style_name, params = _split_style_params("bond", rest)
            style = _get_or_def_style(ff, "bond", style_name)
            if style is None:
                continue
            # Best-effort: we need two AtomTypes; synthesise wildcard types
            at1 = _ensure_atomtype(ff, f"{tname}_i")
            at2 = _ensure_atomtype(ff, f"{tname}_j")
            _call_def_type(style, "bond", at1, at2, params, type_name=tname)

        elif kind == "angle":
            if not toks[0].startswith("@angle:"):
                continue
            tname = _strip_prefix(toks[0])
            rest = toks[1:]
            style_name, params = _split_style_params("angle", rest)
            style = _get_or_def_style(ff, "angle", style_name)
            if style is None:
                continue
            ats = [
                _ensure_atomtype(ff, f"{tname}_{s}")
                for s in ("i", "j", "k")
            ]
            _call_def_type(style, "angle", *ats, params, type_name=tname)

        elif kind == "dihedral":
            if not toks[0].startswith("@dihedral:"):
                continue
            tname = _strip_prefix(toks[0])
            rest = toks[1:]
            style_name, params = _split_style_params("dihedral", rest)
            style = _get_or_def_style(ff, "dihedral", style_name)
            if style is None:
                continue
            ats = [
                _ensure_atomtype(ff, f"{tname}_{s}")
                for s in ("i", "j", "k", "l")
            ]
            _call_def_type(style, "dihedral", *ats, params, type_name=tname)

        elif kind == "improper":
            if not toks[0].startswith("@improper:"):
                continue
            tname = _strip_prefix(toks[0])
            rest = toks[1:]
            style_name, params = _split_style_params("improper", rest)
            style = _get_or_def_style(ff, "improper", style_name)
            if style is None:
                continue
            ats = [
                _ensure_atomtype(ff, f"{tname}_{s}")
                for s in ("i", "j", "k", "l")
            ]
            _call_def_type(style, "improper", *ats, params, type_name=tname)


def _split_style_params(
    kind: str, tokens: list[str]
) -> tuple[str | None, list[float]]:
    """Split a coeff tail into (style_name, numeric_params).

    If the first token is non-numeric treat it as a style hint; otherwise
    return ``(None, all_numeric_params)`` and let the caller fall back to
    the default style (``"harmonic"`` for bond/angle, ``"lj/cut/coul/cut"``
    for pair, etc.).
    """
    if not tokens:
        return None, []
    try:
        float(tokens[0])
        return None, [float(t) for t in tokens if _is_num(t)]
    except ValueError:
        return tokens[0], [float(t) for t in tokens[1:] if _is_num(t)]


def _is_num(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False


_DEFAULTS = {
    "bond": "harmonic",
    "angle": "harmonic",
    "dihedral": "opls",
    "improper": "harmonic",
    "pair": "lj/cut/coul/cut",
}


def _get_or_def_style(ff: ForceField, kind: str, style_name: str | None) -> Style | None:
    name = style_name or _DEFAULTS[kind]
    registry = ForceField._kernel_registry.get(kind, {})
    if name not in registry:
        return None
    # Find a style class: we need the matching Style subclass, not the Potential.
    # Convention: BondHarmonicStyle, BondMorseStyle, etc. — search by registering
    # Style instances lazily.
    style_cls = _STYLE_REGISTRY.get((kind, name))
    if style_cls is None:
        return None
    try:
        return ff.def_style(style_cls())
    except TypeError:
        # Style needs kwargs (e.g. PairLJ126CoulCutStyle)
        return ff.def_style(style_cls())


def _call_def_type(style: Style, kind: str, *args, type_name: str | None = None) -> Type | None:
    """Invoke ``style.def_type`` with a best-effort positional unpack.

    ``args`` ends with a ``list[float]`` of numeric params; everything before
    it is AtomTypes. If ``type_name`` is set it's passed as the ``name=`` kwarg.
    """
    params = args[-1]
    atom_types = args[:-1]
    kwargs: dict[str, Any] = {}
    if type_name is not None:
        kwargs["name"] = type_name
    try:
        return style.def_type(*atom_types, *params, **kwargs)
    except TypeError:
        # Style expects keyword params — map positional to common names.
        # Best-effort: just attach as params kwargs dict by convention.
        return None


# Static lookup table of Style subclasses keyed by (kind, name).
# Populated lazily on first call to `_get_or_def_style`.
_STYLE_REGISTRY: dict[tuple[str, str], type] = {}


def _init_style_registry() -> None:
    if _STYLE_REGISTRY:
        return
    from molpy.potential.angle import AngleHarmonicStyle, AngleClass2Style
    from molpy.potential.bond import (
        BondHarmonicStyle,
        BondMorseStyle,
        BondClass2Style,
    )
    from molpy.potential.dihedral import (
        DihedralOPLSStyle,
        DihedralPeriodicStyle,
        DihedralCharmmStyle,
        DihedralMultiHarmonicStyle,
        DihedralClass2Style,
    )
    from molpy.potential.improper import (
        ImproperPeriodicStyle,
        ImproperHarmonicStyle,
        ImproperCvffStyle,
        ImproperClass2Style,
    )
    from molpy.potential.pair import (
        PairLJ126CoulCutStyle,
        PairLJ126CoulLongStyle,
        PairBuckStyle,
        PairMorseStyle,
        PairLJClass2Style,
    )

    _STYLE_REGISTRY.update(
        {
            ("bond", "harmonic"): BondHarmonicStyle,
            ("bond", "morse"): BondMorseStyle,
            ("bond", "class2"): BondClass2Style,
            ("angle", "harmonic"): AngleHarmonicStyle,
            ("angle", "class2"): AngleClass2Style,
            ("dihedral", "opls"): DihedralOPLSStyle,
            ("dihedral", "periodic"): DihedralPeriodicStyle,
            ("dihedral", "charmm"): DihedralCharmmStyle,
            ("dihedral", "multi/harmonic"): DihedralMultiHarmonicStyle,
            ("dihedral", "class2"): DihedralClass2Style,
            ("improper", "periodic"): ImproperPeriodicStyle,
            ("improper", "harmonic"): ImproperHarmonicStyle,
            ("improper", "cvff"): ImproperCvffStyle,
            ("improper", "class2"): ImproperClass2Style,
            ("pair", "lj/cut/coul/cut"): PairLJ126CoulCutStyle,
            ("pair", "lj/cut/coul/long"): PairLJ126CoulLongStyle,
            ("pair", "buck"): PairBuckStyle,
            ("pair", "morse"): PairMorseStyle,
            ("pair", "lj/class2"): PairLJClass2Style,
        }
    )


def build_forcefield(doc: Document, *, base_dir: Path | None = None) -> ForceField:
    """Build a ``ForceField`` from a parsed MolTemplate document.

    Args:
        doc: Parsed document (from :func:`parse_file` / :func:`parse_string`).
        base_dir: Base directory for resolving ``import`` statements. Defaults
            to the current working directory.

    Returns:
        Populated ``ForceField`` instance.
    """
    _init_style_registry()
    if base_dir is None:
        base_dir = Path.cwd()
    doc = _resolve_imports(doc, base_dir)
    flat = _flatten(doc)

    ff = ForceField(name="moltemplate", units="real")
    ff.metadata = {}  # type: ignore[attr-defined]

    for stmt in flat:
        if isinstance(stmt, (WriteBlock, WriteOnceBlock)):
            section = stmt.section
            if section in ("Data Masses",):
                _parse_data_masses(ff, stmt.body_lines)
            elif section in ("Data Charges", "In Charges"):
                _parse_data_charges(ff, stmt.body_lines)
            elif section == "In Settings":
                _parse_in_settings(ff, stmt.body_lines)
            else:
                # Preserve other sections verbatim for emitters to consume.
                getattr(ff, "metadata", {}).setdefault(section, []).extend(
                    stmt.body_lines
                )
    return ff


# ---------------------------------------------------------------------------
# System (Atomistic) construction
# ---------------------------------------------------------------------------

_DATA_ATOM_RE = re.compile(
    r"""
    ^\s*\$atom(?::(?P<id>\S+))?\s+
    (?:\$mol(?::\S+)?\s+)?
    @atom(?::(?P<type>\S+))?\s+
    (?P<charge>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+
    (?P<x>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+
    (?P<y>-?\d+\.?\d*(?:[eE][+-]?\d+)?)\s+
    (?P<z>-?\d+\.?\d*(?:[eE][+-]?\d+)?)
    """,
    re.VERBOSE,
)

# Bond with type: $bond:id @bond:T $atom:a $atom:b
_DATA_BOND_TYPED_RE = re.compile(
    r"^\s*\$bond:(?P<id>\S+)\s+@bond:(?P<type>\S+)\s+\$atom:(?P<a>\S+)\s+\$atom:(?P<b>\S+)"
)
# Bond list (no type): $bond:id $atom:a $atom:b
_DATA_BOND_LIST_RE = re.compile(
    r"^\s*\$bond:(?P<id>\S+)\s+\$atom:(?P<a>\S+)\s+\$atom:(?P<b>\S+)"
)

_DATA_ANGLE_RE = re.compile(
    r"^\s*\$angle:(?P<id>\S+)\s+@angle:(?P<type>\S+)\s+\$atom:(?P<a>\S+)\s+\$atom:(?P<b>\S+)\s+\$atom:(?P<c>\S+)"
)

_DATA_DIHE_RE = re.compile(
    r"^\s*\$dihedral:(?P<id>\S+)\s+@dihedral:(?P<type>\S+)\s+"
    r"\$atom:(?P<a>\S+)\s+\$atom:(?P<b>\S+)\s+\$atom:(?P<c>\S+)\s+\$atom:(?P<d>\S+)"
)
_DATA_IMPR_RE = re.compile(
    r"^\s*\$improper:(?P<id>\S+)\s+@improper:(?P<type>\S+)\s+"
    r"\$atom:(?P<a>\S+)\s+\$atom:(?P<b>\S+)\s+\$atom:(?P<c>\S+)\s+\$atom:(?P<d>\S+)"
)


def _build_template(
    cls_stmts: list,
    *,
    class_map: dict | None = None,
    class_body_fn=None,
) -> Atomistic:
    """Turn a class body into an ``Atomistic`` template.

    Handles nested ``new`` statements inside classes (e.g. butane's
    ``monomer1 = new CH3``) and scoped ``$atom:submol/atom`` references
    by first recursively materialising sub-instances, then resolving
    paths against the local symbol table ``atoms``.
    """
    atoms: dict[str, Atom] = {}
    tmpl = Atomistic()

    # Pass 1: expand inner `new` statements into a merged sub-Atomistic per
    # instance name. Their atoms are keyed "<inst_name>/<atom_id>" in the
    # local `atoms` map so bonds in "Data Bond List" and similar blocks can
    # refer to them via `$atom:inst/aid` or plain `$atom:inst_aid` paths.
    if class_map is not None and class_body_fn is not None:
        for s in cls_stmts:
            if isinstance(s, NewStmt):
                sub_stmts = class_body_fn(s.class_name)
                if not sub_stmts:
                    continue
                sub_tmpl = _build_template(
                    sub_stmts, class_map=class_map, class_body_fn=class_body_fn
                )
                for t in s.transforms:
                    _apply_transform(sub_tmpl, t)
                # Rename sub atoms as "inst/atom_id" and merge
                for a in list(sub_tmpl.atoms):
                    base_name = a.get("name", "")
                    key = f"{s.instance_name}/{base_name}"
                    atoms[key] = a
                    a["name"] = key
                tmpl += sub_tmpl

    # Pass 2: own Data Atoms / Bonds / Angles / Dihedrals / Impropers
    for s in cls_stmts:
        if not isinstance(s, (WriteBlock, WriteOnceBlock)):
            continue
        if s.section == "Data Atoms":
            for ln in s.body_lines:
                m = _DATA_ATOM_RE.match(ln)
                if m:
                    aid = m.group("id") or f"_anon{len(atoms)}"
                    atom = tmpl.def_atom(
                        name=aid,
                        type=m.group("type") or "",
                        charge=float(m.group("charge")),
                        xyz=[
                            float(m.group("x")),
                            float(m.group("y")),
                            float(m.group("z")),
                        ],
                    )
                    atoms[aid] = atom
        elif s.section in ("Data Bonds", "Data Bond List"):
            for ln in s.body_lines:
                # Try typed form first
                m = _DATA_BOND_TYPED_RE.match(ln)
                if m:
                    a = _resolve_atom(atoms, m.group("a"))
                    b = _resolve_atom(atoms, m.group("b"))
                    if a is not None and b is not None:
                        tmpl.def_bond(a, b, type=m.group("type"))
                        continue
                m = _DATA_BOND_LIST_RE.match(ln)
                if m:
                    a = _resolve_atom(atoms, m.group("a"))
                    b = _resolve_atom(atoms, m.group("b"))
                    if a is not None and b is not None:
                        tmpl.def_bond(a, b, type="")
        elif s.section == "Data Angles":
            for ln in s.body_lines:
                m = _DATA_ANGLE_RE.match(ln)
                if m:
                    ra = _resolve_atom(atoms, m.group("a"))
                    rb = _resolve_atom(atoms, m.group("b"))
                    rc = _resolve_atom(atoms, m.group("c"))
                    if ra and rb and rc:
                        tmpl.def_angle(ra, rb, rc, type=m.group("type"))
        elif s.section == "Data Dihedrals":
            for ln in s.body_lines:
                m = _DATA_DIHE_RE.match(ln)
                if m:
                    endpoints = [
                        _resolve_atom(atoms, m.group(g))
                        for g in ("a", "b", "c", "d")
                    ]
                    if all(endpoints):
                        tmpl.def_dihedral(*endpoints, type=m.group("type"))
        elif s.section == "Data Impropers":
            for ln in s.body_lines:
                m = _DATA_IMPR_RE.match(ln)
                if m:
                    endpoints = [
                        _resolve_atom(atoms, m.group(g))
                        for g in ("a", "b", "c", "d")
                    ]
                    if all(endpoints):
                        # Store improper as a plain Dihedral with a marker
                        # in its type attribute (MolPy has no Improper link
                        # class yet; keeping data so emitters can round-trip).
                        tmpl.def_dihedral(
                            *endpoints, type=f"improper:{m.group('type')}"
                        )
    return tmpl


def _resolve_atom(atoms: dict[str, Atom], ref: str) -> Atom | None:
    """Resolve a scoped ``$atom:...`` reference.

    Tries exact match, then ``a/b`` → ``a_b`` / nested lookups to support
    moltemplate's path syntax (``$atom:monomer1/c``). Falls back to a
    longest-suffix scan if none of the above matches.
    """
    if ref in atoms:
        return atoms[ref]
    # Swap / for _ and vice-versa
    if "/" in ref:
        alt = ref.replace("/", "_")
        if alt in atoms:
            return atoms[alt]
    if "_" in ref:
        alt = ref.replace("_", "/")
        if alt in atoms:
            return atoms[alt]
    # Try matching by tail of the path
    for k, v in atoms.items():
        if k.endswith("/" + ref) or k.endswith("_" + ref) or k == ref:
            return v
    return None


def _apply_transform(mol: Atomistic, t: Transform) -> None:
    import math

    op = t.op
    a = t.args
    if op == "move" and len(a) == 3:
        mol.move([a[0], a[1], a[2]])
    elif op == "rot" and len(a) >= 4:
        angle_deg, ax, ay, az = a[0], a[1], a[2], a[3]
        about = list(a[4:7]) if len(a) >= 7 else None
        mol.rotate(
            [ax, ay, az], math.radians(angle_deg), about=about
        )
    elif op == "scale" and len(a) == 1:
        mol.scale(a[0])


def build_system(
    doc: Document,
    ff: ForceField | None = None,
    *,
    base_dir: Path | None = None,
    auto_topology: bool = True,
) -> tuple[Atomistic, ForceField]:
    """Build an ``Atomistic`` system + ``ForceField`` from a MolTemplate doc.

    ``new`` statements are instantiated by deep-copying the declaring class's
    template atoms/bonds/angles/dihedrals, applying the transform chain, and
    merging into the resulting container Atomistic.

    Args:
        doc: Parsed MolTemplate document.
        ff: Optional pre-built ForceField to populate further. If None a
            new one is built from ``doc``.
        base_dir: Base directory used to resolve ``import`` statements.
        auto_topology: When True (default), angles and dihedrals missing
            from the source are auto-generated from bond connectivity
            after all ``new`` statements are expanded. This approximates
            moltemplate's ``Angles By Type`` / ``Dihedrals By Type`` rules
            for the common case where rules are purely connectivity based.
            Set False to retain only explicitly declared angles/dihedrals.
    """
    _init_style_registry()
    if base_dir is None:
        base_dir = Path.cwd()
    doc = _resolve_imports(doc, base_dir)

    if ff is None:
        ff = build_forcefield(doc, base_dir=base_dir)

    # Collect class templates (resolving inheritance by merging statements)
    class_map: dict[str, ClassDef] = {}

    def collect(stmts):
        for s in stmts:
            if isinstance(s, ClassDef):
                class_map[s.name] = s
                collect(s.statements)

    collect(doc.statements)

    def class_body(name: str, stack=()) -> list:
        if name in stack:
            return []
        cls = class_map.get(name)
        if cls is None:
            return []
        body: list = []
        for base in cls.bases:
            body.extend(class_body(base, stack + (name,)))
        body.extend(s for s in cls.statements if not isinstance(s, ClassDef))
        return body

    system = Atomistic()
    for stmt in doc.statements:
        if isinstance(stmt, NewStmt):
            if stmt.class_name == "random":
                continue  # handled separately (no deterministic expansion)
            cls_stmts = class_body(stmt.class_name)
            base_count = max(stmt.count, 1)
            grid = _enumerate_array_positions(stmt.arrays)
            for _ in range(base_count):
                for offsets in grid:
                    tmpl = _build_template(
                        cls_stmts,
                        class_map=class_map,
                        class_body_fn=class_body,
                    )
                    for t in stmt.transforms:
                        _apply_transform(tmpl, t)
                    for offset in offsets:
                        _apply_transform(tmpl, offset)
                    system += tmpl

    if auto_topology and len(list(system.bonds)) > 0:
        # Use MolPy's existing topology generator to fill in missing
        # angles/dihedrals based on bond connectivity. This approximates
        # moltemplate's `Angles By Type` / `Dihedrals By Type` rules for
        # connectivity-driven FFs. Type assignment remains best-effort.
        existing_angles = len(list(system.angles))
        existing_dihedrals = len(list(system.dihedrals))
        try:
            system = system.get_topo(
                gen_angle=existing_angles == 0,
                gen_dihe=existing_dihedrals == 0,
            )  # type: ignore[assignment]
        except Exception:
            pass

    return system, ff


def _enumerate_array_positions(
    arrays: list[ArrayDim],
) -> list[list[Transform]]:
    """Enumerate every grid cell as a list-of-transforms.

    For ``[N1].move(dx1) [N2].move(dx2)`` returns N1 * N2 lists; cell
    ``(k1, k2)`` is ``[Transform("move", k1*dx1), Transform("move", k2*dx2)]``.
    If ``arrays`` is empty, returns a single empty list (one copy, no shift).
    """
    if not arrays:
        return [[]]
    cells: list[list[Transform]] = [[]]
    for dim in arrays:
        new_cells: list[list[Transform]] = []
        for cell in cells:
            for k in range(dim.count):
                if dim.transform is None or k == 0:
                    new_cells.append(cell + [])
                else:
                    scaled = Transform(
                        op=dim.transform.op,
                        args=[a * k for a in dim.transform.args],
                    )
                    new_cells.append(cell + [scaled])
        cells = new_cells
    return cells
