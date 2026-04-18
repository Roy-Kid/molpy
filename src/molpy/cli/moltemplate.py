"""``molpy moltemplate`` — execute moltemplate .lt scripts natively.

Sub-sub-commands:
  * ``run``      -- parse .lt and emit a full engine input set.
  * ``parse``    -- parse .lt and dump the IR (debug).
  * ``info``     -- print counts (atom types, molecules, ...).
  * ``convert``  -- FF-only conversion to MolPy XML.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict, is_dataclass
from pathlib import Path


def register(sub: argparse._SubParsersAction) -> None:
    p = sub.add_parser(
        "moltemplate",
        help="Execute moltemplate .lt scripts natively.",
    )
    mtsub = p.add_subparsers(dest="mt_cmd", required=True)

    # run
    r = mtsub.add_parser(
        "run",
        help="Parse .lt and emit complete input set for one or more engines.",
    )
    r.add_argument("script", type=Path, help="Path to the .lt script.")
    r.add_argument(
        "--emit",
        "-e",
        action="append",
        default=None,
        help=(
            "Target engine (repeatable). Choices: lammps, openmm, gromacs, "
            "xml, all. Default: lammps."
        ),
    )
    r.add_argument(
        "--out-dir",
        "-o",
        type=Path,
        default=Path("."),
        help="Output directory (default: cwd).",
    )
    r.add_argument(
        "--prefix",
        default="system",
        help="Filename prefix (default: 'system').",
    )
    r.set_defaults(func=_cmd_run)

    # parse
    pr = mtsub.add_parser("parse", help="Parse a .lt file and dump the IR.")
    pr.add_argument("script", type=Path)
    pr.add_argument(
        "--json", type=Path, default=None, help="Write IR as JSON to this path."
    )
    pr.set_defaults(func=_cmd_parse)

    # info
    info = mtsub.add_parser("info", help="Summarise a .lt file.")
    info.add_argument("script", type=Path)
    info.set_defaults(func=_cmd_info)

    # convert
    c = mtsub.add_parser(
        "convert", help="Convert a moltemplate FF to MolPy XML."
    )
    c.add_argument("src", type=Path, help="Input .lt file.")
    c.add_argument("dst", type=Path, help="Output .xml path.")
    c.set_defaults(func=_cmd_convert)


# ---------------------------------------------------------------------------
# Handlers
# ---------------------------------------------------------------------------

def _cmd_run(args: argparse.Namespace) -> int:
    from molpy.io.emit import EMITTERS, emit, emit_all
    from molpy.io.forcefield.moltemplate import read_moltemplate_system

    if not args.script.exists():
        print(f"molpy: error: {args.script} not found", file=sys.stderr)
        return 1

    atomistic, ff = read_moltemplate_system(args.script)
    out_dir = args.out_dir
    out_dir.mkdir(parents=True, exist_ok=True)

    engines = args.emit or ["lammps"]
    if "all" in engines:
        results = emit_all(atomistic, ff, out_dir, prefix=args.prefix)
        for engine, paths in results.items():
            _print_emit_result(engine, paths)
        return 0

    unknown = [e for e in engines if e not in EMITTERS]
    if unknown:
        print(
            f"molpy: error: unknown engine(s): {unknown}. "
            f"Registered: {sorted(EMITTERS)}",
            file=sys.stderr,
        )
        return 2

    for engine in engines:
        paths = emit(engine, atomistic, ff, out_dir, prefix=args.prefix)
        _print_emit_result(engine, paths)
    return 0


def _cmd_parse(args: argparse.Namespace) -> int:
    from molpy.parser.moltemplate import parse_file

    doc = parse_file(args.script)
    if args.json:
        args.json.write_text(json.dumps(_ir_to_jsonable(doc), indent=2))
        print(f"IR written to {args.json}")
        return 0
    # Default: print a summary
    kinds: dict[str, int] = {}
    for s in doc.statements:
        kinds[type(s).__name__] = kinds.get(type(s).__name__, 0) + 1
    print(f"{args.script}:")
    for k, n in sorted(kinds.items()):
        print(f"  {k}: {n}")
    return 0


def _cmd_info(args: argparse.Namespace) -> int:
    from molpy.core.forcefield import AtomStyle, AtomType
    from molpy.io.forcefield.moltemplate import read_moltemplate_system

    atomistic, ff = read_moltemplate_system(args.script)
    astyle = ff.get_style_by_name("full", AtomStyle)
    n_atomtypes = len(astyle.types.bucket(AtomType)) if astyle else 0
    print(f"{args.script}:")
    print(f"  atom types: {n_atomtypes}")
    print(f"  atoms:      {len(list(atomistic.atoms))}")
    print(f"  bonds:      {len(list(atomistic.bonds))}")
    print(f"  angles:     {len(list(atomistic.angles))}")
    print(f"  dihedrals:  {len(list(atomistic.dihedrals))}")
    print(f"  ff styles:  {len(list(ff.styles.bucket(object)))}")
    return 0


def _cmd_convert(args: argparse.Namespace) -> int:
    from molpy.io.forcefield.moltemplate import read_moltemplate
    from molpy.io.forcefield.xml import XMLForceFieldWriter

    ff = read_moltemplate(args.src)
    args.dst.parent.mkdir(parents=True, exist_ok=True)
    XMLForceFieldWriter(args.dst).write(ff)
    print(f"{args.src} -> {args.dst}")
    return 0


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _print_emit_result(engine: str, paths: list[Path]) -> None:
    print(f"[{engine}]")
    for p in paths:
        print(f"  {p}")


def _ir_to_jsonable(obj):
    if is_dataclass(obj):
        d = {"_kind": type(obj).__name__}
        d.update(asdict(obj))
        # Recurse into nested dataclasses in lists
        for k, v in list(d.items()):
            d[k] = _ir_to_jsonable(v)
        return d
    if isinstance(obj, list):
        return [_ir_to_jsonable(x) for x in obj]
    if isinstance(obj, dict):
        return {k: _ir_to_jsonable(v) for k, v in obj.items()}
    return obj
