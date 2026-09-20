"""Top-level CLI entry: builds argparse tree and dispatches subcommands.

Registered subcommands:
  * ``molpy moltemplate`` -- Native execution of moltemplate .lt scripts.
"""

from __future__ import annotations

import argparse

from . import moltemplate


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="molpy", description="MolPy CLI")
    parser.add_argument("--version", action="store_true", help="Show version and exit.")
    sub = parser.add_subparsers(dest="command")
    moltemplate.register(sub)
    return parser


def main(argv: list[str] | None = None) -> int:
    parser = _build_parser()
    args = parser.parse_args(argv)

    if getattr(args, "version", False):
        from molpy.version import version

        print(f"molpy {version}")
        return 0

    cmd = args.command
    handler = getattr(args, "func", None)
    if cmd is None or handler is None:
        parser.print_help()
        return 0
    rc = handler(args)
    return int(rc) if rc is not None else 0
