"""XYZ file I/O — the native core backend with thin molpy column normalization.

Parse/serialize: the native ``io`` module. After read, molpy may merge split multi-
columns (``CS_1``+``CS_2``→``CS``), map ``species``→``element``, and fill
``atomic_number`` when missing.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

import molrs.io
from molrs import Element, Frame

from molpy.core.fields import ATOMIC_NUMBER

from .base import DataReader, DataWriter


def _normalize_xyz_frame(frame: Frame) -> Frame:
    """Apply molpy column conventions on a native XYZ Frame (in place)."""
    for block_name in list(frame.keys()):
        block = frame[block_name]
        keys = set(block.keys())
        # molrs splits an n-wide property into ``base_1`` .. ``base_n``;
        # rejoin every such run, whatever its width.
        merged: list[tuple[str, list[str]]] = []
        for key in sorted(keys):
            if not key.endswith("_1"):
                continue
            base = key[:-2]
            parts = [key]
            while f"{base}_{len(parts) + 1}" in keys:
                parts.append(f"{base}_{len(parts) + 1}")
            if len(parts) > 1:
                merged.append((base, parts))
        for base, parts in merged:
            block[base] = np.column_stack([np.asarray(block[k]) for k in parts])
            for k in parts:
                del block[k]
        if "species" in block and "element" not in block:
            block["element"] = np.asarray(block["species"])
        if "element" in block and ATOMIC_NUMBER not in block:
            z_list = [Element.get_atomic_number(str(s)) for s in block["element"]]
            block[ATOMIC_NUMBER] = np.array(z_list, dtype=np.int64)
    return frame


class XYZReader(DataReader):
    """Read XYZ natively + :func:`_normalize_xyz_frame`."""

    def __init__(self, path: str | Path, **kwargs: object) -> None:
        super().__init__(Path(path), **kwargs)

    def read(self, frame: Frame | None = None) -> Frame:
        del frame
        return _normalize_xyz_frame(molrs.io.read_xyz(str(self._path)))


class XYZWriter(DataWriter):
    """Write XYZ natively."""

    def __init__(self, path: str | Path, **kwargs: object) -> None:
        super().__init__(Path(path), **kwargs)

    def write(self, frame: Frame) -> None:
        molrs.io.write_xyz(str(self._path), frame)
