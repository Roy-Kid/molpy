from __future__ import annotations

from datatree import DataTree
import xarray as xr
import numpy as np
from typing import Mapping, Any


class Frame:
    """Experimental frame implementation based on ``xarray`` and ``DataTree``."""

    def __init__(self, data: Mapping[str, Any] | None = None) -> None:
        data = data or {}
        datasets = {}
        for key, value in data.items():
            datasets[key] = self._to_dataset(value)
        self.tree = DataTree(datasets)

    def _to_dataset(self, value: Any) -> xr.Dataset:
        if isinstance(value, xr.Dataset):
            return value
        if isinstance(value, Mapping):
            arrays = {k: ("index", np.asarray(v)) for k, v in value.items()}
            return xr.Dataset(arrays)
        raise TypeError("Frame values must be mapping or xarray.Dataset")

    def __getitem__(self, key: str) -> xr.Dataset:
        return self.tree[key].ds

    def __setitem__(self, key: str, value: Any) -> None:
        self.tree[key] = self._to_dataset(value)

    def keys(self):
        return self.tree.keys()

    def values(self):
        return (self.tree[k].ds for k in self.tree.keys())

    def items(self):
        for k in self.tree.keys():
            yield k, self.tree[k].ds

    def __iter__(self):
        return iter(self.tree.keys())

    def __contains__(self, key: str) -> bool:
        return key in self.tree

    def __len__(self) -> int:
        return len(self.tree.keys())

    def copy(self) -> "Frame":
        copied = {k: v.copy() for k, v in self.items()}
        return Frame(copied)
