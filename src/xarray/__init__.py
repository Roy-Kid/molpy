import numpy as np
from collections.abc import Mapping, Sequence

class DataArray:
    def __init__(self, data, dims=("index",)):
        self.values = np.asarray(data)
        if isinstance(dims, str):
            self.dims = (dims,)
        else:
            self.dims = tuple(dims)

    @property
    def shape(self):
        return self.values.shape

    def __len__(self):
        return self.values.shape[0]

class Dataset:
    def __init__(self, data=None):
        self.data_vars = {}
        if data:
            for name, value in data.items():
                dims, arr = value
                self.data_vars[name] = DataArray(arr, dims)
        self._update_dims()

    def _update_dims(self):
        dims = {}
        for arr in self.data_vars.values():
            for d, s in zip(arr.dims, arr.values.shape):
                dims[d] = s
        self.dims = dims

    @property
    def sizes(self):
        return self.dims

    @property
    def shape(self):
        return (self.dims.get("index", 0), len(self.data_vars))

    @property
    def array_length(self):
        return self.dims.get("index", 0)

    def __getitem__(self, key):
        return self.data_vars[key]

    def get(self, key, default=None):
        return self.data_vars.get(key, default)

    def __setitem__(self, key, value):
        if isinstance(value, DataArray):
            self.data_vars[key] = value
        elif isinstance(value, tuple) and len(value) == 2:
            dims, arr = value
            self.data_vars[key] = DataArray(arr, dims)
        elif isinstance(value, Mapping):
            self.data_vars[key] = DataArray(value)
        else:
            raise TypeError("Invalid value for Dataset")
        self._update_dims()

    def __delitem__(self, key):
        del self.data_vars[key]
        self._update_dims()

    def __iter__(self):
        return iter(self.data_vars)

    def __len__(self):
        return self.dims.get("index", 0)

    def copy(self):
        data = {k: (v.dims, v.values.copy()) for k, v in self.data_vars.items()}
        return Dataset(data)

    def isel(self, *, index=None):
        if index is None:
            raise ValueError("index selection required")
        index = np.asarray(index)
        data = {k: (v.dims, v.values[index]) for k, v in self.data_vars.items()}
        return Dataset(data)

class DataTreeNode:
    def __init__(self, ds=None):
        self.ds = ds if ds is not None else Dataset()
        self.children = {}

    def __getitem__(self, key):
        return self.children[key]

    def __setitem__(self, key, value):
        if isinstance(value, DataTreeNode):
            self.children[key] = value
        elif isinstance(value, Dataset):
            self.children[key] = DataTreeNode(value)
        else:
            raise TypeError("DataTreeNode value must be Dataset or DataTreeNode")

    def keys(self):
        return self.children.keys()

class DataTree:
    def __init__(self, mapping=None):
        self.root = DataTreeNode()
        mapping = mapping or {}
        for k, v in mapping.items():
            self.root[k] = v

    def __getitem__(self, key):
        return self.root[key]

    def __setitem__(self, key, value):
        self.root[key] = value

    def __delitem__(self, key):
        del self.root.children[key]

    def keys(self):
        return self.root.keys()

__all__ = ["Dataset", "DataArray", "DataTree", "concat"]

def concat(datasets: Sequence[Dataset], dim="index") -> Dataset:
    arrays = {}
    for ds in datasets:
        for k, arr in ds.data_vars.items():
            arrays.setdefault(k, []).append(arr.values)
    new_data = {}
    for k, arrs in arrays.items():
        new_data[k] = (dim, np.concatenate(arrs, axis=0))
    return Dataset(new_data)

