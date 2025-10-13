import csv
from collections.abc import Iterator, MutableMapping
from io import StringIO
from pathlib import Path
from typing import Any, Mapping, Self, TypeAlias, overload

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .selector import Selector
from .topology import Topology

BlockLike: TypeAlias = Mapping[str, ArrayLike]


class Block(MutableMapping[str, np.ndarray]):
    """
    Lightweight container that maps variable names → NumPy arrays.

    • Behaves like a dict but auto-casts any assigned value to ndarray.
    • All built-in `dict`/`MutableMapping` helpers work out of the box.

    Examples
    --------
    >>> blk = Block()
    >>> blk["xyz"] = [[0, 0, 0], [1, 1, 1]]
    >>> "xyz" in blk
    True
    >>> len(blk)
    1
    >>> blk["xyz"].dtype
    dtype('float64')
    """

    __slots__ = ("_vars",)

    def __init__(self, vars_: BlockLike | None = None) -> None:
        self._vars: dict[str, np.ndarray] = {k: np.asarray(v) for k, v in {}.items()}
        if vars_ is not None:
            if not isinstance(vars_, dict):
                raise ValueError(f"vars_ must be a dict, got {type(vars_)}")
            for k, v in vars_.items():
                try:
                    self._vars[k] = np.asarray(v)
                except Exception as e:
                    raise ValueError(
                        f"Value must be a BlockLike, i.e. dict[str, np.ndarray], but got {type(v)} for key {k}"
                    ) from e

    # --- core mapping API

    @overload
    def __getitem__(self, key: str) -> np.ndarray: ...

    @overload
    def __getitem__(self, key: int | slice) -> dict[str, np.ndarray]: ...  # type: ignore[override]

    @overload
    def __getitem__(self, key: list[str]) -> np.ndarray: ...  # type: ignore[override]

    @overload
    def __getitem__(self, key: np.ndarray) -> "Block": ...  # type: ignore[override]

    @overload
    def __getitem__(self, key: Selector) -> "Block": ...  # type: ignore[override]

    def __getitem__(self, key):  # type: ignore[override]
        if isinstance(key, (int, slice)):
            return {
                k: (v[key] if v[key].ndim > 0 else v[key].item())
                for k, v in self._vars.items()
            }
        elif isinstance(key, str):
            return self._vars[key]
        elif isinstance(key, list):
            # Handle list of column names for concatenation
            if not key:
                raise KeyError("Empty list not allowed for indexing")

            # Check if all keys exist
            for k in key:
                if k not in self._vars:
                    raise KeyError(f"Key '{k}' not found in Block")

            # Get the arrays
            arrays = [self._vars[k] for k in key]

            # Check if all arrays have the same shape and dtype
            if not arrays:
                raise ValueError("No arrays to concatenate")

            first_array = arrays[0]
            for i, arr in enumerate(arrays[1:], 1):
                if arr.shape != first_array.shape:
                    raise ValueError(
                        f"Arrays must have the same shape. Array {key[0]} has shape {first_array.shape}, but array {key[i]} has shape {arr.shape}"
                    )
                if arr.dtype != first_array.dtype:
                    raise ValueError(
                        f"Arrays must have the same dtype. Array {key[0]} has dtype {first_array.dtype}, but array {key[i]} has dtype {arr.dtype}"
                    )

            # Concatenate along the last axis
            return np.column_stack(arrays)
        elif isinstance(key, tuple):
            return np.array([self[k] for k in key])
        elif isinstance(key, np.ndarray):
            return Block({k: v[key] for k, v in self._vars.items()})
        elif isinstance(key, Selector):
            return key(self)
        else:
            raise KeyError(
                f"Invalid key type: {type(key)}. Expected str, int, slice, list[str], or np.ndarray."
            )

    def __setitem__(self, key: str, value: Any) -> None:  # type: ignore[override]
        self._vars[key] = np.asarray(value)

    def __delitem__(self, key: str) -> None:  # type: ignore[override]
        del self._vars[key]

    def __iter__(self) -> Iterator[str]:  # type: ignore[override]
        return iter(self._vars)

    def __len__(self) -> int:  # type: ignore[override]
        return len(self._vars)

    def __contains__(self, key: str) -> bool:  # type: ignore[override]
        """Check if a variable exists in this block."""
        return key in self._vars

    # ------------------------------------------------------------------ helpers
    def to_dict(self) -> dict[str, np.ndarray]:
        """Return a JSON-serialisable copy (arrays → Python lists)."""
        return {k: v for k, v in self._vars.items()}

    @classmethod
    def from_dict(cls, data: dict[str, np.ndarray]) -> "Block":
        """Inverse of :meth:`to_dict`."""
        return cls({k: np.asarray(v) for k, v in data.items()})

    @classmethod
    def from_csv(
        cls,
        filepath: str | Path | StringIO,
        *,
        delimiter: str = ",",
        encoding: str = "utf-8",
        header: list[str] | None = None,
        **kwargs,
    ) -> "Block":
        """
        Create a Block from a CSV file or StringIO.

        Parameters
        ----------
        filepath : str, Path, or StringIO
            Path to the CSV file or StringIO object
        delimiter : str, default=","
            CSV delimiter character
        encoding : str, default="utf-8"
            File encoding (ignored for StringIO)
        header : list[str] or None, default=None
            Column names. If None, first row is used as headers.
            If provided, CSV is assumed to have no header row.
        **kwargs
            Additional arguments passed to csv.reader

        Returns
        -------
        Block
            A new Block instance with data from the CSV file

        Examples
        --------
        >>> block = Block.from_csv("data.csv")
        >>> block = Block.from_csv("data.csv", delimiter=";")
        >>> from io import StringIO
        >>> csv_data = StringIO("x,y,z\\n0,1,2\\n3,4,5")
        >>> block = Block.from_csv(csv_data)
        >>> # No header CSV
        >>> csv_no_header = StringIO("0,1,2\\n3,4,5")
        >>> block = Block.from_csv(csv_no_header, header=["x", "y", "z"])
        """
        # 判断类型
        if isinstance(filepath, StringIO):
            csvfile = filepath
            csvfile.seek(0)
            close_file = False
        else:
            filepath = Path(filepath)
            if not filepath.exists():
                raise FileNotFoundError(f"CSV file not found: {filepath}")
            csvfile = open(filepath, "r", encoding=encoding, newline="")
            close_file = True

        try:
            reader = csv.reader(csvfile, delimiter=delimiter, **kwargs)

            # Handle headers
            if header is None:
                # Use first row as headers
                try:
                    headers = next(reader)
                except StopIteration:
                    raise ValueError("CSV file is empty")
            else:
                # Use provided headers, no header row in CSV
                headers = header

            raw_data = {h: [] for h in headers}

            for row in reader:
                for i, header_name in enumerate(headers):
                    raw_data[header_name].append(row[i])

            data = {}
            for k, v in raw_data.items():
                for dtype in (int, float, str):
                    try:
                        data[k] = np.array(v, dtype=dtype)
                        break
                    except ValueError:
                        continue
                else:
                    raise ValueError(f"Failed to convert {k} to any of int, float, str")

            return cls(data)
        finally:
            if close_file:
                csvfile.close()

    def copy(self) -> "Block":
        """Shallow copy (arrays are **not** copied)."""
        return Block(self._vars.copy())  # type: ignore[arg-type]

    def _sort(self, key: str, *, reverse: bool = False) -> dict[str, NDArray[Any]]:
        """Sort variables by a specific key and return sorted data.

        This is a private helper method that performs the actual sorting logic.

        Args:
            key: The variable name to sort by.
            reverse: If True, sort in descending order. Defaults to False.

        Returns:
            Dictionary with sorted variable data.

        Raises:
            KeyError: If the key variable doesn't exist in the block.
            ValueError: If any variable has different length than the key variable.
        """
        if not self._vars:
            return {}

        if key not in self._vars:
            raise KeyError(f"Variable '{key}' not found in block")

        # Get the sorting indices
        sort_indices = np.argsort(self._vars[key])
        if reverse:
            sort_indices = sort_indices[::-1]

        # Create sorted data
        sorted_vars: dict[str, np.ndarray] = {}
        for var_name, var_data in self._vars.items():
            if len(var_data) != len(self._vars[key]):
                raise ValueError(
                    f"Variable '{var_name}' has different length than '{key}'"
                )
            sorted_vars[var_name] = var_data[sort_indices]

        return sorted_vars

    def sort(self, key: str, *, reverse: bool = False) -> "Block":
        """Sort the block by a specific variable and return a new sorted Block.

        This method creates a new Block instance with sorted data, leaving the
        original Block unchanged.

        Args:
            key: The variable name to sort by.
            reverse: If True, sort in descending order. Defaults to False.

        Returns:
            A new Block with sorted data.

        Raises:
            KeyError: If the key variable doesn't exist in the block.
            ValueError: If any variable has different length than the key variable.

        Example:
            >>> blk = Block({"x": [3, 1, 2], "y": [30, 10, 20]})
            >>> sorted_blk = blk.sort("x")
            >>> sorted_blk["x"]
            array([1, 2, 3])
            >>> sorted_blk["y"]
            array([10, 20, 30])
            >>> # Original block is unchanged
            >>> blk["x"]
            array([3, 1, 2])
        """
        sorted_vars = self._sort(key, reverse=reverse)
        return Block(sorted_vars)

    def sort_(self, key: str, *, reverse: bool = False) -> "Self":
        """Sort the block in-place by a specific variable.

        This method modifies the current Block instance by sorting all variables
        according to the specified key. The original data is overwritten.

        Args:
            key: The variable name to sort by.
            reverse: If True, sort in descending order. Defaults to False.

        Returns:
            Self (for method chaining).

        Raises:
            KeyError: If the key variable doesn't exist in the block.
            ValueError: If any variable has different length than the key variable.

        Example:
            >>> blk = Block({"x": [3, 1, 2], "y": [30, 10, 20]})
            >>> blk.sort_("x")
            >>> blk["x"]
            array([1, 2, 3])
            >>> blk["y"]
            array([10, 20, 30])
            >>> # Original data is now sorted
        """
        sorted_vars = self._sort(key, reverse=reverse)
        if sorted_vars:  # Only update if we have data to sort
            self._vars.update(sorted_vars)
        return self

    # ------------------------------------------------------------------ repr / str
    def __repr__(self) -> str:
        contents = ", ".join(f"{k}: shape={v.shape}" for k, v in self._vars.items())
        return f"Block({contents})"

    @property
    def nrows(self) -> int:
        """Return the number of rows in the first variable (if any)."""
        if not self._vars:
            return 0
        return len(next(iter(self._vars.values())))

    @property
    def shape(self) -> tuple[int, ...]:
        """Return the shape of the first variable (if any)."""
        if not self._vars:
            return ()
        return self.nrows, len(self)

    def iterrows(self, n: int | None = None) -> Iterator[tuple[int, dict[str, Any]]]:
        """
        Iterate over rows of the block.

        Returns
        -------
        Iterator[tuple[int, dict[str, Any]]]
            An iterator yielding (index, row_data) pairs where:
            - index: int, the row index
            - row_data: dict, mapping variable names to their values for this row

        Examples
        --------
        >>> blk = Block({
        ...     "id": [1, 2, 3],
        ...     "type": ["C", "O", "N"],
        ...     "x": [0.0, 1.0, 2.0],
        ...     "y": [0.0, 0.0, 1.0],
        ...     "z": [0.0, 0.0, 0.0]
        ... })
        >>> for index, row in blk.iterrows():
        ...     print(f"Row {index}: {row}")
        Row 0: {'id': 1, 'type': 'C', 'x': 0.0, 'y': 0.0, 'z': 0.0}
        Row 1: {'id': 2, 'type': 'O', 'x': 1.0, 'y': 0.0, 'z': 0.0}
        Row 2: {'id': 3, 'type': 'N', 'x': 2.0, 'y': 1.0, 'z': 0.0}

        Notes
        -----
        This method is similar to pandas DataFrame.iterrows() but returns
        a dictionary for each row instead of a pandas Series.
        """
        if not self._vars:
            return

        # Get the number of rows from the first variable
        nrows = self.nrows if n is None else n
        if nrows == 0:
            return

        # Get all variable names
        var_names = list(self._vars.keys())

        for i in range(nrows):
            row_data = {}
            for var_name in var_names:
                var_data = self._vars[var_name]
                if i < len(var_data):
                    # Handle scalar values
                    if var_data.ndim == 0:
                        row_data[var_name] = var_data.item()
                    else:
                        row_data[var_name] = var_data[i]
                else:
                    # Handle case where variable has fewer rows
                    row_data[var_name] = None

            yield i, row_data

    def itertuples(self, index: bool = True, name: str = "Row") -> Iterator[Any]:
        """
        Iterate over rows of the block as named tuples.

        Parameters
        ----------
        index : bool, default=True
            If True, include the row index as the first element
        name : str, default="Row"
            The name of the named tuple class

        Returns
        -------
        Iterator[Any]
            An iterator yielding named tuples for each row

        Examples
        --------
        >>> blk = Block({
        ...     "id": [1, 2, 3],
        ...     "type": ["C", "O", "N"],
        ...     "x": [0.0, 1.0, 2.0]
        ... })
        >>> for row in blk.itertuples():
        ...     print(f"Index: {row.Index}, ID: {row.id}, Type: {row.type}")
        Index: 0, ID: 1, Type: C
        Index: 1, ID: 2, Type: O
        Index: 2, ID: 3, Type: N

        Notes
        -----
        This method is similar to pandas DataFrame.itertuples().
        """
        from collections import namedtuple

        if not self._vars:
            return

        # Get the number of rows from the first variable
        nrows = self.nrows
        if nrows == 0:
            return

        # Get all variable names
        var_names = list(self._vars.keys())

        # Create field names for the named tuple
        if index:
            field_names = ["Index"] + var_names
        else:
            field_names = var_names

        # Create the named tuple class
        RowTuple = namedtuple(name, field_names)

        for i in range(nrows):
            row_values = []
            if index:
                row_values.append(i)

            for var_name in var_names:
                var_data = self._vars[var_name]
                if i < len(var_data):
                    # Handle scalar values
                    if var_data.ndim == 0:
                        row_values.append(var_data.item())
                    else:
                        row_values.append(var_data[i])
                else:
                    # Handle case where variable has fewer rows
                    row_values.append(None)

            yield RowTuple(*row_values)


class Frame:
    """
    Hierarchical numerical data container.

        Frame
        ├- blocks (dict[str, Block])
        ├- metadata   (dict[str, Any])  # metadata
    """

    def __init__(self, blocks=None, **props) -> None:
        # guarantee a root block even if none supplied
        self._blocks: dict[str, Block] = {}
        if blocks is not None:
            self._blocks = self._validate_and_convert_blocks(blocks)
        self.metadata = props

    def _validate_and_convert_blocks(self, blocks) -> dict[str, Block]:
        """
        Validate and convert input blocks to ensure all values are Block instances.

        This method recursively processes nested dictionaries and converts
        all leaf values to Block instances.

        Parameters
        ----------
        blocks : dict[str, Block] | dict[str, dict] | dict[str, Any]
            Input blocks. Can be:
            - dict[str, Block]: Already correct format
            - dict[str, dict]: Nested dictionaries that will be converted to Block
            - dict[str, Any]: Mixed format that will be validated and converted

        Returns
        -------
        dict[str, Block]
            Validated blocks where all values are Block instances

        Raises
        ------
        ValueError
            If any leaf value cannot be converted to Block
        """
        if not isinstance(blocks, dict):
            raise ValueError(f"blocks must be a dict, got {type(blocks)}")

        validated_blocks = {}

        for key, value in blocks.items():
            if not isinstance(key, str):
                raise ValueError(
                    f"Block keys must be strings, got {type(key)} for key {key}"
                )

            if isinstance(value, Block):
                # Already a Block, use as is
                validated_blocks[key] = value
            elif isinstance(value, dict):
                # Nested dict, convert to Block
                try:
                    validated_blocks[key] = Block(value)
                except Exception as e:
                    raise ValueError(
                        f"Failed to convert nested dict to Block for key '{key}': {e}"
                    )
            else:
                # Try to convert to Block (e.g., list, array, etc.)
                try:
                    validated_blocks[key] = Block(value)
                except Exception as e:
                    raise ValueError(
                        f"Failed to convert value to Block for key '{key}' (type {type(value)}): {e}"
                    )

        return validated_blocks

    # ---------- main get/set --------------------------------------------

    @overload
    def __getitem__(self, key: str) -> Block: ...  # str  → Block

    @overload
    def __getitem__(self, key: tuple[str, str]) -> np.ndarray: ...  # tuple→ ndarray

    def __getitem__(self, key: str | tuple[str, str]) -> np.ndarray | Block:
        if isinstance(key, tuple):
            grp, var = key
            return self._blocks[grp][var]
        return self._blocks[key]

    def __setitem__(self, key: str | tuple[str, str], value: BlockLike | Block):

        if isinstance(key, tuple):
            grp, var = key
            self._blocks.setdefault(grp, Block())[var] = value
        elif isinstance(key, str):
            if not isinstance(value, Block):
                value = Block(value)
            self._blocks[key] = value
        else:
            raise KeyError(
                f"Invalid key type: {type(key)}. Expected str, tuple[str, str]."
            )

    def __delitem__(self, key: str | tuple[str, str]) -> None:
        del self[key]

    def __contains__(self, key: str | tuple[str, str]) -> bool:
        if isinstance(key, tuple):
            grp, var = key
            return grp in self._blocks and var in self._blocks[grp]
        return key in self._blocks

    # ---------- helpers -------------------------------------------------
    def blocks(self) -> Iterator[str]:
        return iter(self._blocks)

    def variables(self, block: str) -> Iterator[str]:
        return iter(self._blocks[block])

    # ---------- (de)serialization --------------------------------------
    def to_dict(self) -> dict:
        block_dict = {g: grp.to_dict() for g, grp in self._blocks.items()}
        meta_dict = {k: v for k, v in self.metadata.items()}
        return {"blocks": block_dict, "metadata": meta_dict}

    @classmethod
    def from_dict(cls, data: dict) -> "Frame":
        blocks = {g: Block.from_dict(grp) for g, grp in data["blocks"].items()}
        frame = cls(blocks=blocks)
        frame.metadata = data.get("metadata", {})
        return frame

    # ---------- repr ----------------------------------------------------
    def __repr__(self) -> str:
        txt = ["Frame("]
        for g, grp in self._blocks.items():
            for k, v in grp._vars.items():
                txt.append(f"  [{g}] {k}: shape={v.shape}")
        return "\n".join(txt) + "\n)"

    def get_topology(self) -> Topology:
        """Get the topology of the frame."""
        bonds_block = self["bonds"]
        i = bonds_block["i"]
        j = bonds_block["j"]
        bonds_list = [(int(ii), int(jj)) for ii, jj in zip(i.tolist(), j.tolist())]
        n_atoms = self["atoms"].nrows
        topo = Topology()
        topo.add_atoms(n_atoms)
        topo.add_bonds(bonds_list)
        return topo
