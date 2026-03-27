# Adding an I/O Format

This page shows how to add readers and writers for new file formats and force field backends.

## Data file readers and writers

Subclass `DataReader` or `DataWriter` from `molpy.io.data.base`.

### Reader

```python
from pathlib import Path
from molpy.core.frame import Frame, Block
from molpy.io.data.base import DataReader

class MyFormatReader(DataReader):
    """Read .myformat files into a Frame."""

    def __init__(self, file: Path, **kwargs):
        super().__init__(path=file, **kwargs)

    def read(self, frame: Frame | None = None) -> Frame:
        if frame is None:
            frame = Frame()

        # Parse the file (self._path is set by FileBase)
        with open(self._path) as f:
            lines = f.readlines()

        # Populate blocks
        frame["atoms"] = Block({
            "element": [...],
            "x": [...],
            "y": [...],
            "z": [...],
        })
        return frame
```

### Writer

```python
from molpy.io.data.base import DataWriter

class MyFormatWriter(DataWriter):
    """Write a Frame to .myformat."""

    def __init__(self, file: Path, **kwargs):
        super().__init__(path=file, **kwargs)

    def write(self, frame: Frame) -> None:
        atoms = frame["atoms"]
        with open(self._path, "w") as f:
            for i in range(atoms.nrows):
                f.write(f"{atoms['element'][i]} {atoms['x'][i]} ...\n")
```

### Register in factory functions

Add your reader/writer to `molpy/io/readers.py` and `molpy/io/writers.py` so they are accessible via `mp.io.read_myformat()` and `mp.io.write_myformat()`.


## Force field writers with custom formatters

The force field export system uses a **formatter registry**. Each `ForceFieldWriter` subclass has its own isolated set of formatters. When you add a custom `Style` (e.g., `MorseBondStyle`), you must register formatters for each backend that should support it.

### The formatter contract

A **style formatter** receives a `Style` and returns `FormattedParams`:

```python
def format_my_style(style) -> FormattedParams:
    return FormattedParams(keyword={"source": "custom"})
```

A **type formatter** receives a `Type` and a `Style` and returns `FormattedParams`:

```python
def format_my_type(typ, style) -> FormattedParams:
    return FormattedParams(positional=[typ["k0"], typ["r0"]])
```

`FormattedParams` carries `positional` (list) and `keyword` (dict). Each writer interprets them according to its own output rules (LAMMPS: space-separated tokens; GROMACS: column values; XML: attributes).

### Registering formatters

```python
from molpy.io.forcefield import LAMMPSForceFieldWriter
from molpy.io.forcefield.base import FormattedParams

LAMMPSForceFieldWriter.register_style_formatter(
    MorseBondStyle,
    format_my_style,
    style_name="morse",
)
LAMMPSForceFieldWriter.register_formatter(
    MorseBondStyle,
    format_my_type,
    style_name="morse",
)
```

Registrations are **isolated per writer subclass** — adding a formatter to one writer does not affect another. This isolation is enforced by `__init_subclass__` copying the registry.


## Trajectory readers and writers

Trajectory readers use memory-mapped files and a persistent frame index for efficient random access. Subclass `BaseTrajectoryReader` and implement `_scan_frames` (build byte-offset index) and `_parse_frame_bytes` (parse one frame):

```python
import mmap
from molpy.io.trajectory.base import BaseTrajectoryReader
from molpy.io.trajectory.index import FrameEntry
from molpy.core.frame import Frame

class MyTrajectoryReader(BaseTrajectoryReader):
    _format_id = "myformat"

    def _scan_frames(self, file_idx: int, mm: mmap.mmap) -> list[FrameEntry]:
        entries = []
        # scan file for frame boundaries, record byte offsets
        return entries

    def _parse_frame_bytes(self, mm: mmap.mmap, entry: FrameEntry) -> Frame:
        # parse one frame from mm[entry.offset:entry.offset+entry.length]
        return frame
```

The persistent index (`.tridx`) is built automatically on first read and cached for subsequent accesses. Subclass `TrajectoryWriter` and implement `write_frame()` for writing.


## Checklist

- [ ] Subclass `DataReader`/`DataWriter` or `BaseTrajectoryReader`/`TrajectoryWriter`
- [ ] Implement `read()` or `write()` with correct Frame/Trajectory signatures
- [ ] Add factory function in `readers.py` / `writers.py`
- [ ] Register force field formatters if adding a custom Style
- [ ] Write round-trip tests (`write → read → compare`) in `tests/test_io/`
