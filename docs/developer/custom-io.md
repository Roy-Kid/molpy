# Custom I/O Systems: Extending MolPy's File Handling

This chapter explains how to create custom I/O systems in MolPy. I/O systems are responsible for reading and writing molecular data in various formats, and MolPy provides a flexible framework for extending these capabilities.

## Understanding MolPy's I/O Architecture

### Base Classes and Hierarchy

MolPy's I/O system is built around several base classes that provide common functionality:

```python
# Data I/O base classes
class FileBase(ABC):
    """Common logic for Context-manager + lazy file handle."""

    def __init__(self, path: str | Path, mode: str, **open_kwargs):
        self._path = Path(path)
        self._mode = mode
        self._open_kwargs = open_kwargs
        self._fh: IO[Any] | None = None

    def __enter__(self):
        self._fh = self._path.open(self._mode, **self._open_kwargs)
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self._fh:
            self._fh.close()
            self._fh = None

    @property
    def fh(self) -> IO[Any]:
        if self._fh is None:
            self._fh = self._path.open(self._mode, **self._open_kwargs)
        return self._fh

class DataReader(FileBase, ABC):
    """Text reader that filters out blank lines."""

    def __init__(self, path: str | Path, **open_kwargs):
        super().__init__(path, mode="r", **open_kwargs)

    def _iter_nonblank(self) -> Iterator[str]:
        self.fh.seek(0)
        for raw in self.fh:
            line = raw.strip()
            if line:
                yield line

    @abstractmethod
    def read(self, frame: mp.Frame) -> mp.Frame:
        """Populate / update a Frame from the underlying text."""
        ...

class DataWriter(FileBase, ABC):
    """Text writer that guarantees the file is open in write-mode."""

    def __init__(self, path: str | Path, **open_kwargs):
        super().__init__(path, mode="w", **open_kwargs)

    @abstractmethod
    def write(self, frame: mp.Frame) -> None:
        """Serialize *frame* into the underlying text file."""
        ...
```

**Key Features**:
- **Context manager support**: Automatic file handling with `with` statements
- **Lazy file access**: Files are only opened when needed
- **Abstract methods**: Force implementation of core functionality
- **Common utilities**: Built-in blank line filtering and file management

### Trajectory I/O with Memory Mapping

MolPy's trajectory system uses memory mapping for efficient handling of large files:

```python
class TrajectoryReader(ABC):
    """Base class for trajectory file readers with memory mapping."""

    def __init__(self, fpath: Union[Path, str, List[Path], List[str]]):
        # Handle both single file and multiple files
        if isinstance(fpath, (str, Path)):
            self.fpaths = [Path(fpath)]
        else:
            self.fpaths = [Path(p) for p in fpath]

        self._frame_locations: List[FrameLocation] = []
        self._mms: List[mmap.mmap] = []
        self._total_frames = 0

        self._open_files()

    def _open_files(self):
        """Open trajectory files with memory mapping and build global index."""
        self._mms = []

        for file_index, fpath in enumerate(self.fpaths):
            # Open file
            fp = open(fpath, "rb")
            # Check if empty
            fp.seek(0, 2)
            if fp.tell() == 0:
                raise ValueError(f"File is empty: {fpath}")
            fp.seek(0)  # Seek back to beginning

            # Create memory mapping
            mm = mmap.mmap(fp.fileno(), 0, access=mmap.ACCESS_READ)
            self._mms.append(mm)

            # Parse this file to get frame locations
            self._parse_trajectory(file_index)

    def read_frame(self, index: int) -> "Frame":
        """Read a specific frame from the trajectory file(s)."""
        if index < 0:
            index = self._total_frames + index

        if index < 0 or index >= self._total_frames:
            raise IndexError(f"Frame index {index} out of range [0, {self._total_frames})")

        # Get location info for this frame
        location = self._get_frame_location(index)

        # Calculate frame end position
        if index + 1 < len(self._frame_locations):
            next_location = self._frame_locations[index + 1]
            if next_location.file_index == location.file_index:
                frame_end = next_location.byte_offset
            else:
                frame_end = None  # End of file
        else:
            frame_end = None  # Last frame

        # Get the memory-mapped file and read frame data
        mm = self._get_mmap(location.file_index)
        frame_bytes = mm[location.byte_offset:frame_end]
        frame_lines = frame_bytes.decode().splitlines()

        # Parse the frame lines using the derived class implementation
        return self._parse_frame(frame_lines)
```

**Key Features**:
- **Memory mapping**: Efficient access to large files without loading everything into memory
- **Frame indexing**: Pre-computed frame locations for fast random access
- **Multi-file support**: Handle trajectories split across multiple files
- **Lazy loading**: Only read frames when requested

## Creating Custom Data Readers

### When to Create Custom Readers

Create custom data readers when you need to:

1. **Support new file formats** not covered by built-in readers
2. **Handle domain-specific data** with custom parsing logic
3. **Optimize performance** for specific file types
4. **Add validation** or preprocessing to existing formats
5. **Integrate with external tools** or databases

### Basic Reader Structure

Every custom reader must inherit from `DataReader` and implement the `read` method:

```python
class CustomReader(DataReader):
    """Template for custom data readers."""

    def __init__(self, path: str | Path, **kwargs):
        super().__init__(path)
        # Store any additional parameters
        self.param1 = kwargs.get('param1', 'default')
        self.param2 = kwargs.get('param2', 0)

    def read(self, frame: mp.Frame) -> mp.Frame:
        """
        Read data from file and populate the frame.

        Args:
            frame: Existing frame to populate (or create new one)

        Returns:
            Populated frame
        """
        # Your parsing logic here
        # Must return a Frame object
        pass
```

### Reader Design Patterns

#### Pattern 1: Simple Text Parser

For simple text-based formats:

```python
class SimpleTextReader(DataReader):
    """Reader for simple text-based molecular formats."""

    def __init__(self, path: str | Path, comment_char: str = "#"):
        super().__init__(path)
        self.comment_char = comment_char

    def read(self, frame: mp.Frame) -> mp.Frame:
        """Parse simple text format into frame."""
        # Initialize data containers
        elements = []
        x_coords = []
        y_coords = []
        z_coords = []

        # Parse file line by line
        for line in self._iter_nonblank():
            # Skip comment lines
            if line.startswith(self.comment_char):
                continue

            # Parse data line (assuming: element x y z format)
            try:
                parts = line.split()
                if len(parts) >= 4:
                    element = parts[0]
                    x = float(parts[1])
                    y = float(parts[2])
                    z = float(parts[3])

                    elements.append(element)
                    x_coords.append(x)
                    y_coords.append(y)
                    z_coords.append(z)
            except (ValueError, IndexError):
                # Skip malformed lines
                continue

        # Create atoms block
        if elements:
            atoms = mp.Block({
                'element': elements,
                'x': x_coords,
                'y': y_coords,
                'z': z_coords
            })

            # Add to frame
            frame['atoms'] = atoms

        return frame

# Usage
reader = SimpleTextReader("molecule.txt")
frame = reader.read(mp.Frame())
print(f"Loaded {len(frame['atoms'])} atoms")
```

**Key Design Principles**:
- **Error handling**: Skip malformed lines gracefully
- **Flexible parsing**: Handle variations in format
- **Efficient iteration**: Use built-in `_iter_nonblank()` method

#### Pattern 2: Structured Format Parser

For more complex, structured formats:

```python
class StructuredFormatReader(DataReader):
    """Reader for structured molecular formats with headers and sections."""

    def __init__(self, path: str | Path):
        super().__init__(path)
        self.sections = {}

    def _parse_header(self) -> dict:
        """Parse file header section."""
        header = {}

        for line in self._iter_nonblank():
            if line.startswith('END_HEADER'):
                break

            if '=' in line:
                key, value = line.split('=', 1)
                header[key.strip()] = value.strip()

        return header

    def _parse_atoms_section(self) -> mp.Block:
        """Parse atoms section."""
        atoms_data = {
            'element': [],
            'x': [],
            'y': [],
            'z': [],
            'type': [],
            'charge': []
        }

        in_atoms = False
        for line in self._iter_nonblank():
            if line.startswith('BEGIN_ATOMS'):
                in_atoms = True
                continue
            elif line.startswith('END_ATOMS'):
                break

            if in_atoms:
                try:
                    parts = line.split()
                    if len(parts) >= 7:
                        atoms_data['element'].append(parts[0])
                        atoms_data['x'].append(float(parts[1]))
                        atoms_data['y'].append(float(parts[2]))
                        atoms_data['z'].append(float(parts[3]))
                        atoms_data['type'].append(int(parts[4]))
                        atoms_data['charge'].append(float(parts[5]))
                except (ValueError, IndexError):
                    continue

        # Convert to Block, handling missing data
        block_data = {}
        for key, values in atoms_data.items():
            if values:  # Only include non-empty arrays
                block_data[key] = values

        return mp.Block(block_data)

    def read(self, frame: mp.Frame) -> mp.Frame:
        """Parse structured format into frame."""
        # Parse header
        header = self._parse_header()

        # Parse atoms
        atoms = self._parse_atoms_section()

        # Add to frame
        frame['atoms'] = atoms

        # Store header information
        frame['metadata'] = header

        return frame

# Usage
reader = StructuredFormatReader("molecule.struc")
frame = reader.read(mp.Frame())
print(f"Loaded {len(frame['atoms'])} atoms")
print(f"Metadata: {frame['metadata']}")
```

**Key Design Principles**:
- **Section-based parsing**: Handle different parts of the file separately
- **Metadata preservation**: Store format-specific information
- **Flexible data handling**: Only include fields that have data

#### Pattern 3: Binary Format Reader

For binary formats:

```python
import struct

class BinaryFormatReader(DataReader):
    """Reader for binary molecular formats."""

    def __init__(self, path: str | Path, byte_order: str = '<'):
        super().__init__(path, mode='rb')  # Binary mode
        self.byte_order = byte_order

    def _read_header(self) -> dict:
        """Read binary header."""
        header = {}

        # Read magic number
        magic = self.fh.read(4)
        if magic != b'MOLP':
            raise ValueError("Invalid file format")

        # Read version
        version = struct.unpack(f'{self.byte_order}I', self.fh.read(4))[0]
        header['version'] = version

        # Read number of atoms
        n_atoms = struct.unpack(f'{self.byte_order}I', self.fh.read(4))[0]
        header['n_atoms'] = n_atoms

        return header

    def _read_atoms(self, n_atoms: int) -> mp.Block:
        """Read binary atom data."""
        # Pre-allocate arrays
        elements = []
        x_coords = np.zeros(n_atoms, dtype=np.float32)
        y_coords = np.zeros(n_atoms, dtype=np.float32)
        z_coords = np.zeros(n_atoms, dtype=np.float32)

        # Read atom data
        for i in range(n_atoms):
            # Read element (4 bytes, null-padded)
            element_bytes = self.fh.read(4)
            element = element_bytes.decode('ascii').rstrip('\x00')
            elements.append(element)

            # Read coordinates (3 floats)
            coords = struct.unpack(f'{self.byte_order}3f', self.fh.read(12))
            x_coords[i] = coords[0]
            y_coords[i] = coords[1]
            z_coords[i] = coords[2]

        return mp.Block({
            'element': elements,
            'x': x_coords,
            'y': y_coords,
            'z': z_coords
        })

    def read(self, frame: mp.Frame) -> mp.Frame:
        """Read binary format into frame."""
        # Read header
        header = self._read_header()

        # Read atoms
        atoms = self._read_atoms(header['n_atoms'])

        # Add to frame
        frame['atoms'] = atoms
        frame['metadata'] = header

        return frame

# Usage
reader = BinaryFormatReader("molecule.bin")
frame = reader.read(mp.Frame())
print(f"Loaded {len(frame['atoms'])} atoms from binary file")
```

**Key Design Principles**:
- **Binary mode**: Use `mode='rb'` for binary files
- **Structured unpacking**: Use `struct` module for binary data
- **Memory efficiency**: Pre-allocate arrays when possible

## Creating Custom Data Writers

### When to Create Custom Writers

Create custom data writers when you need to:

1. **Support new output formats** for external tools
2. **Customize data formatting** for specific applications
3. **Add metadata** or annotations to output files
4. **Optimize performance** for specific use cases
5. **Integrate with databases** or other systems

### Basic Writer Structure

Every custom writer must inherit from `DataWriter` and implement the `write` method:

```python
class CustomWriter(DataWriter):
    """Template for custom data writers."""

    def __init__(self, path: str | Path, **kwargs):
        super().__init__(path)
        # Store any additional parameters
        self.format_version = kwargs.get('version', '1.0')
        self.precision = kwargs.get('precision', 6)

    def write(self, frame: mp.Frame) -> None:
        """
        Write frame data to file.

        Args:
            frame: Frame to write
        """
        # Your writing logic here
        pass
```

### Writer Design Patterns

#### Pattern 1: Simple Text Writer

For simple text-based formats:

```python
class SimpleTextWriter(DataWriter):
    """Writer for simple text-based molecular formats."""

    def __init__(self, path: str | Path, precision: int = 6):
        super().__init__(path)
        self.precision = precision

    def write(self, frame: mp.Frame) -> None:
        """Write frame to simple text format."""
        if 'atoms' not in frame:
            raise ValueError("Frame must contain atoms")

        atoms = frame['atoms']

        # Write header
        self.fh.write(f"# MolPy Simple Format\n")
        self.fh.write(f"# {len(atoms)} atoms\n")
        self.fh.write(f"# element x y z\n")

        # Write atom data
        for i in range(len(atoms)):
            element = atoms['element'][i]
            x = atoms['x'][i]
            y = atoms['y'][i]
            z = atoms['z'][i]

            line = f"{element:<4} {x:>{self.precision+6}.{self.precision}f} "
            line += f"{y:>{self.precision+6}.{self.precision}f} "
            line += f"{z:>{self.precision+6}.{self.precision}f}\n"

            self.fh.write(line)

# Usage
writer = SimpleTextWriter("output.txt", precision=8)
writer.write(frame)
```

**Key Design Principles**:
- **Formatted output**: Use proper spacing and precision
- **Header information**: Include metadata and format information
- **Error checking**: Validate frame contents before writing

#### Pattern 2: Structured Format Writer

For more complex, structured formats:

```python
class StructuredFormatWriter(DataWriter):
    """Writer for structured molecular formats with headers and sections."""

    def __init__(self, path: str | Path, include_metadata: bool = True):
        super().__init__(path)
        self.include_metadata = include_metadata

    def _write_header(self, frame: mp.Frame) -> None:
        """Write file header section."""
        self.fh.write("BEGIN_HEADER\n")
        self.fh.write(f"format=MolPy_Structured\n")
        self.fh.write(f"version=1.0\n")
        self.fh.write(f"n_atoms={len(frame['atoms'])}\n")

        if self.include_metadata and 'metadata' in frame:
            for key, value in frame['metadata'].items():
                self.fh.write(f"{key}={value}\n")

        self.fh.write("END_HEADER\n\n")

    def _write_atoms_section(self, frame: mp.Frame) -> None:
        """Write atoms section."""
        atoms = frame['atoms']

        self.fh.write("BEGIN_ATOMS\n")

        # Determine available fields
        fields = ['element', 'x', 'y', 'z']
        if 'type' in atoms:
            fields.append('type')
        if 'charge' in atoms:
            fields.append('charge')

        # Write field header
        self.fh.write(f"# {' '.join(fields)}\n")

        # Write atom data
        for i in range(len(atoms)):
            line_parts = []
            for field in fields:
                if field in ['x', 'y', 'z']:
                    line_parts.append(f"{atoms[field][i]:12.6f}")
                elif field == 'type':
                    line_parts.append(f"{atoms[field][i]:6d}")
                elif field == 'charge':
                    line_parts.append(f"{atoms[field][i]:10.4f}")
                else:
                    line_parts.append(f"{atoms[field][i]:<6}")

            self.fh.write(f"{' '.join(line_parts)}\n")

        self.fh.write("END_ATOMS\n")

    def write(self, frame: mp.Frame) -> None:
        """Write frame to structured format."""
        if 'atoms' not in frame:
            raise ValueError("Frame must contain atoms")

        # Write header
        self._write_header(frame)

        # Write atoms
        self._write_atoms_section(frame)

# Usage
writer = StructuredFormatWriter("output.struc", include_metadata=True)
writer.write(frame)
```

**Key Design Principles**:
- **Section-based writing**: Organize output into logical sections
- **Flexible field handling**: Only write fields that exist
- **Consistent formatting**: Use consistent spacing and precision

#### Pattern 3: Binary Format Writer

For binary formats:

```python
class BinaryFormatWriter(DataWriter):
    """Writer for binary molecular formats."""

    def __init__(self, path: str | Path, byte_order: str = '<'):
        super().__init__(path, mode='wb')  # Binary mode
        self.byte_order = byte_order

    def _write_header(self, frame: mp.Frame) -> None:
        """Write binary header."""
        atoms = frame['atoms']

        # Write magic number
        self.fh.write(b'MOLP')

        # Write version
        version = 1
        self.fh.write(struct.pack(f'{self.byte_order}I', version))

        # Write number of atoms
        n_atoms = len(atoms)
        self.fh.write(struct.pack(f'{self.byte_order}I', n_atoms))

    def _write_atoms(self, frame: mp.Frame) -> None:
        """Write binary atom data."""
        atoms = frame['atoms']

        for i in range(len(atoms)):
            # Write element (4 bytes, null-padded)
            element = atoms['element'][i]
            element_bytes = element.ljust(4, '\x00').encode('ascii')
            self.fh.write(element_bytes)

            # Write coordinates (3 floats)
            x = float(atoms['x'][i])
            y = float(atoms['y'][i])
            z = float(atoms['z'][i])

            coords = struct.pack(f'{self.byte_order}3f', x, y, z)
            self.fh.write(coords)

    def write(self, frame: mp.Frame) -> None:
        """Write frame to binary format."""
        if 'atoms' not in frame:
            raise ValueError("Frame must contain atoms")

        # Write header
        self._write_header(frame)

        # Write atoms
        self._write_atoms(frame)

# Usage
writer = BinaryFormatWriter("output.bin")
writer.write(frame)
```

**Key Design Principles**:
- **Binary mode**: Use `mode='wb'` for binary files
- **Structured packing**: Use `struct` module for binary data
- **Data validation**: Ensure data types are correct before writing

## Creating Custom Trajectory Readers

### When to Create Custom Trajectory Readers

Create custom trajectory readers when you need to:

1. **Support new trajectory formats** not covered by built-in readers
2. **Handle domain-specific trajectory data** with custom parsing
3. **Optimize performance** for specific trajectory types
4. **Add preprocessing** or filtering to trajectory data
5. **Integrate with external trajectory tools**

### Basic Trajectory Reader Structure

Every custom trajectory reader must inherit from `TrajectoryReader` and implement abstract methods:

```python
class CustomTrajectoryReader(TrajectoryReader):
    """Template for custom trajectory readers."""

    def _parse_trajectory(self, file_index: int):
        """Parse trajectory file at given index, storing frame locations."""
        mm = self._get_mmap(file_index)

        # Your parsing logic here
        # Must populate self._frame_locations and update self._total_frames
        pass

    def _parse_frame(self, frame_lines: List[str]) -> "Frame":
        """Parse frame lines into a Frame object."""
        # Your frame parsing logic here
        # Must return a Frame object
        pass
```

### Trajectory Reader Design Patterns

#### Pattern 1: Line-Based Trajectory Parser

For trajectory formats where frames are separated by specific markers:

```python
class LineBasedTrajectoryReader(TrajectoryReader):
    """Reader for line-based trajectory formats."""

    def __init__(self, fpath: Union[Path, str, List[Path], List[str]],
                 frame_marker: str = "FRAME"):
        super().__init__(fpath)
        self.frame_marker = frame_marker.encode('ascii')

    def _parse_trajectory(self, file_index: int):
        """Parse trajectory file to find frame locations."""
        mm = self._get_mmap(file_index)

        # Find all frame markers
        frame_offsets = []
        start = 0

        while True:
            pos = mm.find(self.frame_marker, start)
            if pos == -1:
                break

            frame_offsets.append(pos)
            start = pos + 1

        # Create frame locations
        from .base import FrameLocation
        for offset in frame_offsets:
            location = FrameLocation(
                file_index=file_index,
                byte_offset=offset,
                file_path=self.fpaths[file_index]
            )
            self._frame_locations.append(location)

        self._total_frames += len(frame_offsets)

    def _parse_frame(self, frame_lines: List[str]) -> "Frame":
        """Parse frame lines into a Frame object."""
        frame = mp.Frame()

        # Skip frame marker line
        if frame_lines and frame_lines[0].startswith("FRAME"):
            frame_lines = frame_lines[1:]

        # Parse atom data
        elements = []
        x_coords = []
        y_coords = []
        z_coords = []

        for line in frame_lines:
            if line.strip() and not line.startswith('#'):
                try:
                    parts = line.split()
                    if len(parts) >= 4:
                        elements.append(parts[0])
                        x_coords.append(float(parts[1]))
                        y_coords.append(float(parts[2]))
                        z_coords.append(float(parts[3]))
                except (ValueError, IndexError):
                    continue

        if elements:
            atoms = mp.Block({
                'element': elements,
                'x': x_coords,
                'y': y_coords,
                'z': z_coords
            })
            frame['atoms'] = atoms

        return frame

# Usage
reader = LineBasedTrajectoryReader("trajectory.txt", frame_marker="FRAME")
print(f"Found {len(reader)} frames")
frame = reader.read_frame(0)
```

**Key Design Principles**:
- **Efficient parsing**: Use memory mapping for large files
- **Frame indexing**: Pre-compute frame locations for fast access
- **Flexible parsing**: Handle variations in frame format

#### Pattern 2: Binary Trajectory Parser

For binary trajectory formats:

```python
class BinaryTrajectoryReader(TrajectoryReader):
    """Reader for binary trajectory formats."""

    def __init__(self, fpath: Union[Path, str, List[Path], List[str]],
                 frame_size: int = 1024):
        super().__init__(fpath)
        self.frame_size = frame_size

    def _parse_trajectory(self, file_index: int):
        """Parse binary trajectory file to find frame locations."""
        mm = self._get_mmap(file_index)
        file_size = len(mm)

        # Calculate frame locations based on fixed frame size
        frame_offsets = []
        offset = 0

        while offset < file_size:
            frame_offsets.append(offset)
            offset += self.frame_size

        # Create frame locations
        from .base import FrameLocation
        for offset in frame_offsets:
            location = FrameLocation(
                file_index=file_index,
                byte_offset=offset,
                file_path=self.fpaths[file_index]
            )
            self._frame_locations.append(location)

        self._total_frames += len(frame_offsets)

    def _parse_frame(self, frame_lines: List[str]) -> "Frame":
        """Parse binary frame data into a Frame object."""
        # Convert frame_lines (which are actually bytes) to frame
        frame_bytes = b''.join(line.encode('latin1') for line in frame_lines)

        frame = mp.Frame()

        # Parse binary data
        if len(frame_bytes) >= 12:  # Minimum size for header
            # Read header
            n_atoms = struct.unpack('<I', frame_bytes[0:4])[0]
            timestamp = struct.unpack('<f', frame_bytes[4:8])[0]

            # Read atom data
            if len(frame_bytes) >= 12 + n_atoms * 16:  # 16 bytes per atom
                elements = []
                x_coords = np.zeros(n_atoms, dtype=np.float32)
                y_coords = np.zeros(n_atoms, dtype=np.float32)
                z_coords = np.zeros(n_atoms, dtype=np.float32)

                for i in range(n_atoms):
                    offset = 12 + i * 16

                    # Read element (4 bytes)
                    element_bytes = frame_bytes[offset:offset+4]
                    element = element_bytes.decode('ascii').rstrip('\x00')
                    elements.append(element)

                    # Read coordinates (12 bytes, 3 floats)
                    coords = struct.unpack('<3f', frame_bytes[offset+4:offset+16])
                    x_coords[i] = coords[0]
                    y_coords[i] = coords[1]
                    z_coords[i] = coords[2]

                atoms = mp.Block({
                    'element': elements,
                    'x': x_coords,
                    'y': y_coords,
                    'z': z_coords
                })
                frame['atoms'] = atoms

                # Store metadata
                frame['timestamp'] = timestamp

        return frame

# Usage
reader = BinaryTrajectoryReader("trajectory.bin", frame_size=2048)
print(f"Found {len(reader)} frames")
frame = reader.read_frame(0)
```

**Key Design Principles**:
- **Fixed frame size**: Use known frame sizes for efficient indexing
- **Binary parsing**: Handle binary data carefully with proper struct unpacking
- **Error handling**: Validate frame data sizes and contents

## Best Practices for I/O Development

### 1. **Always Use Context Managers**

```python
# ✅ Good: Automatic file handling
with CustomReader("file.txt") as reader:
    frame = reader.read(mp.Frame())

# ❌ Avoid: Manual file management
reader = CustomReader("file.txt")
try:
    frame = reader.read(mp.Frame())
finally:
    reader.fh.close()
```

### 2. **Handle Errors Gracefully**

```python
def read(self, frame: mp.Frame) -> mp.Frame:
    try:
        # Your parsing logic here
        pass
    except Exception as e:
        raise IOError(f"Failed to read file: {e}") from e
```

### 3. **Validate Input Data**

```python
def write(self, frame: mp.Frame) -> None:
    if 'atoms' not in frame:
        raise ValueError("Frame must contain atoms")

    atoms = frame['atoms']
    if len(atoms) == 0:
        raise ValueError("Frame must contain at least one atom")

    # Proceed with writing
```

### 4. **Use Efficient Data Structures**

```python
# ✅ Good: Pre-allocate arrays
x_coords = np.zeros(n_atoms, dtype=np.float32)
for i, x in enumerate(x_values):
    x_coords[i] = x

# ❌ Avoid: Growing lists
x_coords = []
for x in x_values:
    x_coords.append(x)
x_coords = np.array(x_coords)
```

### 5. **Document Your Formats**

```python
class CustomFormatReader(DataReader):
    """
    Reader for Custom Molecular Format (CMF).

    Format specification:
    - Header: # Custom Format v1.0
    - Atoms: element x y z [type] [charge]
    - Comments: Lines starting with #

    Examples:
        >>> reader = CustomFormatReader("molecule.cmf")
        >>> frame = reader.read(mp.Frame())
    """
```

## Summary

MolPy's I/O system provides a flexible and efficient framework for handling molecular data:

- **Base classes**: Common functionality for file handling and context management
- **Memory mapping**: Efficient access to large trajectory files
- **Extensible design**: Easy to create new readers and writers for any format
- **Performance optimization**: Built-in support for efficient data structures
- **Error handling**: Robust error handling and validation

The key insight is that I/O systems should be **format-agnostic** and **performance-focused**, providing a clean interface for reading and writing molecular data in any format.

### Next Steps

Continue exploring MolPy's extension capabilities by learning about custom analysis pipelines, integration with external tools, and advanced development patterns.
