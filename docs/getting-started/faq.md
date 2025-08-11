# Frequently Asked Questions

This FAQ addresses common questions about MolPy. If you don't find your answer here, please create a [GitHub issue](https://github.com/MolCrafts/molpy/issues).

## Getting Started

### What is MolPy?

MolPy is a modern, high-performance Python framework for molecular simulation and analysis. It provides elegant data structures and comprehensive toolkits for working with atomic structures, molecular dynamics trajectories, and force field data.

### What makes MolPy different from other molecular simulation packages?

MolPy is designed with several key principles:

- **LLM-friendly**: Clean APIs with explicit typing and minimal magic
- **Composable**: Build complex systems from simple, reusable parts
- **High-performance**: Optimized for large-scale molecular systems
- **Extensible**: Plugin-ready architecture for custom functionality
- **Unified interface**: Consistent patterns across all components

### What Python versions does MolPy support?

MolPy requires **Python 3.10 or higher**. This ensures compatibility with modern Python features and performance improvements.

### Is MolPy free to use?

Yes! MolPy is open-source software released under the BSD-3-Clause license. You can use it freely for academic, commercial, or personal projects.

## Installation & Setup

### How do I install MolPy?

The simplest way is using pip:

```bash
pip install molcrafts-molpy
```

For development and documentation:
```bash
pip install molcrafts-molpy[dev,doc]
```

### I'm getting installation errors. What should I do?

1. **Check Python version**: Ensure you have Python 3.10+
2. **Use virtual environment**: Avoid dependency conflicts
3. **Install build tools**: Some systems need additional compilers
4. **Try conda**: `conda install -c conda-forge molcrafts-molpy`

### Can I install MolPy on Windows/macOS/Linux?

Yes! MolPy supports all major operating systems:
- **Windows**: Use pip or conda (conda recommended)
- **macOS**: Works on both Intel and Apple Silicon
- **Linux**: Standard pip installation works on most distributions

### What are the system requirements?

- **Python**: 3.10 or higher
- **Memory**: 4GB RAM minimum, 8GB+ recommended
- **Storage**: 1GB free space for installation
- **OS**: Windows 10+, macOS 10.14+, or Linux (kernel 3.0+)

## Core Concepts

### What is a Block in MolPy?

A `Block` is a lightweight container that maps variable names to NumPy arrays. It's the atomic unit for storing data like atomic coordinates, properties, and connectivity.

```python
import molpy as mp

# Create atomic data
atoms = mp.Block({
    'x': [0.0, 1.0, 2.0],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C']
})

print(f"Block has {atoms.nrows} rows")
print(f"Available variables: {list(atoms.variables())}")
```

### What is a Frame in MolPy?

A `Frame` is a hierarchical container that holds multiple blocks, simulation box, and metadata. It represents a complete molecular system at a single point in time.

```python
# Create a frame with atoms and box
frame = mp.Frame(
    data={'atoms': atoms_data},
    box=mp.Box.cubic(10.0),
    metadata={'temperature': 300.0}
)

print(f"Frame has {len(frame['atoms'])} atoms")
print(f"Box volume: {frame.box.volume:.1f} Å³")
```

### How does MolPy handle periodic boundary conditions?

MolPy uses `Box` objects to define simulation boundaries and periodic conditions:

```python
# Different box types
cubic = mp.Box.cubic(10.0)                    # 10×10×10 Å³
orth = mp.Box.orth([8.0, 10.0, 12.0])        # 8×10×12 Å³
tric = mp.Box.tric([10.0, 10.0, 10.0],       # Triclinic with tilts
                    angles=[90, 90, 90])

# Check properties
print(f"Box volume: {cubic.volume:.1f} Å³")
print(f"Periodic boundaries: {cubic.pbc}")
```

## Data & File Formats

### What file formats does MolPy support?

MolPy supports many common formats:

- **Trajectory files**: LAMMPS dump, AMBER netCDF, GROMACS XTC
- **Structure files**: PDB, XYZ, LAMMPS data
- **Force field files**: OPLS-AA, AMBER, custom formats
- **Data files**: HDF5, CSV, JSON

### Can I read my existing simulation data?

Yes! MolPy is designed to work with existing data. If your format isn't directly supported, you can:

1. **Convert to supported format** using external tools
2. **Write custom reader** using MolPy's I/O framework
3. **Use intermediate format** like HDF5 or CSV

### How do I handle large trajectory files?

MolPy provides several strategies for large files:

- **Lazy loading**: Only load data when needed
- **Chunked processing**: Process frames in batches
- **Memory mapping**: Efficient access to disk-based data
- **Parallel processing**: Distribute work across cores

```python
# Process large trajectory in chunks
trajectory = mp.Trajectory('large_file.h5')
for i in range(0, len(trajectory), 100):
    chunk = trajectory[i:i+100]
    process_chunk(chunk)
```

## Performance & Optimization

### How do I optimize performance for large systems?

Several strategies can help:

1. **Use NumPy operations**: Vectorize calculations when possible
2. **Process in chunks**: Avoid loading entire datasets into memory
3. **Use efficient data types**: Choose appropriate NumPy dtypes
4. **Leverage parallel processing**: Use multiprocessing for independent operations

```python
# Vectorized distance calculation
positions = frame['atoms'][['x', 'y', 'z']]
distances = np.linalg.norm(positions[1:] - positions[:-1], axis=1)
```

### Does MolPy support GPU acceleration?

Currently, MolPy focuses on CPU-based computation with NumPy optimization. For GPU acceleration, you can:

1. **Use CuPy**: Convert MolPy data to CuPy arrays
2. **Integrate with PyTorch/TensorFlow**: Export data for GPU computation
3. **Use external GPU libraries**: Interface with CUDA-enabled packages

## Common Issues & Solutions

### Import Error: No module named 'molpy'

**Problem**: MolPy not found in your Python environment.

**Solutions**:
```python
# Check your Python environment
import sys
print(f"Python executable: {sys.executable}")
print(f"Python version: {sys.version}")

# Verify MolPy installation
import molpy
print(f"MolPy location: {molpy.__file__}")
```

**Common causes**:
- Wrong Python environment activated
- MolPy not installed in current environment
- PATH issues

### Memory Error when loading large files

**Problem**: System runs out of memory loading large trajectories.

**Solutions**:
```python
# Use lazy loading
trajectory = mp.Trajectory('large_file.h5')

# Process in chunks
for i in range(0, len(trajectory), 1000):
    chunk = trajectory[i:i+1000]
    process_chunk(chunk)
    del chunk  # Free memory
```

### Data Type Issues

**Problem**: Errors related to data types or array shapes.

**Solutions**:
```python
# Ensure numeric data for coordinates
atoms_data = {
    'x': [0.0, 1.0, 2.0],  # Must be numeric
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C']  # Can be strings
}

# Check data consistency
print(f"X coordinates: {len(atoms_data['x'])}")
print(f"Y coordinates: {len(atoms_data['y'])}")
print(f"Z coordinates: {len(atoms_data['z'])}")
```

## Advanced Usage

### How do I create custom potential functions?

MolPy provides a flexible framework for custom potentials:

```python
class CustomPotential(mp.Potential):
    def __init__(self, k=100.0, r0=1.0):
        self.k = k
        self.r0 = r0

    def calc_energy(self, frame):
        # Your energy calculation here
        return energy

    def calc_forces(self, frame):
        # Your force calculation here
        return forces

# Use custom potential
potential = CustomPotential(k=100.0, r0=1.0)
energy = potential.calc_energy(frame)
```

### How do I extend MolPy with custom functionality?

MolPy is designed for extensibility:

1. **Create custom classes**: Inherit from MolPy base classes
2. **Use composition**: Combine existing components
3. **Plugin system**: Register custom functionality
4. **Custom I/O**: Implement readers/writers for new formats

### Can I use MolPy with other molecular simulation packages?

Yes! MolPy is designed for interoperability:

- **LAMMPS**: Read/write data files and trajectories
- **AMBER**: Import/export prmtop and coordinate files
- **GROMACS**: Work with gro and top files
- **OpenMM**: Convert between data structures
- **Custom formats**: Extend I/O capabilities

## Getting Help

### Where can I find more documentation?

- **[Tutorials](../tutorials/block-frame.md)**: Step-by-step guides
- **[How-to Guides](../how-to/io-workflows.md)**: Practical examples
- **[API Reference](../reference/index.md)**: Complete function reference
- **[GitHub Repository](https://github.com/MolCrafts/molpy)**: Source code and issues

### How do I report bugs or request features?

1. **Check existing issues**: Search [GitHub Issues](https://github.com/MolCrafts/molpy/issues)
2. **Create new issue**: Provide detailed description and minimal example
3. **Include system info**: Python version, OS, MolPy version
4. **Share code**: Minimal reproducible example

### Can I contribute to MolPy?

Absolutely! MolPy welcomes contributions:

- **Bug reports**: Help improve stability
- **Feature requests**: Suggest new functionality
- **Code contributions**: Submit pull requests
- **Documentation**: Improve guides and examples
- **Testing**: Help ensure quality

## Summary

This FAQ covered the most common questions about MolPy. Remember:

- **Start simple**: Use the [Quickstart Guide](quickstart.md) to get familiar
- **Check tutorials**: Learn advanced features step by step
- **Use examples**: Copy and modify working code
- **Ask for help**: Create GitHub issues for problems

MolPy is designed to be user-friendly and extensible. If you have questions not covered here, the community is ready to help!
