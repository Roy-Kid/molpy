# Installation

MolPy is designed to work across multiple platforms and Python versions. This guide covers all installation methods and requirements.

## Requirements

- **Python**: 3.10 or higher
- **Operating Systems**: Windows, macOS, Linux
- **Architecture**: x86_64, ARM64 (Apple Silicon)

## Quick Install

The simplest way to install MolPy:

```bash
pip install molcrafts-molpy
```

## Installation with Extras

MolPy provides optional dependency groups for different use cases:

### Development Dependencies
```bash
pip install molcrafts-molpy[dev]
```
Includes: pytest, black, isort, pytest-cov

### Documentation Dependencies
```bash
pip install molcrafts-molpy[doc]
```
Includes: mkdocs, mkdocs-material, mkdocstrings

### All Extras
```bash
pip install molcrafts-molpy[all]
```

## Core Dependencies

MolPy automatically installs these required packages:

- **numpy** (≥2.0) - Numerical computing
- **molcrafts-molq** - Core quantum chemistry backend
- **python-igraph** (≥0.9.0) - Graph analysis
- **lark** (≥1.1.0) - SMILES/SMARTS parsing
- **pint** - Unit handling
- **freud-analysis** - Analysis algorithms

## Platform-Specific Notes

### Linux

Standard pip installation works on most Linux distributions:

```bash
pip install molcrafts-molpy
```

For HPC clusters, consider using conda:
```bash
conda install -c conda-forge molcrafts-molpy
```

### macOS

#### Intel Macs
```bash
pip install molcrafts-molpy
```

#### Apple Silicon (M1/M2)
```bash
pip install molcrafts-molpy
```
*Note: All dependencies are available as native ARM64 packages*

### Windows

#### Using pip
```bash
pip install molcrafts-molpy
```

#### Using conda (recommended for Windows)
```bash
conda install -c conda-forge molcrafts-molpy
```

## Development Installation

For contributors and developers:

```bash
# Clone the repository
git clone https://github.com/MolCrafts/molpy.git
cd molpy

# Install in editable mode with dev dependencies
pip install -e ".[dev]"

# Verify installation
python -c "import molpy; print(molpy.__version__)"
```

## Virtual Environments

We recommend using virtual environments:

```bash
# Create virtual environment
python -m venv molpy-env

# Activate (Linux/macOS)
source molpy-env/bin/activate

# Activate (Windows)
molpy-env\Scripts\activate

# Install MolPy
pip install molcrafts-molpy
```

## Conda Installation

```bash
# Create new environment
conda create -n molpy python=3.10

# Activate environment
conda activate molpy

# Install MolPy
conda install -c conda-forge molcrafts-molpy
```

## Verification

Test your installation:

```python
import molpy as mp
import numpy as np

# Create a simple box
box = mp.Box.cubic(10.0)
print(f"Box volume: {box.volume}")

# Create a frame
frame = mp.Frame(data={'atoms': {'x': [0], 'y': [0], 'z': [0]}}, box=box)
print(f"Frame created with {len(frame['atoms'])} atoms")
```

## Troubleshooting

### Common Issues

**Import Error**: If you get import errors, ensure you're using Python 3.10+:
```bash
python --version
```

**Permission Errors**: Use virtual environments or add `--user` flag:
```bash
pip install --user molcrafts-molpy
```

**Build Errors**: On some systems, you may need build tools:
```bash
# Ubuntu/Debian
sudo apt-get install build-essential python3-dev

# macOS
xcode-select --install

# Windows
# Install Visual Studio Build Tools
```

### Getting Help

- **GitHub Issues**: [Report installation problems](https://github.com/MolCrafts/molpy/issues)
- **Discussions**: [Community support](https://github.com/MolCrafts/molpy/discussions)
- **Documentation**: Check the [FAQ](faq.md) for common issues

## Next Steps

Once installed, proceed to:
- **[Quickstart Guide](quickstart.md)** - Your first MolPy simulation
- **[Quickstart](quickstart.md)** - Getting started with MolPy
- **[Tutorials](../tutorials/block-frame.md)** - Learning core concepts
