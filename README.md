# MolPy 🧬

[![Python](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/downloads/)
[![License](https://img.shields.io/badge/license-BSD-green.svg)](LICENSE)
[![Documentation](https://img.shields.io/badge/docs-mkdocs-blue.svg)](https://molcrafts.github.io/molpy)

**MolPy** is a modular Python toolkit for molecular modeling, simulation setup, and structural analysis. It provides a unified framework for manipulating atomic structures, generating force field inputs, and analyzing simulation data.

## ✨ Key Features

- **🏗️ Modular Architecture**: Clean, composable design with composition-over-inheritance philosophy
- **📊 Flexible Data Structures**: Frame/Block system for atomic data, Trajectory for time series
- **🎯 Advanced Selection System**: Boolean algebra for atom selection with custom predicates
- **🔧 Force Field Management**: Comprehensive typing and parameter management
- **📁 Multi-format I/O**: Support for PDB, XYZ, LAMMPS, AMBER, GROMACS formats
- **🧮 Analysis Tools**: Built-in analysis for diffusion, clustering, and molecular properties
- **🏗️ Molecular Building**: Polymer construction, bulk system generation, reaction modeling
- **📦 Packing Algorithms**: Molecular packing with Packmol integration and optimization
- **🔌 Extensible Design**: Custom wrappers, selections, and I/O components
- **⚡ High Performance**: Memory-mapped trajectory reading, efficient data structures

## 🚀 Quick Start

### Installation

```bash
# Clone the repository
git clone https://github.com/MolCrafts/molpy.git
cd molpy

# Install in development mode
pip install -e .
```

### Basic Usage

```python
import molpy as mp
import numpy as np

# Create a simulation box
box = mp.Box.cubic(10.0)

# Create atoms with coordinates and properties
atoms_data = {
    'x': [0.0, 1.0, 2.0],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C'],
    'type': [1, 1, 1]
}

# Create a frame with atoms and box
frame = mp.Frame(data={'atoms': atoms_data}, box=box)

# Select specific atoms using the selection system
carbon_atoms = frame['atoms'][mp.AtomTypeSelection(1)]
print(f"Selected {len(carbon_atoms)} carbon atoms")

# Save to file
frame.save('system.pdb')
```

## 📚 Documentation

Comprehensive documentation is available at [https://molcrafts.github.io/molpy](https://molcrafts.github.io/molpy)

### 📖 Getting Started
- **Quickstart Guide**: Basic concepts and first steps
- **Installation**: Setup and configuration
- **FAQ**: Common questions and solutions

### 🎓 Tutorials
- **Core Concepts**: Frame/Block system, Selection algebra, Systems
- **Data Structures**: Atoms, Bonds, Topology, Force Fields
- **Analysis**: Trajectory analysis, MSD calculations, Clustering
- **Building**: Molecular construction, Polymer building, Packing

### 🔧 How-to Guides
- **I/O Workflows**: Reading/writing various formats
- **Molecular Building**: Creating complex systems
- **Potential Calculations**: Energy and force computations
- **Molecular Packing**: System generation and optimization

### 👨‍💻 Developer Guides
- **Custom Wrappers**: Extending MolPy with custom functionality
- **Custom Selections**: Building new selection predicates
- **Custom I/O**: Adding new file format support

### 📖 API Reference
- **Core Modules**: Fundamental data structures and classes
- **I/O Modules**: File format support and trajectory handling
- **Analysis Modules**: Property calculation tools
- **Builder Modules**: Molecular construction tools
- **Operations Modules**: Transformation and manipulation
- **Potential Modules**: Energy and force calculations
- **Packing Modules**: System generation algorithms
- **Typifier Modules**: Automatic atom typing
- **Engine Modules**: Simulation engine interfaces

## 🏗️ Architecture

MolPy is built around several core design principles:

- **Composition over Inheritance**: Flexible wrapper system for extending functionality
- **Selection Algebra**: Boolean logic for atom selection with custom predicates
- **Memory Efficiency**: Lazy loading and memory mapping for large trajectories
- **Type Safety**: Comprehensive type hints and validation
- **Modular Design**: Clean separation of concerns across modules

### Core Components

- **Frame/Block**: Container system for atomic data with flexible indexing
- **Selection System**: Boolean algebra for atom selection with custom predicates
- **Wrapper System**: Composition-based extension mechanism
- **Trajectory I/O**: Memory-mapped reading with lazy loading
- **Force Field System**: Comprehensive typing and parameter management

## 📋 Dependencies

**Core Dependencies:**
- [numpy](https://github.com/numpy/numpy) - Numerical computing
- [python-igraph](https://github.com/igraph/python-igraph) - Graph analysis
- [lark](https://github.com/lark-parser/lark) - SMARTS/SMILES parsing
- [pint](https://github.com/hgrecco/pint) - Physical quantities and units
- [freud-analysis](https://github.com/glotzerlab/freud) - Analysis algorithms

**Optional Dependencies:**
- [mkdocs](https://github.com/mkdocs/mkdocs) - Documentation generation
- [mkdocs-material](https://github.com/squidfunk/mkdocs-material) - Documentation theme
- [mkdocstrings](https://github.com/mkdocstrings/mkdocstrings) - API documentation

## 🌟 Ecosystem

### 🧠 Machine Learning Integration
**[MolNex](https://github.com/MolCrafts/molnex)** - Universal potential training platform
- Neural network potential development
- Transfer learning for molecular systems
- Integration with popular ML frameworks

### 🎨 Interactive Visualization
**[MolVis](https://github.com/Roy-Kid/molvis)** - Production-level visualization
- WebGL-accelerated rendering
- Real-time molecular manipulation
- Interactive debugging and analysis tools

### 🚀 High Performance Computing
**[MolCPP](https://github.com/MolCrafts/molcpp)** - C++ backend for performance-critical operations

## 🤝 Contributing

We welcome contributions from the community! Here's how you can help:

1. **🐛 Bug Reports**: Use [GitHub Issues](https://github.com/MolCrafts/molpy/issues) to report bugs
2. **💡 Feature Requests**: Suggest new features and improvements
3. **📖 Documentation**: Help improve documentation and examples
4. **🔧 Code Contributions**: Submit pull requests with new features or fixes

### Development Setup
```bash
# Clone and setup development environment
git clone https://github.com/MolCrafts/molpy.git
cd molpy

# Install development dependencies
pip install -e ".[dev]"

# Install pre-commit hooks (automatically formats code)
pre-commit install

# Run tests
pytest tests/

# Manual formatting (if needed)
black src/
isort src/
```

### Code Quality

This project uses pre-commit hooks to ensure code quality:

- **Black**: Automatic code formatting (line length: 88)
- **isort**: Import sorting and organization
- **Pre-commit hooks**: Run automatically on every commit

See [DEVELOPMENT.md](DEVELOPMENT.md) for detailed development guidelines.

## 📄 License

This project is licensed under the BSD-3-Clause License - see the [LICENSE](LICENSE) file for details.

---

**Built with ❤️ by the MolCrafts team**
