# MolPy 🧬

**MolPy** is a modern, high-performance Python framework for molecular simulation and analysis. Built with flexibility and extensibility in mind, MolPy provides elegant data structures and comprehensive toolkits that let researchers focus on science rather than infrastructure.

## ✨ Key Features

- **🏗️ Modern Architecture**: Built on hierarchical data management with xarray-like interfaces
- **🔄 Unified Serialization**: Consistent `to_dict`/`from_dict` interface across all components
- **📦 Flexible Data Structures**: Support for atoms, bonds, angles, dihedrals, and complex molecular systems
- **🎯 Force Field Integration**: Comprehensive force field management and type assignment
- **📊 Trajectory Analysis**: Efficient handling of molecular dynamics trajectories
- **🧮 Advanced Algorithms**: Built-in optimization, packing, and reaction modeling
- **🔌 Extensible Design**: Plugin-ready architecture for custom functionality
- **⚡ High Performance**: Optimized for large-scale molecular systems

## 🚀 Quick Start

```python
import molpy as mp
import numpy as np

# Create a simulation box
box = mp.Box.cubic(10.0)

# Create static data
atoms_data = {
    'xyz': [[0.0, -0.06, 0.0], [0.75, 0.52, 0.0], [-0.75, 0.52, 0.0]],
    'element': ['O', 'H', 'H'],
    'type': [1, 2, 2]
}
h2o = mp.Frame(data={'atoms': atoms_data}, box=box)

# Create dynamic data
struct = mp.Atomistic()
struct.def_atom(xyz=[1.0, 0.0, 0.0], element='C', type=3)
struct.def_atom(xyz=[0.0, 0.0, 0.0], element='O', type=1)
struct.def_atom(xyz=[2.0, 0.0, 0.0], element='O', type=1)
struct.def_bond(i=0, j=1)

co2 = struct.to_frame()

# Save and load
mp.io.write_pdb("frame.pdb", co2)
frame = mp.io.read_pdb("frame.pdb")
```

## 📦 Installation

```bash
pip install molcrafts-molpy
```

For development and documentation:
```bash
pip install molcrafts-molpy[dev,doc]
```

## 🎯 What's Next?

- **[Getting Started](getting-started/installation.md)** - Set up your environment
- **[Quickstart Guide](getting-started/quickstart.md)** - Your first MolPy simulation
- **[Tutorials](tutorials/block-frame.md)** - Learn core concepts step-by-step
- **[API Reference](reference/index.md)** - Complete API documentation

---

**Built with ❤️ by the MolCrafts team**
