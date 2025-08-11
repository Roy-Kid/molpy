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

# Create atoms data
atoms_data = {
    'x': [0.0, 1.0, 2.0],
    'y': [0.0, 0.0, 0.0],
    'z': [0.0, 0.0, 0.0],
    'element': ['C', 'C', 'C'],
    'type': [1, 1, 1]
}

# Create a frame
frame = mp.Frame(data={'atoms': atoms_data}, box=box)

# Save and load
frame_dict = frame.to_dict()  # Serialize to dictionary
restored_frame = mp.Frame.from_dict(frame_dict)  # Restore from dictionary
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
