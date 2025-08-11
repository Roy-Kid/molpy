# API Reference

Complete API documentation for MolPy, organized by module functionality.

## Quick Navigation

::: molpy
    handler: python
    options:
      show_source: true
      separate_signature: true
      heading_level: 2

## Module Categories

### [Core Modules](core/index.md) - Fundamental Data Structures
- Frame & Block containers
- Selection system
- Topology management
- Force field system
- Atomistic structures
- System organization

### [I/O Modules](io/index.md) - File Format Support
- Data readers/writers (PDB, XYZ, LAMMPS, AMBER, GROMACS)
- Trajectory I/O with memory mapping
- Force field parameter files

### [Analysis Modules](analysis/index.md) - Property Calculations
- Diffusion analysis
- Density analysis
- Clustering algorithms
- Molecular properties

### [Builder Modules](builder/index.md) - Molecular Construction
- Polymer building
- Bulk system generation
- AmberTools integration
- Reaction building

### [Operations Modules](op/index.md) - Transformations
- Geometric operations
- Spatial transformations
- Topological manipulations

### [Potential Modules](potential/index.md) - Energy & Forces
- Bond potentials
- Angle potentials
- Pair potentials (LJ, Coulomb)
- Custom potential functions

### [Packing Modules](pack/index.md) - System Generation
- Molecular packing algorithms
- Packmol integration
- Optimization tools
- Constraint handling

### [Typifier Modules](typifier/index.md) - Automatic Typing
- SMARTS-based typing
- Graph-based algorithms
- Force field assignment

### [Engine Modules](engine/index.md) - Simulation Interfaces
- LAMMPS integration
- CP2K integration
- Engine abstraction layer

## Getting Started

For beginners, start with:
1. **[Core Modules](core/index.md)** - Understand basic data structures
2. **[I/O Modules](io/index.md)** - Learn how to read/write files
3. **[Analysis Modules](analysis/index.md)** - See what analysis tools are available

For advanced users:
- **[Builder Modules](builder/index.md)** - Create complex molecular systems
- **[Operations Modules](op/index.md)** - Transform and manipulate structures
- **[Potential Modules](potential/index.md)** - Define custom force fields
