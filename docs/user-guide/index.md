# Guides

These pages solve concrete tasks from input to output. Each guide assumes you understand the [core concepts](../tutorials/index.md).

## Foundations

- [Tool Layer](tools.md) — packaged recipes (`PrepareMonomer`, `polymer`, `polymer_system`) for common multi-step jobs
- [I/O](io.md) — reading, writing, and extending file formats for data, trajectories, and force fields
- [Parsing Chemistry](01_parsing_chemistry.md) — SMILES, SMARTS, BigSMILES, CGSmiles notation into MolPy structures

## Polymer Workflows

- [Stepwise Polymer Construction](02_polymer_stepwise.md) — manual reaction, PolymerBuilder, and facade paths
- [Topology-Driven Assembly](03_polymer_topology.md) — CGSmiles expressions for linear, ring, and branched architectures
- [Crosslinked Networks](04_crosslinking.md) — template generation for LAMMPS `fix bond/react`
- [Polydisperse Systems](05_polydisperse_systems.md) — distribution sampling, chain building, and packing

## Parameterization and Integration

- [Force Field Typification](06_typifier.md) — SMARTS-based atom typing and parameter assignment
- [PEO-LiTFSI with AmberTools](07_ambertools_integration.md) — external tool integration for polymer electrolytes
