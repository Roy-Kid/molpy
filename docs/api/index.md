# API Reference

## How to use this reference

Each page documents one MolPy package. Symbols are auto-generated from source docstrings. Use the tables below to locate the right page by task.

### Task → Symbol lookup

| I want to... | Symbol | Page |
|---|---|---|
| Build a molecule from atoms and bonds | `Atomistic`, `def_atom`, `def_bond` | [Core](core.md) |
| Store tabular molecular data | `Block`, `Frame` | [Core](core.md) |
| Define a periodic simulation cell | `Box` | [Core](core.md) |
| Handle time-ordered frame sequences | `Trajectory` | [Core](core.md) |
| Query bond graph (angles, paths, rings) | `Topology` | [Core](core.md) |
| Define force field parameters | `AtomisticForcefield`, `Style`, `Type` | [Core](core.md) |
| Parse SMILES / BigSMILES / SMARTS / CGSmiles | `parse_molecule`, `parse_monomer`, `parse_smarts`, `parse_cgsmiles` | [Parser](parser.md) |
| Run a chemical reaction (bond formation) | `Reacter`, `find_port`, `select_neighbor` | [Reacter](reacter.md) |
| Generate `fix bond/react` templates | `TemplateReacter` | [Reacter](reacter.md) |
| Assemble polymer chains from CGSmiles | `PolymerBuilder`, `Connector`, `Placer` | [Builder](builder.md) |
| Pack molecules into a periodic box | `Molpack`, `InsideBoxConstraint` | [Pack](pack.md) |
| Assign force field types via SMARTS matching | `OplsAtomisticTypifier`, `GaffTypifier` | [Typifier](typifier.md) |
| Compute bond/angle/pair energies and forces | `BondHarmonic`, `LJ126` | [Potential](potential.md) |
| Read / write molecular files (PDB, LAMMPS, GRO, ...) | `read_pdb`, `write_lammps_data`, `read_xml_forcefield` | [I/O](io.md) |
| Bridge to RDKit or OpenBabel objects | `RDKitAdapter`, `OpenBabelAdapter` | [Adapter](adapter.md) |
| Run external CLI tools (antechamber, tleap) | `Wrapper`, `AntechamberWrapper` | [Wrapper](wrapper.md) |
| Use packaged multi-step recipes | `PrepareMonomer`, `polymer`, `generate_3d` | [Tool](tool.md) |
| Compute MSD or displacement correlations | `MSD`, `DisplacementCorrelation` | [Tool](tool.md) |
| Run LAMMPS or CP2K simulations | `LAMMPSEngine`, `CP2KEngine` | [Engine](engine.md) |

### Package map

| Package | Responsibility |
|---------|---------------|
| [Core](core.md) | Data structures: `Atomistic`, `Frame`, `Block`, `Box`, `Trajectory`, `Topology`, `ForceField` |
| [Parser](parser.md) | Grammar-based parsing: SMILES, SMARTS, BigSMILES, CGSmiles, GBigSMILES |
| [Reacter](reacter.md) | Reaction framework: site/leaving selectors, bond formers, templates |
| [Builder](builder.md) | System assembly: polymer builders, connectors, placers |
| [Pack](pack.md) | Spatial packing via Packmol |
| [Typifier](typifier.md) | SMARTS-based atom typing (OPLS, GAFF) |
| [Potential](potential.md) | Numerical potential kernels (bond, angle, dihedral, pair) |
| [I/O](io.md) | File readers/writers for data, force fields, and trajectories |
| [Adapter](adapter.md) | In-memory bridges to RDKit, OpenBabel |
| [Wrapper](wrapper.md) | Subprocess wrappers for AmberTools CLIs |
| [Tool](tool.md) | High-level recipes and analysis operations |
| [Engine](engine.md) | MD engine abstractions (LAMMPS, CP2K) |
| [Optimization](optimize.md) | Potential wrappers for geometry optimization |
