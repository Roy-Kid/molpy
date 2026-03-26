# Concepts

MolPy represents a molecular system differently at each stage of work. These pages introduce the core data structures in the order you will encounter them: editable chemistry first, then system snapshots, then periodic geometry, force fields, trajectories, selections, and external tool integration.

Read them in order if you are new. Jump to a specific page if you already know what you need.


## How the pieces fit together

The diagram below shows the typical data flow through a MolPy workflow. Each box is a core data structure; each arrow is a transformation step.

```text
                    ┌─────────────────────────┐
  SMILES / file     │  Atomistic              │
  ────parser────>   │  (editable molecular    │
                    │   graph: atoms + bonds)  │
                    └───────────┬─────────────┘
                                │
                  typifier + ForceField
                                │
                    ┌───────────▼─────────────┐
                    │  Typed Atomistic         │
                    │  (atoms carry type,      │
                    │   charge, ff params)     │
                    └───────────┬─────────────┘
                                │
                          .to_frame()
                                │
                    ┌───────────▼─────────────┐
                    │  Frame                   │
                    │  (Block tables +         │
                    │   Box + metadata)        │
                    └───────────┬─────────────┘
                                │
                          io.write_*
                                │
                    ┌───────────▼─────────────┐
                    │  LAMMPS / GROMACS /      │
                    │  PDB / HDF5 files        │
                    └─────────────────────────┘
```

- **Atomistic** is where you edit chemistry — add atoms, remove bonds, run reactions, build polymers.
- **Frame** is where you do numerical work — vectorized distances, file I/O, engine export.
- **ForceField** is a separate data structure that travels alongside the system — it is not stored inside the molecule.
- **Box** defines the periodic cell and attaches to a Frame as metadata.
- **Trajectory** is a time-ordered sequence of Frames.


## Pages

- [Atomistic and Topology](01_atomistic_and_topology.md) — editable molecular graph and derived connectivity
- [Block and Frame](02_block_and_frame.md) — columnar tables and system snapshots
- [Box and Periodicity](03_box_and_periodicity.md) — simulation cells and minimum-image distances
- [Force Field](04_force_field.md) — parameter data, potentials, and multi-format export
- [Trajectory](05_trajectory.md) — time-ordered frame sequences with lazy access
- [Selector](06_selector.md) — composable atom filters over Block columns
- [Wrapper and Adapter](07_wrapper_and_adapter.md) — execution boundaries and representation boundaries
