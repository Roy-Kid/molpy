# Getting Started

This section takes you from installation to a working mental model in four steps. Read them in order.

1. **[Installation](installation.md)** — install MolPy and verify your environment
2. **[Quickstart](quickstart.md)** — build a water box, type it, export LAMMPS inputs
3. **[Core Concepts](core-concepts.md)** — understand the `Atomistic → Topology → Frame` pipeline
4. **[FAQ](faq.md)** — troubleshooting and comparisons with other tools


## Five-minute smoke test

Before diving into the full quickstart, verify that MolPy works by building a water molecule and writing a PDB file. This example requires no external dependencies — not even RDKit.

```python
import molpy as mp

# Build a water molecule from atoms and bonds
water = mp.Atomistic(name="water")
o  = water.def_atom(symbol="O", x=0.000, y=0.000, z=0.000)
h1 = water.def_atom(symbol="H", x=0.957, y=0.000, z=0.000)
h2 = water.def_atom(symbol="H", x=-0.239, y=0.927, z=0.000)
water.def_bond(o, h1)
water.def_bond(o, h2)

# Convert to a Frame and export
frame = water.to_frame()
mp.io.write_pdb("water.pdb", frame)
print(f"Wrote {frame['atoms'].nrows} atoms to water.pdb")
```

If this prints `Wrote 3 atoms to water.pdb`, your installation is working.


## Built-in data files

MolPy ships with commonly used force field files so you do not need to download them separately. When a guide or example references `"oplsaa.xml"`, the file is resolved from MolPy's built-in data directory automatically:

```python
# This loads the bundled OPLS-AA force field — no extra download needed
ff = mp.io.read_xml_forcefield("oplsaa.xml")
```

To see all available built-in files:

```python
from molpy.data import list_files
print(list(list_files("forcefield")))   # e.g. ['oplsaa.xml', 'tip3p.xml']
```


## Workflow at a glance

Most MolPy workflows follow the same pipeline:

```text
SMILES / file               Atomistic                   Frame
  input       ──parser──>   (editable graph)  ──────>   (columnar arrays)
                                  │                          │
                            typifier + ff              io.write_*
                                  │                          │
                            Typed Atomistic          LAMMPS / GROMACS
                                                     simulation files
```

1. **Parse or build** — create an `Atomistic` structure from SMILES, a file, or manually
2. **Edit** — add/remove atoms, run reactions, build polymers
3. **Typify** — assign force field types via SMARTS pattern matching
4. **Convert** — `atomistic.to_frame()` produces columnar arrays for numerical work
5. **Export** — write to LAMMPS, GROMACS, PDB, or other formats


After completing these pages, you should be able to answer:

- When should I edit a system as an `Atomistic` graph instead of as a `Frame`?
- Why is `Topology` derived from bonds instead of stored independently?
- Where do force field typing and export fit into the workflow?

If any of those still feel unclear, continue to [Concepts](../tutorials/index.md) for deeper explanations. If you want to jump straight to a concrete task, go to [Guides](../user-guide/index.md).
