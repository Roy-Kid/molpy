# Moltemplate CLI

MolPy ships a native moltemplate engine. The `molpy moltemplate` subcommand parses `.lt` scripts and emits ready-to-run inputs for **LAMMPS**, **OpenMM**, **GROMACS** — or MolPy's canonical XML force-field format.

## Quick start

```bash
# Summarise a script (atom types, molecule count, styles)
molpy moltemplate info water.lt

# Generate LAMMPS inputs (data + in.settings + in.init + starter in)
molpy moltemplate run water.lt --emit lammps --out-dir out/

# Generate every engine at once
molpy moltemplate run water.lt --emit all --out-dir out/

# FF-only: convert .lt force field to MolPy XML
molpy moltemplate convert gaff2.lt gaff2.xml

# Dump parsed IR (debug)
molpy moltemplate parse water.lt --json ir.json
```

## Subcommands

### `run` — emit engine inputs

```
molpy moltemplate run SCRIPT [--emit ENGINE ...] [--out-dir DIR] [--prefix NAME]
```

| Engine    | Files produced (given `--prefix system`)                               |
|-----------|------------------------------------------------------------------------|
| `lammps`  | `system.data`, `system.in.settings`, `system.in.init`, `system.in`     |
| `openmm`  | `system.xml`, `system.pdb`, `system.py`                                |
| `gromacs` | `system.gro`, `system.top`, `em.mdp`, `nvt.mdp`                        |
| `xml`     | `system.xml`, `system.pdb` (MolPy canonical)                           |
| `all`     | every engine above                                                     |

`--emit` may be repeated: `--emit lammps --emit openmm`.

### `parse` — debug the IR

```
molpy moltemplate parse SCRIPT [--json OUT]
```

Without `--json`, prints a one-line-per-kind statement summary.

### `info` — one-liner summary

```
molpy moltemplate info SCRIPT
```

Prints atom-type/atom/bond/angle/dihedral counts after all `new` instances are expanded.

### `convert` — FF-only to XML

```
molpy moltemplate convert SRC.lt DST.xml
```

## Supported moltemplate features

| Feature                                              | Status |
|------------------------------------------------------|--------|
| `ClassName { ... }`, nested classes                  | ✔      |
| `inherits Parent1, Parent2`                          | ✔      |
| `import "file.lt"` (recursive)                       | ✔      |
| `write("...")`, `write_once("...")`                  | ✔      |
| `Data Masses`, `Data Charges`, `In Charges`          | ✔      |
| `In Settings` coeff lines (pair/bond/angle/…)        | ✔      |
| `Data Atoms`, `Data Bonds`, `Data Angles`            | ✔      |
| `Data Dihedrals`, `Data Impropers`                   | ✔      |
| `inst = new Cls`                                     | ✔      |
| `.move(x,y,z)`, `.rot(θ,ax,ay,az)`, `.scale(s)`      | ✔      |
| `.rotvv(v1,v2)`                                      | partial (axis-chain only) |
| `Data Bonds By Type` auto-bond rules                 | planned |
| `[N]` array replication                              | planned |

## Python API

Everything the CLI does is available programmatically.

```python
from molpy.io.forcefield.moltemplate import read_moltemplate_system
from molpy.io.emit import emit, emit_all

atomistic, ff = read_moltemplate_system("water.lt")

# Single engine
emit("lammps", atomistic, ff, "out/", prefix="w")

# All engines
emit_all(atomistic, ff, "out/", prefix="w")
```

Core editing primitives on `Atomistic` / `CoarseGrain` / `ForceField` (see `core.atomistic`, `core.cg`, `core.forcefield`) include:

- `del_atom`, `del_bond`, `del_angle`, `del_dihedral`
- `rename_type(old, new, *, kind=Atom)`
- `set_property(selector, key, value, *, kind=Atom)`
- `select(predicate) -> Atomistic`
- `ForceField.rename_type / remove_type / remove_style`

These are kernel-level operations usable outside the moltemplate pipeline.
