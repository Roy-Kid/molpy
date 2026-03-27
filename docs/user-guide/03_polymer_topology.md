[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/03_polymer_topology.ipynb)

# Topology-Driven Assembly with CGSmiles

After reading this page you will be able to translate CGSmiles expressions into linear, ring, and branched polymer products using a single builder configuration.

!!! note "Prerequisites"
    This guide requires RDKit (for `generate_3d`), the `oplsaa.xml` force field, and familiarity with [Stepwise Polymer Construction](02_polymer_stepwise.md).

## One builder, three topologies

In the stepwise guide, the reaction kernel was always the same — only the loop structure changed. CGSmiles pushes that idea further: the same builder configuration produces different architectures depending solely on the topology expression.

```text
linear:  {[#EO2]|4[#PS]}
ring:    {[#EO2]1[#PS][#EO2][#PS][#EO2]1}
branch:  {[#PS][#EO3]([#PS])[#PS]}
```

The builder does not change between these three products. Only the string changes.


## Step 1: inspect the topology graph

Always validate the CGSmiles parse before building. This catches label typos and structural mistakes cheaply.

```python
from molpy.parser import parse_cgsmiles

expressions = {
    "linear": "{[#EO2]|4[#PS]}",
    "ring":   "{[#EO2]1[#PS][#EO2][#PS][#EO2]1}",
    "branch": "{[#PS][#EO3]([#PS])[#PS]}",
}

for name, expr in expressions.items():
    ir = parse_cgsmiles(expr)
    labels = [node.label for node in ir.base_graph.nodes]
    print(f"{name}: nodes={len(ir.base_graph.nodes)}, labels={labels}")
```


## Step 2: build a typed monomer library

Every label in the CGSmiles expression must have a corresponding template in the library. All templates use `$` as the reactive port marker.

```python
import molpy as mp
from molpy.typifier import OplsAtomisticTypifier

ff = mp.io.read_xml_forcefield("oplsaa.xml")
typifier = OplsAtomisticTypifier(ff, strict_typing=True)

BIGSMILES = {
    "EO2": "{[][$]OCCO[$][]}",
    "EO3": "{[][$]OCC(CO[$])(CO[$])[]}",
    "PS":  "{[][$]OCC(c1ccccc1)CO[$][]}",
}

def build_monomer(bigsmiles, typifier):
    monomer = mp.parser.parse_monomer(bigsmiles)
    monomer = mp.tool.generate_3d(monomer, add_hydrogens=True, optimize=False)
    monomer.get_topo(gen_angle=True, gen_dihe=True)
    typifier.typify(monomer)
    return monomer

library = {label: build_monomer(bs, typifier) for label, bs in BIGSMILES.items()}

for label, mon in library.items():
    ports = [a.get("port") for a in mon.atoms if a.get("port")]
    print(f"{label}: atoms={len(mon.atoms)}, ports={ports}")
```


## Step 3: define reaction and assembly policies

The reaction kernel is the same dehydration condensation from the [Stepwise Polymer Construction](02_polymer_stepwise.md) guide. If you have already worked through that page, this code will be familiar — the only difference is that the leaving selector on the right side is written as an inline function here instead of using `select_hydrogens(1)`.

??? note "Reaction setup (identical to the stepwise guide)"
    The `select_hydroxyl_group` function finds the -OH leaving group on the left monomer's reaction site. The `select_one_hydrogen` function picks one H from the right monomer's port atom. Together they implement dehydration condensation: -OH + -H are removed, forming a new C-O bond and releasing water.

```python
from molpy.core.atomistic import Atom, Atomistic
from molpy.builder.polymer import (
    Connector, CovalentSeparator, LinearOrienter,
    Placer, PolymerBuilder,
)
from molpy.reacter import (
    Reacter, find_neighbors, form_single_bond,
    select_neighbor, select_self,
)

def select_hydroxyl_group(struct: Atomistic, site: Atom) -> list[Atom]:
    for nb in find_neighbors(struct, site):
        if nb.get("symbol") != "O":
            continue
        hs = [a for a in find_neighbors(struct, nb, element="H")]
        if hs:
            return [nb, hs[0]]
    raise ValueError("No hydroxyl group found")

def select_one_hydrogen(struct: Atomistic, site: Atom) -> list[Atom]:
    hs = [a for a in find_neighbors(struct, site, element="H")]
    if not hs:
        raise ValueError("No hydrogen found")
    return [hs[0]]

rxn = Reacter(
    name="dehydration",
    anchor_selector_left=select_neighbor("C"),
    anchor_selector_right=select_self,
    leaving_selector_left=select_hydroxyl_group,
    leaving_selector_right=select_one_hydrogen,
    bond_former=form_single_bond,
)

rules = {(l, r): ("$", "$") for l in library for r in library}
connector = Connector(port_map=rules, reacter=rxn)
placer = Placer(
    separator=CovalentSeparator(buffer=-0.1),
    orienter=LinearOrienter(),
)

builder = PolymerBuilder(
    library=library,
    connector=connector,
    placer=placer,
    typifier=typifier,
)
```


## Step 4: build all three topologies

```python
for name, expr in expressions.items():
    result = builder.build(expr)
    polymer = result.polymer
    print(f"{name}: atoms={len(polymer.atoms)}, steps={result.total_steps}")
```

The same builder, the same reaction, the same library — only the CGSmiles string changes. This is the key advantage of topology-driven assembly: new architectures do not require new code.


## Export to LAMMPS

Each product can be exported using the same pattern.

```python
import numpy as np
from pathlib import Path
from molpy.io.writers import write_lammps_data

output_dir = Path("03_output")
output_dir.mkdir(exist_ok=True)

for name, expr in expressions.items():
    result = builder.build(expr)
    typed_polymer = typifier.typify(result.polymer)
    frame = typed_polymer.to_frame()
    atoms = frame["atoms"]
    if "q" not in atoms:
        atoms["q"] = np.zeros(atoms.nrows, dtype=float)
    if "mol" not in atoms:
        atoms["mol"] = np.ones(atoms.nrows, dtype=int)
    write_lammps_data(output_dir / f"{name}.data", frame)
```


## Troubleshooting

Debug in this order:

1. Parse CGSmiles and verify node/bond counts first
2. Confirm each label exists in the library
3. Confirm connector rules exist for each reacting label pair
4. Print selected site/leaving atoms if reaction fails
5. Tune `CovalentSeparator(buffer=...)` if geometry overlaps

See also: [Stepwise Polymer Construction](02_polymer_stepwise.md), [Crosslinked Networks](04_crosslinking.md).
