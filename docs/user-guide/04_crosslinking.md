[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/04_crosslinking.ipynb)

# Crosslinked Networks

After reading this page you will be able to generate `fix bond/react` templates, pack a crosslinkable monomer mixture, and export all files LAMMPS needs for reactive MD.

!!! note "Prerequisites"
    This guide requires RDKit, Packmol, and the `oplsaa.xml` force field. Familiarity with [Stepwise Polymer Construction](02_polymer_stepwise.md) is assumed.

## Workflow overview

A crosslinking simulation in LAMMPS needs four deliverables:

1. **Reaction templates** — `pre.mol`, `post.mol`, and `.map` files for `fix bond/react`
2. **Packed configuration** — `.data` file with non-overlapping monomer placement
3. **Force field coefficients** — `.ff` file covering all types in the system and templates
4. **LAMMPS input script** — wiring templates, data, and force field together

MolPy provides `TemplateReacter` to generate the template files, and `Molpack` to handle packing.


## Build monomers

We use two monomers: EO2 (linear, two `$` ports) and EO3 (branched, three `$` ports). Both share identical arm chemistry (`...COCCO[$]`), which means a single reaction template covers all pairings.

```python
import molpy as mp
from molpy.core.atomistic import Atom, Atomistic
from molpy.typifier import OplsTypifier

ff = mp.io.read_xml_forcefield("oplsaa.xml")
typifier = OplsTypifier(ff, strict_typing=False)

def build_monomer(bigsmiles):
    monomer = mp.parser.parse_monomer(bigsmiles)
    monomer = mp.tool.generate_3d(monomer, add_hydrogens=True, optimize=True)
    monomer.get_topo(gen_angle=True, gen_dihe=True)
    for idx, atom in enumerate(monomer.atoms, start=1):
        atom["id"] = idx
    typifier.typify(monomer)
    return monomer

eo2 = build_monomer("{[][$]OCCOCCOCCO[$][]}")
eo3 = build_monomer("{[]C(COCCO[$])(COCCO[$])COCCO[$][]}")

for label, mon in [("EO2", eo2), ("EO3", eo3)]:
    ports = [a.get("port") for a in mon.atoms if a.get("port")]
    print(f"{label}: atoms={len(mon.atoms)}, ports={ports}")
```


## Generate fix bond/react templates

`TemplateReacter` extends `Reacter` with subgraph extraction. It captures the local environment around each reaction site (controlled by `radius`) and writes pre-reaction and post-reaction molecule files plus an atom map.

With `radius=3`, the BFS stops before reaching the EO3 branching carbon. This makes the extracted subgraph identical for EO2+EO2, EO2+EO3, and EO3+EO3 pairings — one template covers all three.

The `select_hydroxyl_group` function is the same dehydration selector used in [Stepwise Polymer Construction](02_polymer_stepwise.md).

```python
from pathlib import Path
from molpy.reacter import (
    find_neighbors, find_port, form_single_bond,
    select_hydrogens, select_neighbor, select_self,
)
from molpy.reacter.template import TemplateReacter

output_dir = Path("crosslinking_output")
output_dir.mkdir(exist_ok=True)

def select_hydroxyl_group(struct: Atomistic, site: Atom) -> list[Atom]:
    """Remove the -OH leaving group. Same as in the stepwise guide."""
    for nb in find_neighbors(struct, site):
        if nb.get("symbol") != "O":
            continue
        hs = [a for a in find_neighbors(struct, nb, element="H")]
        if hs:
            return [nb, hs[0]]
    raise ValueError("No hydroxyl group found")

template_reacter = TemplateReacter(
    name="rxn1",
    site_selector_left=select_neighbor("C"),
    site_selector_right=select_self,
    leaving_selector_left=select_hydroxyl_group,
    leaving_selector_right=select_hydrogens(1),
    bond_former=form_single_bond,
    radius=3,
)

# Use EO2+EO2 as representative pair
left = eo2.copy()
right = eo2.copy()
_result, template = template_reacter.run(
    left=left, right=right,
    port_atom_L=find_port(left, "$"),
    port_atom_R=find_port(right, "$"),
    compute_topology=True,
)

template.write(base_path=output_dir / "rxn1", typifier=typifier)
print(f"pre atoms: {len(template.pre.atoms)}")
print(f"post atoms: {len(template.post.atoms)}")
```


## Pack the monomer mixture

Compute box size from total mass and target density, then use Packmol for non-overlapping placement.

```python
import numpy as np
from molpy.core.element import Element
from molpy.pack import InsideBoxConstraint, Molpack

n_eo2, n_eo3 = 27, 9
target_density = 1.2  # g/cm^3

def molecular_weight(monomer):
    return sum(Element(a.get("symbol")).mass for a in monomer.atoms)

total_mass_g = (n_eo2 * molecular_weight(eo2) + n_eo3 * molecular_weight(eo3)) / 6.022e23
box_length = ((total_mass_g / target_density) * 1e24) ** (1/3)

eo2_frame = eo2.to_frame()
eo3_frame = eo3.to_frame()

packer = Molpack(workdir=output_dir / "packmol")
constraint = InsideBoxConstraint(
    length=np.array([box_length] * 3),
    origin=np.zeros(3),
)
packer.add_target(eo2_frame, number=n_eo2, constraint=constraint)
packer.add_target(eo3_frame, number=n_eo3, constraint=constraint)

packed = packer.optimize(max_steps=10000, seed=42)
packed.metadata["box"] = mp.Box.cubic(length=box_length)

print(f"packed atoms: {packed['atoms'].nrows}")
print(f"box length: {box_length:.1f} A")
```


## Export LAMMPS files

The force field file must include types from both the packed configuration and the reaction templates, since `fix bond/react` may create atom types that do not exist in the initial state.

```python
from molpy.core.frame import Frame
from molpy.io.data.lammps import LammpsDataWriter
from molpy.io.forcefield.lammps import LAMMPSForceFieldWriter

# Collect all types across packed frame + template frames
all_frames = [packed, template.pre.to_frame(), template.post.to_frame()]
atom_types = set()
for frame in all_frames:
    if "atoms" in frame and "type" in frame["atoms"]:
        for t in frame["atoms"]["type"]:
            if t and not str(t).isdigit():
                atom_types.add(str(t))

LAMMPSForceFieldWriter(output_dir / "system.ff").write(
    ff, atom_types=atom_types,
)
LammpsDataWriter(output_dir / "system.data", atom_style="full").write(packed)
```


## LAMMPS input script

The script references the template files and runs minimization, equilibration, and reactive MD.

```python
script = """
units           real
atom_style      full
read_data       system.data
include         system.ff

minimize        1.0e-4 1.0e-4 1000 10000

molecule rxn1_pre rxn1_pre.mol
molecule rxn1_post rxn1_post.mol

fix rxns all bond/react stabilization yes npt_grp .03 &
  react rxn1 all 1 0.0 5 rxn1_pre rxn1_post rxn1.map

fix npt_grp_react all npt temp 300.0 300.0 100.0 iso 1.5 1.5 1000.0
run 5000

write_data final.data
"""

(output_dir / "run.in").write_text(script)
```


## Verifying the output

Before running the LAMMPS simulation, check that the generated files are internally consistent:

```python
# Verify template atom counts
from molpy.io.readers import read_lammps_molecule

pre_check = read_lammps_molecule(output_dir / "rxn1_pre.mol")
post_check = read_lammps_molecule(output_dir / "rxn1_post.mol")
print(f"pre template:  {pre_check['atoms'].nrows} atoms")
print(f"post template: {post_check['atoms'].nrows} atoms")

# The post template should have fewer atoms (leaving groups removed)
assert post_check["atoms"].nrows < pre_check["atoms"].nrows

# Verify the map file exists and is non-empty
map_file = output_dir / "rxn1.map"
assert map_file.exists() and map_file.stat().st_size > 0
print(f"map file: {map_file.stat().st_size} bytes")

# Verify data file has expected atom count
packed_check = mp.io.read_lammps_data(output_dir / "system.data", atom_style="full")
expected = n_eo2 * len(eo2.atoms) + n_eo3 * len(eo3.atoms)
print(f"packed: {packed_check['atoms'].nrows} atoms (expected ~{expected})")
```


## Troubleshooting

**Template generation fails**

Print the selected site and leaving-group atoms to see what the selectors found:

```python
site = find_port(left, "$")
carbon = [a for a in find_neighbors(left, site) if a.get("symbol") == "C"][0]
print(f"site: {carbon.get('symbol')} name={carbon.get('name')}")
for nb in find_neighbors(left, carbon):
    print(f"  neighbor: {nb.get('symbol')} name={nb.get('name')}")
```

**Packing fails**

Reduce the target density or monomer count. Packmol needs enough room to place molecules without overlap.

**Force field writing fails**

Inspect the collected type labels — if any type name is purely numeric or `None`, the typifier did not assign types to all atoms:

```python
for frame in all_frames:
    if "atoms" in frame and "type" in frame["atoms"]:
        types = frame["atoms"]["type"]
        missing = [i for i, t in enumerate(types) if not t or str(t).isdigit()]
        if missing:
            print(f"Missing types at indices: {missing}")
```

**LAMMPS rejects templates**

Check that the `.map` file contains valid equivalence IDs by reading the first few lines:

```python
print((output_dir / "rxn1.map").read_text()[:500])
```

See also: [Topology-Driven Assembly](03_polymer_topology.md), [Polydisperse Systems](05_polydisperse_systems.md).
