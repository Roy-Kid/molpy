[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/05_polydisperse_systems.ipynb)

# Polydisperse Systems

After reading this page you will be able to sample chain lengths from a target distribution, build atomistic chains, pack them into a periodic box, and export LAMMPS files.

!!! note "Prerequisites"
    This guide requires RDKit, Packmol, and the `oplsaa.xml` force field. Familiarity with [Stepwise Polymer Construction](02_polymer_stepwise.md) is assumed.

## From distribution to simulation box

Real polymer systems are polydisperse — chains have varying lengths. MolPy provides distribution samplers, a system planner, and packing tools to go from a target molecular weight distribution to a simulation-ready box in one continuous workflow.

The three summary statistics that describe polydispersity are:

- **Mn** — number-average molecular weight
- **Mw** — weight-average molecular weight
- **PDI = Mw / Mn** — dispersity (1.0 = monodisperse)


## Step 1: define the distribution

MolPy provides four built-in distributions. Parameters are chosen so that mean chain length falls in a comparable range.

```python
from molpy.builder.polymer import (
    SchulzZimmPolydisperse,
    UniformPolydisperse,
    PoissonPolydisperse,
    FlorySchulzPolydisperse,
)

distributions = {
    "Schulz-Zimm": SchulzZimmPolydisperse(Mn=2000, Mw=2500, random_seed=42),
    "Uniform":     UniformPolydisperse(min_dp=10, max_dp=55, random_seed=42),
    "Poisson":     PoissonPolydisperse(lambda_param=32, random_seed=42),
    "Flory-Schulz": FlorySchulzPolydisperse(a=0.06, random_seed=42),
}
```


## Step 2: sample chain ensembles

Sampling uses a three-level architecture:

1. `WeightedSequenceGenerator` — monomer identity sequence
2. `PolydisperseChainGenerator` — DP/mass sampling from a distribution
3. `SystemPlanner` — accumulate chains until a target total mass is reached

```python
import numpy as np
import molpy as mp
from molpy.core.element import Element
from molpy.builder.polymer import (
    WeightedSequenceGenerator,
    PolydisperseChainGenerator,
    SystemPlanner,
)

# Compute monomer mass
eo_ref = mp.parser.parse_monomer("{[][<]OCCO[>][]}")
eo_ref = mp.tool.generate_3d(eo_ref, add_hydrogens=True, optimize=True)
monomer_mass = sum(Element(a.get("symbol")).mass for a in eo_ref.atoms)

seq_gen = WeightedSequenceGenerator(monomer_weights={"EO": 1.0})
target_total_mass = 5e5

results = {}
for name, dist in distributions.items():
    chain_gen = PolydisperseChainGenerator(
        seq_generator=seq_gen,
        monomer_mass={"EO": monomer_mass},
        end_group_mass=0.0,
        distribution=dist,
    )
    planner = SystemPlanner(
        chain_generator=chain_gen,
        target_total_mass=target_total_mass,
        max_rel_error=0.02,
    )
    plan = planner.plan_system(np.random.default_rng(42))
    results[name] = plan.chains

for name, chains in results.items():
    mw = np.array([c.mass for c in chains])
    Mn = float(np.mean(mw))
    Mw = float(np.sum(mw**2) / np.sum(mw))
    print(f"{name:15s}: {len(chains):4d} chains, Mn={Mn:.0f}, PDI={Mw/Mn:.3f}")
```


## Step 3: build atomistic chains

To keep runtime manageable, build a subset of the sampled chains. The statistical properties are preserved as long as the subset is reasonably large.

The reaction setup below is the same dehydration condensation used in [Stepwise Polymer Construction](02_polymer_stepwise.md) and [Topology-Driven Assembly](03_polymer_topology.md).

```python
from molpy.builder.polymer import (
    Connector, CovalentSeparator, LinearOrienter,
    Placer, PolymerBuilder,
)
from molpy.reacter import (
    Reacter, find_neighbors, form_single_bond,
    select_hydrogens, select_neighbor, select_self,
)
from molpy.core.atomistic import Atom, Atomistic
from molpy.typifier import OplsAtomisticTypifier

def select_hydroxyl_group(struct: Atomistic, site: Atom) -> list[Atom]:
    """Remove the -OH leaving group. Same as in the stepwise guide."""
    for nb in find_neighbors(struct, site):
        if nb.get("symbol") != "O":
            continue
        hs = [a for a in find_neighbors(struct, nb, element="H")]
        if hs:
            return [nb, hs[0]]
    raise ValueError("No hydroxyl group found")

rxn = Reacter(
    name="dehydration",
    anchor_selector_left=select_neighbor("C"),
    anchor_selector_right=select_self,
    leaving_selector_left=select_hydroxyl_group,
    leaving_selector_right=select_hydrogens(1),
    bond_former=form_single_bond,
)

eo_template = mp.parser.parse_monomer("{[][<]OCCO[>][]}")
eo_template = mp.tool.generate_3d(eo_template, add_hydrogens=True, optimize=True)
eo_template.get_topo(gen_angle=True, gen_dihe=True)

connector = Connector(port_map={("EO", "EO"): (">", "<")}, reacter=rxn)
placer = Placer(separator=CovalentSeparator(buffer=-0.1), orienter=LinearOrienter())
builder = PolymerBuilder(library={"EO": eo_template}, connector=connector, placer=placer)

n_build = 10
sz_chains = results["Schulz-Zimm"][:n_build]
atomistic_chains = []

for chain in sz_chains:
    cgsmiles = "{" + f"[#EO]|{chain.dp}" + "}"
    result = builder.build(cgsmiles)
    polymer = result.polymer
    polymer.get_topo(gen_angle=True, gen_dihe=True)
    atomistic_chains.append(polymer)

total_atoms = sum(len(c.atoms) for c in atomistic_chains)
print(f"built {len(atomistic_chains)} chains, total atoms: {total_atoms}")
```


## Step 4: typify and pack

Assign force field types, then pack chains into a periodic box with Packmol.

```python
ff = mp.io.read_xml_forcefield("oplsaa.xml")
typifier = OplsAtomisticTypifier(ff, strict_typing=False)

atomistic_chains = [typifier.typify(chain) for chain in atomistic_chains]

# Compute box size
total_mw = sum(
    sum(Element(a.get("symbol")).mass for a in c.atoms)
    for c in atomistic_chains
)
target_density = 0.05  # g/cm^3 (use ~1.0 for production)
volume = (total_mw / 6.022e23) / target_density * 1e24
box_length = volume ** (1/3)

from pathlib import Path
from molpy.pack import InsideBoxConstraint, Molpack

packer = Molpack(workdir=Path("05_output/packmol"))
constraint = InsideBoxConstraint(
    length=np.array([box_length] * 3),
    origin=np.zeros(3),
)
for chain in atomistic_chains:
    packer.add_target(chain.to_frame(), number=1, constraint=constraint)

packed = packer.optimize(max_steps=10000, seed=42)
packed.box = mp.Box.cubic(length=box_length)

print(f"packed: {packed['atoms'].nrows} atoms, box: {box_length:.1f} A")
```


## Step 5: export LAMMPS files

```python
atoms = packed["atoms"]
if "mol_id" not in atoms:
    atoms["mol_id"] = np.ones(atoms.nrows, dtype=int)

mp.io.write_lammps_system("05_output/lammps", packed, ff)
```


## GBigSMILES shortcut

If you already know the repeat unit, distribution, and target mass, encode them in one GBigSMILES string:

```python
# Pattern: {[<]monomer[>]}|distribution(params)|endgroup.|total_mass|
gbigsmiles = '{[<]OCCO[>]}|schulz_zimm(2000,2500)|[H].|5e5|'
ir = mp.parser.parse_gbigsmiles(gbigsmiles)
```


## Troubleshooting

| Step | Check |
|------|-------|
| Monomer mass wrong | Verify monomer has explicit hydrogens before mass calculation |
| SystemPlanner total mass off | Check `max_rel_error` setting |
| Chain topology missing | Call `get_topo(gen_angle=True, gen_dihe=True)` before typification |
| Packing fails | Lower target density or increase `max_steps` |

See also: [Topology-Driven Assembly](03_polymer_topology.md), [Crosslinked Networks](04_crosslinking.md).
