[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/02_polymer_stepwise.ipynb)

# Stepwise Polymer Construction

After reading this page you will be able to construct a poly(ethylene oxide) (PEO) chain by repeating a single coupling reaction, and then reproduce the same result using MolPy's automated builder.

!!! note "Prerequisites"
    This guide requires RDKit (for `generate_3d`) and the `oplsaa.xml` force field file.


## What we are building

The monomer is **ethylene oxide (EO)**, represented in BigSMILES as `{[][<]OCCOCCOCCO[>][]}`. The `<` and `>` markers are reactive ports — they tell the builder where two monomers connect.

The polymerization follows **dehydration condensation**. At each step, the right port (`>`) of the growing chain reacts with the left port (`<`) of an incoming monomer. The reaction removes a hydroxyl group (OH) from one side and a hydrogen (H) from the other, forming a new C–O ether bond and releasing water.

```text
Monomer:    HO–[–O–CH₂–CH₂–O–CH₂–CH₂–O–CH₂–CH₂–O–]–H
            port <                                port >

Step 1 (dimer):

    HO–[EO]–OH  +  HO–[EO]–H
           ↑ port >    ↑ port <
           └──react──→ new C–O bond, remove OH + H (= H₂O)
    ─────────────────────────────
    Product:  HO–[EO]–O–[EO]–H
              port <         port >   (still available for next step)

Step 2 (trimer):

    HO–[EO]–O–[EO]–OH  +  HO–[EO]–H
                    ↑ port >   ↑ port <
    ─────────────────────────────
    Product:  HO–[EO]–O–[EO]–O–[EO]–H

Step n:
    Repeat until the chain reaches the desired length.
```

The key insight is that the product of each step still has a `>` port on its right end, so the same reaction can be applied again. This is how a loop builds an n-mer from n monomers.


## Shared setup

All three paths below (manual, builder, facade) use the same monomer, force field, and reaction definition.

```python
import molpy as mp
from molpy.core.atomistic import Atom, Atomistic
from molpy.reacter import (
    Reacter, find_neighbors, find_port,
    form_single_bond, select_hydrogens, select_neighbor, select_self,
)
from molpy.typifier import OplsAtomisticTypifier

MONOMER_BIGSMILES = "{[][<]OCCOCCOCCO[>][]}"

def make_eo_monomer() -> mp.Atomistic:
    monomer = mp.parser.parse_monomer(MONOMER_BIGSMILES)
    monomer = mp.tool.generate_3d(monomer, add_hydrogens=True, optimize=True)
    monomer.get_topo(gen_angle=True, gen_dihe=True)
    return monomer
```

The reaction needs two custom decisions: what to select as the reaction site, and what to remove as the leaving group.

On the **left side** (the growing chain's `>` port), the site selector picks the carbon neighbor of the port atom, and the leaving selector removes the hydroxyl group (O + H) attached to that carbon.

On the **right side** (the incoming monomer's `<` port), the site selector uses the port atom itself, and the leaving selector removes one hydrogen.

```python
def select_hydroxyl_group(struct: Atomistic, reaction_site: Atom) -> list[Atom]:
    """Find and return [O, H] — the hydroxyl leaving group."""
    for neighbor in find_neighbors(struct, reaction_site):
        if neighbor.get("symbol") != "O":
            continue
        h_neighbors = [a for a in find_neighbors(struct, neighbor, element="H")]
        if h_neighbors:
            return [neighbor, h_neighbors[0]]
    raise ValueError("No hydroxyl group found near reaction site")

rxn = Reacter(
    name="dehydration",
    anchor_selector_left=select_neighbor("C"),
    anchor_selector_right=select_self,
    leaving_selector_left=select_hydroxyl_group,
    leaving_selector_right=select_hydrogens(1),
    bond_former=form_single_bond,
)

ff = mp.io.read_xml_forcefield("oplsaa.xml")
typifier = OplsAtomisticTypifier(ff, strict_typing=True)
```

Every selector follows the same signature: `(struct: Atomistic, atom: Atom) -> list[Atom]`. Site selectors must return exactly one atom; leaving selectors return zero or more atoms to remove.


## Path 1: manual reaction loop

The manual path makes every step visible. First, one coupling reaction produces a dimer. Then a loop repeats the same kernel to grow the chain to any length.

### One coupling step (monomer + monomer → dimer)

```python
mono_a = make_eo_monomer()
typifier.typify(mono_a)

mono_b = make_eo_monomer()
typifier.typify(mono_b)
mono_b.move([10.0, 0.0, 0.0])   # separate in space to avoid overlap

port_a = find_port(mono_a, ">")  # right end of chain
port_b = find_port(mono_b, "<")  # left end of incoming monomer

atoms_before = len(mono_a.atoms) + len(mono_b.atoms)

reaction = rxn.run(
    left=mono_a, right=mono_b,
    port_atom_L=port_a, port_atom_R=port_b,
    compute_topology=True,
)

dimer = reaction.product_info.product
dimer.get_topo(gen_angle=True, gen_dihe=True)
typifier.typify(dimer)

removed = len(reaction.topology_changes.removed_atoms)
assert atoms_before - len(dimer.atoms) == removed
print(f"dimer: {len(dimer.atoms)} atoms, removed {removed} (OH + H = water)")
```

The dimer still has a `>` port on its right end and a `<` port on its left end. This is what makes iteration possible.

### Growth loop (dimer → n-mer)

The loop applies the same reaction repeatedly. At each step, the growing chain's `>` port reacts with a fresh monomer's `<` port. Conservation checks inside the loop catch errors immediately at the step where they occur.

```python
def build_chain_manual(n_units, typifier, rxn):
    """Build a PEO chain of n_units by repeating the coupling reaction."""
    chain = make_eo_monomer()
    typifier.typify(chain)

    for i in range(1, n_units):
        # 1. Create and position a fresh monomer
        unit = make_eo_monomer()
        unit.move([10.0 * i, 0.0, 0.0])

        atoms_before = len(chain.atoms) + len(unit.atoms)

        # 2. React: chain's > port + monomer's < port
        reaction = rxn.run(
            left=chain, right=unit,
            port_atom_L=find_port(chain, ">"),
            port_atom_R=find_port(unit, "<"),
            compute_topology=True,
        )

        # 3. The product becomes the new chain
        chain = reaction.product_info.product
        chain.get_topo(gen_angle=True, gen_dihe=True)
        typifier.typify(chain)

        # 4. Conservation check
        removed = len(reaction.topology_changes.removed_atoms)
        assert atoms_before - len(chain.atoms) == removed, f"step {i} failed"

    return chain

pentamer = build_chain_manual(5, typifier, rxn)
print(f"5-mer: {len(pentamer.atoms)} atoms")
```

This is the complete manual workflow. Every decision — which port to use, what to remove, how to place monomers — is explicit in the code.


## Path 2: PolymerBuilder

The builder encodes the same decisions as reusable policies, so you don't repeat them in a loop.

| Manual step | Builder policy |
|-------------|---------------|
| `find_port(chain, ">")` + `find_port(unit, "<")` | `Connector(port_map={("EO","EO"): (">","<")})` |
| `unit.move([10*i, 0, 0])` | `Placer(separator=..., orienter=...)` |
| The `for` loop itself | `builder.build("{[#EO]|5}")` parses CGSmiles and traverses |

```python
from molpy.builder.polymer import (
    Connector, CovalentSeparator, LinearOrienter,
    Placer, PolymerBuilder,
)

eo_template = make_eo_monomer()

connector = Connector(
    port_map={("EO", "EO"): (">", "<")},
    reacter=rxn,
)
placer = Placer(
    separator=CovalentSeparator(buffer=-0.1),
    orienter=LinearOrienter(),
)

builder = PolymerBuilder(
    library={"EO": eo_template},
    connector=connector,
    placer=placer,
    typifier=typifier,
)

result = builder.build("{[#EO]|5}")
print(f"builder 5-mer: {len(result.polymer.atoms)} atoms, {result.total_steps} steps")
```

The CGSmiles string `{[#EO]|5}` means "5 copies of the EO monomer in a linear chain." The builder handles port lookup, monomer placement, and reaction execution internally.


## Path 3: polymer() tool function

When the reaction and placement policies are already defined, the `polymer()` tool function collapses everything into one call.

```python
from molpy.tool import polymer

chain = polymer(
    "{[#EO]|5}",
    library={"EO": eo_template},
    reaction_preset="dehydration",
    use_placer=True,
)
print(f"tool 5-mer: {len(chain.atoms)} atoms")
```

All three paths produce the same result — a 5-unit PEO chain with identical chemistry. The difference is how much control you retain over each step.


## Export

Save the typed polymer from the builder path to LAMMPS files:

```python
from pathlib import Path
import numpy as np

output_dir = Path("02_output")
output_dir.mkdir(exist_ok=True)

typed_polymer = typifier.typify(result.polymer)
frame = typed_polymer.to_frame()
atoms = frame["atoms"]
if "mol_id" not in atoms:
    atoms["mol_id"] = np.ones(atoms.nrows, dtype=int)

coords = np.column_stack([atoms["x"], atoms["y"], atoms["z"]])
lo = coords.min(axis=0) - 5.0
hi = coords.max(axis=0) + 5.0
frame.box = mp.Box(matrix=hi - lo, origin=lo)

mp.io.write_lammps_system(output_dir / "peo5", frame, ff)

print(f"exported to {output_dir}")
```


## Troubleshooting

**Port not found**

Inspect the monomer's port markers to verify they exist and are on the expected atoms:

```python
for a in monomer.atoms:
    if a.get("port"):
        print(f"  {a.get('symbol')} port={a.get('port')} name={a.get('name')}")
```

If no output appears, the BigSMILES string is missing port markers (`<`, `>`) or the parser did not assign them.

**Reaction fails with "No hydroxyl group found"**

Print the neighbors of the reaction site to see what the selector is working with:

```python
site = find_port(monomer, ">")
carbon = [a for a in find_neighbors(monomer, site) if a.get("symbol") == "C"][0]
for nb in find_neighbors(monomer, carbon):
    print(f"  {nb.get('symbol')} name={nb.get('name')}")
```

**Geometry overlap after building**

Increase the placement buffer. Negative values bring monomers closer together; positive values push them apart:

```python
placer = Placer(separator=CovalentSeparator(buffer=1.0), orienter=LinearOrienter())
```

**Atom count drift**

Add conservation assertions at each loop step as shown in the growth loop above. If the assertion fails, print `reaction.topology_changes` to see what was removed and what was added.

See also: [Topology-Driven Assembly](03_polymer_topology.md), [Crosslinked Networks](04_crosslinking.md).
