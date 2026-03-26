[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/01_parsing_chemistry.ipynb)

# Parsing Chemistry

After reading this page you will be able to convert SMILES, SMARTS, BigSMILES, and CGSmiles strings into MolPy structures, and know which parser to use for each task.

## Four notations, four purposes

MolPy's parser module handles a stack of chemical notations. Each one answers a different question:

- **SMILES** — what is this exact molecule?
- **SMARTS** — what pattern should I match?
- **BigSMILES** — what is this monomer, and where does it connect?
- **CGSmiles** — how are building blocks arranged into an architecture?

The parser module lives at `mp.parser`. Convenience functions return `Atomistic` objects directly; IR-level functions return intermediate representations for advanced inspection.


## SMILES: one molecule, one structure

`parse_molecule` converts a SMILES string into an `Atomistic` object with atoms and bonds. Hydrogens are implicit in SMILES and can be added later during coordinate generation.

```python
import molpy as mp

mol = mp.parser.parse_molecule("CC(=O)OCC")   # ethyl acetate
print(f"atoms: {len(mol.atoms)}, bonds: {len(mol.bonds)}")

symbols = [atom.get("symbol") for atom in mol.atoms]
print(symbols)   # ['C', 'C', 'O', 'O', 'C', 'C']
```

For dot-separated mixtures like `[Li+].[F-]`, use `parse_mixture` instead. It always returns a list.

```python
ions = mp.parser.parse_mixture("[Li+].[F-]")
print(len(ions))   # 2

# parse_mixture also works for single molecules
mols = mp.parser.parse_mixture("CCO")
print(len(mols))   # 1
```

Aromatic atoms are lowercase in SMILES. Benzene is `c1ccccc1`, where the matching digits close the ring.

```python
benzene = mp.parser.parse_molecule("c1ccccc1")
for atom in benzene.atoms:
    print(f"{atom.get('symbol')}, aromatic={atom.get('aromatic')}")
```

Stereochemical information defined in SMILES/SMARTS (e.g., tetrahedral and double-bond stereochemistry) is parsed and preserved during topology construction when explicitly specified.


## SMARTS: pattern matching, not structure building

SMILES describes one concrete molecule. SMARTS describes a matching rule that can apply to many molecules. The parser produces a query graph with logical constraints on each atom, not a physical structure.

```python
query = mp.parser.parse_smarts("[C;X4][O;H1]")

print(f"query atoms: {len(query.atoms)}")
print(f"query bonds: {len(query.bonds)}")

for i, atom in enumerate(query.atoms):
    print(f"  atom {i}: {atom.expression}")
```

SMARTS returns `SmartsIR`, not `Atomistic`. This distinction matters: you are building matching rules for typification, not molecular structures.


## BigSMILES: monomers with connection intent

Polymer modeling needs explicit connection points between repeat units. BigSMILES adds port markers (`<`, `>`, `$`) to standard SMILES.

`parse_monomer` produces an `Atomistic` with port metadata on the relevant atoms.

```python
monomer = mp.parser.parse_monomer("{[][<]CC(c1ccccc1)[>][]}")

print(f"atoms: {len(monomer.atoms)}")
ports = [a for a in monomer.atoms if a.get("port")]
print(f"ports: {len(ports)}")
for p in ports:
    print(f"  port '{p.get('port')}' on {p.get('symbol')}")
```

For multi-monomer specifications, `parse_polymer` preserves segment-level organization.

```python
spec = mp.parser.parse_polymer("{[<]CC[>],[<]CC(C)[>]}")

print(f"topology:  {spec.topology}")
print(f"monomers:  {len(spec.all_monomers())}")
```


## CGSmiles: topology-level architecture

CGSmiles describes how labeled building blocks connect, without specifying their internal chemistry. Nodes are labeled blocks; edges are connections.

```python
from molpy.parser import parse_cgsmiles

# Linear chain: 5 copies of PEO with a fragment definition
cg = parse_cgsmiles("{[#PEO]|5}.{#PEO=[$]COC[$]}")

print(f"nodes: {len(cg.base_graph.nodes)}")
print(f"bonds: {len(cg.base_graph.bonds)}")
print(f"fragments: {len(cg.fragments)}")
```

Key CGSmiles syntax:

- `[#LABEL]` — one node referring to a monomer label
- `|n` — repeat operator
- `(...)` — branch
- matched digits — ring closure

CGSmiles is the input to `PolymerBuilder`. BigSMILES defines *what* each block is; CGSmiles defines *how* blocks are arranged.


## Two-step parsing for advanced inspection

The convenience functions (`parse_molecule`, `parse_monomer`) do parsing and conversion in one call. For diagnostics, you can split these into two steps: parse to IR, then convert.

```python
from molpy.parser import parse_smiles, smilesir_to_atomistic

ir = parse_smiles("CCO")
print(f"IR atoms: {len(ir.atoms)}")

for atom_ir in ir.atoms:
    print(f"  element={atom_ir.element}, aromatic={atom_ir.aromatic}")

mol = smilesir_to_atomistic(ir)
print(f"Atomistic atoms: {len(mol.atoms)}")
```

The same pattern works for BigSMILES:

```python
from molpy.parser import parse_bigsmiles, bigsmilesir_to_polymerspec

ir = parse_bigsmiles("{[<]CC[>],[<]CC(C)[>]}")
spec = bigsmilesir_to_polymerspec(ir)
print(f"topology: {spec.topology}")
```


## Choosing the right parser

| Notation | Function | Returns | Use when |
|----------|----------|---------|----------|
| SMILES | `parse_molecule(s)` | `Atomistic` | One specific molecule |
| SMILES (mixture) | `parse_mixture(s)` | `list[Atomistic]` | Dot-separated components |
| BigSMILES | `parse_monomer(s)` | `Atomistic` (with ports) | One repeat unit |
| BigSMILES | `parse_polymer(s)` | `PolymerSpec` | Multi-monomer specification |
| SMARTS | `parse_smarts(s)` | `SmartsIR` | Matching rules for typification |
| CGSmiles | `parse_cgsmiles(s)` | `CGSmilesIR` | Topology-level architecture |

For 3D coordinates after parsing, use `mp.tool.generate_3d(mol)` (requires RDKit).

See also: [Stepwise Polymer Construction](02_polymer_stepwise.md), [Atomistic and Topology](../tutorials/01_atomistic_and_topology.md).
