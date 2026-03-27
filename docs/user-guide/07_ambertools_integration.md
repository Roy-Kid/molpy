[![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/molcrafts/molpy/blob/master/docs/user-guide/07_ambertools_integration.ipynb)

# PEO-LiTFSI with AmberTools

After reading this page you will be able to parameterize molecules with the Amber workflow (antechamber, parmchk2, tleap), build polymer chains with `AmberPolymerBuilder`, and assemble a multi-component electrolyte system.

!!! warning "External dependencies"
    This guide requires **AmberTools** (via conda), **RDKit**, and **Packmol**. All three must be installed and accessible. Without AmberTools, no code on this page will run.

??? note "Setting up AmberTools"
    Install AmberTools in a dedicated conda environment:

    ```bash
    conda create -n AmberTools25 -c conda-forge ambertools=25
    conda activate AmberTools25
    # Verify installation
    which antechamber   # should print a path
    which tleap         # should print a path
    ```

    MolPy's wrapper classes activate the conda environment automatically when running commands, so you do not need to keep it active in your shell. The `env="AmberTools25"` parameter in the code below tells the wrapper which environment to activate.

    If you use a different environment name, replace `"AmberTools25"` throughout this guide.

## System overview

A PEO-LiTFSI polymer electrolyte consists of:

- **PEO chains** — poly(ethylene oxide) backbone
- **Li+** — lithium cation
- **TFSI-** — bis(trifluoromethylsulfonyl)imide anion

The workflow has five stages:

1. Create ion structures from SMILES and parameterize with Amber
2. Create PEO monomer with port markers
3. Build polymer chains with `AmberPolymerBuilder`
4. Pack all components into a periodic box
5. Export to LAMMPS


## Stage 1: parameterize TFSI with Amber

The Amber workflow for small molecules is: antechamber (assign types + charges) → parmchk2 (missing parameters) → tleap (topology + coordinates).

```python
from pathlib import Path
import molpy as mp
from molpy.adapter import RDKitAdapter
from molpy.tool import Generate3D
from molpy.io.writers import write_pdb
from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper

output_dir = Path("07_output")
ions_dir = output_dir / "ions"
ions_dir.mkdir(parents=True, exist_ok=True)

# Create TFSI from SMILES and generate 3D coordinates
tfsi = mp.parser.parse_molecule("O=S(=O)(C(F)(F)F)[N-]S(=O)(=O)C(F)(F)F")
adapter = RDKitAdapter(internal=tfsi)
adapter = Generate3D(
    add_hydrogens=False, embed=True, optimize=True, update_internal=True,
)(adapter)
tfsi = adapter.get_internal()

# Write PDB for antechamber input
write_pdb(ions_dir / "tfsi.pdb", tfsi.to_frame())
```

```python
conda_env = "AmberTools25"

# Step 1: antechamber — assign GAFF types and BCC charges
ac = AntechamberWrapper(name="antechamber", workdir=ions_dir,
                        env=conda_env, env_manager="conda")
ac.atomtype_assign(
    input_file=(ions_dir / "tfsi.pdb").absolute(),
    output_file=(ions_dir / "tfsi.mol2").absolute(),
    input_format="pdb", output_format="mol2",
    charge_method="bcc", atom_type="gaff2", net_charge=-1,
)

# Step 2: parmchk2 — generate missing parameters
parmchk2 = Parmchk2Wrapper(name="parmchk2", workdir=ions_dir,
                           env=conda_env, env_manager="conda")
parmchk2.run(args=["-i", "tfsi.mol2", "-o", "tfsi.frcmod", "-f", "mol2", "-s", "gaff2"])

# Step 3: tleap — generate prmtop and inpcrd
leap_script = """source leaprc.gaff2
TFSI = loadmol2 tfsi.mol2
loadamberparams tfsi.frcmod
saveamberparm TFSI tfsi.prmtop tfsi.inpcrd
quit
"""
(ions_dir / "tfsi_leap.in").write_text(leap_script)

tleap = TLeapWrapper(name="tleap", workdir=ions_dir,
                     env=conda_env, env_manager="conda")
tleap.run(args=["-f", "tfsi_leap.in"])
```


## Stage 2: create PEO monomer with port markers

Port markers in BigSMILES (`[>]` and `[<]`) define connection points. The monomer `{[][<]COC[>][]}` places port atoms on the CH₂ groups flanking the ether oxygen, ensuring the oxygen gets typed as ether (not hydroxyl).

```python
monomer = mp.parser.parse_monomer("{[][<]COC[>][]}")
adapter = RDKitAdapter(internal=monomer)
adapter = Generate3D(
    add_hydrogens=True, embed=True, optimize=True, update_internal=True,
)(adapter)
eo_monomer = adapter.get_internal()
```

End caps use single-port methane monomers.

```python
def make_cap(bigsmiles):
    cap = mp.parser.parse_monomer(bigsmiles)
    adapter = RDKitAdapter(internal=cap)
    adapter = Generate3D(
        add_hydrogens=True, embed=True, optimize=True, update_internal=True,
    )(adapter)
    return adapter.get_internal()

me_left = make_cap("{[][<]C[]}")
me_right = make_cap("{[]C[>][]}")
```


## Stage 3: build PEO with AmberPolymerBuilder

`AmberPolymerBuilder` wraps the monomer library, connector rules, and Amber tool chain (prepgen + tleap) into one builder that produces fully parameterized chains.

```python
from molpy.builder.polymer.ambertools import AmberPolymerBuilder, AmberPolymerBuilderConfig

config = AmberPolymerBuilderConfig(
    force_field="gaff2",
    charge_method="bcc",
    env="AmberTools25",
    env_manager="conda",
)

builder = AmberPolymerBuilder(
    library={"MeL": me_left, "EO": eo_monomer, "MeR": me_right},
    config=config,
)

result = builder.build("{[#MeL][#EO]|10[#MeR]}")
```

`AmberPolymerBuilder.build()` internally runs prepgen to create residue templates, then tleap to assemble the chain and produce `polymer.prmtop` and `polymer.inpcrd` in the output directory. These are standard Amber topology and coordinate files that Stage 4 reads back.


## Stage 4: assemble and pack

Read Amber topology files, merge force fields, add Li+ parameters, and pack with Packmol.

```python
import numpy as np
from molpy.io import read_amber
from molpy.core.frame import Frame
from molpy.pack import Molpack, InsideBoxConstraint

# Read parameterized components
tfsi_frame, tfsi_ff = read_amber(
    ions_dir / "tfsi.prmtop", ions_dir / "tfsi.inpcrd", frame=Frame(),
)
peo_frame, peo_ff = read_amber(
    output_dir / "polymer.prmtop", output_dir / "polymer.inpcrd", frame=Frame(),
)

# Merge force fields and add Li+ manually
combined_ff = peo_ff.merge(tfsi_ff)

# Create Li+ frame
li_frame = mp.parser.parse_molecule("[Li+]").to_frame()
li_frame["atoms"]["q"] = np.array([1.0])
li_frame["atoms"]["mass"] = np.array([6.94])
li_frame["atoms"]["type"] = np.array(["Li+"])

# Pack system
box_size = 60.0
packer = Molpack(workdir=output_dir)
constraint = InsideBoxConstraint(length=[box_size]*3, origin=[0.0]*3)
packer.add_target(peo_frame, number=3, constraint=constraint)
packer.add_target(li_frame, number=10, constraint=constraint)
packer.add_target(tfsi_frame, number=10, constraint=constraint)

system = packer.optimize(max_steps=20000, seed=12345)
system.box = mp.Box.cubic(box_size)
```


## Stage 5: export to LAMMPS

```python
from molpy.io.writers import write_lammps_data, write_lammps_forcefield

write_lammps_data(output_dir / "system.data", system, atom_style="full")
write_lammps_forcefield(output_dir / "system.ff", combined_ff)
```


## Troubleshooting

| Symptom | Check |
|---------|-------|
| Antechamber fails | Verify PDB has correct atom names and no duplicate IDs |
| TFSI charge wrong | Use `charge_method="bcc"` and verify `-nc -1` |
| Polymer build fails | Check port markers in monomer SMILES |
| Force field merge conflict | Inspect atom type names for collisions between PEO and TFSI |
| Packing fails | Increase box size or reduce molecule count |

See also: [Force Field Typification](06_typifier.md), [Wrapper and Adapter](../tutorials/07_wrapper_and_adapter.md).
