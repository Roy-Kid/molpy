# How to Use AmberTools API

This guide explains how to use MolPy's AmberTools API to automate molecular force field parameterization workflows.

## Overview

The AmberTools API provides modular components for each step of the AmberTools workflow:
- **Antechamber**: Atom type assignment and charge calculation
- **Prepgen**: Generate prepi files from AnteChamber output
- **Parmchk**: Generate frcmod files with missing force field parameters
- **TLeap**: Generate final parameter files and coordinates
- **AmberToolsTypifier**: Complete workflow orchestrator

All components inherit from `AmberToolsComponent` base class, providing consistent working directory management and conda environment configuration.

## Why This API Design?

The new API addresses several key problems with traditional AmberTools workflows:

1. **Repetitive Command Line Calls**: Instead of manually typing commands, you call Python methods
2. **File Management**: Automatic working directory creation and file organization
3. **Error Handling**: Built-in validation and error checking
4. **Workflow Composition**: Mix and match components as needed
5. **Integration**: Seamless integration with MolPy data structures

## Basic Setup

```python
import molpy as mp
from molpy.builder.ambertools import (
    Antechamber,
    Prepgen,
    Parmchk,
    TLeap,
    AmberToolsTypifier
)

# Create component instances
antech = Antechamber(conda_env="AmberTools25")
prepgen = Prepgen(conda_env="AmberTools25")
parmchk = Parmchk(conda_env="AmberTools25")
tleap = TLeap(conda_env="AmberTools25")
```

## Understanding the Workflow

The AmberTools workflow follows this sequence:

```
Molecular Structure → AnteChamber → PrepGen → ParmChk → TLeap → Final Files
```

Each step transforms the output of the previous step:

1. **AnteChamber** takes a molecular structure and assigns atom types and charges
2. **PrepGen** converts the .ac file to .prepi format for TLeap
3. **ParmChk** checks for missing force field parameters and generates .frcmod
4. **TLeap** combines everything into final .prmtop and .inpcrd files

## Component-by-Component Guide

### Antechamber Component

**What it does**: Automatically assigns atom types and calculates charges using quantum chemistry methods.

**When to use**: Always the first step for any new molecule that needs force field parameters.

```python
antech = Antechamber(conda_env="AmberTools25")

# Basic usage
ac_path, frame = antech(
    struct=struct,
    net_charge=0.0,           # Total system charge
    forcefield="gaff",        # Force field type
    charge_type="bcc",        # Charge calculation method
    output_format="ac"        # Output format
)

# With custom working directory
ac_path, frame = antech(
    struct=struct,
    net_charge=0.0,
    forcefield="gaff",
    workdir="/scratch/custom_work"  # Override default working directory
)
```

**Key Parameters Explained**:
- `net_charge`: Total charge of your system. Use 0.0 for neutral molecules
- `forcefield`: Choose from "gaff" (general), "gaff2" (improved), or "amber14" (protein-specific)
- `charge_type`: "bcc" (bond charge correction) is most common, "resp" for high accuracy
- `output_format`: "ac" is standard, but you can also get "mol2" or "pdb"
- `workdir`: Optional custom working directory (overrides the component's default)

**What you get**:
- `.ac` file: Contains atom types, charges, and connectivity
- `.mc` file: Monomer connectivity info (if input is a Monomer object)

**Additional Methods**:
```python
# Read an existing .ac file without running AnteChamber
frame = antech.read_ac("/path/to/existing/file.ac")
print(f"Read {len(frame['atoms'])} atoms from existing file")
```

### Prepgen Component

**What it does**: Converts AnteChamber's .ac output to .prepi format that TLeap can read.

**When to use**: After AnteChamber, before ParmChk.

```python
prepgen = Prepgen(conda_env="AmberTools25")

# Simple usage - just pass the .ac file path
prepi_path = prepgen(ac_path)
```

**What happens internally**:
1. Copies your .ac file to PrepGen's working directory
2. Runs `prepgen` command to convert format
3. Returns path to the new .prepi file

### Parmchk Component

**What it does**: Checks if all force field parameters exist and generates missing ones.

**When to use**: After PrepGen, before TLeap.

```python
parmchk = Parmchk(conda_env="AmberTools25")

# Pass the .prepi file from PrepGen
frcmod_path = parmchk(prepi_path)
```

**Why this step matters**: Even common molecules might have missing parameters for:
- Unusual bond lengths
- Non-standard dihedral angles
- Special atom types

### TLeap Component

**What it does**: Generates final Amber parameter files (.prmtop, .inpcrd) and PDB structure.

**When to use**: Final step to combine all components into simulation-ready files.

```python
tleap = TLeap(conda_env="AmberTools25")

# Combine multiple molecules
prmtop, inpcrd = tleap(
    name="system_name",
    prepi_files=[prepi_path1, prepi_path2],    # List of .prepi files
    frcmod_files=[frcmod_path1, frcmod_path2]  # List of .frcmod files
)
```

**What TLeap does automatically**:
1. Generates a TLeap script with all necessary commands
2. Sources the appropriate force field libraries
3. Loads all your molecule files
4. Combines molecules if you have multiple
5. Saves the final parameter files

## Complete Workflow Example

Let's walk through parameterizing ethanol step by step:

```python
import molpy as mp
from molpy.builder.ambertools import Antechamber, Prepgen, Parmchk, TLeap

# Step 1: Create molecular structure
ethanol = mp.Atomistic()
ethanol["name"] = "ethanol"

# Add atoms (simplified ethanol structure)
ethanol.def_atom(xyz=[0.0, 0.0, 0.0], element='C', name='C1')
ethanol.def_atom(xyz=[1.4, 0.0, 0.0], element='C', name='C2')
ethanol.def_atom(xyz=[2.1, 1.2, 0.0], element='O', name='O1')
ethanol.def_atom(xyz=[0.0, 1.0, 0.0], element='H', name='H1')
ethanol.def_atom(xyz=[0.0, -1.0, 0.0], element='H', name='H2')
ethanol.def_atom(xyz=[1.4, 1.0, 0.0], element='H', name='H3')
ethanol.def_atom(xyz=[1.4, -1.0, 0.0], element='H', name='H4')
ethanol.def_atom(xyz=[2.1, 1.2, 1.0], element='H', name='H5')

# Add bonds
ethanol.def_bond(i=0, j=1)  # C1-C2
ethanol.def_bond(i=1, j=2)  # C2-O1
ethanol.def_bond(i=0, j=3)  # C1-H1
ethanol.def_bond(i=0, j=4)  # C1-H2
ethanol.def_bond(i=1, j=5)  # C2-H3
ethanol.def_bond(i=1, j=6)  # C2-H4
ethanol.def_bond(i=2, j=7)  # O1-H5

# Step 2: Run AnteChamber
antech = Antechamber(conda_env="AmberTools25")

# Option 1: Use default working directory
ac_path, frame = antech(ethanol, net_charge=0.0, forcefield="gaff")
print(f"✓ AnteChamber completed: {ac_path}")

# Option 2: Specify custom working directory
ac_path, frame = antech(
    ethanol,
    net_charge=0.0,
    forcefield="gaff",
    workdir="/scratch/custom_work"
)
print(f"✓ AnteChamber completed in custom directory: {ac_path}")

# Option 3: Read existing .ac file without running AnteChamber
if ac_path.exists():
    frame = antech.read_ac(ac_path)
    print(f"✓ Read existing .ac file: {len(frame['atoms'])} atoms")

# Step 3: Generate prepi file
prepgen = Prepgen(conda_env="AmberTools25")
prepi_path = prepgen(ac_path)
print(f"✓ PrepGen completed: {prepi_path}")

# Step 4: Check parameters
parmchk = Parmchk(conda_env="AmberTools25")
frcmod_path = parmchk(prepi_path)
print(f"✓ ParmChk completed: {frcmod_path}")

# Step 5: Generate final files
tleap = TLeap(conda_env="AmberTools25")
prmtop, inpcrd = tleap("ethanol", [prepi_path], [frcmod_path])
print(f"✓ TLeap completed: {prmtop}, {inpcrd}")

# Read the final structure
frame, forcefield = mp.io.read_amber(prmtop, inpcrd)
print(f"✓ Final structure: {len(frame['atoms'])} atoms")
```

## Advanced Usage Patterns

### Custom Working Directories

By default, each component creates its own working directory. You can customize this:

```python
from pathlib import Path

# Specify custom working directories
antech = Antechamber(
    conda_env="AmberTools25",
    workdir=Path("/scratch/antechamber_work")
)

prepgen = Prepgen(
    conda_env="AmberTools25",
    workdir=Path("/scratch/prepgen_work")
)
```

### Batch Processing Multiple Molecules

```python
# Process multiple molecules efficiently
molecules = [mol1, mol2, mol3]
ac_paths = []
prepi_paths = []
frcmod_paths = []

for i, mol in enumerate(molecules):
    mol["name"] = f"molecule_{i}"

    # Run AnteChamber
    ac_path, _ = antech(mol, net_charge=0.0)
    ac_paths.append(ac_path)

    # Generate prepi file
    prepi_path = prepgen(ac_path)
    prepi_paths.append(prepi_path)

    # Generate frcmod file
    frcmod_path = parmchk(prepi_path)
    frcmod_paths.append(frcmod_path)

# Combine all molecules with TLeap
prmtop, inpcrd = tleap("combined_system", prepi_paths, frcmod_paths)
```

### Using AmberToolsTypifier for Complete Workflow

The `AmberToolsTypifier` class provides a convenient way to run the entire workflow:

```python
typifier = AmberToolsTypifier(conda_env="AmberTools25")

# Complete parameterization including TLeap
typed_struct = typifier(
    struct=struct,
    forcefield="gaff",
    charge_type="bcc",
    net_charge=0.0,
    is_frcmod=True,    # Generate frcmod file
    is_prepi=True,     # Generate prepi file
    is_tleap=True      # Run TLeap
)

# Only type assignment (no TLeap)
typed_struct = typifier(
    struct=struct,
    forcefield="gaff",
    is_tleap=False
)
```

## Troubleshooting Common Issues

### Conda Environment Problems

```bash
# Make sure AmberTools is installed
conda install -c conda-forge ambertools

# Or use mamba for faster installation
mamba install -c conda-forge ambertools

# Activate the environment
conda activate AmberTools25
```

### Command Not Found Errors

```bash
# Check if commands are available
which antechamber
which prepgen
which parmchk
which tleap

# If not found, check your PATH
echo $PATH
```

### Permission Issues

```bash
# Ensure working directory has write permissions
chmod 755 /path/to/workdir

# Check disk space
df -h
```

### Debugging Tips

```python
import logging

# Enable detailed logging
logging.basicConfig(level=logging.DEBUG)

# Check working directories
print(f"AnteChamber workdir: {antech.workdir}")
print(f"PrepGen workdir: {prepgen.workdir}")

# List generated files
import pathlib
for file in pathlib.Path(antech.workdir).glob("**/*"):
    print(f"  {file}")
```

## Performance Optimization

### Parallel Processing

```python
import concurrent.futures

def process_molecule(mol, i):
    mol["name"] = f"molecule_{i}"
    ac_path, _ = antech(mol, net_charge=0.0)
    prepi_path = prepgen(ac_path)
    frcmod_path = parmchk(prepi_path)
    return ac_path, prepi_path, frcmod_path

# Process multiple molecules in parallel
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as executor:
    futures = [executor.submit(process_molecule, mol, i)
               for i, mol in enumerate(molecules)]
    results = [future.result() for future in futures]
```

### Memory Management

```python
# Clean up intermediate files if needed
import shutil

# After successful completion, you can clean up
shutil.rmtree(antech.workdir)
shutil.rmtree(prepgen.workdir)
```

## Integration with Other MolPy Components

### Working with Polymers

```python
from molpy.builder.polymer import Monomer

# Create a monomer
monomer = Monomer()
monomer["name"] = "ethylene_glycol"
# ... add atoms and bonds ...

# Parameterize the monomer
typifier = AmberToolsTypifier(conda_env="AmberTools25")
typed_monomer = typifier(
    monomer,
    forcefield="gaff",
    charge_type="bcc",
    net_charge=0.0
)
```

### Reading Existing Structures

```python
# Read a PDB file and parameterize it
frame = mp.io.read_pdb("molecule.pdb")
struct = mp.Atomistic.from_frame(frame)
struct["name"] = "molecule"

# Parameterize
typed_struct = typifier(struct, forcefield="gaff")
```
