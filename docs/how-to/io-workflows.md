# How to Work with Molecular File Formats

This guide shows practical examples of reading and writing molecular data in various formats.

## Reading Molecular Structures

### Loading a Protein from PDB

```python
import molpy as mp

# Load a protein structure
protein = mp.read_pdb("protein.pdb")
print(f"Protein has {protein['atoms'].nrows} atoms")

# Check elements present
elements = protein['atoms']['element']
unique_elements = set(elements)
for element in unique_elements:
    count = list(elements).count(element)
    print(f"{element}: {count} atoms")
```

### Converting PDB to LAMMPS Format

```python
# Load PDB and create system
protein = mp.read_pdb("protein.pdb")
system = mp.FrameSystem(frame=protein)

# Write to LAMMPS format
workdir = Path("lammps_output")
mp.write_lammps(workdir, system)
print(f"LAMMPS files written to {workdir}")
```

## Working with Force Fields

### Loading AMBER Force Fields

```python
# Load AMBER files
system, forcefield = mp.read_amber("protein.prmtop", "protein.inpcrd")
print(f"System: {system['atoms'].nrows} atoms")
print(f"Force field: {forcefield.name}")
print(f"Atom types: {forcefield.n_atomtypes}")
```

### Converting Force Fields

```python
# Convert AMBER to LAMMPS
system, amber_ff = mp.read_amber("protein.prmtop", "protein.inpcrd")
lammps_system = mp.FrameSystem(frame=system, forcefield=amber_ff)
mp.write_lammps(Path("converted"), lammps_system)
```

## Handling Trajectories

### Reading LAMMPS Trajectories

```python
# Read trajectory
traj_reader = mp.read_lammps_trajectory("trajectory.lammpstrj")

# Process frames
frames = []
for frame in traj_reader:
    frames.append(frame)
    if len(frames) >= 100:
        break

print(f"Loaded {len(frames)} frames")

# Write subset
mp.write_lammps_trajectory("subset.lammpstrj", frames[:50])
```

### Converting Trajectory Formats

```python
# Convert LAMMPS to PDB for visualization
traj_reader = mp.read_lammps_trajectory("trajectory.lammpstrj")
frames = list(traj_reader)

for i, frame in enumerate(frames[:10]):
    pdb_filename = f"frame_{i:04d}.pdb"
    mp.write_pdb(pdb_filename, frame)
```

## Advanced Workflows

### Building Complete Systems

```python
def build_complete_system():
    # Load protein
    protein = mp.read_pdb("protein.pdb")

    # Add force field
    system = mp.read_amber("protein.prmtop", "protein.inpcrd")

    # Add solvent if available
    try:
        solvent = mp.read_pdb("solvent.pdb")
        # Combine protein and solvent
        combined_atoms = mp.Block({
            'x': np.concatenate([protein['atoms']['x'], solvent['atoms']['x']]),
            'y': np.concatenate([protein['atoms']['y'], solvent['atoms']['y']]),
            'z': np.concatenate([protein['atoms']['z'], solvent['atoms']['z']]),
            'element': np.concatenate([protein['atoms']['element'], solvent['atoms']['element']])
        })
        complete_frame = mp.Frame({'atoms': combined_atoms})
        system = mp.FrameSystem(frame=complete_frame, forcefield=system.forcefield)
    except FileNotFoundError:
        print("No solvent file found")

    return system

# Build and export
complete_system = build_complete_system()
mp.write_pdb("complete.pdb", complete_system.frame)
mp.write_lammps(Path("complete"), complete_system)
```

### Batch Processing

```python
def batch_convert_pdb_to_lammps(input_dir, output_dir):
    input_path = Path(input_dir)
    output_path = Path(output_dir)
    output_path.mkdir(exist_ok=True)

    pdb_files = list(input_path.glob("*.pdb"))

    for pdb_file in pdb_files:
        try:
            structure = mp.read_pdb(pdb_file)
            system = mp.FrameSystem(frame=structure)
            output_subdir = output_path / pdb_file.stem
            mp.write_lammps(output_subdir, system)
            print(f"Converted {pdb_file.name}")
        except Exception as e:
            print(f"Error converting {pdb_file.name}: {e}")

# batch_convert_pdb_to_lammps("pdb_structures", "lammps_structures")
```

## Error Handling

### Safe File Reading

```python
def safe_read_molecular_file(file_path):
    file_path = Path(file_path)

    if not file_path.exists():
        raise FileNotFoundError(f"File not found: {file_path}")

    try:
        if file_path.suffix.lower() == ".pdb":
            return mp.read_pdb(file_path)
        elif file_path.suffix.lower() == ".gro":
            return mp.read_gro(file_path)
        elif file_path.suffix.lower() == ".mol2":
            return mp.read_mol2(file_path)
        else:
            raise ValueError(f"Unknown file type: {file_path.suffix}")
    except Exception as e:
        print(f"Error reading {file_path}: {e}")
        raise

# Usage
try:
    structure = safe_read_molecular_file("protein.pdb")
    print(f"Successfully loaded structure")
except Exception as e:
    print(f"Failed to load structure: {e}")
```

## Summary

This guide covered practical workflows for molecular file formats:

- Read/write different formats (PDB, LAMMPS, AMBER, etc.)
- Convert between simulation package formats
- Handle force fields and trajectories
- Build complete molecular systems
- Automate batch processing
- Handle errors safely

MolPy's IO system provides a unified interface for working with different file types and converting between them as needed.
