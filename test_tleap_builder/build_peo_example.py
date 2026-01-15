"""
Example: Build PEO polymer using AmberPolymerBuilder.

Demonstrates:
1. Loading monomer from PDB as Frame
2. Converting Frame to Atomistic with proper port markers
3. Building polymer from CGSmiles
4. Exporting to LAMMPS format
"""

from pathlib import Path

from molpy.builder.ambertools import (
    AmberPolymerBuilder,
    AmberPolymerBuilderConfig,
)
from molpy.core.atomistic import Atomistic, Atom, Bond
from molpy.core.frame import Frame
from molpy.io import read_pdb

WORK_DIR = Path(__file__).parent / "peo_build_output"
PEO_PDB = Path(__file__).parent / "PEO_initial.pdb"


def frame_to_atomistic(frame: Frame) -> Atomistic:
    """Convert Frame to Atomistic.
    
    Args:
        frame: Frame object with atoms and bonds blocks.
        
    Returns:
        Atomistic object with atoms and bonds.
    """
    atomistic = Atomistic()
    
    atoms_block = frame["atoms"]
    n_atoms = atoms_block.nrows
    
    # Create atoms - build index -> Atom mapping
    atom_map: dict[int, Atom] = {}
    for i in range(n_atoms):
        atom = Atom()
        for key in atoms_block.keys():
            val = atoms_block[key][i]
            # Convert numpy types to Python types
            if hasattr(val, 'item'):
                val = val.item()
            atom[key] = val
        atomistic.atoms.add(atom)
        atom_map[i] = atom
    
    # Create bonds if present
    if "bonds" in frame:
        bonds_block = frame["bonds"]
        for i in range(bonds_block.nrows):
            idx_i = int(bonds_block["atomi"][i])
            idx_j = int(bonds_block["atomj"][i])
            bond = Bond(atom_map[idx_i], atom_map[idx_j])
            atomistic.bonds.add(bond)
    
    return atomistic


def main():
    """Build PEO polymer and export to LAMMPS format."""
    
    print("=== Build PEO Polymer with AmberPolymerBuilder ===\n")
    
    # Step 1: Load monomer from PDB
    print("Step 1: Load monomer from PDB")
    frame = read_pdb(PEO_PDB)
    print(f"  Loaded Frame with {frame['atoms'].nrows} atoms, {frame['bonds'].nrows} bonds")
    
    # Step 2: Convert to Atomistic
    print("\nStep 2: Convert Frame to Atomistic")
    peo_monomer = frame_to_atomistic(frame)
    atoms = list(peo_monomer.atoms)
    print(f"  Created Atomistic with {len(atoms)} atoms, {len(list(peo_monomer.bonds))} bonds")
    
    # Step 3: Mark ports
    # PEO_initial.pdb: atom 1 and atom 2 are terminal carbons
    # (1-indexed in PDB, so atoms[0] and atoms[1] in our 0-indexed list)
    print("\nStep 3: Mark ports on monomer")
    
    head_atom = atoms[0]  # C1 - terminal carbon
    tail_atom = atoms[1]  # C2 - terminal carbon
    
    head_atom["port"] = "<"
    tail_atom["port"] = ">"
    
    print(f"  Head port (<): atom {head_atom['name']} (element: {head_atom['element']})")
    print(f"  Tail port (>): atom {tail_atom['name']} (element: {tail_atom['element']})")
    
    # Step 4: Create builder
    print("\nStep 4: Create AmberPolymerBuilder")
    
    config = AmberPolymerBuilderConfig(
        force_field="gaff2",
        charge_method="bcc",
        work_dir=WORK_DIR / "amber_work",
        keep_intermediates=True,
        env="AmberTools25",
        env_manager="conda",
    )
    
    builder = AmberPolymerBuilder(
        library={"EO": peo_monomer},
        config=config,
    )
    print(f"  Work directory: {config.work_dir}")
    print(f"  Force field: {config.force_field}")
    
    # Step 5: Build polymer
    print("\nStep 5: Build polymer from CGSmiles")
    
    cgsmiles = "{[#EO]|5}"
    print(f"  CGSmiles: {cgsmiles}")
    
    result = builder.build(cgsmiles)
    
    print(f"\n  ✓ Built polymer with {result.monomer_count} monomers")
    print(f"  AMBER topology: {result.prmtop_path}")
    print(f"  AMBER coords: {result.inpcrd_path}")
    if result.pdb_path:
        print(f"  PDB output: {result.pdb_path}")
    
    # Step 6: Export to LAMMPS
    print("\nStep 6: Export to LAMMPS format")
    
    from molpy.io.writers import write_lammps_data, write_lammps_forcefield
    
    lammps_dir = WORK_DIR / "lammps"
    lammps_dir.mkdir(parents=True, exist_ok=True)
    
    lammps_data = lammps_dir / "peo.data"
    lammps_ff = lammps_dir / "peo.ff"
    
    write_lammps_data(lammps_data, result.frame, atom_style="full")
    print(f"  ✓ LAMMPS data: {lammps_data}")
    
    write_lammps_forcefield(lammps_ff, result.forcefield)
    print(f"  ✓ LAMMPS forcefield: {lammps_ff}")
    
    print("\n✓ Done!")


if __name__ == "__main__":
    main()
