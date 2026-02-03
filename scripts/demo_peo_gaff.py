"""Demo: 3-mer PEO from SMILES, add H with RDKit, type with GAFF.

Workflow:
  1. Parse PEO trimer from SMILES: OCCOCCOCCO
  2. RDKit: add explicit H, generate 3D coordinates
  3. Generate topology (angles, dihedrals)
  4. Assign GAFF force field types
"""

from molpy.adapter import RDKitAdapter
from molpy.compute import Generate3D
from molpy.io.forcefield.xml import read_xml_forcefield
from molpy.parser.smiles import parse_smiles, smilesir_to_atomistic
from molpy.typifier.gaff import GaffAtomisticTypifier


def main():
    # 1. Build PEO trimer from SMILES: HO-CH2CH2-O-CH2CH2-O-CH2CH2-OH
    print("[1] Parsing SMILES ...")
    ir = parse_smiles("OCCOCCOCCO")
    peo = smilesir_to_atomistic(ir)
    print(f"    heavy atoms: {len(list(peo.atoms))}, bonds: {len(list(peo.bonds))}")

    # 2. Add hydrogens + 3D coordinates via RDKit
    print("[2] RDKit: add H, embed 3D, optimize ...")
    adapter = RDKitAdapter(internal=peo)
    adapter = Generate3D(
        add_hydrogens=True,
        embed=True,
        optimize=True,
        update_internal=True,
    )(adapter)
    peo = adapter.get_internal()
    print(f"    total atoms: {len(list(peo.atoms))}, bonds: {len(list(peo.bonds))}")

    # 3. Generate topology
    print("[3] Generating angles & dihedrals ...")
    peo.get_topo(gen_angle=True, gen_dihe=True)
    n_bonds = len(list(peo.bonds))
    n_angles = len(list(peo.angles))
    n_dihedrals = len(list(peo.dihedrals))
    print(f"    bonds: {n_bonds}, angles: {n_angles}, dihedrals: {n_dihedrals}")

    # 4. GAFF typing
    print("[4] GAFF force field typing ...")
    gaff_ff = read_xml_forcefield("gaff.xml")
    typifier = GaffAtomisticTypifier(gaff_ff, strict_typing=False)
    typifier.typify(peo)

    # 5. Report
    type_counts: dict[str, int] = {}
    untyped = 0
    for atom in peo.atoms:
        t = atom.get("type")
        if t is None:
            untyped += 1
        else:
            type_counts[t] = type_counts.get(t, 0) + 1

    print(f"\n    Atom types ({sum(type_counts.values())} typed, {untyped} untyped):")
    for t in sorted(type_counts):
        print(f"      {t:6s}  x{type_counts[t]}")

    bond_typed = sum(1 for b in peo.bonds if b.get("type") is not None)
    angle_typed = sum(1 for a in peo.angles if a.get("type") is not None)
    dihe_typed = sum(1 for d in peo.dihedrals if d.get("type") is not None)
    print(f"\n    Bonds:     {bond_typed}/{n_bonds}")
    print(f"    Angles:    {angle_typed}/{n_angles}")
    print(f"    Dihedrals: {dihe_typed}/{n_dihedrals}")

    print("\n    First 12 atoms:")
    for i, atom in enumerate(peo.atoms):
        if i >= 12:
            print("    ...")
            break
        sym = atom.get("symbol", "?")
        t = atom.get("type", "???")
        x, y, z = atom.get("x", 0.0), atom.get("y", 0.0), atom.get("z", 0.0)
        print(f"      {sym:2s}  type={t:6s}  ({x:7.3f}, {y:7.3f}, {z:7.3f})")

    print("\nDone.")


if __name__ == "__main__":
    main()
