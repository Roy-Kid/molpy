"""Compare GAFF typing: MolPy implementation vs OpenBabel.

This script:
1. Builds a PEO trimer from SMILES
2. Adds hydrogens with RDKit and generates 3D coords
3. Types it with our GAFF implementation
4. Types it with OpenBabel's GAFF
5. Compares results
"""

from molpy.adapter import RDKitAdapter, OpenBabelAdapter
from molpy.compute import Generate3D
from molpy.io.forcefield.xml import read_xml_forcefield
from molpy.parser.smiles import parse_smiles, smilesir_to_atomistic
from molpy.typifier.gaff import GaffAtomisticTypifier


def main():
    print("=" * 70)
    print("  GAFF Typing Comparison: MolPy vs OpenBabel")
    print("=" * 70)

    # 1. Build PEO trimer
    print("\n[1] Building PEO trimer (SMILES: OCCOCCOCCO) ...")
    ir = parse_smiles("OCCOCCOCCO")
    peo = smilesir_to_atomistic(ir)

    # 2. Add H + 3D with RDKit
    print("[2] Adding hydrogens and 3D coordinates with RDKit ...")
    adapter_rdkit = RDKitAdapter(internal=peo)
    adapter_rdkit = Generate3D(
        add_hydrogens=True,
        embed=True,
        optimize=True,
        update_internal=True,
    )(adapter_rdkit)
    peo = adapter_rdkit.get_internal()
    print(f"    Total atoms: {len(list(peo.atoms))}")

    # 3. Generate topology
    print("[3] Generating topology ...")
    peo.get_topo(gen_angle=True, gen_dihe=True)

    # 4. Type with our GAFF implementation
    print("[4] Typing with MolPy GAFF implementation ...")
    gaff_ff = read_xml_forcefield("gaff.xml")
    typifier = GaffAtomisticTypifier(gaff_ff, strict_typing=False)
    typifier.typify(peo)

    # 5. Type with OpenBabel
    print("[5] Typing with OpenBabel GAFF ...")
    try:
        adapter_ob = OpenBabelAdapter(internal=peo)
        ob_comparison = adapter_ob.get_type_comparison()
    except Exception as e:
        print(f"    ERROR: {e}")
        print("    OpenBabel may not be installed or GAFF support missing.")
        ob_comparison = None

    # 6. Report
    print("\n" + "=" * 70)
    print("  TYPING RESULTS")
    print("=" * 70)

    # MolPy types
    molpy_types: dict[str, int] = {}
    for atom in peo.atoms:
        t = atom.get("type", "???")
        molpy_types[t] = molpy_types.get(t, 0) + 1

    print("\nMolPy GAFF types:")
    for t in sorted(molpy_types.keys()):
        print(f"  {t:6s}  x{molpy_types[t]}")

    if ob_comparison:
        print("\n" + "-" * 70)
        print(f"Comparison with OpenBabel:")
        print(f"  Total atoms:        {ob_comparison['total_atoms']}")
        print(f"  Agreement:          {ob_comparison['agreement']}")
        print(f"  Disagreement:       {ob_comparison['disagreement']}")
        print(f"  MolPy missing type: {ob_comparison['molpy_missing']}")
        print(f"  OpenBabel missing:  {ob_comparison['openbabel_missing']}")

        if ob_comparison["disagreement"] > 0:
            print("\n  Disagreements:")
            for detail in ob_comparison["details"]:
                if not detail["match"] and detail["match"] is not None:
                    print(
                        f"    Atom {detail['index']:2d} ({detail['symbol']}): "
                        f"MolPy={detail['molpy']:6s}, OpenBabel={detail['openbabel']:6s}"
                    )

        print("\n  Detailed comparison (first 20 atoms):")
        print(f"  {'Idx':3s} {'El':2s} {'MolPy':8s} {'OpenBabel':10s} {'Match':6s}")
        print(f"  {'-' * 50}")
        for detail in ob_comparison["details"][:20]:
            match_str = (
                "✓" if detail["match"] is True else
                "✗" if detail["match"] is False else
                "?"
            )
            print(
                f"  {detail['index']:3d} {detail['symbol']:2s} "
                f"{detail['molpy']:8s} {detail['openbabel']:10s} {match_str:6s}"
            )

        if len(ob_comparison["details"]) > 20:
            print(f"  ... ({len(ob_comparison['details']) - 20} more atoms)")

    print("\n" + "=" * 70)
    print("  Done.")
    print("=" * 70)


if __name__ == "__main__":
    main()
