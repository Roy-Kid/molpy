"""Example: Using the SMARTS matcher for atomtyping.

This example demonstrates:
1. Building patterns programmatically
2. Creating molecules
3. Finding and resolving atom type candidates
4. Using explain mode for debugging
"""

import json
from molpy.core.atomistic import Atomistic
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import (
    is_element, is_aromatic, charge, degree, bond_order
)
from molpy.typifier.matcher import SmartsMatcher
from molpy.typifier.adapter import build_mol_graph


def create_ethanol():
    """Create ethanol molecule: CH3-CH2-OH"""
    mol = Atomistic()
    
    # Atoms
    c1 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    c2 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    o = mol.add_atom(element="O", is_aromatic=False, charge=0)
    
    # Bonds
    mol.add_bond(c1, c2, order=1)
    mol.add_bond(c2, o, order=1)
    
    return mol, {"C1": c1, "C2": c2, "O": o}


def create_patterns():
    """Create atomtype patterns."""
    patterns = []
    
    # Generic carbon
    p1 = quick_pattern("CT", "C", priority=0)
    patterns.append(p1)
    
    # Hydroxyl oxygen
    pb2 = PatternBuilder("OH", priority=5, source="hydroxyl_pattern")
    v_c = pb2.add_vertex(preds=[is_element("C")])
    v_o = pb2.add_vertex(preds=[is_element("O")])
    pb2.add_edge(v_c, v_o, preds=[bond_order(1)])
    pb2.set_target_vertices([v_o])  # Only type the oxygen
    p2 = pb2.build()
    patterns.append(p2)
    
    # Alcohol carbon (C-O)
    pb3 = PatternBuilder("C_OH", priority=2, source="alcohol_pattern")
    v_c3 = pb3.add_vertex(preds=[is_element("C")])
    v_o3 = pb3.add_vertex(preds=[is_element("O")])
    pb3.add_edge(v_c3, v_o3, preds=[bond_order(1)])
    pb3.set_target_vertices([v_c3])  # Only type the carbon
    p3 = pb3.build()
    patterns.append(p3)
    
    return patterns


def main():
    print("=" * 70)
    print("SMARTS Matcher Example: Ethanol Atomtyping")
    print("=" * 70)
    
    # Create molecule
    print("\n1. Creating ethanol molecule (CH3-CH2-OH)...")
    mol, atoms = create_ethanol()
    print(f"   Created molecule with {len(mol.atoms)} atoms, {len(mol.bonds)} bonds")
    
    # Create patterns
    print("\n2. Creating atomtype patterns...")
    patterns = create_patterns()
    for i, p in enumerate(patterns):
        priority = p.get_priority() if hasattr(p, 'get_priority') else p._priority
        print(f"   P{i+1}: {p.atomtype_name:10s} (priority={priority}, "
              f"vertices={p.vcount()}, edges={p.ecount()})")
    
    # Build graph
    print("\n3. Converting molecule to graph...")
    graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)
    print(f"   Graph: {graph.vcount()} vertices, {graph.ecount()} edges")
    
    # Find candidates
    print("\n4. Finding atom type candidates...")
    matcher = SmartsMatcher(patterns)
    candidates = matcher.find_candidates(graph, vs_to_atomid)
    print(f"   Found {len(candidates)} candidate assignments")
    
    # Show candidates grouped by atom
    from collections import defaultdict
    by_atom = defaultdict(list)
    for c in candidates:
        by_atom[c.atom_id].append(c)
    
    print("\n   Candidates by atom:")
    atom_names = {id(v): k for k, v in atoms.items()}
    for atom_id, atom_candidates in by_atom.items():
        atom_name = atom_names.get(atom_id, f"atom_{atom_id}")
        print(f"   - {atom_name}:")
        for c in sorted(atom_candidates):
            print(f"      {c.atomtype:10s} (priority={c.priority}, score={c.score})")
    
    # Resolve conflicts
    print("\n5. Resolving conflicts...")
    result = matcher.resolve(candidates)
    
    print("\n   Final atom types:")
    for atom_name, atom in atoms.items():
        atomtype = result.get(id(atom), "UNTYPED")
        print(f"   - {atom_name}: {atomtype}")
    
    # Explain mode
    print("\n6. Explain mode (detailed debugging)...")
    explain = matcher.explain(candidates)
    
    print("\n   Detailed explanation:")
    for atom_name, atom in atoms.items():
        atom_id = id(atom)
        if atom_id in explain:
            data = explain[atom_id]
            print(f"\n   {atom_name}:")
            print(f"     Winner: {data['winner']}")
            print(f"     All candidates:")
            for cand in data['candidates']:
                print(f"       {cand['rank']}. {cand['atomtype']:10s} "
                      f"(prio={cand['priority']}, score={cand['score']}, "
                      f"size={cand['pattern_size']})")
    
    # Export explain as JSON
    print("\n7. Exporting explain data as JSON...")
    explain_json = {}
    for atom_name, atom in atoms.items():
        atom_id = id(atom)
        if atom_id in explain:
            explain_json[atom_name] = explain[atom_id]
    
    json_output = json.dumps(explain_json, indent=2)
    print("\n" + json_output)
    
    print("\n" + "=" * 70)
    print("Example completed successfully!")
    print("=" * 70)


if __name__ == "__main__":
    main()
