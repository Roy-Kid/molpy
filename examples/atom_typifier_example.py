"""Example: Using OplsAtomTypifier with SMARTS matcher.

This example demonstrates:
1. Creating patterns for atom typing
2. Building a simple molecule
3. Using OplsAtomTypifier to assign atom types
4. Inspecting the results
"""

from molpy.core.atomistic import Atomistic
from molpy.core.forcefield import ForceField, AtomType
from molpy.typifier.atomistic import OplsAtomTypifier
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import is_element, bond_order


def create_forcefield_with_patterns():
    """Create a forcefield with atom types and their SMARTS patterns."""
    ff = ForceField(name="example_ff")
    
    # 定义 atom styles（这里简化，实际可能需要更复杂的设置）
    # 我们需要为 AtomType 添加 SMARTS pattern 或足够的信息来构造 pattern
    
    # 方式1: 在 AtomType 中存储元素信息（适用于简单情况）
    at_ct = AtomType("CT", element="C", priority=0, mass=12.01, charge=0.0)
    at_oh = AtomType("OH", element="O", priority=5, mass=16.00, charge=-0.5)
    at_c_oh = AtomType("C_OH", element="C", priority=2, mass=12.01, charge=0.2)
    
    # 将 AtomType 添加到 forcefield
    # 由于 ForceField.get_types() 从 styles 中获取 types，我们需要先创建一个 style
    from molpy.core.forcefield import Style
    atom_style = Style("atom")
    atom_style.types.add(at_ct)
    atom_style.types.add(at_oh)
    atom_style.types.add(at_c_oh)
    ff.styles.add(atom_style)
    
    return ff


def create_forcefield_with_custom_matcher():
    """Create a forcefield and manually create patterns for the typifier."""
    ff = ForceField(name="custom_ff")
    
    # 创建 AtomType 对象
    at_ct = AtomType("CT", element="C", priority=0, mass=12.01, charge=0.0)
    at_oh = AtomType("OH", element="O", priority=5, mass=16.00, charge=-0.5)
    at_c_oh = AtomType("C_OH", element="C", priority=2, mass=12.01, charge=0.2)
    
    # 添加到 forcefield
    from molpy.core.forcefield import Style
    atom_style = Style("atom")
    atom_style.types.add(at_ct)
    atom_style.types.add(at_oh)
    atom_style.types.add(at_c_oh)
    ff.styles.add(atom_style)
    
    # 手动创建 patterns（这些将在 OplsAtomTypifier._extract_patterns 中使用）
    # 如果需要更复杂的 patterns，可以在这里定义
    patterns = []
    
    # Generic carbon (lowest priority)
    p1 = quick_pattern("CT", "C", priority=0)
    patterns.append(p1)
    
    # Hydroxyl oxygen (high priority)
    pb2 = PatternBuilder("OH", priority=5, source="hydroxyl_pattern")
    v_c = pb2.add_vertex(preds=[is_element("C")])
    v_o = pb2.add_vertex(preds=[is_element("O")])
    pb2.add_edge(v_c, v_o, preds=[bond_order(1)])
    pb2.set_target_vertices([v_o])  # Only type the oxygen
    p2 = pb2.build()
    patterns.append(p2)
    
    # Alcohol carbon (medium priority)
    pb3 = PatternBuilder("C_OH", priority=2, source="alcohol_pattern")
    v_c3 = pb3.add_vertex(preds=[is_element("C")])
    v_o3 = pb3.add_vertex(preds=[is_element("O")])
    pb3.add_edge(v_c3, v_o3, preds=[bond_order(1)])
    pb3.set_target_vertices([v_c3])  # Only type the carbon
    p3 = pb3.build()
    patterns.append(p3)
    
    # 将 patterns 附加到 forcefield（自定义扩展）
    ff._atom_patterns = patterns
    
    return ff


class CustomOplsAtomTypifier(OplsAtomTypifier):
    """Custom atom typifier that uses pre-defined patterns."""
    
    def _extract_patterns(self):
        """Override to use custom patterns if available."""
        # 如果 forcefield 有预定义的 patterns，使用它们
        if hasattr(self.ff, '_atom_patterns'):
            return self.ff._atom_patterns
        
        # 否则使用默认逻辑
        return super()._extract_patterns()


def create_ethanol():
    """Create ethanol molecule: CH3-CH2-OH"""
    mol = Atomistic()
    
    # Atoms (需要设置足够的属性供 matcher 使用)
    c1 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    c2 = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
    o = mol.add_atom(element="O", is_aromatic=False, charge=0)
    
    # Bonds
    mol.add_bond(c1, c2, order=1)
    mol.add_bond(c2, o, order=1)
    
    return mol, {"C1": c1, "C2": c2, "O": o}


def main():
    print("=" * 70)
    print("OplsAtomTypifier Example: Ethanol Atom Typing")
    print("=" * 70)
    
    # Create forcefield with patterns
    print("\n1. Creating forcefield with custom patterns...")
    ff = create_forcefield_with_custom_matcher()
    atom_types = ff.get_types(AtomType)
    print(f"   Defined {len(atom_types)} atom types:")
    for at in atom_types:
        print(f"     - {at.name}: element={at.get('element')}, priority={at.get('priority')}")
    
    # Create typifier
    print("\n2. Creating OplsAtomTypifier...")
    typifier = CustomOplsAtomTypifier(ff)
    print(f"   Loaded {len(typifier.patterns)} patterns")
    
    # Create molecule
    print("\n3. Creating ethanol molecule (CH3-CH2-OH)...")
    mol, atoms = create_ethanol()
    print(f"   Created molecule with {len(mol.atoms)} atoms, {len(mol.bonds)} bonds")
    
    # Show atoms before typing
    print("\n4. Atoms before typing:")
    atom_names = {id(v): k for k, v in atoms.items()}
    for atom_name, atom in atoms.items():
        elem = atom.get("element")
        print(f"   - {atom_name}: element={elem}, type={atom.get('type', 'UNTYPED')}")
    
    # Apply typifier
    print("\n5. Applying typifier...")
    typifier.typify(mol)
    
    # Show atoms after typing
    print("\n6. Atoms after typing:")
    for atom_name, atom in atoms.items():
        elem = atom.get("element")
        atype = atom.get("type", "UNTYPED")
        charge = atom.get("charge", 0.0)
        mass = atom.get("mass", 0.0)
        print(f"   - {atom_name}: element={elem}, type={atype}, "
              f"charge={charge:.2f}, mass={mass:.2f}")
    
    # Verify expected results
    print("\n7. Verification:")
    expected = {
        "C1": "CT",      # Generic carbon
        "C2": "C_OH",    # Alcohol carbon (C-O)
        "O": "OH"        # Hydroxyl oxygen
    }
    
    all_correct = True
    for atom_name, expected_type in expected.items():
        atom = atoms[atom_name]
        actual_type = atom.get("type", "UNTYPED")
        status = "✓" if actual_type == expected_type else "✗"
        if actual_type != expected_type:
            all_correct = False
        print(f"   {status} {atom_name}: expected={expected_type}, actual={actual_type}")
    
    print("\n" + "=" * 70)
    if all_correct:
        print("SUCCESS: All atoms typed correctly!")
    else:
        print("FAILED: Some atoms were not typed as expected")
    print("=" * 70)


if __name__ == "__main__":
    main()
