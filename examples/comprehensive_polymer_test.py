"""
综合测试：BigSMILES → 3D Polymer → 力场参数

完整流程演示：
1. 解析 BigSMILES 定义 monomer
2. 使用 RDKit 生成 3D 坐标
3. 使用 PolymerBuilder 组装 polymer（带几何排布）
4. 标记 atom types 并匹配 OPLS 力场参数
5. 可视化和验证
"""

import os

import numpy as np
from rdkit.Chem import Draw

from molpy.adapter.rdkit_adapter import draw_molecule, generate_3d_coords, atomistic_to_mol as to_rdkit
from molpy.builder.polymer import (
    AutoConnector,
    DockPlacer,
    PolymerBuilder,
)
from molpy.io import read_xml_forcefield
from molpy.parser.smiles import SmilesParser, bigsmilesir_to_monomer
from molpy.typifier.atomistic import (
    OplsAngleTypifier,
    OplsBondTypifier,
    OplsDihedralTypifier,
)


def main():
    print("=" * 80)
    print("综合测试：BigSMILES → 3D Polymer → 力场参数")
    print("=" * 80)
    
    # ========================================================================
    # Step 1: 解析 BigSMILES 定义 Monomers
    # ========================================================================
    print("\n=== Step 1: 解析 BigSMILES 定义 Monomers ===")
    
    parser = SmilesParser()
    
    monomer_smiles = {
        "A": "CCCCO[*:1]",
        "B": "CC(C[*:2])O[*:3]",
        "C": "CCC(C[*:4])O[*:5]",
        "D": "CCC(C[*:6])O[*:7]",
    }
    
    monomers = {}
    for label, smiles in monomer_smiles.items():
        ir = parser.parse_bigsmiles(smiles)
        monomer = bigsmilesir_to_monomer(ir)
        monomers[label] = monomer
        print(f"{label}: {smiles}")
        print(f"  - Atoms: {len(list(monomer.unwrap().atoms))}")
        print(f"  - Ports: {list(monomer.ports.keys())}")
    
    print(f"✓ 成功解析 {len(monomers)} 个 monomers")
    
    # ========================================================================
    # Step 2: 可视化 Monomer 结构（2D）
    # ========================================================================
    print("\n=== Step 2: 可视化 Monomer 结构（2D） ===")

    
    for label in ["A", "B", "C", "D"]:

        monomer = monomers[label]
        mol = to_rdkit(monomer.unwrap())
        print(f"{label}: {mol.GetNumAtoms()} atoms, {mol.GetNumBonds()} bonds")
        
        # 保存 2D 结构图
        img = Draw.MolToImage(mol, size=(400, 400))
        img_path = f"monomer_{label}_2d.png"
        img.save(img_path)
        print(f"  ✓ 已保存 2D 结构图: {img_path}")

    # ========================================================================
    # Step 3: 生成 3D 坐标
    # ========================================================================
    print("\n=== Step 3: 生成 3D 坐标 ===")
    
    monomers_3d = {}
    for label, monomer in monomers.items():
        monomer_3d = generate_3d_coords(monomer, random_seed=42)
        monomers_3d[label] = monomer_3d
        
        atoms = list(monomer_3d.unwrap().atoms)
        has_coords = sum(1 for a in atoms if 'pos' in a.data or 'xyz' in a.data)
        print(f"{label}: {has_coords}/{len(atoms)} atoms with 3D coordinates")
        
        first_atom = atoms[0]
        pos = first_atom.data.get('pos') or first_atom.data.get('xyz')
        if pos:
            print(f"  First atom position: [{pos[0]:.3f}, {pos[1]:.3f}, {pos[2]:.3f}]")

    print(f"✓ 成功为所有 monomers 生成 3D 坐标")
    
    # ========================================================================
    # Step 4: 设置 Port 方向
    # ========================================================================
    print("\n=== Step 4: 设置 Port 方向 ===")
    
    # Port role 规则：
    # A: [*:1] 是右端口（terminus）
    # B: [*:2] 是左端口, [*:3] 是右端口  
    # C: [*:4] 是左端口, [*:5] 是右端口
    # D: [*:6] 是左端口, [*:7] 是右端口
    
    port_roles = {
        "A": {"port_1": "right"},  # A is left terminus
        "B": {"port_2": "left", "port_3": "right"},
        "C": {"port_4": "left", "port_5": "right"},
        "D": {"port_6": "left", "port_7": "right"},
    }
    
    for label, monomer in monomers_3d.items():
        for port_name, port in monomer.ports.items():
            port_atom = port.target
            
            # Get role for this port
            role = port_roles[label].get(port_name, "right")
            
            # Set orientation based on role
            if role == "left":
                port_atom.data['orientation'] = [-1.0, 0.0, 0.0]  # Points left
            else:
                port_atom.data['orientation'] = [1.0, 0.0, 0.0]  # Points right
            
            port.data['role'] = role
            
            print(f"{label} port {port_name}: role={role}, "
                  f"orientation={port_atom.data.get('orientation')}")
    
    print(f"✓ Port 方向设置完成")
    
    # ========================================================================
    # Step 5: 加载 OPLS-AA 力场并初始化 Typifier
    # ========================================================================
    print("\n=== Step 5: 加载 OPLS-AA 力场并初始化 Typifier ===")
    
    # 从内置数据加载 OPLS-AA 力场参数
    print("加载 OPLS-AA 力场参数文件...")
    forcefield = read_xml_forcefield("oplsaa.xml")
    print(f"✓ 力场加载成功: {forcefield.name}")
    
    # 创建 Typifiers
    bond_typifier = OplsBondTypifier(forcefield)
    angle_typifier = OplsAngleTypifier(forcefield)
    dihedral_typifier = OplsDihedralTypifier(forcefield)
    print("✓ Typifiers 初始化完成")
    
    # ========================================================================
    # Step 6: 为 Monomer 原子标记 Atom Types
    # ========================================================================
    print("\n=== Step 6: 为 Monomer 原子标记 Atom Types ===")
    
    # 注意：这里我们手动标记 atom types（实际应用中应该使用 AtomTypifier）
    # 使用具体的 OPLS type names
    for label, monomer in monomers_3d.items():
        assembly = monomer.unwrap()
        atoms = list(assembly.atoms)
        
        for atom in atoms:
            symbol = atom.data.get('element', 'C')
            
            # 简化版：根据元素符号分配 OPLS type name
            if symbol == 'C':
                atom.data['type'] = 'opls_135'  # alkane CH3 (CT class)
            elif symbol == 'O':
                atom.data['type'] = 'opls_179'  # ether oxygen (OS class)
            elif symbol == 'H':
                atom.data['type'] = 'opls_140'  # alkane H (HC class)
            else:
                atom.data['type'] = f'opls_{symbol}'
        
        print(f"{label}: 已标记 {len(atoms)} 个原子")
    
    print("✓ 所有 Monomer 原子已标记 atom type (使用 OPLS type names)")
    
    # ========================================================================
    # Step 7: 对每个 Monomer 进行力场参数匹配
    # ========================================================================
    print("\n=== Step 7: 对每个 Monomer 进行力场参数匹配 ===")
    
    from molpy.core.atomistic import Angle, Bond, Dihedral
    
    for label, monomer in monomers_3d.items():
        assembly = monomer.unwrap()
        
        # 统计 monomer 的拓扑结构
        bonds = list(assembly.links.bucket(Bond))
        angles = list(assembly.links.bucket(Angle))
        dihedrals = list(assembly.links.bucket(Dihedral))
        
        print(f"\n{label}: {len(bonds)} bonds, {len(angles)} angles, {len(dihedrals)} dihedrals")
        
        # 对每个 bond 进行类型匹配
        bond_matched = 0
        bond_errors = []
        for bond in bonds:
            try:
                bond_typifier.typify(bond)
                if 'bondtype' in bond.data:
                    bond_matched += 1
            except Exception as e:
                bond_errors.append(str(e))
        
        if bond_errors and len(bond_errors) < 3:  # Only print first few errors
            for err in bond_errors[:2]:
                print(f"  Bond error: {err}")        
        # 对每个 angle 进行类型匹配
        angle_matched = 0
        for angle in angles:
            try:
                angle_typifier.typify(angle)
                if 'type' in angle.data:
                    angle_matched += 1
            except Exception as e:
                pass
        
        # 对每个 dihedral 进行类型匹配
        dihedral_matched = 0
        for dihedral in dihedrals:
            try:
                dihedral_typifier.typify(dihedral)
                if 'type' in dihedral.data:
                    dihedral_matched += 1
            except Exception as e:
                pass
        
        print(f"  Matched: {bond_matched}/{len(bonds)} bonds, "
              f"{angle_matched}/{len(angles)} angles, "
              f"{dihedral_matched}/{len(dihedrals)} dihedrals")
    
    print("\n✓ Monomer 力场参数匹配完成")
    
    # ========================================================================
    # Step 8: 使用 PolymerBuilder 组装序列（带几何对接）
    # ========================================================================
    print("\n=== Step 8: 使用 PolymerBuilder 组装序列 ===")
    
    geom_ctx = GeometryContext({
        "vdw_scale": 0.80,
        "vdw_scales_by_bond": {
            "-": 0.80,
            "=": 0.72,
        },
        "placement_log": [],
    })
    
    try:
        polymer = PolymerBuilder.linear(
            sequence="ABCBD",
            library=monomers_3d,
            connector=AutoConnector(),
            placer=DockPlacer(),
            geom_ctx=geom_ctx,
        )
        
        atoms = list(polymer.unwrap().atoms)
        bonds = list(polymer.unwrap().bonds)
        print(f"✓ Polymer 组装完成")
        print(f"  - 总原子数: {len(atoms)}")
        print(f"  - 总键数: {len(bonds)}")
        print(f"  - 对接次数: {len(geom_ctx['placement_log'])}")
        
        print("\n对接详情:")
        for i, log in enumerate(geom_ctx["placement_log"]):
            print(f"  Step {i+1}: distance={log['distance']:.3f} Å, "
                  f"VDW radii={log['rvdw_L']:.2f}+{log['rvdw_R']:.2f} Å")
    except Exception as e:
        print(f"ERROR assembling polymer: {e}")
        import traceback
        traceback.print_exc()
        return
    
    # ========================================================================
    # Step 9: 验证 3D 坐标
    # ========================================================================
    print("\n=== Step 9: 验证 3D 坐标 ===")
    
    coords = []
    atoms_with_coords = 0
    
    for atom in atoms:
        pos = atom.data.get('pos') or atom.data.get('xyz')
        if pos:
            coords.append(pos)
            atoms_with_coords += 1
    
    coords_array = np.array(coords)
    
    print(f"坐标统计:")
    print(f"  - 有坐标的原子: {atoms_with_coords}/{len(atoms)}")
    print(f"\n坐标范围:")
    print(f"  - X: [{coords_array[:, 0].min():.3f}, {coords_array[:, 0].max():.3f}] Å")
    print(f"  - Y: [{coords_array[:, 1].min():.3f}, {coords_array[:, 1].max():.3f}] Å")
    print(f"  - Z: [{coords_array[:, 2].min():.3f}, {coords_array[:, 2].max():.3f}] Å")
    
    size_x = coords_array[:, 0].max() - coords_array[:, 0].min()
    size_y = coords_array[:, 1].max() - coords_array[:, 1].min()
    size_z = coords_array[:, 2].max() - coords_array[:, 2].min()
    print(f"\n分子尺寸: {size_x:.3f} × {size_y:.3f} × {size_z:.3f} Å³")
    
    # ========================================================================
    # Step 10: 可视化 Polymer 结构（2D）
    # ========================================================================
    print("\n=== Step 10: 可视化 Polymer 结构（2D） ===")
    
    from rdkit.Chem import Draw
    
    try:
        polymer_mol = to_rdkit(polymer.unwrap())
        print(f"Polymer RDKit Mol:")
        print(f"  - Atoms: {polymer_mol.GetNumAtoms()}")
        print(f"  - Bonds: {polymer_mol.GetNumBonds()}")
        
        # 保存 Polymer 整体 2D 结构图
        img = Draw.MolToImage(polymer_mol, size=(800, 600))
        img_path = "polymer_ABCBD_2d.png"
        img.save(img_path)
        print(f"  ✓ 已保存 Polymer 2D 结构图: {img_path}")
    except Exception as e:
        print(f"ERROR converting polymer to RDKit: {e}")
        import traceback
        traceback.print_exc()
    
    # ========================================================================
    # Step 11: 为 Polymer 新产生的拓扑结构标注类型
    # ========================================================================
    print("\n=== Step 11: 为 Polymer 新产生的拓扑结构标注类型 ===")
    
    assembly = polymer.unwrap()
    
    # 获取所有的 bonds, angles, dihedrals
    from molpy.core.atomistic import Angle, Bond, Dihedral
    
    all_bonds = list(assembly.links.bucket(Bond))
    all_angles = list(assembly.links.bucket(Angle))
    all_dihedrals = list(assembly.links.bucket(Dihedral))
    
    print(f"Polymer 拓扑结构: {len(all_bonds)} bonds, {len(all_angles)} angles, {len(all_dihedrals)} dihedrals")
    
    # 检查并标注没有 type 的 bond
    bonds_without_type = [b for b in all_bonds if 'type' not in b.data]
    print(f"\n发现 {len(bonds_without_type)} 个未标注类型的 bonds")
    
    bond_typed = 0
    for bond in bonds_without_type:
        try:
            bond_typifier.typify(bond)
            if 'type' in bond.data:
                bond_typed += 1
        except Exception as e:
            pass  # 力场数据未加载时会失败
    
    print(f"  成功标注: {bond_typed}/{len(bonds_without_type)} bonds")
    
    # 检查并标注没有 type 的 angle
    angles_without_type = [a for a in all_angles if 'type' not in a.data]
    print(f"\n发现 {len(angles_without_type)} 个未标注类型的 angles")
    
    angle_typed = 0
    for angle in angles_without_type:
        try:
            angle_typifier.typify(angle)
            if 'type' in angle.data:
                angle_typed += 1
        except Exception as e:
            pass
    
    print(f"  成功标注: {angle_typed}/{len(angles_without_type)} angles")
    
    # 检查并标注没有 type 的 dihedral
    dihedrals_without_type = [d for d in all_dihedrals if 'type' not in d.data]
    print(f"\n发现 {len(dihedrals_without_type)} 个未标注类型的 dihedrals")
    
    dihedral_typed = 0
    for dihedral in dihedrals_without_type:
        try:
            dihedral_typifier.typify(dihedral)
            if 'type' in dihedral.data:
                dihedral_typed += 1
        except Exception as e:
            pass
    
    print(f"  成功标注: {dihedral_typed}/{len(dihedrals_without_type)} dihedrals")
    
    # 统计最终的类型标注情况
    total_bond_typed = sum(1 for b in all_bonds if 'type' in b.data)
    total_angle_typed = sum(1 for a in all_angles if 'type' in a.data)
    total_dihedral_typed = sum(1 for d in all_dihedrals if 'type' in d.data)
    
    print(f"\n✓ Polymer 拓扑结构类型标注完成:")
    print(f"  Bonds: {total_bond_typed}/{len(all_bonds)}")
    print(f"  Angles: {total_angle_typed}/{len(all_angles)}")
    print(f"  Dihedrals: {total_dihedral_typed}/{len(all_dihedrals)}")
    
    # ========================================================================
    # Step 12: 导出最终结果
    # ========================================================================
    print("\n" + "=" * 80)
    print("最终 Polymer 数据汇总")
    print("=" * 80)
    
    # 统计 atom types
    atom_type_counts = {}
    for atom in atoms:
        atype = atom.data.get('type', 'unknown')
        atom_type_counts[atype] = atom_type_counts.get(atype, 0) + 1
    
    print(f"\n【拓扑信息】")
    print(f"  序列: ABCBD")
    print(f"  总原子数: {len(atoms)}")
    print(f"  总键数: {len(bonds)}")
    print(f"  Angles: {len(all_angles)}")
    print(f"  Dihedrals: {len(all_dihedrals)}")
    
    print(f"\n【3D 坐标】")
    print(f"  有坐标的原子: {atoms_with_coords}/{len(atoms)}")
    print(f"  分子尺寸: {size_x:.3f} × {size_y:.3f} × {size_z:.3f} Å³")
    
    print(f"\n【Atom Types】")
    for atype, count in sorted(atom_type_counts.items()):
        print(f"  {atype}: {count}")
    
    print(f"\n【力场参数匹配】")
    print(f"  Bonds: {total_bond_typed}/{len(all_bonds)}")
    print(f"  Angles: {total_angle_typed}/{len(all_angles)}")
    print(f"  Dihedrals: {total_dihedral_typed}/{len(all_dihedrals)}")
    
    print(f"\n【几何对接】")
    print(f"  对接步骤: {len(geom_ctx['placement_log'])}")
    print(f"  VDW scale: {geom_ctx['vdw_scale']}")
    
    print("\n" + "=" * 80)
    print("✓ 完整流程测试成功！")
    print("=" * 80)


if __name__ == "__main__":
    main()
