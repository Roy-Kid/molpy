"""
示例：使用 OplsAtomisticTypifier 为结构分配 bond/angle/dihedral 类型

前提：
1. 已有力场文件（如 OPLS 力场）
2. Atomistic 结构中的原子已经有 'type' 属性
"""

from pathlib import Path
from molpy.io.forcefield import read_xml_forcefield
from molpy.core.atomistic import Atom, Bond, Angle, Dihedral, Atomistic
from molpy.typifier.atomistic import (
    OplsBondTypifier,
    OplsAngleTypifier,
    OplsDihedralTypifier,
    OplsAtomisticTypifier,
)


def example_bond_typification():
    """示例：为单个 bond 分配类型"""
    # 1. 读取力场（需要提供实际的文件路径）
    ff = read_xml_forcefield(Path("path/to/opls.xml"))
    
    # 2. 创建 bond typifier
    bond_typifier = OplsBondTypifier(ff)
    
    # 3. 创建一个测试 bond（假设两个原子已经有类型）
    atom1 = Atom(type="CT", symbol="C")
    atom2 = Atom(type="HC", symbol="H")
    bond = Bond(atom1, atom2)
    
    # 4. 为 bond 分配类型
    try:
        bond_typifier.typify(bond)
        print(f"Bond type: {bond.get('bondtype')}")
        print(f"Bond parameters: {bond.get('bondtype_params')}")
    except ValueError as e:
        print(f"Error: {e}")


def example_atomistic_typification():
    """示例：为整个 Atomistic 结构分配所有类型"""
    # 1. 读取力场（需要提供实际的文件路径）
    ff = read_xml_forcefield(Path("path/to/opls.xml"))
    
    # 2. 创建 atomistic typifier
    typifier = OplsAtomisticTypifier(ff)
    
    # 3. 创建或读取一个 Atomistic 结构
    # 假设你已经有一个结构，并且原子已经有了 'type' 属性
    struct = Atomistic()
    
    # 添加一些原子（示例）
    c1 = Atom(type="CT", symbol="C", pos=[0, 0, 0])
    h1 = Atom(type="HC", symbol="H", pos=[1.09, 0, 0])
    h2 = Atom(type="HC", symbol="H", pos=[-0.545, 0.943, 0])
    
    struct.entities.add(c1)
    struct.entities.add(h1)
    struct.entities.add(h2)
    
    # 添加键
    bond1 = Bond(c1, h1)
    bond2 = Bond(c1, h2)
    struct.links.add(bond1)
    struct.links.add(bond2)
    
    # 添加角（如果需要）
    angle = Angle(h1, c1, h2)
    struct.links.add(angle)
    
    # 4. 为整个结构分配类型
    try:
        typifier.typify(struct)
        print("Typification successful!")
        
        # 查看结果
        for bond in struct.bonds:
            print(f"{bond}: type={bond.get('bondtype')}")
        
        for angle in struct.links.bucket(Angle):
            print(f"{angle}: type={angle.get('angletype')}")
            
    except ValueError as e:
        print(f"Error during typification: {e}")


def example_with_real_structure():
    """
    示例：从文件读取结构并分配类型
    
    工作流程：
    1. 读取结构文件（如 PDB, XYZ 等）
    2. 使用 atom typifier 为原子分配类型（基于 SMARTS 匹配等）
    3. 使用 bond/angle/dihedral typifier 为连接分配类型
    """
    # 1. 读取力场
    ff = read_xml_forcefield(Path("path/to/opls.xml"))
    
    # 2. 读取结构（示例）
    # from molpy.io import read_pdb
    # struct = read_pdb("molecule.pdb")
    
    # 3. 首先为原子分配类型（需要使用 atom typifier）
    # atom_typifier = OplsAtomTypifier(ff)
    # atom_typifier.typify(struct)
    
    # 4. 然后为 bonds/angles/dihedrals 分配类型
    # typifier = OplsAtomisticTypifier(ff)
    # typifier.typify(struct)
    
    # 5. 输出或保存结果
    # from molpy.io import write_lammps_data
    # write_lammps_data("output.data", struct, ff)
    
    pass


if __name__ == "__main__":
    print("=== Bond Typification Example ===")
    example_bond_typification()
    
    print("\n=== Atomistic Typification Example ===")
    example_atomistic_typification()
    
    print("\n=== Real Structure Example ===")
    print("See example_with_real_structure() for workflow")
