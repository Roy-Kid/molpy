# Typifier 使用说明

## 概述

Typifier 模块用于为分子结构中的原子、键、角和二面角分配力场类型。主要功能包括：

- **OplsBondTypifier**: 根据键两端原子的类型匹配键类型
- **OplsAngleTypifier**: 根据角三个原子的类型匹配角类型
- **OplsDihedralTypifier**: 根据二面角四个原子的类型匹配二面角类型
- **OplsAtomisticTypifier**: 为整个分子结构分配所有类型

## 前提条件

1. **力场文件**: 需要有包含原子类型、键类型、角类型和二面角类型定义的力场文件（如 OPLS-AA）
2. **原子类型**: 结构中的所有原子必须已经分配了 `type` 属性

## 核心机制

### 1. 键类型匹配（Bond Typification）

根据键两端原子的类型查找对应的键类型。系统会构建一个查找表：

```
(atom_type_1, atom_type_2) -> BondType
```

由于键是无向的，所以会双向存储：
- `(CT, HC)` -> BondType_1
- `(HC, CT)` -> BondType_1

### 2. 角类型匹配（Angle Typification）

根据角三个原子的类型查找对应的角类型：

```
(atom_type_1, atom_type_2, atom_type_3) -> AngleType
```

角也可以翻转，所以会双向存储：
- `(HC, CT, HC)` -> AngleType_1
- `(HC, CT, HC)` -> AngleType_1 （中心原子相同）

### 3. 二面角类型匹配（Dihedral Typification）

根据二面角四个原子的类型查找对应的二面角类型：

```
(atom_type_1, atom_type_2, atom_type_3, atom_type_4) -> DihedralType
```

## 使用示例

### 示例 1: 为单个键分配类型

```python
from pathlib import Path
from molpy.io.forcefield import read_xml_forcefield
from molpy.core.atomistic import Atom, Bond
from molpy.typifier.atomistic import OplsBondTypifier

# 1. 读取力场
ff = read_xml_forcefield(Path("path/to/opls.xml"))

# 2. 创建 bond typifier
bond_typifier = OplsBondTypifier(ff)

# 3. 创建一个键（原子已有类型）
atom1 = Atom(type="CT", symbol="C")
atom2 = Atom(type="HC", symbol="H")
bond = Bond(atom1, atom2)

# 4. 分配类型
bond_typifier.typify(bond)

# 5. 查看结果
print(f"Bond type: {bond.get('bondtype')}")
print(f"Parameters: {bond.get('bondtype_params')}")
```

### 示例 2: 为整个分子结构分配类型

```python
from molpy.typifier.atomistic import OplsAtomisticTypifier
from molpy.core.atomistic import Atomistic, Atom, Bond, Angle

# 1. 读取力场
ff = read_xml_forcefield(Path("path/to/opls.xml"))

# 2. 创建 typifier
typifier = OplsAtomisticTypifier(ff)

# 3. 创建或读取结构（原子必须已有 type 属性）
struct = Atomistic()

# 添加原子
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

# 添加角
angle = Angle(h1, c1, h2)
struct.links.add(angle)

# 4. 为整个结构分配类型
typifier.typify(struct)

# 5. 查看结果
for bond in struct.bonds:
    print(f"{bond}: type={bond.get('bondtype')}")

for angle in struct.links.bucket(Angle):
    print(f"{angle}: type={angle.get('angletype')}")
```

### 示例 3: 完整工作流程

```python
# 1. 读取力场
ff = read_xml_forcefield(Path("opls_aa.xml"))

# 2. 读取分子结构
from molpy.io import read_pdb
struct = read_pdb("molecule.pdb")

# 3. 首先为原子分配类型（使用 SMARTS 匹配或其他方法）
# from molpy.typifier import OplsAtomTypifier
# atom_typifier = OplsAtomTypifier(ff)
# atom_typifier.typify(struct)

# 4. 然后为 bonds/angles/dihedrals 分配类型
typifier = OplsAtomisticTypifier(ff)
typifier.typify(struct)

# 5. 输出结果
from molpy.io import write_lammps_data
write_lammps_data("output.data", struct, ff)
```

## 力场文件格式

力场文件中需要包含类型定义，例如：

```xml
<ForceField>
  <AtomTypes>
    <Type name="CT" class="CT" element="C" mass="12.011"/>
    <Type name="HC" class="HC" element="H" mass="1.008"/>
  </AtomTypes>
  
  <BondForce>
    <Bond type1="CT" type2="HC" length="0.1090" k="284512.0"/>
  </BondForce>
  
  <AngleForce>
    <Angle type1="HC" type2="CT" type3="HC" angle="1.88146" k="276.144"/>
  </AngleForce>
  
  <TorsionForce>
    <Proper type1="HC" type2="CT" type3="CT" type4="HC" 
            k1="0.0" k2="0.0" k3="0.6276" k4="0.0"/>
  </TorsionForce>
</ForceField>
```

## 数据存储

分配类型后，相关信息会存储在对象的 `data` 字典中：

- **Bond**: 
  - `bondtype`: 键类型名称
  - `bondtype_params`: 键参数（如键长、力常数）
  
- **Angle**:
  - `angletype`: 角类型名称
  - `angletype_params`: 角参数（如平衡角、力常数）
  
- **Dihedral**:
  - `dihedraltype`: 二面角类型名称
  - `dihedraltype_params`: 二面角参数

## 错误处理

如果找不到匹配的类型，系统会抛出 `ValueError`：

```python
try:
    typifier.typify(bond)
except ValueError as e:
    print(f"No matching bond type found: {e}")
```

## 性能优化

- 查找表在初始化时构建，避免重复计算
- 使用字典查找，时间复杂度 O(1)
- 支持批量处理整个结构

## 扩展

可以通过继承 `TypifierBase` 创建自定义 typifier：

```python
class CustomBondTypifier(TypifierBase[Bond]):
    def __init__(self, forcefield: ForceField):
        super().__init__(forcefield)
        # 自定义初始化逻辑
        
    def typify(self, bond: Bond) -> Bond:
        # 自定义类型分配逻辑
        return bond
```
