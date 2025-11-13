# IO Module Refactoring

## 概述

IO模块已经重构，将原本混乱的383行`__init__.py`重新组织为清晰的模块结构，同时保持100%向后兼容性。

## 设计原则

### 核心设计模式

1. **Reader/Writer模式**
   - 每种文件格式都有专门的Reader和Writer类
   - 所有Reader都实现`read()`方法
   - 所有Writer都实现`write()`方法

2. **工厂函数模式**
   - 提供便利的`read_xxx()`和`write_xxx()`函数
   - 隐藏Reader/Writer类的实例化细节
   - 统一的参数处理逻辑

3. **延迟导入**
   - 只在函数被调用时才导入具体的Reader/Writer类
   - 避免启动时加载所有依赖
   - 减少内存占用

4. **统一接口**
   - 所有reader函数接受可选的`frame`参数
   - 所有writer函数接受数据对象作为参数
   - 一致的参数命名和顺序

## 新模块结构

```
molpy/io/
├── __init__.py              # 简化的公共API，只负责导入和导出
├── readers.py               # 所有读取器工厂函数
├── writers.py               # 所有写入器工厂函数
├── data/                    # 数据文件Reader/Writer实现
│   ├── base.py              # DataReader, DataWriter基类
│   ├── lammps.py            # LAMMPS数据文件
│   ├── pdb.py               # PDB文件
│   ├── gro.py               # GROMACS文件
│   └── ...
├── forcefield/              # 力场文件Reader/Writer实现
│   ├── base.py              # ForceFieldReader, ForceFieldWriter基类
│   ├── xml.py               # XML力场文件
│   ├── lammps.py            # LAMMPS力场文件
│   └── ...
└── trajectory/              # 轨迹文件Reader/Writer实现
    ├── base.py              # BaseTrajectoryReader, TrajectoryWriter基类
    ├── lammps.py            # LAMMPS轨迹文件
    └── xyz.py               # XYZ轨迹文件
```

## 改进点

### 1. 代码组织

**之前**：
- `__init__.py` 383行，包含所有工厂函数
- 逻辑混乱，难以维护
- 重复的frame初始化代码

**之后**：
- `__init__.py` 仅负责导入和导出
- `readers.py` 集中所有读取器
- `writers.py` 集中所有写入器
- 统一的辅助函数`_ensure_frame()`

### 2. 参数一致性

**之前**：
```python
# 有些自动创建frame
read_pdb(file, frame=None)  
if frame is None:
    frame = Frame()

# 有些不创建
read_lammps_data(file, atom_style, frame=None)
# 直接传给reader
```

**之后**：
```python
# 统一使用_ensure_frame()
frame = _ensure_frame(frame)
reader = SomeReader(Path(file))
return reader.read(frame)
```

### 3. 文档完善

**之前**：
- 缺少模块级文档
- 函数文档不完整

**之后**：
- 详细的模块级文档字符串
- 每个函数都有完整的参数说明
- 清晰的使用示例

### 4. 向后兼容性

所有现有代码无需修改，以下导入方式全部支持：

```python
# 公共API导入
from molpy.io import read_xml_forcefield, write_lammps_data

# 直接导入Reader/Writer类
from molpy.io.data.lammps import LammpsDataWriter
from molpy.io.forcefield.xml import XMLForceFieldReader

# 旧的函数名称（通过别名支持）
from molpy.io import read_amber  # read_amber_prmtop的别名
```

## 使用示例

### 基本用法

```python
from molpy.io import (
    read_pdb,
    write_pdb,
    read_xml_forcefield,
    write_lammps_data
)

# 读取数据文件
frame = read_pdb("structure.pdb")

# 写入数据文件
write_pdb("output.pdb", frame)

# 读取力场
ff = read_xml_forcefield("oplsaa.xml")

# 写入LAMMPS数据
write_lammps_data("system.data", frame, atom_style="full")
```

### 高级用法（直接使用Reader/Writer类）

```python
from molpy.io.data.lammps import LammpsDataReader, LammpsDataWriter
from molpy.core.frame import Frame

# 使用Reader类（支持上下文管理器）
with LammpsDataReader("data.lammps", "full") as reader:
    frame = reader.read()

# 使用Writer类
frame = Frame()
writer = LammpsDataWriter("output.data", atom_style="full")
writer.write(frame)
```

### 轨迹文件

```python
from molpy.io import read_lammps_trajectory, write_xyz_trajectory

# 读取轨迹（返回迭代器）
traj = read_lammps_trajectory("dump.lammpstrj")

# 迭代处理每一帧
frames = []
for frame in traj:
    # 处理frame
    frames.append(frame)

# 写入新轨迹
write_xyz_trajectory("output.xyz", frames)
```

## 迁移指南

对于现有代码，**无需任何修改**！所有旧的API都保持不变。

但建议新代码遵循以下最佳实践：

### ✅ 推荐做法

```python
# 使用清晰的工厂函数
from molpy.io import read_pdb, write_pdb

frame = read_pdb("input.pdb")
write_pdb("output.pdb", frame)
```

### ⚠️ 不推荐但支持

```python
# 旧的复合函数（保留用于向后兼容）
from molpy.io import read_lammps, read_amber_system

# 建议拆分为独立调用
from molpy.io import read_lammps_data, read_lammps_forcefield

frame = read_lammps_data("data.lammps", "full")
ff = read_lammps_forcefield("force.field")
```

## 测试验证

重构后所有测试通过：
```
396 passed, 3 skipped in 3.33s
```

向后兼容性验证：
- ✅ 所有公共API导入正常
- ✅ 直接导入Reader/Writer类正常
- ✅ notebook中的代码无需修改
- ✅ 旧的函数名称（别名）正常工作

## 未来改进

虽然当前重构已经大幅改善了代码质量，但还有一些可以进一步优化的地方：

1. **类型注解**：为所有函数添加完整的类型注解
2. **异常处理**：统一异常处理策略
3. **日志系统**：添加结构化日志
4. **配置系统**：支持全局IO配置（如默认编码、缓冲大小等）
5. **插件机制**：支持用户注册自定义文件格式

## 总结

此次重构：
- ✅ **大幅简化**了`__init__.py`（从383行到约220行）
- ✅ **提高了可维护性**：代码按功能清晰组织
- ✅ **增强了可读性**：每个模块职责单一
- ✅ **保持了兼容性**：现有代码无需修改
- ✅ **改进了文档**：完整的文档字符串和使用示例
- ✅ **统一了接口**：一致的参数处理逻辑

**最重要的是**：这次重构遵循了"不破坏业务代码"的原则，所有现有功能都得到保留！

