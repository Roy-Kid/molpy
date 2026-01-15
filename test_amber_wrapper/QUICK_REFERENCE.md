# AMBER Wrapper 快速参考

## 安装和导入

```python
from pathlib import Path
from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper
from molpy.io.readers import read_amber_prmtop, read_amber_frcmod
from molpy.io.writers import write_amber_frcmod
```

## 完整工作流（3 步）

```python
workdir = Path("./amber_work")

# 1️⃣ Antechamber: PDB → MOL2（原子类型 + 电荷）
ante = AntechamberWrapper(name="ante", workdir=workdir / "ante")
ante.atomtype_assign(
    input_file="ligand.pdb",
    output_file="ligand.mol2",
    charge_method="bcc",      # 推荐
    atom_type="gaff2",        # 推荐
    net_charge=-1,
)

# 2️⃣ Parmchk2: MOL2 → FRCMOD（缺失参数）
parmchk = Parmchk2Wrapper(name="parmchk", workdir=workdir / "parmchk")
parmchk.generate_parameters(
    input_file=workdir / "ante" / "ligand.mol2",
    output_file="ligand.frcmod",
    parameter_level=2,        # 推荐
)

# 3️⃣ Tleap: 构建系统 → PRMTOP + INPCRD
tleap = TLeapWrapper(name="tleap", workdir=workdir / "tleap")
tleap.run_from_script(
    script_text=f"""source leaprc.gaff2
loadAmberParams {workdir / 'parmchk' / 'ligand.frcmod'}
LIG = loadMol2 {workdir / 'ante' / 'ligand.mol2'}
saveAmberParm LIG ligand.prmtop ligand.inpcrd
quit
"""
)

# 4️⃣ 读取为 Frame
frame, ff = read_amber_prmtop(
    workdir / "tleap" / "ligand.prmtop",
    workdir / "tleap" / "ligand.inpcrd",
)
```

## API 速查

### AntechamberWrapper.atomtype_assign()

| 参数 | 类型 | 可选值 | 默认 |
|------|------|--------|------|
| `input_format` | str | `"pdb"`, `"mol2"`, `"ac"` | `"pdb"` |
| `output_format` | str | `"mol2"`, `"ac"` | `"mol2"` |
| `charge_method` | str | `"gas"`, `"bcc"`, `"be3"`, `"cm2"`, `"esp"` | `"bcc"` |
| `atom_type` | str | `"gaff"`, `"gaff2"`, `"amber"`, `"sybyl"` | `"gaff2"` |
| `net_charge` | int | 任意整数 | `0` |
| `formal_charges` | bool | `True`, `False` | `False` |

**推荐设置**：
- 小分子：`charge_method="bcc"`, `atom_type="gaff2"`
- 蛋白质：`charge_method="be3"`, `atom_type="amber"`

### Parmchk2Wrapper.generate_parameters()

| 参数 | 类型 | 可选值 | 默认 |
|------|------|--------|------|
| `input_format` | str | `"mol2"`, `"ac"`, `"mol"`, `"pdb"` | `"mol2"` |
| `parameter_level` | int | `1` (仅缺失), `2` (全部) | `2` |

**推荐设置**：`parameter_level=2`

### TLeapWrapper.run_from_script()

| 参数 | 类型 | 默认 |
|------|------|------|
| `script_text` | str | 必需 |
| `script_name` | str | `"tleap.in"` |

**常用 tleap 命令**：
```tcl
source leaprc.gaff2           # 加载力场
loadAmberParams file.frcmod   # 加载自定义参数
LIG = loadMol2 lig.mol2       # 加载分子
saveAmberParm LIG out.prmtop out.inpcrd  # 保存
quit                          # 退出
```

## 常见场景

### 场景 1: 单个小分子配体

```python
ante = AntechamberWrapper(workdir=Path("work/ante"))
ante.atomtype_assign("ligand.pdb", "ligand.mol2", net_charge=0)

parmchk = Parmchk2Wrapper(workdir=Path("work/parmchk"))
parmchk.generate_parameters("work/ante/ligand.mol2", "ligand.frcmod")

tleap = TLeapWrapper(workdir=Path("work/tleap"))
tleap.run_from_script("""
source leaprc.gaff2
loadAmberParams work/parmchk/ligand.frcmod
LIG = loadMol2 work/ante/ligand.mol2
saveAmberParm LIG ligand.prmtop ligand.inpcrd
quit
""")
```

### 场景 2: 配体 + 反离子（如 Li-TFSI）

```python
# TFSI 阴离子
ante.atomtype_assign("tfsi.pdb", "tfsi.mol2", net_charge=-1)
parmchk.generate_parameters("tfsi.mol2", "tfsi.frcmod")

# 构建复合物
tleap.run_from_script("""
source leaprc.gaff2
loadAmberParams tfsi.frcmod
TFSI = loadMol2 tfsi.mol2
Li = Li+
complex = combine { Li TFSI }
saveAmberParm complex litfsi.prmtop litfsi.inpcrd
quit
""")
```

### 场景 3: 读取和修改 FRCMOD

```python
# 读取
data = read_amber_frcmod("ligand.frcmod")
print(f"Bond params: {data['bond']}")

# 修改并写回
write_amber_frcmod(
    "modified.frcmod",
    remark="Modified parameters",
    bond=data['bond'] + "\nc3-oh   320.0   1.42",
    angle=data['angle'],
    dihe=data['dihe'],
)
```

## 环境配置（可选）

### 使用 Conda 环境

```python
ante = AntechamberWrapper(
    workdir=Path("work"),
    env="AmberTools25",      # Conda 环境名
    env_manager="conda",
)
```

### 使用虚拟环境

```python
ante = AntechamberWrapper(
    workdir=Path("work"),
    env=Path("/path/to/venv"),
    env_manager="venv",
)
```

### 系统级别（默认）

```python
ante = AntechamberWrapper(workdir=Path("work"))
# 工具需在 PATH 中
```

## 错误处理

```python
# 检查工具是否可用
if not ante.is_available():
    print("Antechamber 不可用")
    
# 或者直接检查（会抛出异常）
ante.check()  # FileNotFoundError if not found

# 检查返回码
result = ante.atomtype_assign(...)
if result.returncode != 0:
    print(f"失败: {result.stderr}")
```

## 文件输出位置

所有输出文件都在 `wrapper.workdir` 中：

```python
workdir = Path("my_work")
ante = AntechamberWrapper(workdir=workdir)

# 输出文件会在: my_work/ligand.mol2
ante.atomtype_assign("ligand.pdb", "ligand.mol2")
```

## 提示和技巧

✅ **推荐命名规范**：
- 分子名 + 后缀：`tfsi.pdb` → `tfsi.mol2` → `tfsi.frcmod`

✅ **路径处理**：
- 使用 `Path` 对象处理路径
- 相对路径相对于 `workdir`

✅ **电荷计算**：
- BCC (`"bcc"`)：最常用，速度快，精度好
- ESP (`"esp"`)：需要量化结果，精度高但慢

✅ **力场选择**：
- GAFF2 (`"gaff2"`)：推荐用于小分子
- GAFF (`"gaff"`)：旧版本，向后兼容

❌ **避免**：
- 不要混用不同版本的力场（如 gaff + gaff2）
- 不要忘记设置正确的 `net_charge`
