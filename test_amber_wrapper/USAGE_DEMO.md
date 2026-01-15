# AMBER 工具链完整工作流演示

从 `tfsi.pdb` → `inpcrd` + `prmtop` → `Frame` 的四个层级方案

---

## 方案 A: 低层次 - 手工调用各个 Wrapper

**特点**：完全手动控制，代码较长但灵活

```python
from pathlib import Path
from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper
from molpy.io.readers import read_amber_prmtop

workdir = Path("./amber_work")
workdir.mkdir(exist_ok=True)

# 1. 初始化三个 wrapper
ante = AntechamberWrapper(name="ante", workdir=workdir / "ante")
parmchk = Parmchk2Wrapper(name="parmchk", workdir=workdir / "parmchk")
tleap = TLeapWrapper(name="tleap", workdir=workdir / "tleap")

# 验证工具可用
ante.check()
parmchk.check()
tleap.check()

# 2. Stage 1: Antechamber - PDB → MOL2（原子类型和电荷计算）
result = ante.run_raw(
    args=[
        "-i", "tfsi.pdb",
        "-fi", "pdb",
        "-o", "tfsi_gaff2.mol2",
        "-fo", "mol2",
        "-c", "bcc",
        "-nc", "-1",  # 净电荷 -1
        "-at", "gaff2",
    ],
    cwd=workdir / "ante",
)
assert result.returncode == 0, f"Antechamber failed: {result.stderr}"

mol2_file = workdir / "ante" / "tfsi_gaff2.mol2"
print(f"✓ Generated: {mol2_file}")

# 3. Stage 2: Parmchk2 - MOL2 → FRCMOD（参数检查和补充）
result = parmchk.run_raw(
    args=[
        "-i", str(mol2_file),
        "-f", "mol2",
        "-o", "tfsi.frcmod",
    ],
    cwd=workdir / "parmchk",
)
assert result.returncode == 0, f"Parmchk2 failed: {result.stderr}"

frcmod_file = workdir / "parmchk" / "tfsi.frcmod"
print(f"✓ Generated: {frcmod_file}")

# 4. Stage 3: Tleap - 构建系统（MOL2 + FRCMOD → PRMTOP + INPCRD）
tleap_script = f"""source leaprc.gaff2

loadAmberParams {frcmod_file}
TFSI = loadMol2 {mol2_file}

Li = Li+

complex = combine {{ Li TFSI }}

saveAmberParm complex litfsi.prmtop litfsi.inpcrd

quit
"""

result = tleap.run_script(
    script_text=tleap_script,
    script_name="tleap.in",
    cwd=workdir / "tleap",
)
assert result.returncode == 0, f"Tleap failed: {result.stderr}"

prmtop_file = workdir / "tleap" / "litfsi.prmtop"
inpcrd_file = workdir / "tleap" / "litfsi.inpcrd"
print(f"✓ Generated: {prmtop_file}, {inpcrd_file}")

# 5. Stage 4: 读取结果到 Frame
frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
print(f"✓ Frame loaded with {len(frame['atoms'])} atoms")
print(f"✓ ForceField loaded with {len(forcefield.atom_types)} atom types")
```

---

## 方案 B: 中层次 - 使用 Wrapper 辅助方法

**特点**：简洁易用，减少重复代码

```python
from pathlib import Path
from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper
from molpy.io.readers import read_amber_prmtop

workdir = Path("./amber_work")

# 1. Antechamber：一个方法完成 PDB → MOL2
ante = AntechamberWrapper(name="ante", workdir=workdir / "ante")
mol2_file = ante.prepare_pdb(
    input_pdb="tfsi.pdb",
    output_format="mol2",
    charge_method="bcc",
    atom_type="gaff2",
    net_charge=-1,
)
print(f"✓ Antechamber: {mol2_file}")

# 2. Parmchk2：一个方法完成 MOL2 → FRCMOD
parmchk = Parmchk2Wrapper(name="parmchk", workdir=workdir / "parmchk")
frcmod_file = parmchk.generate_frcmod(
    input_mol2=mol2_file,
    output_frcmod="tfsi.frcmod",
)
print(f"✓ Parmchk2: {frcmod_file}")

# 3. Tleap：使用模板方法，避免手写脚本
tleap = TLeapWrapper(name="tleap", workdir=workdir / "tleap")
prmtop_file, inpcrd_file = tleap.run_with_template(
    template="""source leaprc.gaff2

loadAmberParams {frcmod_path}
TFSI = loadMol2 {mol2_path}

Li = Li+

complex = combine {{ Li TFSI }}

saveAmberParm complex litfsi.prmtop litfsi.inpcrd

quit
""",
    context={
        "frcmod_path": str(frcmod_file),
        "mol2_path": str(mol2_file),
    },
    prmtop_output="litfsi.prmtop",
    inpcrd_output="litfsi.inpcrd",
)
print(f"✓ Tleap: {prmtop_file}, {inpcrd_file}")

# 4. 读取结果
frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
print(f"✓ Frame: {len(frame['atoms'])} atoms")
```

---

## 方案 C: 高层次 - 使用 Compute 节点

**特点**：一行代码完成全流程（最终目标）

```python
from pathlib import Path
from molpy.core.atomistic import Atomistic
from molpy.compute.amber import PrepareLigandAMBER

# 加载初始 PDB
atomistic = Atomistic.from_pdb("tfsi.pdb")

# 一个节点，一个方法，完成所有步骤
prepare = PrepareLigandAMBER(
    charge_method="bcc",
    atom_type="gaff2",
    net_charge=-1,
    workdir=Path("./amber_work"),
)

result = prepare.compute(atomistic)

frame = result.frame
forcefield = result.forcefield
print(f"✓ Frame: {len(frame['atoms'])} atoms")
print(f"✓ ForceField: {len(forcefield.atom_types)} atom types")
```

---

## 方案 D: 带缓存和批处理

**特点**：如果中间结果已存在，直接使用缓存；支持并行处理

```python
from pathlib import Path
from molpy.compute.amber import CachedAMBERPipeline

# 创建带缓存的管道
pipeline = CachedAMBERPipeline(
    cache_dir=Path("./cache"),
    keep_intermediates=True,  # 保留 mol2, frcmod 等中间文件
)

# 第一次运行：执行所有步骤
frame1, ff1 = pipeline.prepare_ligand(
    pdb_file="tfsi.pdb",
    charge_method="bcc",
    atom_type="gaff2",
    net_charge=-1,
)

# 第二次相同参数：从缓存读取，秒级完成
frame2, ff2 = pipeline.prepare_ligand(
    pdb_file="tfsi.pdb",
    charge_method="bcc",
    atom_type="gaff2",
    net_charge=-1,
)

# 批量处理多个分子
pdb_files = ["mol1.pdb", "mol2.pdb", "mol3.pdb"]
results = pipeline.prepare_ligands_batch(
    pdb_files=pdb_files,
    charge_method="bcc",
    atom_type="gaff2",
    num_workers=4,  # 并行 4 个
)

for pdb_file, (frame, ff) in zip(pdb_files, results):
    print(f"{pdb_file}: {len(frame['atoms'])} atoms")
```

---

## 工作流对比

| 特性 | 方案 A | 方案 B | 方案 C | 方案 D |
|------|--------|--------|--------|--------|
| **代码行数** | ~40 | ~25 | ~15 | ~20 |
| **灵活性** | ★★★★★ | ★★★★☆ | ★★☆☆☆ | ★★★☆☆ |
| **易用性** | ★☆☆☆☆ | ★★★☆☆ | ★★★★★ | ★★★★☆ |
| **性能** | 基准 | 基准 | 基准+验证 | 基准+缓存 |
| **推荐场景** | 定制工作流 | 一般使用 | 常规任务 | 批量处理 |

---

## 实现路线图

### 第一阶段（必需）
- [ ] `Parmchk2Wrapper` - 基本 `run_raw()` 支持
- [ ] `Frcmod` I/O - `read_frcmod()` 和 `write_frcmod()`

### 第二阶段（推荐）
- [ ] `AntechamberWrapper.prepare_pdb()` - 方案 B 支持
- [ ] `Parmchk2Wrapper.generate_frcmod()` - 方案 B 支持
- [ ] `TLeapWrapper.run_with_template()` - 方案 B 支持

### 第三阶段（可选）
- [ ] `PrepareLigandAMBER` Compute 节点 - 方案 C 支持
- [ ] `CachedAMBERPipeline` - 方案 D 支持
