"""
演示：从 tfsi.pdb 到 inpcrd+prmtop 再到 Frame 的完整 AMBER 工作流

这展示了理想的 API 设计和用法。
"""

from pathlib import Path
import tempfile

from molpy.wrapper import (
    AntechamberWrapper,
    Parmchk2Wrapper,  # 新增
    TLeapWrapper,
)
from molpy.io.readers import read_amber_prmtop
from molpy.core.frame import Frame

# ============================================================================
# 方案 A: 低层次 - 逐个工具调用（用户手动编排）
# ============================================================================
def workflow_manual():
    """手工调用各个 wrapper，完全控制流程。"""
    
    workdir = Path("./amber_work")
    workdir.mkdir(exist_ok=True)
    
    # 1. 准备环境（可选，如果工具在 PATH 中则不需要）
    ante_wrapper = AntechamberWrapper(
        name="antechamber",
        workdir=workdir / "antechamber_step",
        # env="AmberTools25",
        # env_manager="conda",
    )
    
    parmchk_wrapper = Parmchk2Wrapper(
        name="parmchk2",
        workdir=workdir / "parmchk2_step",
    )
    
    tleap_wrapper = TLeapWrapper(
        name="tleap",
        workdir=workdir / "tleap_step",
    )
    
    # 2. 验证工具可用
    ante_wrapper.check()
    parmchk_wrapper.check()
    tleap_wrapper.check()
    
    # ========== Stage 1: Antechamber ==========
    # 运行 antechamber 进行原子类型分配和电荷计算
    ante_result = ante_wrapper.run_raw(
        args=[
            "-i", "tfsi.pdb",
            "-fi", "pdb",
            "-o", "tfsi_gaff2.mol2",
            "-fo", "mol2",
            "-c", "bcc",
            "-nc", "-1",
            "-at", "gaff2",
        ],
        cwd=workdir / "antechamber_step",
    )
    
    if ante_result.returncode != 0:
        print(f"Antechamber failed:\n{ante_result.stderr}")
        return None
    
    mol2_file = workdir / "antechamber_step" / "tfsi_gaff2.mol2"
    print(f"✓ Antechamber generated: {mol2_file}")
    
    # ========== Stage 2: Parmchk2 ==========
    # 运行 parmchk2 检查并填补缺失的参数
    parmchk_result = parmchk_wrapper.run_raw(
        args=[
            "-i", str(mol2_file),
            "-f", "mol2",
            "-o", "tfsi.frcmod",
        ],
        cwd=workdir / "parmchk2_step",
    )
    
    if parmchk_result.returncode != 0:
        print(f"Parmchk2 failed:\n{parmchk_result.stderr}")
        return None
    
    frcmod_file = workdir / "parmchk2_step" / "tfsi.frcmod"
    print(f"✓ Parmchk2 generated: {frcmod_file}")
    
    # ========== Stage 3: Tleap ==========
    # 生成 tleap 脚本并执行
    tleap_script = f"""source leaprc.gaff2

loadAmberParams {frcmod_file}
TFSI = loadMol2 {mol2_file}

Li = Li+

complex = combine {{ Li TFSI }}

saveAmberParm complex litfsi.prmtop litfsi.inpcrd

quit
"""
    
    tleap_result = tleap_wrapper.run_script(
        script_text=tleap_script,
        script_name="tleap.in",
        cwd=workdir / "tleap_step",
    )
    
    if tleap_result.returncode != 0:
        print(f"Tleap failed:\n{tleap_result.stderr}")
        return None
    
    prmtop_file = workdir / "tleap_step" / "litfsi.prmtop"
    inpcrd_file = workdir / "tleap_step" / "litfsi.inpcrd"
    print(f"✓ Tleap generated: {prmtop_file}, {inpcrd_file}")
    
    # ========== Stage 4: Read to Frame ==========
    frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
    print(f"✓ Loaded Frame with {len(frame['atoms'])} atoms")
    
    return frame, forcefield


# ============================================================================
# 方案 B: 中层次 - Wrapper 辅助方法（简化常见操作）
# ============================================================================
def workflow_with_helpers():
    """使用 wrapper 的辅助方法，减少重复代码。"""
    
    workdir = Path("./amber_work_v2")
    
    # 1. Antechamber：使用 prepare_pdb 方法（新增）
    ante = AntechamberWrapper(name="ante", workdir=workdir / "ante")
    mol2_file = ante.prepare_pdb(
        input_pdb="tfsi.pdb",
        output_format="mol2",
        charge_method="bcc",
        atom_type="gaff2",
        net_charge=-1,
    )
    print(f"✓ Antechamber: {mol2_file}")
    
    # 2. Parmchk2：使用 generate_frcmod 方法（新增）
    parmchk = Parmchk2Wrapper(name="parmchk", workdir=workdir / "parmchk")
    frcmod_file = parmchk.generate_frcmod(
        input_mol2=mol2_file,
        output_frcmod="tfsi.frcmod",
    )
    print(f"✓ Parmchk2: {frcmod_file}")
    
    # 3. Tleap：使用 run_with_template 方法（新增）
    tleap = TLeapWrapper(name="tleap", workdir=workdir / "tleap")
    prmtop_file, inpcrd_file = tleap.run_with_template(
        template="""source leaprc.gaff2

loadAmberParams {frcmod_file}
TFSI = loadMol2 {mol2_file}

Li = Li+

complex = combine {{ Li TFSI }}

saveAmberParm complex litfsi.prmtop litfsi.inpcrd

quit
""",
        context={
            "frcmod_file": str(frcmod_file),
            "mol2_file": str(mol2_file),
        },
        prmtop_output="litfsi.prmtop",
        inpcrd_output="litfsi.inpcrd",
    )
    print(f"✓ Tleap: {prmtop_file}, {inpcrd_file}")
    
    # 4. 读取结果
    frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
    print(f"✓ Frame: {len(frame['atoms'])} atoms")
    
    return frame, forcefield


# ============================================================================
# 方案 C: 高层次 - Compute 节点（自动编排）
# ============================================================================
def workflow_with_compute():
    """使用高层 Compute 节点，一行代码完成整个工作流。"""
    
    from molpy.core.atomistic import Atomistic
    from molpy.compute.amber import PrepareLigandAMBER
    
    # 加载初始 PDB
    atomistic = Atomistic.from_pdb("tfsi.pdb")
    
    # 创建 compute 节点，配置参数
    prepare = PrepareLigandAMBER(
        charge_method="bcc",
        atom_type="gaff2",
        net_charge=-1,
        workdir=Path("./amber_work_v3"),
    )
    
    # 一个调用完成整个工作流
    result = prepare.compute(atomistic)
    
    # 结果包含 Frame 和 ForceField
    frame = result.frame
    forcefield = result.forcefield
    
    print(f"✓ Frame: {len(frame['atoms'])} atoms")
    print(f"✓ ForceField loaded")
    
    return frame, forcefield


# ============================================================================
# 方案 D: 高层次 + 缓存（可重复调用）
# ============================================================================
def workflow_with_caching():
    """如果缓存了中间输出，可以跳过已完成的步骤。"""
    
    from molpy.compute.amber import PrepareLigandAMBER, CachedAMBERWorkflow
    
    # 单步骤执行（带缓存）
    prepare = CachedAMBERWorkflow(
        cache_dir=Path("./cache"),
        num_workers=1,
    )
    
    frame, forcefield = prepare.from_pdb(
        pdb_file="tfsi.pdb",
        charge_method="bcc",
        atom_type="gaff2",
        net_charge=-1,
    )
    
    # 第二次调用相同参数 → 从缓存读取，无需重新运行
    frame2, forcefield2 = prepare.from_pdb(
        pdb_file="tfsi.pdb",
        charge_method="bcc",
        atom_type="gaff2",
        net_charge=-1,
    )
    
    return frame, forcefield


# ============================================================================
# 使用示例
# ============================================================================
if __name__ == "__main__":
    print("\n" + "="*70)
    print("方案 A: 手工调用（完全控制，代码冗长）")
    print("="*70)
    # result = workflow_manual()
    
    print("\n" + "="*70)
    print("方案 B: 使用 Wrapper 辅助方法（简洁，推荐）")
    print("="*70)
    # result = workflow_with_helpers()
    
    print("\n" + "="*70)
    print("方案 C: 高层 Compute 节点（最简洁）")
    print("="*70)
    # result = workflow_with_compute()
    
    print("\n" + "="*70)
    print("方案 D: 带缓存的工作流")
    print("="*70)
    # result = workflow_with_caching()
