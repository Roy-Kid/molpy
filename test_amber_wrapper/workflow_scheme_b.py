"""方案 B: 使用 Wrapper 辅助方法进行 AMBER 工作流

这个脚本演示了如何使用改进的 Wrapper API（方案 B），
从 PDB 文件一步步到 Frame 和 ForceField。
"""

from pathlib import Path

from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper
from molpy.io.readers import read_amber_prmtop


def workflow_scheme_b():
    """方案 B：使用 Wrapper 辅助方法，参数完整、易于使用。"""
    
    # 准备工作目录
    workdir = Path("./amber_work_scheme_b")
    workdir.mkdir(exist_ok=True)
    
    print("=" * 70)
    print("AMBER 工作流 - 方案 B（中层 API）")
    print("=" * 70)
    
    # ========= Stage 1: Antechamber - PDB → MOL2 =========
    print("\n[Stage 1] 原子类型分配和电荷计算 (antechamber)")
    print("-" * 70)
    
    ante = AntechamberWrapper(
        name="antechamber",
        workdir=workdir / "stage1_antechamber",
    )
    
    # 检查工具可用性
    ante.check()
    print(f"✓ Antechamber 可用")
    
    # 运行原子类型分配
    result = ante.atomtype_assign(
        input_file="tfsi.pdb",
        output_file="tfsi_gaff2.mol2",
        input_format="pdb",
        output_format="mol2",
        charge_method="bcc",  # Bond charge correction
        atom_type="gaff2",
        net_charge=-1,
    )
    
    if result.returncode != 0:
        print(f"✗ Antechamber 失败:")
        print(result.stderr)
        return None
    
    mol2_file = workdir / "stage1_antechamber" / "tfsi_gaff2.mol2"
    print(f"✓ 生成: {mol2_file}")
    print(f"  包含: 原子类型、部分电荷、键信息")
    
    # ========= Stage 2: Parmchk2 - MOL2 → FRCMOD =========
    print("\n[Stage 2] 参数检查和补充 (parmchk2)")
    print("-" * 70)
    
    parmchk = Parmchk2Wrapper(
        name="parmchk2",
        workdir=workdir / "stage2_parmchk",
    )
    
    parmchk.check()
    print(f"✓ Parmchk2 可用")
    
    # 运行参数检查和生成
    result = parmchk.generate_parameters(
        input_file=mol2_file,
        output_file="tfsi.frcmod",
        input_format="mol2",
        parameter_level=2,  # level 2: all parameters (recommended)
    )
    
    if result.returncode != 0:
        print(f"✗ Parmchk2 失败:")
        print(result.stderr)
        return None
    
    frcmod_file = workdir / "stage2_parmchk" / "tfsi.frcmod"
    print(f"✓ 生成: {frcmod_file}")
    print(f"  包含: 缺失的键、角、二面角参数")
    
    # ========= Stage 3: Tleap - 系统构建 =========
    print("\n[Stage 3] 系统构建和参数化 (tleap)")
    print("-" * 70)
    
    tleap = TLeapWrapper(
        name="tleap",
        workdir=workdir / "stage3_tleap",
    )
    
    tleap.check()
    print(f"✓ Tleap 可用")
    
    # 构建 tleap 脚本
    # 注意：路径需要相对于 tleap 的工作目录调整
    tleap_script = f"""# 加载力场参数
source leaprc.gaff2

# 加载自定义参数
loadAmberParams {frcmod_file}

# 加载配体分子
TFSI = loadMol2 {mol2_file}

# 定义反离子
Li = Li+

# 构建复合物
complex = combine {{ Li TFSI }}

# 保存参数化的系统
saveAmberParm complex litfsi.prmtop litfsi.inpcrd

quit
"""
    
    print("执行 tleap...")
    result = tleap.run_from_script(
        script_text=tleap_script,
        script_name="tleap.in",
    )
    
    if result.returncode != 0:
        print(f"✗ Tleap 失败:")
        print(result.stderr)
        return None
    
    prmtop_file = workdir / "stage3_tleap" / "litfsi.prmtop"
    inpcrd_file = workdir / "stage3_tleap" / "litfsi.inpcrd"
    print(f"✓ 生成: {prmtop_file}, {inpcrd_file}")
    print(f"  包含: 完整的参数化系统和初始坐标")
    
    # ========= Stage 4: 读取为 Frame =========
    print("\n[Stage 4] 读取到 Frame")
    print("-" * 70)
    
    frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
    
    print(f"✓ 成功加载:")
    print(f"  - 原子数: {len(frame['atoms'])}")
    print(f"  - 坐标范围: {frame['atoms']['x'].min():.2f} - {frame['atoms']['x'].max():.2f} Å")
    print(f"  - ForceField: {len(forcefield.atom_types)} 种原子类型")
    
    print("\n" + "=" * 70)
    print("工作流完成！")
    print("=" * 70)
    
    return frame, forcefield


if __name__ == "__main__":
    try:
        result = workflow_scheme_b()
        if result:
            frame, ff = result
            print(f"\nFrame 对象: {type(frame)}")
            print(f"ForceField 对象: {type(ff)}")
    except Exception as e:
        print(f"\n错误: {e}")
        import traceback
        traceback.print_exc()
