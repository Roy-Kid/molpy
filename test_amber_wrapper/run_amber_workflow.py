#!/usr/bin/env python
"""AMBER 工作流模板 - 从 PDB 到 Frame

使用方法：
    python run_amber_workflow.py input.pdb --charge -1 --output output_dir
"""

import argparse
import sys
from pathlib import Path

from molpy.wrapper import AntechamberWrapper, Parmchk2Wrapper, TLeapWrapper
from molpy.io.readers import read_amber_prmtop


def run_amber_workflow(
    input_pdb: Path,
    output_dir: Path,
    net_charge: int = 0,
    charge_method: str = "bcc",
    atom_type: str = "gaff2",
    add_counterion: str | None = None,
):
    """运行完整的 AMBER 参数化工作流。
    
    Args:
        input_pdb: 输入 PDB 文件路径
        output_dir: 输出目录
        net_charge: 分子净电荷
        charge_method: 电荷计算方法（bcc, gas, esp 等）
        atom_type: 原子类型方案（gaff2, gaff 等）
        add_counterion: 可选的反离子（如 "Li+", "Na+", "Cl-"）
    
    Returns:
        (Frame, ForceField) 元组
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    mol_name = input_pdb.stem
    
    print("=" * 70)
    print(f"AMBER 参数化工作流: {input_pdb.name}")
    print("=" * 70)
    print(f"  净电荷: {net_charge}")
    print(f"  电荷方法: {charge_method}")
    print(f"  原子类型: {atom_type}")
    if add_counterion:
        print(f"  反离子: {add_counterion}")
    print()
    
    # ========== Stage 1: Antechamber ==========
    print("[1/4] 运行 Antechamber...")
    print("-" * 70)
    
    ante_dir = output_dir / "step1_antechamber"
    ante = AntechamberWrapper(name="antechamber", workdir=ante_dir)
    
    # 检查可用性
    try:
        ante.check()
    except FileNotFoundError as e:
        print(f"✗ 错误: {e}")
        print("  请确保 antechamber 已安装且在 PATH 中")
        sys.exit(1)
    
    # 运行原子类型分配
    result = ante.atomtype_assign(
        input_file=input_pdb,
        output_file=f"{mol_name}.mol2",
        input_format="pdb",
        output_format="mol2",
        charge_method=charge_method,
        atom_type=atom_type,
        net_charge=net_charge,
    )
    
    if result.returncode != 0:
        print(f"✗ Antechamber 失败:")
        print(result.stderr)
        sys.exit(1)
    
    mol2_file = ante_dir / f"{mol_name}.mol2"
    print(f"✓ 完成: {mol2_file.relative_to(output_dir)}")
    print()
    
    # ========== Stage 2: Parmchk2 ==========
    print("[2/4] 运行 Parmchk2...")
    print("-" * 70)
    
    parmchk_dir = output_dir / "step2_parmchk2"
    parmchk = Parmchk2Wrapper(name="parmchk2", workdir=parmchk_dir)
    
    try:
        parmchk.check()
    except FileNotFoundError as e:
        print(f"✗ 错误: {e}")
        print("  请确保 parmchk2 已安装且在 PATH 中")
        sys.exit(1)
    
    # 生成缺失参数
    result = parmchk.generate_parameters(
        input_file=mol2_file,
        output_file=f"{mol_name}.frcmod",
        input_format="mol2",
        parameter_level=2,
    )
    
    if result.returncode != 0:
        print(f"✗ Parmchk2 失败:")
        print(result.stderr)
        sys.exit(1)
    
    frcmod_file = parmchk_dir / f"{mol_name}.frcmod"
    print(f"✓ 完成: {frcmod_file.relative_to(output_dir)}")
    print()
    
    # ========== Stage 3: Tleap ==========
    print("[3/4] 运行 Tleap...")
    print("-" * 70)
    
    tleap_dir = output_dir / "step3_tleap"
    tleap = TLeapWrapper(name="tleap", workdir=tleap_dir)
    
    try:
        tleap.check()
    except FileNotFoundError as e:
        print(f"✗ 错误: {e}")
        print("  请确保 tleap 已安装且在 PATH 中")
        sys.exit(1)
    
    # 构建 tleap 脚本
    system_name = mol_name
    if add_counterion:
        system_name = f"{mol_name}_ion"
    
    if add_counterion:
        tleap_script = f"""# 加载力场
source leaprc.gaff2

# 加载自定义参数
loadAmberParams {frcmod_file}

# 加载分子
MOL = loadMol2 {mol2_file}

# 添加反离子
ION = {add_counterion}

# 构建复合物
system = combine {{ MOL ION }}

# 保存参数化系统
saveAmberParm system {system_name}.prmtop {system_name}.inpcrd

quit
"""
    else:
        tleap_script = f"""# 加载力场
source leaprc.gaff2

# 加载自定义参数
loadAmberParams {frcmod_file}

# 加载分子
MOL = loadMol2 {mol2_file}

# 保存参数化系统
saveAmberParm MOL {system_name}.prmtop {system_name}.inpcrd

quit
"""
    
    result = tleap.run_from_script(
        script_text=tleap_script,
        script_name="tleap.in",
    )
    
    if result.returncode != 0:
        print(f"✗ Tleap 失败:")
        print(result.stderr)
        sys.exit(1)
    
    prmtop_file = tleap_dir / f"{system_name}.prmtop"
    inpcrd_file = tleap_dir / f"{system_name}.inpcrd"
    print(f"✓ 完成: {prmtop_file.relative_to(output_dir)}")
    print(f"        {inpcrd_file.relative_to(output_dir)}")
    print()
    
    # ========== Stage 4: 读取为 Frame ==========
    print("[4/4] 读取为 Frame...")
    print("-" * 70)
    
    frame, forcefield = read_amber_prmtop(prmtop_file, inpcrd_file)
    
    print(f"✓ Frame: {len(frame['atoms'])} 个原子")
    print(f"✓ ForceField: {len(forcefield.atom_types)} 种原子类型")
    print()
    
    print("=" * 70)
    print("✓ 工作流完成！")
    print("=" * 70)
    print(f"\n输出文件位于: {output_dir.absolute()}")
    
    return frame, forcefield


def main():
    parser = argparse.ArgumentParser(
        description="AMBER 参数化工作流：从 PDB 到 Frame",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例:
  # 中性小分子
  python run_amber_workflow.py ligand.pdb -o output/
  
  # 带电荷的离子
  python run_amber_workflow.py tfsi.pdb -c -1 -o output/
  
  # 添加反离子
  python run_amber_workflow.py tfsi.pdb -c -1 --counterion "Li+" -o output/
        """,
    )
    
    parser.add_argument(
        "input_pdb",
        type=Path,
        help="输入 PDB 文件",
    )
    
    parser.add_argument(
        "-o", "--output",
        type=Path,
        default=Path("amber_output"),
        help="输出目录（默认: amber_output）",
    )
    
    parser.add_argument(
        "-c", "--charge",
        type=int,
        default=0,
        help="分子净电荷（默认: 0）",
    )
    
    parser.add_argument(
        "--charge-method",
        choices=["gas", "bcc", "be3", "cm2", "esp"],
        default="bcc",
        help="电荷计算方法（默认: bcc）",
    )
    
    parser.add_argument(
        "--atom-type",
        choices=["gaff", "gaff2", "amber", "sybyl"],
        default="gaff2",
        help="原子类型方案（默认: gaff2）",
    )
    
    parser.add_argument(
        "--counterion",
        type=str,
        help='添加反离子（如 "Li+", "Na+", "Cl-"）',
    )
    
    args = parser.parse_args()
    
    # 检查输入文件
    if not args.input_pdb.exists():
        print(f"✗ 错误: 输入文件不存在: {args.input_pdb}")
        sys.exit(1)
    
    # 运行工作流
    try:
        frame, forcefield = run_amber_workflow(
            input_pdb=args.input_pdb,
            output_dir=args.output,
            net_charge=args.charge,
            charge_method=args.charge_method,
            atom_type=args.atom_type,
            add_counterion=args.counterion,
        )
    except Exception as e:
        print(f"\n✗ 工作流失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
