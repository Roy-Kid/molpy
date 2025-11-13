"""
Reacter + PolymerBuilder 简单使用示例

演示如何使用新的 ReacterConnector API 构建聚合物。
"""

from molpy.parser.smiles import SmilesParser, bigsmilesir_to_monomer
from molpy.reacter import Reacter
from molpy.reacter.selectors import port_anchor_selector, remove_one_H
from molpy.reacter.transformers import make_single_bond
from molpy.builder.polymer import PolymerBuilder, ReacterConnector

print("=" * 70)
print("Reacter + PolymerBuilder 简单示例")
print("=" * 70)

# Step 1: 解析 Monomers
print("\n[Step 1] 解析 Monomers...")
parser = SmilesParser()

monomer_smiles = {
    "A": "CC[*:1]",      # 简单的乙基
    "B": "CCC[*:2]",     # 简单的丙基
}

monomers = {}
for label, smiles in monomer_smiles.items():
    ir = parser.parse_bigsmiles(smiles)
    monomer = bigsmilesir_to_monomer(ir)
    monomers[label] = monomer
    print(f"  Monomer {label}: {smiles} → {len(list(monomer.unwrap().atoms))} atoms")

# Step 2: 定义 Reacter
print("\n[Step 2] 定义化学反应...")
default_reacter = Reacter(
    name="C-C_coupling_with_H_loss",
    anchor_left=port_anchor_selector,
    anchor_right=port_anchor_selector,
    leaving_left=remove_one_H,
    leaving_right=remove_one_H,
    bond_maker=make_single_bond,
)
print(f"  Reacter: {default_reacter.name}")

# Step 3: 创建 ReacterConnector（使用显式端口选择策略）
print("\n[Step 3] 创建 ReacterConnector...")


# 定义端口选择函数
def first_available_ports(left, right, left_ports, right_ports, ctx):
    """选择每个单体的第一个可用端口"""
    if not left_ports or not right_ports:
        raise ValueError("No available ports")
    return list(left_ports.keys())[0], list(right_ports.keys())[0]


connector = ReacterConnector(
    default=default_reacter,
    port_strategy=first_available_ports,  # 使用自定义函数显式选择端口
)
print("  ✓ ReacterConnector 创建成功")
print("  ✓ 端口策略: first_available_ports (显式选择第一个可用端口)")

# Step 4: 使用 PolymerBuilder 组装
print("\n[Step 4] 组装 Polymer (序列: ABA)...")
try:
    polymer = PolymerBuilder.linear(
        sequence="ABA",
        library=monomers,
        connector=connector,  # 一个参数包含所有反应逻辑！
    )
    
    atoms = list(polymer.unwrap().atoms)
    bonds = list(polymer.unwrap().bonds)
    print(f"  ✓ 组装成功")
    print(f"    - 总原子数: {len(atoms)}")
    print(f"    - 总键数: {len(bonds)}")
    
    # 查看反应历史
    history = connector.get_history()
    print(f"\n[Step 5] 反应历史:")
    print(f"  总反应数: {len(history)}")
    for i, product in enumerate(history):
        removed_count = product.notes.get('removed_count', 0)
        new_bonds_count = len(product.notes.get('new_bonds', []))
        print(f"  Step {i+1}: 移除 {removed_count} 个原子, 创建 {new_bonds_count} 个新键")
    
except Exception as e:
    print(f"  ❌ 组装失败: {e}")
    import traceback
    traceback.print_exc()

print("\n" + "=" * 70)
print("✅ 示例完成")
print("=" * 70)
print("\n💡 关键优势:")
print("  1. 一个 connector 参数 = 端口选择 + 化学反应")
print("  2. Reacter 自动指导端口选择")
print("  3. 支持 default + overrides 模式")
print("  4. 反应历史自动记录")
