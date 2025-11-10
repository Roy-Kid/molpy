"""测试 XML 力场文件读取功能"""

import pytest
from pathlib import Path
from molpy.io.forcefield.xml import read_xml_forcefield
from molpy.core.forcefield import AtomType, BondType, AngleType


class TestTIP3PForceField:
    """使用 TIP3P 水模型测试力场读取功能"""
    
    def test_atomtype_stores_class_info(self):
        """测试 AtomType 是否正确存储 class 信息"""
        ff = read_xml_forcefield("tip3p.xml")
        
        # 获取所有 atom types
        atom_types = ff.get_types(AtomType)
        assert len(atom_types) == 2, "TIP3P should have 2 atom types"
        
        # 查找 tip3p-O
        o_type = None
        h_type = None
        for at in atom_types:
            if at.name == "tip3p-O":
                o_type = at
            elif at.name == "tip3p-H":
                h_type = at
        
        assert o_type is not None, "Should find tip3p-O atom type"
        assert h_type is not None, "Should find tip3p-H atom type"
        
        # 检查 class_name 属性（现在是 AtomType 的属性，而不是 params.kwargs）
        assert o_type.class_name == 'tip3p-O', "AtomType should have class_name property"
        assert h_type.class_name == 'tip3p-H'
        
        # 检查 element 和 mass 信息仍在 params.kwargs 中
        assert 'element' in o_type.params.kwargs
        assert o_type.params.kwargs['element'] == 'O'
        assert 'mass' in o_type.params.kwargs
        assert abs(o_type.params.kwargs['mass'] - 15.99943) < 0.0001
    
    def test_xmlreader_caches_by_name_and_class(self):
        """测试 XMLReader 是否同时用 name 和 class 缓存 AtomType"""
        from molpy.io.forcefield.xml import XMLForceFieldReader
        
        reader = XMLForceFieldReader("tip3p.xml")
        ff = reader.read()
        
        # 检查缓存中是否同时有 name 和 class 作为键
        assert 'tip3p-O' in reader._atomtype_cache, "Should cache by name"
        assert 'tip3p-H' in reader._atomtype_cache, "Should cache by name"
        
        # 对于 TIP3P，name 和 class 相同，所以应该指向同一个对象
        assert reader._atomtype_cache['tip3p-O'].name == 'tip3p-O'
        assert reader._atomtype_cache['tip3p-H'].name == 'tip3p-H'
    
    def test_bond_types_reference_correct_atomtypes(self):
        """测试 BondType 是否正确引用 AtomType"""
        ff = read_xml_forcefield("tip3p.xml")
        
        bond_types = ff.get_types(BondType)
        assert len(bond_types) == 1, "TIP3P should have 1 bond type"
        
        bond_type = bond_types[0]
        
        # 检查 bond type 的参数
        assert len(bond_type.params.args) >= 2, "Bond should have at least 2 atom types"
        
        at1, at2 = bond_type.params.args[0], bond_type.params.args[1]
        assert isinstance(at1, AtomType), "First param should be AtomType"
        assert isinstance(at2, AtomType), "Second param should be AtomType"
        
        # 检查是否是正确的 atom types (O-H bond)
        atom_names = {at1.name, at2.name}
        assert atom_names == {'tip3p-O', 'tip3p-H'}, f"Bond should be O-H, got {atom_names}"
        
        # 检查键参数
        assert 'r0' in bond_type.params.kwargs or 'length' in bond_type.params.kwargs
        assert 'k' in bond_type.params.kwargs


class TestTIP3PTypifier:
    """测试 Typifier 基于 class 匹配参数"""
    
    def test_bond_typifier_matches_by_class(self):
        """测试 BondTypifier 能够基于 atom class 匹配键类型"""
        from molpy.core.atomistic import Atom, Bond
        from molpy.typifier.atomistic import OplsBondTypifier
        
        # 加载 TIP3P 力场
        ff = read_xml_forcefield("tip3p.xml")
        
        # 创建 typifier
        bond_typifier = OplsBondTypifier(ff)
        
        # 检查 bond_type_map 是否正确构建
        assert len(bond_typifier.bond_type_map) > 0, "Bond type map should not be empty"
        
        # 对于 TIP3P，应该有 O-H 键的映射
        # 检查是否可以用 name 查找
        has_oh_bond = False
        for key in bond_typifier.bond_type_map.keys():
            if set(key) == {'tip3p-O', 'tip3p-H'}:
                has_oh_bond = True
                break
        
        assert has_oh_bond, "Should have O-H bond in type map"
        
        # 创建模拟的原子和键
        o_atom = Atom(element='O', type='tip3p-O')
        h_atom = Atom(element='H', type='tip3p-H')
        bond = Bond(o_atom, h_atom)
        
        # 测试 typify
        bond_typifier.typify(bond)
        
        # 检查键是否被正确标注类型
        assert 'bondtype' in bond.data, "Bond should have bondtype after typify"
        print(f"Bond type assigned: {bond.data['bondtype']}")


class TestOPLSAAClassHandling:
    """测试 OPLS-AA 力场的 class 处理"""
    
    def test_oplsaa_atomtype_class_storage(self):
        """测试 OPLS-AA atom types 是否正确存储 class 信息"""
        ff = read_xml_forcefield("oplsaa.xml")
        
        atom_types = ff.get_types(AtomType)
        assert len(atom_types) > 0, "OPLS-AA should have atom types"
        
        # 查找 opls_135 (alkane CH3, class=CT)
        opls_135 = None
        for at in atom_types:
            if at.name == "opls_135":
                opls_135 = at
                break
        
        assert opls_135 is not None, "Should find opls_135"
        assert 'class_name' in opls_135.params.kwargs, "Should store class_name"
        assert opls_135.params.kwargs['class_name'] == 'CT', f"opls_135 class should be CT, got {opls_135.params.kwargs.get('class_name')}"
    
    def test_oplsaa_bond_uses_class_for_matching(self):
        """测试 OPLS-AA bond 定义使用 class 进行匹配"""
        from molpy.io.forcefield.xml import XMLForceFieldReader
        
        reader = XMLForceFieldReader("oplsaa.xml")
        ff = reader.read()
        
        # 检查缓存中是否有 CT class
        assert 'CT' in reader._atomtype_cache, "Should cache atom types by class name"
        
        # CT 应该指向某个 opls_XXX type
        ct_atomtype = reader._atomtype_cache['CT']
        assert ct_atomtype.name.startswith('opls_'), f"CT should map to opls_XXX, got {ct_atomtype.name}"
        
        # 检查是否有 CT-CT 键
        bond_types = ff.get_types(BondType)
        
        # 查找使用 CT class 的键
        ct_bonds = []
        for bt in bond_types:
            if len(bt.params.args) >= 2:
                at1, at2 = bt.params.args[0], bt.params.args[1]
                if isinstance(at1, AtomType) and isinstance(at2, AtomType):
                    at1_class = at1.params.kwargs.get('class_name', at1.name)
                    at2_class = at2.params.kwargs.get('class_name', at2.name)
                    if at1_class == 'CT' and at2_class == 'CT':
                        ct_bonds.append(bt)
        
        assert len(ct_bonds) > 0, "Should find CT-CT bonds in OPLS-AA"


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
