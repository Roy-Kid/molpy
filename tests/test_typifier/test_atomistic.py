import pytest

import molpy as mp
from molpy.core.atomistic import Angle, Atom, Atomistic, Bond, Dihedral
from molpy.io import read_xml_forcefield
from molpy.typifier.atomistic import (
    OplsAngleTypifier,
    OplsBondTypifier,
    OplsDihedralTypifier,
    OplsAtomisticTypifier,
)


@pytest.fixture
def oplsaa_forcefield():
    """加载 OPLS-AA 力场"""
    return read_xml_forcefield("oplsaa.xml")

@pytest.fixture
def tip3p_forcefield():
    """加载 TIP3P 水模型力场"""
    return read_xml_forcefield("tip3p.xml")

@pytest.fixture
def h2o():
    _h2o = Atomistic()
    o = Atom({'element': 'O', 'type': 'tip3p-O'})
    h1 = Atom({'element': 'H', 'type': 'tip3p-H'})
    h2 = Atom({'element': 'H', 'type': 'tip3p-H'})
    _h2o.atoms.extend([o, h1, h2])
    _h2o.bonds.append(Bond(o, h1))
    _h2o.angles.append(Angle(h1, o, h2))
    return _h2o


class TestTip3pTypifying:
    """测试 TIP3P 水模型参数匹配"""

    def test_tip3p_typifying(self, tip3p_forcefield, h2o):
        """测试 TIP3P 水模型参数匹配"""
        bond_typifier = OplsBondTypifier(tip3p_forcefield)
        angle_typifier = OplsAngleTypifier(tip3p_forcefield)

        # Typify bonds
        for bond in h2o.bonds:
            bond_typifier.typify(bond)
            assert 'type' in bond, "Bond type should be matched"
            assert bond['type'] == 'tip3p-O-tip3p-H', f"Expected bond type 'tip3p-O-tip3p-H', got {bond['type']}"

        # Typify angles
        for angle in h2o.angles:
            angle_typifier.typify(angle)
            assert 'type' in angle, "Angle type should be matched"
            assert angle['type'] == 'tip3p-H-tip3p-O-tip3p-H', f"Expected angle type 'tip3p-H-tip3p-O-tip3p-H', got {angle['type']}"


@pytest.fixture
def simple_ether():
    """构建简单的醚分子：CH3-O-CH3
    
    结构：
        H       H
        |       |
    H - C - O - C - H
        |       |
        H       H
    
    Atom types (OPLS):
    - C: opls_135 (CT class, alkane CH3)
    - O: opls_180 (OS class, ether oxygen)
    - H: opls_140 (HC class, alkane H)
    """
    struct = Atomistic()
    
    # 创建原子
    # Left CH3
    c1 = Atom({'element': 'C', 'type': 'opls_135'})
    h1 = Atom({'element': 'H', 'type': 'opls_140'})
    h2 = Atom({'element': 'H', 'type': 'opls_140'})
    h3 = Atom({'element': 'H', 'type': 'opls_140'})
    
    # Oxygen
    o = Atom({'element': 'O', 'type': 'opls_180'})
    
    # Right CH3
    c2 = Atom({'element': 'C', 'type': 'opls_135'})
    h4 = Atom({'element': 'H', 'type': 'opls_140'})
    h5 = Atom({'element': 'H', 'type': 'opls_140'})
    h6 = Atom({'element': 'H', 'type': 'opls_140'})
    
    # 添加到 assembly
    atoms = [c1, h1, h2, h3, o, c2, h4, h5, h6]
    for atom in atoms:
        struct.atoms.append(atom)
    
    # 创建键
    # C1-H bonds
    struct.bonds.append(Bond(c1, h1))
    struct.bonds.append(Bond(c1, h2))
    struct.bonds.append(Bond(c1, h3))
    
    # C1-O bond (关键！)
    struct.bonds.append(Bond(c1, o))
    
    # O-C2 bond (关键！)
    struct.bonds.append(Bond(o, c2))
    
    # C2-H bonds
    struct.bonds.append(Bond(c2, h4))
    struct.bonds.append(Bond(c2, h5))
    struct.bonds.append(Bond(c2, h6))
    
    return struct

class TestOplsTypifying:

    def test_bond(self, oplsaa_forcefield):

        bond = mp.Bond(Atom({'type': 'opls_135', 'name': 'C'}), Atom({'type': 'opls_140', 'name': 'H'}))
        bond_typifier = OplsBondTypifier(oplsaa_forcefield)
        bond_typifier.typify(bond)
        assert bond['type'] == 'CT-HC'

    def test_angle(self, oplsaa_forcefield):

        angle = mp.Angle(
            mp.Atom({'type': 'opls_140', 'name': 'H'}),
            mp.Atom({'type': 'opls_135', 'name': 'C'}),
            mp.Atom({'type': 'opls_180', 'name': 'O'})
        )
        angle_typifier = OplsAngleTypifier(oplsaa_forcefield)
        angle_typifier.typify(angle)
        assert angle['type'] == 'HC-CT-OS'

    def test_monomer(self, oplsaa_forcefield, simple_ether):

        typifier = OplsAtomisticTypifier(oplsaa_forcefield)
        typifier.typify(simple_ether)
        for bond in simple_ether.bonds:
            assert 'type' in bond, "Bond should be typed"
        for angle in simple_ether.angles:
            assert 'type' in angle, "Angle should be typed"