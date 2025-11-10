
# pyright: reportIncompatibleMethodOverride=false
from abc import ABC, abstractmethod
from typing import Iterable, TypeVar, override
from molpy.core.forcefield import (
    ForceField,
    AtomType,
    BondType,
    AngleType,
    DihedralType,
)
from molpy.core.entity import AssemblyLike
from molpy.core.atomistic import Bond, Angle, Dihedral, Atomistic
from collections import defaultdict


def atomtype_matches(atomtype: AtomType, type_str: str) -> bool:
    """
    判断 AtomType 是否匹配给定的 type 字符串
    
    匹配规则：
    1. 如果 atomtype 有具体的 type（不是 "*"），则比较 type
    2. 如果 type不匹配，再比较 class
    
    Args:
        atomtype: AtomType 实例
        type_str: 要匹配的类型字符串（来自 Atom.data["type"] 或 class 名）
        
    Returns:
        是否匹配
    """
    at_type = atomtype.params.kwargs.get("type_", "*")
    at_class = atomtype.params.kwargs.get("class_", "*")
    
    # 优先匹配 type
    if at_type != "*" and at_type == type_str:
        return True
    
    # 其次匹配 class
    if at_class != "*" and at_class == type_str:
        return True
    
    # 如果两者都是通配符，也匹配
    if at_type == "*" and at_class == "*":
        return True
    
    return False


class TypifierBase[T](ABC):
    def __init__(self, forcefield: ForceField) -> None:
        self.ff = forcefield

    @abstractmethod
    def typify(self, elem: T) -> T:
        ...


class OplsBondTypifier(TypifierBase[Bond]):
    """根据 Bond 两端原子的 atom type 匹配 bond type
    
    策略：建立 class_to_types 表，手动比对每个 AtomType
    """

    def __init__(self, forcefield: ForceField) -> None:
        super().__init__(forcefield)
        self._build_table()
    
    def _build_table(self):
        """构建 class_to_types 表和 bond 表"""
        # 建立 class -> types 映射
        self.class_to_types: dict[str, list[str]] = defaultdict(list)
        for at in self.ff.get_types(AtomType):
            at_type = at.params.kwargs.get("type_", "*")
            at_class = at.params.kwargs.get("class_", "*")
            if at_class != "*":
                if at_type != "*":
                    self.class_to_types[at_class].append(at_type)
                else:
                    self.class_to_types[at_class].append(at_class)
        
        # 建立 bond 表
        self._bond_table = {}
        for bond in self.ff.get_types(BondType):
            self._bond_table[(bond.itom, bond.jtom)] = bond

    @override
    def typify(self, bond: Bond) -> Bond:
        """为键分配类型"""
        itom_type = bond.itom.get("type", None)
        jtom_type = bond.jtom.get("type", None)
        
        if itom_type is None or jtom_type is None:
            raise ValueError(f"Bond atoms must have 'type' attribute: {bond}")
        
        # 遍历所有 bond types，手动匹配
        for (at1, at2), bond_type in self._bond_table.items():
            # 尝试正向和反向匹配
            if (atomtype_matches(at1, itom_type) and atomtype_matches(at2, jtom_type)) or \
               (atomtype_matches(at1, jtom_type) and atomtype_matches(at2, itom_type)):
                bond.data["type"] = bond_type.name
                bond.data.update(**bond_type.params.kwargs)
                return bond
        
        # 没找到，尝试 class 匹配
        # 找到 itom_type 和 jtom_type 对应的 class
        itom_class = None
        jtom_class = None
        for cls, types in self.class_to_types.items():
            if itom_type in types:
                itom_class = cls
            if jtom_type in types:
                jtom_class = cls
        
        if itom_class and jtom_class:
            for (at1, at2), bond_type in self._bond_table.items():
                if (atomtype_matches(at1, itom_class) and atomtype_matches(at2, jtom_class)) or \
                   (atomtype_matches(at1, jtom_class) and atomtype_matches(at2, itom_class)):
                    bond.data["type"] = bond_type.name
                    bond.data.update(**bond_type.params.kwargs)
                    return bond
        
        raise ValueError(f"No bond type found for atom types: {itom_type} - {jtom_type}")


class OplsAngleTypifier(TypifierBase[Angle]):
    """根据 Angle 三个原子的 atom type 匹配 angle type"""

    def __init__(self, forcefield: ForceField) -> None:
        super().__init__(forcefield)
        self._build_table()

    def _build_table(self) -> None:
        """构建 class_to_types 表和 angle 表"""
        # 建立 class -> types 映射
        self.class_to_types: dict[str, list[str]] = defaultdict(list)
        for at in self.ff.get_types(AtomType):
            at_type = at.params.kwargs.get("type_", "*")
            at_class = at.params.kwargs.get("class_", "*")
            if at_class != "*":
                if at_type != "*":
                    self.class_to_types[at_class].append(at_type)
                else:
                    self.class_to_types[at_class].append(at_class)
        
        # 建立 angle 表
        self._angle_table = {}
        for angle in self.ff.get_types(AngleType):
            self._angle_table[(angle.itom, angle.jtom, angle.ktom)] = angle

    @override
    def typify(self, angle: Angle) -> Angle:
        """为角分配类型"""
        itom_type = angle.itom.get("type", None)
        jtom_type = angle.jtom.get("type", None)
        ktom_type = angle.ktom.get("type", None)
        
        if None in (itom_type, jtom_type, ktom_type):
            raise ValueError(f"Angle atoms must have 'type' attribute: {angle}")
        
        assert isinstance(itom_type, str)
        assert isinstance(jtom_type, str)
        assert isinstance(ktom_type, str)
        
        # 遍历所有 angle types，手动匹配
        for (at1, at2, at3), angle_type in self._angle_table.items():
            # 尝试正向和反向匹配（中心原子 at2 不变）
            if (atomtype_matches(at1, itom_type) and 
                atomtype_matches(at2, jtom_type) and 
                atomtype_matches(at3, ktom_type)) or \
               (atomtype_matches(at1, ktom_type) and 
                atomtype_matches(at2, jtom_type) and 
                atomtype_matches(at3, itom_type)):
                angle.data["type"] = angle_type.name
                angle.data.update(**angle_type.params.kwargs)
                return angle
        
        # 没找到，尝试 class 匹配
        itom_class = None
        jtom_class = None
        ktom_class = None
        for cls, types in self.class_to_types.items():
            if itom_type in types:
                itom_class = cls
            if jtom_type in types:
                jtom_class = cls
            if ktom_type in types:
                ktom_class = cls
        
        if itom_class and jtom_class and ktom_class:
            for (at1, at2, at3), angle_type in self._angle_table.items():
                if (atomtype_matches(at1, itom_class) and 
                    atomtype_matches(at2, jtom_class) and 
                    atomtype_matches(at3, ktom_class)) or \
                   (atomtype_matches(at1, ktom_class) and 
                    atomtype_matches(at2, jtom_class) and 
                    atomtype_matches(at3, itom_class)):
                    angle.data["type"] = angle_type.name
                    angle.data.update(**angle_type.params.kwargs)
                    return angle
        
        raise ValueError(
            f"No angle type found for atom types: {itom_type} - {jtom_type} - {ktom_type}"
        )


class OplsDihedralTypifier(TypifierBase[Dihedral]):
    """根据 Dihedral 四个原子的 atom type 匹配 dihedral type"""

    def __init__(self, forcefield: ForceField) -> None:
        super().__init__(forcefield)
        self._build_table()

    def _build_table(self) -> None:
        """构建 class_to_types 表和 dihedral 列表"""
        # 建立 class -> types 映射
        self.class_to_types: dict[str, list[str]] = defaultdict(list)
        for at in self.ff.get_types(AtomType):
            at_type = at.params.kwargs.get("type_", "*")
            at_class = at.params.kwargs.get("class_", "*")
            if at_class != "*":
                if at_type != "*":
                    self.class_to_types[at_class].append(at_type)
                else:
                    self.class_to_types[at_class].append(at_class)
        
        # 建立 dihedral 列表（不是字典！因为多个 dihedral 可能有相同的 AtomType 组合）
        self._dihedral_list: list[DihedralType] = list(self.ff.get_types(DihedralType))

    @override
    def typify(self, dihedral: Dihedral) -> Dihedral:
        """为二面角分配类型"""
        itom_type = dihedral.itom.get("type", None)
        jtom_type = dihedral.jtom.get("type", None)
        ktom_type = dihedral.ktom.get("type", None)
        ltom_type = dihedral.ltom.get("type", None)
        
        if None in (itom_type, jtom_type, ktom_type, ltom_type):
            raise ValueError(f"Dihedral atoms must have 'type' attribute: {dihedral}")
        
        assert isinstance(itom_type, str)
        assert isinstance(jtom_type, str)
        assert isinstance(ktom_type, str)
        assert isinstance(ltom_type, str)
        
        # 遍历所有 dihedral types，手动匹配
        for dihedral_type in self._dihedral_list:
            at1, at2, at3, at4 = dihedral_type.itom, dihedral_type.jtom, dihedral_type.ktom, dihedral_type.ltom
            # 尝试正向和反向匹配
            if (atomtype_matches(at1, itom_type) and 
                atomtype_matches(at2, jtom_type) and 
                atomtype_matches(at3, ktom_type) and 
                atomtype_matches(at4, ltom_type)) or \
               (atomtype_matches(at1, ltom_type) and 
                atomtype_matches(at2, ktom_type) and 
                atomtype_matches(at3, jtom_type) and 
                atomtype_matches(at4, itom_type)):
                dihedral.data["type"] = dihedral_type.name
                dihedral.data.update(**dihedral_type.params.kwargs)
                return dihedral
        
        # 没找到，尝试 class 匹配
        itom_class = None
        jtom_class = None
        ktom_class = None
        ltom_class = None
        for cls, types in self.class_to_types.items():
            if itom_type in types:
                itom_class = cls
            if jtom_type in types:
                jtom_class = cls
            if ktom_type in types:
                ktom_class = cls
            if ltom_type in types:
                ltom_class = cls
        
        if itom_class and jtom_class and ktom_class and ltom_class:
            for dihedral_type in self._dihedral_list:
                at1, at2, at3, at4 = dihedral_type.itom, dihedral_type.jtom, dihedral_type.ktom, dihedral_type.ltom
                if (atomtype_matches(at1, itom_class) and 
                    atomtype_matches(at2, jtom_class) and 
                    atomtype_matches(at3, ktom_class) and 
                    atomtype_matches(at4, ltom_class)) or \
                   (atomtype_matches(at1, ltom_class) and 
                    atomtype_matches(at2, ktom_class) and 
                    atomtype_matches(at3, jtom_class) and 
                    atomtype_matches(at4, itom_class)):
                    dihedral.data["type"] = dihedral_type.name
                    dihedral.data.update(**dihedral_type.params.kwargs)
                    return dihedral
        
        raise ValueError(
            f"No dihedral type found for atom types: {itom_type} - {jtom_type} - {ktom_type} - {ltom_type}"
        )


class OplsAtomisticTypifier(TypifierBase[Atomistic]):
    """为整个 Atomistic 结构分配所有类型（bond, angle, dihedral）"""

    def __init__(self, forcefield: ForceField) -> None:
        super().__init__(forcefield)
        self.bond_typifier = OplsBondTypifier(forcefield)
        self.angle_typifier = OplsAngleTypifier(forcefield)
        self.dihedral_typifier = OplsDihedralTypifier(forcefield)

    @override
    def typify(self, struct: Atomistic) -> Atomistic:
        """
        为 Atomistic 结构中的所有 bonds, angles, dihedrals 分配类型
        
        前提：所有 atoms 已经有 'type' 属性
        """
        # 为所有键分配类型
        for bond in struct.bonds:
            self.bond_typifier.typify(bond)
        
        # 为所有角分配类型（如果存在）
        angles = struct.links.bucket(Angle)
        for angle in angles:
            self.angle_typifier.typify(angle)
        
        # 为所有二面角分配类型（如果存在）
        dihedrals = struct.links.bucket(Dihedral)
        for dihedral in dihedrals:
            self.dihedral_typifier.typify(dihedral)
        
        return struct