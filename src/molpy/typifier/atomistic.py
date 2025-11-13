
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
from molpy.core.entity import StructLike
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


class OplsAtomTypifier(TypifierBase["Atomistic"]):
    """使用 SMARTS matcher 为原子分配类型（支持类型引用和依赖解析）"""

    def __init__(self, forcefield: ForceField) -> None:
        super().__init__(forcefield)
        from .adapter import build_mol_graph
        
        # 从 forcefield 中提取 patterns
        self.pattern_dict = self._extract_patterns()
        self._build_mol_graph = build_mol_graph
        
        # 使用 LayeredTypingEngine
        from .layered_engine import LayeredTypingEngine
        self.engine = LayeredTypingEngine(self.pattern_dict)

    def _extract_patterns(self):
        """从 forcefield 提取或构造 SMARTS patterns
        
        从 OPLS 力场的 AtomType 中提取 SMARTS 定义（def 属性）并转换为 SMARTSGraph 对象。
        支持 overrides 属性来控制优先级和类型引用（%opls_XXX）。
        
        Returns:
            Dictionary mapping atom type name to SMARTSGraph
        """
        from .graph import SMARTSGraph
        from molpy.parser.smarts import SmartsParser
        
        pattern_dict = {}
        atom_types = list(self.ff.get_types(AtomType))
        parser = SmartsParser()
        
        # 构建 overrides 映射
        overrides_map = {}
        for at in atom_types:
            overrides_str = at.params.kwargs.get("overrides")
            if overrides_str:
                overrides_map[at.name] = {s.strip() for s in overrides_str.split(",")}
        
        # 计算优先级：基于 overrides 和 priority 属性
        type_priority = {}
        for at in atom_types:
            # 首先使用显式 priority（如果有）
            explicit_priority = at.params.kwargs.get("priority")
            if explicit_priority is not None:
                try:
                    type_priority[at.name] = int(explicit_priority)
                    continue
                except (ValueError, TypeError):
                    pass
            
            # 否则基于 overrides 计算
            priority = 0
            # 如果这个类型被其他类型 override，降低优先级
            for overrider, overridden_set in overrides_map.items():
                if at.name in overridden_set:
                    priority -= 1
            # 如果这个类型 override 其他类型，提高优先级
            if at.name in overrides_map:
                priority += len(overrides_map[at.name])
            type_priority[at.name] = priority
        
        # 提取 SMARTS patterns
        for at in atom_types:
            smarts_str = at.params.kwargs.get("def_")
            
            if smarts_str:
                # 使用 SMARTSGraph 解析 SMARTS 字符串
                try:
                    priority = type_priority.get(at.name, 0)
                    overrides = overrides_map.get(at.name, set())
                    
                    pattern = SMARTSGraph(
                        smarts_string=smarts_str,
                        parser=parser,
                        atomtype_name=at.name,
                        priority=priority,
                        source=f"oplsaa:{at.name}",
                        overrides=overrides
                    )
                    pattern_dict[at.name] = pattern
                except Exception as e:
                    # 如果解析失败，记录警告但继续
                    import warnings
                    warnings.warn(f"Failed to parse SMARTS for {at.name}: {smarts_str}, error: {e}")
        
        return pattern_dict

    @override
    def typify(self, struct: "Atomistic") -> "Atomistic":
        """为 Atomistic 结构中的所有原子分配类型（使用依赖感知的分层匹配）"""
        # 将分子转换为图
        graph, vs_to_atomid, atomid_to_vs = self._build_mol_graph(struct)
        
        # 使用 LayeredTypingEngine 进行分层匹配
        result = self.engine.typify(graph, vs_to_atomid)
        
        # 将结果应用到原子上
        for atom in struct.atoms:
            atom_id = id(atom)
            if atom_id in result:
                atomtype = result[atom_id]
                atom.data["type"] = atomtype
                
                # 同时从 forcefield 中获取其他参数
                atom_type_obj = self._find_atomtype_by_name(atomtype)
                if atom_type_obj:
                    atom.data.update(**atom_type_obj.params.kwargs)
        
        return struct

    def _find_atomtype_by_name(self, name: str) -> AtomType | None:
        """根据名称查找 AtomType 对象"""
        for at in self.ff.get_types(AtomType):
            if at.name == name:
                return at
        return None


class OplsAtomisticTypifier(TypifierBase[Atomistic]):
    """为整个 Atomistic 结构分配所有类型（bond, angle, dihedral）
    
    注意：此类假设原子已经被分配了类型。如果需要同时分配原子类型，
    请先使用 OplsAtomTypifier，或使用 skip_atom_typing=False 参数。
    """

    def __init__(self, forcefield: ForceField, skip_atom_typing: bool = True) -> None:
        super().__init__(forcefield)
        self.skip_atom_typing = skip_atom_typing
        if not skip_atom_typing:
            self.atom_typifier = OplsAtomTypifier(forcefield)
        self.bond_typifier = OplsBondTypifier(forcefield)
        self.angle_typifier = OplsAngleTypifier(forcefield)
        self.dihedral_typifier = OplsDihedralTypifier(forcefield)

    @override
    def typify(self, struct: Atomistic) -> Atomistic:
        """
        为 Atomistic 结构中的所有 bonds, angles, dihedrals 分配类型
        
        参数：
            struct: Atomistic 结构
        
        前提：
            - 如果 skip_atom_typing=True（默认），所有 atoms 必须已经有 'type' 属性
            - 如果 skip_atom_typing=False，将先为原子分配类型
        """
        # 可选：首先为原子分配类型
        if not self.skip_atom_typing:
            self.atom_typifier.typify(struct)


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