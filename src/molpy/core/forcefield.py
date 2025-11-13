from typing import Any, TypeVar, Self
from .entity import Entity
from .utils import TypeBucket

# --- 泛型变量 ---

S = TypeVar("S", bound="Style")
Ty = TypeVar("Ty", bound="Type")

# --- 基础组件 ---

class Parameters:
    def __init__(self, *args: Any, **kwargs: Any):
        self.args = list(args)
        self.kwargs = kwargs

    def __repr__(self) -> str:
        return f"Parameters(args={self.args}, kwargs={self.kwargs})"
    
    def __getitem__(self, key: int | slice | str):
        if isinstance(key, str):
            return self.kwargs[key]
        else:
            return self.args[key]

class Type:
    def __init__(self, name: str, *args: Any, **kwargs: Any):
        self._name = name
        self.params = Parameters(*args, **kwargs)

    @property
    def name(self) -> str:
        return self._name

    def __hash__(self) -> int:
        return hash((self.__class__, self.name))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.name}>"

    def __eq__(self, other: Self | Entity) -> bool:
        return isinstance(other, self.__class__) and self.name == other.name

    def __gt__(self, other: Self | Entity) -> bool:
        return isinstance(other, self.__class__) and self.name > other.name
    
    def __getitem__(self, key: str) -> Any:
        return self.params.kwargs.get(key, None)
    
    def __setitem__(self, key: str, value: Any) -> None:
        self.params.kwargs[key] = value

    def get(self, key: str, default: Any = None) -> Any:
        return self.params.kwargs.get(key, default)

class Style:
    def __init__(self, name: str, *args: Any, **kwargs: Any):
        self.name = name
        self.params = Parameters(*args, **kwargs)
        self.types: TypeBucket[Type] = TypeBucket(container_type=set)

    def __hash__(self) -> int:
        return hash((self.__class__, self.name))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.name}>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, self.__class__) and self.name == other.name

    def merge(self, other: "Style"):
        self.params.args.extend(other.params.args)
        self.params.kwargs.update(other.params.kwargs)
        for t in other.types.bucket(Type):
            self.types.add(t)
        return self

# ===================================================================
#                    ForceField 基类
# ===================================================================

class ForceField:
    # Kernel registry for potential functions
    _kernel_registry: dict[str, dict[str, type]] = {}
    
    def __init__(self, name: str = "", units: str = "real"):
        self.name = name
        self.units = units
        self.styles: TypeBucket[Style] = TypeBucket(container_type=set)

    def def_style(self, style_class: type[S], name: str, *args: Any, **kwargs: Any) -> S:
        existing_styles = self.styles.bucket(style_class)
        for style in existing_styles:
            if style.name == name:
                return style
        
        new_style = style_class(name, *args, **kwargs)
        self.styles.add(new_style)
        return new_style

    def get_styles(self, style_class: type[S]) -> list[S]:
        return self.styles.bucket(style_class)

    def get_types(self, type_class: type[Ty]) -> list[Ty]:
        all_types = set()
        for style in self.styles.bucket(Style):
            all_types.update(style.types.bucket(type_class))
        return list(all_types)

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.name}>"

    def merge(self, other: "ForceField"):
        for other_style in other.styles.bucket(Style):
            style_bucket = self.styles.bucket(type(other_style))
            found = False
            for style in style_bucket:
                if style.name == other_style.name:
                    style.merge(other_style)
                    found = True
                    break
            if not found:
                self.styles.add(other_style)
        return self

# ===================================================================
#               扩展的 AtomisticForcefield 类
# ===================================================================

class AtomType(Type):
    
    @property
    def is_wildcard(self) -> bool:
        """Check if this is a wildcard atom type (* or empty string)."""
        return self.name in ("*", "")
    
    @property
    def class_name(self) -> str:
        """Get the class name for this atom type."""
        # Check both 'class_name' and 'class_' for compatibility
        return self.params.kwargs.get("class_name", 
                                      self.params.kwargs.get("class_", self.name))
    
    def __eq__(self, other: object) -> bool:
        """Wildcard matching: wildcards match any atom type."""
        if not isinstance(other, AtomType):
            return False
        # If either is wildcard, they match
        if self.is_wildcard or other.is_wildcard:
            return True
        # Otherwise check name equality
        return self.name == other.name
    
    def __hash__(self) -> int:
        """Wildcards all hash to the same value."""
        if self.is_wildcard:
            return hash("__WILDCARD__")
        return hash(self.name)

    def __repr__(self) -> str:
        if self.is_wildcard:
            return f"<{self.__class__.__name__}: WILDCARD({self.name!r})>"
        return f"<{self.__class__.__name__}: {self.name}>"

class BondType(Type):
    """键类型，由两个原子类型定义"""
    
    def __init__(self, name: str, itom: "AtomType", jtom: "AtomType", **kwargs: Any):
        # Store atom types in args for backward compatibility
        super().__init__(name, itom, jtom, **kwargs)
        self.itom = itom
        self.jtom = jtom
    
    def matches(self, at1: "AtomType", at2: "AtomType") -> bool:
        """检查是否匹配给定的原子类型对（支持通配符和顺序无关）"""
        return (self.itom == at1 and self.jtom == at2) or \
               (self.itom == at2 and self.jtom == at1)
    
    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.itom.name}-{self.jtom.name}>"


class AngleType(Type):
    """角类型，由三个原子类型定义"""
    
    def __init__(self, name: str, itom: "AtomType", jtom: "AtomType", ktom: "AtomType", **kwargs: Any):
        # Store atom types in args for backward compatibility
        super().__init__(name, itom, jtom, ktom, **kwargs)
        self.itom = itom
        self.jtom = jtom
        self.ktom = ktom
    
    def matches(self, at1: "AtomType", at2: "AtomType", at3: "AtomType") -> bool:
        """检查是否匹配给定的原子类型三元组（支持通配符和顺序反转）"""
        # 正向匹配
        if self.itom == at1 and self.jtom == at2 and self.ktom == at3:
            return True
        # 反向匹配
        if self.itom == at3 and self.jtom == at2 and self.ktom == at1:
            return True
        return False
    
    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.itom.name}-{self.jtom.name}-{self.ktom.name}>"


class DihedralType(Type):
    """二面角类型，由四个原子类型定义"""
    
    def __init__(self, name: str, itom: "AtomType", jtom: "AtomType", 
                 ktom: "AtomType", ltom: "AtomType", **kwargs: Any):
        # Store atom types in args for backward compatibility
        super().__init__(name, itom, jtom, ktom, ltom, **kwargs)
        self.itom = itom
        self.jtom = jtom
        self.ktom = ktom
        self.ltom = ltom
    
    def matches(self, at1: "AtomType", at2: "AtomType", 
                at3: "AtomType", at4: "AtomType") -> bool:
        """检查是否匹配给定的原子类型四元组（支持通配符和顺序反转）"""
        # 正向匹配
        if (self.itom == at1 and self.jtom == at2 and 
            self.ktom == at3 and self.ltom == at4):
            return True
        # 反向匹配
        if (self.itom == at4 and self.jtom == at3 and 
            self.ktom == at2 and self.ltom == at1):
            return True
        return False
    
    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.itom.name}-{self.jtom.name}-{self.ktom.name}-{self.ltom.name}>"


class ImproperType(Type):
    """不正常二面角类型，由四个原子类型定义"""
    
    def __init__(self, name: str, itom: "AtomType", jtom: "AtomType", 
                 ktom: "AtomType", ltom: "AtomType", **kwargs: Any):
        super().__init__(name, **kwargs)
        self.itom = itom
        self.jtom = jtom
        self.ktom = ktom
        self.ltom = ltom
    
    def matches(self, at1: "AtomType", at2: "AtomType", 
                at3: "AtomType", at4: "AtomType") -> bool:
        """检查是否匹配给定的原子类型四元组（支持通配符）"""
        # Improper通常有特定的中心原子，所以匹配规则可能不同
        # 这里先实现简单的精确匹配
        return (self.itom == at1 and self.jtom == at2 and 
                self.ktom == at3 and self.ltom == at4)
    
    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.itom.name}-{self.jtom.name}-{self.ktom.name}-{self.ltom.name}>"


class PairType(Type):
    """非键相互作用类型，由一个或两个原子类型定义"""
    
    def __init__(self, name: str, *atom_types: "AtomType", **kwargs: Any):
        super().__init__(name, **kwargs)
        if len(atom_types) == 1:
            self.itom = atom_types[0]
            self.jtom = atom_types[0]  # 自相互作用
        elif len(atom_types) == 2:
            self.itom = atom_types[0]
            self.jtom = atom_types[1]
        else:
            raise ValueError("PairType requires 1 or 2 atom types")
    
    def matches(self, at1: "AtomType", at2: "AtomType | None" = None) -> bool:
        """检查是否匹配给定的原子类型对（支持通配符和顺序无关）"""
        if at2 is None:
            at2 = at1  # 自相互作用
        
        return (self.itom == at1 and self.jtom == at2) or \
               (self.itom == at2 and self.jtom == at1)
    
    def __repr__(self) -> str:
        if self.itom == self.jtom:
            return f"<{self.__class__.__name__}: {self.itom.name}>"
        return f"<{self.__class__.__name__}: {self.itom.name}-{self.jtom.name}>"

class AtomStyle(Style):
    def def_type(self, name: str, **kwargs: Any) -> AtomType:
        """定义原子类型
        
        Args:
            type_: 具体类型标识符（如 opls_135）
            class_: 类别标识符（如 CT）
            **kwargs: 其他参数（element, mass 等）
            
        Returns:
            创建的 AtomType 实例
        """
        at = AtomType(name=name, **kwargs)
        self.types.add(at)
        return at

class BondStyle(Style):
    def def_type(self, itom: AtomType, jtom: AtomType, name: str = "", **kwargs: Any) -> BondType:
        """定义键类型
        
        Args:
            itom: 第一个原子类型
            jtom: 第二个原子类型
            name: 可选的名称（默认为 itom-jtom）
            **kwargs: 键参数（如 k, r0 等）
        """
        if not name:
            name = f"{itom.name}-{jtom.name}"
        bt = BondType(name, itom, jtom, **kwargs)
        self.types.add(bt)
        return bt


class AngleStyle(Style):
    def def_type(self, itom: AtomType, jtom: AtomType, ktom: AtomType, 
                 name: str = "", **kwargs: Any) -> AngleType:
        """定义角类型
        
        Args:
            itom: 第一个原子类型
            jtom: 中心原子类型
            ktom: 第三个原子类型
            name: 可选的名称（默认为 itom-jtom-ktom）
            **kwargs: 角参数（如 k, theta0 等）
        """
        if not name:
            name = f"{itom.name}-{jtom.name}-{ktom.name}"
        at = AngleType(name, itom, jtom, ktom, **kwargs)
        self.types.add(at)
        return at


class DihedralStyle(Style):
    def def_type(self, itom: AtomType, jtom: AtomType, ktom: AtomType, 
                 ltom: AtomType, name: str = "", **kwargs: Any) -> DihedralType:
        """定义二面角类型
        
        Args:
            itom: 第一个原子类型
            jtom: 第二个原子类型
            ktom: 第三个原子类型
            ltom: 第四个原子类型
            name: 可选的名称（默认为 itom-jtom-ktom-ltom）
            **kwargs: 二面角参数
        """
        if not name:
            name = f"{itom.name}-{jtom.name}-{ktom.name}-{ltom.name}"
        dt = DihedralType(name, itom, jtom, ktom, ltom, **kwargs)
        self.types.add(dt)
        return dt


class ImproperStyle(Style):
    def def_type(self, itom: AtomType, jtom: AtomType, ktom: AtomType, 
                 ltom: AtomType, name: str = "", **kwargs: Any) -> ImproperType:
        """定义不正常二面角类型
        
        Args:
            itom: 第一个原子类型
            jtom: 第二个原子类型（通常是中心原子）
            ktom: 第三个原子类型
            ltom: 第四个原子类型
            name: 可选的名称（默认为 itom-jtom-ktom-ltom）
            **kwargs: 不正常二面角参数
        """
        if not name:
            name = f"{itom.name}-{jtom.name}-{ktom.name}-{ltom.name}"
        it = ImproperType(name, itom, jtom, ktom, ltom, **kwargs)
        self.types.add(it)
        return it


class PairStyle(Style):
    def def_type(self, itom: AtomType, jtom: AtomType | None = None, 
                 name: str = "", **kwargs: Any) -> PairType:
        """定义非键相互作用类型
        
        Args:
            itom: 第一个原子类型
            jtom: 第二个原子类型（可选，默认为与 itom 相同，即自相互作用）
            name: 可选的名称
            **kwargs: 非键参数（如 sigma, epsilon, charge 等）
        """
        if jtom is None:
            jtom = itom
        
        if not name:
            if itom == jtom:
                name = itom.name
            else:
                name = f"{itom.name}-{jtom.name}"
        
        pt = PairType(name, itom, jtom, **kwargs)
        self.types.add(pt)
        return pt

class AtomisticForcefield(ForceField):
    def def_atomstyle(self, name: str, *args: Any, **kwargs: Any) -> AtomStyle:
        return self.def_style(AtomStyle, name, *args, **kwargs)

    def def_bondstyle(self, name: str, *args: Any, **kwargs: Any) -> BondStyle:
        return self.def_style(BondStyle, name, *args, **kwargs)

    def def_anglestyle(self, name: str, *args: Any, **kwargs: Any) -> AngleStyle:
        return self.def_style(AngleStyle, name, *args, **kwargs)

    def def_dihedralstyle(self, name: str, *args: Any, **kwargs: Any) -> DihedralStyle:
        return self.def_style(DihedralStyle, name, *args, **kwargs)

    def def_improperstyle(self, name: str, *args: Any, **kwargs: Any) -> ImproperStyle:
        return self.def_style(ImproperStyle, name, *args, **kwargs)

    def def_pairstyle(self, name: str, *args: Any, **kwargs: Any) -> PairStyle:
        return self.def_style(PairStyle, name, *args, **kwargs)

    def get_atomtypes(self) -> list[AtomType]:
        return self.get_types(AtomType)

    def get_bondtypes(self) -> list[BondType]:
        return self.get_types(BondType)

    def get_angletypes(self) -> list[AngleType]:
        return self.get_types(AngleType)
    
    def get_dihedraltypes(self) -> list[DihedralType]:
        return self.get_types(DihedralType)