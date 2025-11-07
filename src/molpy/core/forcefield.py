from collections import defaultdict
from typing import Any, TypeVar, Generic, Iterator

# --- 泛型变量和 TypeBucket 实现 ---

T = TypeVar("T")
S = TypeVar("S", bound="Style")
Ty = TypeVar("Ty", bound="Type")

def get_nearest_type(item: T) -> type[T]:
    return type(item)

class TypeBucket(Generic[T]):
    def __init__(self) -> None:
        self._items: dict[type[T], set[T]] = defaultdict(set)

    def add(self, item: T) -> None:
        cls = get_nearest_type(item)
        self._items[cls].add(item)

    def remove(self, item: T) -> None:
        cls = get_nearest_type(item)
        self._items[cls].discard(item)

    def bucket(self, cls: type[T]) -> list[T]:
        result: list[T] = []
        for k, items in self._items.items():
            if issubclass(k, cls):
                result.extend(items)
        return result

    def classes(self) -> Iterator[type[T]]:
        return iter(self._items.keys())

# --- 基础组件 ---

class Parameters:
    def __init__(self, *args: Any, **kwargs: Any):
        self.args = list(args)
        self.kwargs = kwargs

    def __repr__(self) -> str:
        return f"Parameters(args={self.args}, kwargs={self.kwargs})"

class Type:
    def __init__(self, name: str, *args: Any, **kwargs: Any):
        self._name = name
        self.params = Parameters(*args, **kwargs)

    def __hash__(self) -> int:
        return hash((self.__class__, self.name))

    def __repr__(self) -> str:
        return f"<{self.__class__.__name__}: {self.name}>"

    def __eq__(self, other: object) -> bool:
        return isinstance(other, self.__class__) and self.name == other.name

class Style:
    def __init__(self, name: str, *args: Any, **kwargs: Any):
        self.name = name
        self.params = Parameters(*args, **kwargs)
        self.types = TypeBucket[Type]()

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
    def __init__(self, name: str = "", units: str = "real"):
        self.name = name
        self.units = units
        self.styles = TypeBucket[Style]()

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

class AtomType(Type): pass
class BondType(Type): pass
class AngleType(Type): pass
class DihedralType(Type): pass
class ImproperType(Type): pass
class PairType(Type): pass

class AtomStyle(Style):
    def def_type(self, name: str, *args: Any, **kwargs: Any) -> AtomType:
        at = AtomType(name, *args, **kwargs)
        self.types.add(at)
        return at

class BondStyle(Style):
    def def_type(self, *atom_types: AtomType, name: str = "", **kwargs: Any) -> BondType:
        type_names = '-'.join(at.name for at in atom_types)
        bt = BondType(name or type_names, *atom_types, **kwargs)
        self.types.add(bt)
        return bt

class AngleStyle(Style):
    def def_type(self, *atom_types: AtomType, name: str = "", **kwargs: Any) -> AngleType:
        type_names = '-'.join(at.name for at in atom_types)
        at = AngleType(name or type_names, *atom_types, **kwargs)
        self.types.add(at)
        return at

class DihedralStyle(Style): pass
class ImproperStyle(Style): pass
class PairStyle(Style): pass

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