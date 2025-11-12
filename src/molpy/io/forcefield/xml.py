"""XML force field parser for atomistic force fields."""

import logging
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Any

from molpy.core.forcefield import (
    AtomisticForcefield,
    AtomType,
    DihedralType,
)

logger = logging.getLogger(__name__)


def _normalize_to_wildcard(value: str | None) -> str:
    """
    Normalize None or empty string to wildcard "*".
    
    Args:
        value: Input value (can be None, "", or any string)
        
    Returns:
        "*" if value is None or "", otherwise the original value
    """
    if value is None or value == "":
        return "*"
    return value


def _get_canonical_class(class_: str, aliases: dict[str, str]) -> str:
    """
    获取类别的规范形式（解析别名）
    
    例如：CT_2 -> CT, CT_3 -> CT
    
    Args:
        class_: 原始类别名
        aliases: 别名映射字典
        
    Returns:
        规范类别名
    """
    # 如果是通配符或者不在别名表中，直接返回
    if class_ == "*" or class_ not in aliases:
        return class_
    # 递归解析（防止多层别名）
    canonical = aliases[class_]
    if canonical in aliases:
        return _get_canonical_class(canonical, aliases)
    return canonical


def _resolve_forcefield_path(filepath: str | Path) -> Path:
    """
    Resolve force field file path, checking built-in data directory first.
    
    Args:
        filepath: Path to force field file, or just filename for built-in files
        
    Returns:
        Resolved Path object
        
    Raises:
        FileNotFoundError: If file not found in any location
    """
    filepath = Path(filepath)
    
    # If it's just a filename (e.g., "oplsaa.xml"), check built-in data
    if filepath.name == str(filepath):
        # Try built-in data directory
        import molpy
        molpy_root = Path(molpy.__file__).parent
        builtin_path = molpy_root / "data" / "forcefield" / filepath.name
        
        if builtin_path.exists():
            logger.info(f"Using built-in force field: {builtin_path}")
            return builtin_path
    
    # Otherwise use the provided path
    if filepath.exists():
        return filepath
    
    raise FileNotFoundError(f"Force field file not found: {filepath}")


class XMLForceFieldReader:
    """
    XML force field parser for atomistic force fields.

    Parses XML-formatted force field files (e.g., OPLS-AA) and populates
    an AtomisticForcefield object with atom types, bond parameters, angle
    parameters, dihedral parameters, and nonbonded interactions.

    The parser handles:
    - AtomTypes section: atom type definitions
    - HarmonicBondForce: harmonic bond parameters
    - HarmonicAngleForce: harmonic angle parameters
    - RBTorsionForce: Ryckaert-Bellemans dihedral parameters
    - NonbondedForce: LJ and Coulomb parameters
    """

    def __init__(self, filepath: str | Path):
        """
        Initialize the XML force field reader.

        Args:
            filepath: Path to the XML force field file, or filename for built-in files
                     (e.g., "oplsaa.xml" will load from molpy/data/forcefield/)
        """
        self._file = _resolve_forcefield_path(filepath)
        self._type_to_atomtype: dict[str, AtomType] = {}  # type -> AtomType mapping
        self._class_to_atomtype: dict[str, AtomType] = {}  # class -> AtomType mapping
        self._any_atomtype: AtomType | None = None  # 全通配符 ("*", "*")
        # Alias for backward compatibility
        self._atomtype_cache = self._type_to_atomtype
        # 类别别名映射：派生类 -> 规范基类
        # 例如 CT_2 -> CT, CT_3 -> CT
        self._class_aliases: dict[str, str] = {}
        # overrides 映射：type -> overridden_type
        # 例如 opls_961 overrides opls_962
        self._overrides: dict[str, str] = {}
        
        self._ff: AtomisticForcefield | None = None
    
    @property
    def _wildcard_atomtype(self) -> AtomType | None:
        """Alias for _any_atomtype for backward compatibility."""
        return self._any_atomtype

    def read(self, forcefield: AtomisticForcefield | None = None) -> AtomisticForcefield:
        """
        Read and parse the XML force field file.

        Args:
            forcefield: Optional existing force field to populate. If None, creates new one.

        Returns:
            Populated AtomisticForcefield object
        """
        if not self._file.exists():
            raise FileNotFoundError(f"Force field file not found: {self._file}")

        tree = ET.parse(self._file)
        root = tree.getroot()

        # Get force field metadata from root element
        ff_name = root.get("name", "Unknown")
        ff_version = root.get("version", "0.0.0")
        combining_rule = root.get("combining_rule", "geometric")

        # Create or use provided force field
        if forcefield is None:
            self._ff = AtomisticForcefield(name=ff_name, units="real")
        else:
            self._ff = forcefield

        logger.info(f"Parsing force field: {ff_name} v{ff_version}")
        logger.info(f"Combining rule: {combining_rule}")

        # Parse all sections
        for child in root:
            tag = child.tag
            if tag == "AtomTypes":
                self._parse_atomtypes(child)
            elif tag == "HarmonicBondForce":
                self._parse_bonds(child)
            elif tag == "HarmonicAngleForce":
                self._parse_angles(child)
            elif tag == "RBTorsionForce":
                self._parse_dihedrals(child)
            elif tag == "NonbondedForce":
                self._parse_nonbonded(child)
            else:
                logger.debug(f"Skipping unknown section: {tag}")

        logger.info(f"Parsed {len(self._type_to_atomtype)} atom types (by type)")
        return self._ff
    
    def _resolve_atomtype_with_alias(self, class_str: str) -> AtomType:
        """
        解析类别字符串，如果不存在则尝试建立别名映射。
        
        策略：如果 class_str 形如 "XX_N"（N为数字），且 XX 存在于 class_to_atomtype，
        则建立 XX_N -> XX 的别名映射，并返回 XX 对应的 AtomType。
        
        Args:
            class_str: 类别字符串（来自 BondType/AngleType/DihedralType 的 class 属性）
            
        Returns:
            对应的 AtomType
            
        Raises:
            KeyError: 如果找不到对应的 AtomType
        """
        # 先尝试直接查找
        if class_str in self._class_to_atomtype:
            return self._class_to_atomtype[class_str]
        
        # 如果找不到，检查是否是派生类（XX_N 形式）
        if "_" in class_str:
            parts = class_str.rsplit("_", 1)
            if len(parts) == 2 and parts[1].isdigit():
                base_class = parts[0]
                # 检查基类是否存在
                if base_class in self._class_to_atomtype:
                    # 建立别名映射
                    self._class_aliases[class_str] = base_class
                    logger.debug(f"Auto-registered class alias: {class_str} -> {base_class}")
                    return self._class_to_atomtype[base_class]
        
        # 都找不到，抛出异常
        raise KeyError(f"AtomType with class '{class_str}' not found")

    def _get_or_create_atomtype(self, type_: str, class_: str, **kwargs: Any) -> AtomType:
        """
        Get or create an AtomType with exact type_ and class_.
        
        使用别名映射：派生类（如 CT_2, CT_3）会被规范化为基类（CT）
        
        维护三个映射表：
        1. type_to_atomtype: type -> AtomType (仅当 type != "*" 时)
        2. class_to_atomtype: class -> AtomType (仅当 class != "*" 时)
        3. any_atomtype: 全通配符 ("*", "*")
        
        Args:
            type_: 类型标识符（已规范化，None或""已转为"*"）
            class_: 类别标识符（已规范化，None或""已转为"*"）
            **kwargs: 其他参数（element, mass 等）
            
        Returns:
            AtomType 实例
        """
        # 规范化输入
        type_ = _normalize_to_wildcard(type_)
        class_ = _normalize_to_wildcard(class_)
        
        # 规范化 class（解析别名）
        # 如果 class 不在 _class_to_atomtype 中，尝试别名解析
        canonical_class = class_
        if class_ != "*" and class_ not in self._class_to_atomtype:
            # 尝试别名解析
            if class_ in self._class_aliases:
                canonical_class = _get_canonical_class(class_, self._class_aliases)
            elif "_" in class_:
                # 自动检测并注册别名（XX_N -> XX）
                parts = class_.rsplit("_", 1)
                if len(parts) == 2 and parts[1].isdigit():
                    base_class = parts[0]
                    if base_class in self._class_to_atomtype:
                        self._class_aliases[class_] = base_class
                        canonical_class = base_class
                        logger.debug(f"Auto-registered class alias: {class_} -> {base_class}")
        elif class_ in self._class_aliases:
            canonical_class = _get_canonical_class(class_, self._class_aliases)
        
        # 确定 name：type 优先，否则 canonical_class，否则 "*"
        if type_ != "*":
            name = type_
        elif canonical_class != "*":
            name = canonical_class
        else:
            name = "*"
        
        # 将原始 type_ 和 canonical class_ 添加到 kwargs
        kwargs['type_'] = type_
        kwargs['class_'] = canonical_class  # 使用规范化的 class
        # also expose a full_name (type-class) to preserve both identifiers
        if type_ != "*":
            kwargs.setdefault('full_name', f"{type_}-{canonical_class}")
        else:
            kwargs.setdefault('full_name', canonical_class)
        
        # 情况1: 全通配符 ("*", "*")
        if type_ == "*" and canonical_class == "*":
            if self._any_atomtype is None:
                atomstyle = self._ff.def_atomstyle("full")
                self._any_atomtype = atomstyle.def_type(name=name, **kwargs)
                # record overrides if present for a specific type (not wildcard)
                ov = kwargs.get("overrides")
                if ov and type_ != "*":
                    self._overrides[type_] = ov
                logger.debug("Created global wildcard atom type (*,*)")
            return self._any_atomtype
        
        # 情况2: 具体 type，任意 class (type, "*")
        if type_ != "*" and canonical_class == "*":
            if type_ in self._type_to_atomtype:
                return self._type_to_atomtype[type_]
            atomstyle = self._ff.def_atomstyle("full")
            atomtype = atomstyle.def_type(name=name, **kwargs)
            ov = kwargs.get("overrides")
            if ov and type_ != "*":
                self._overrides[type_] = ov
            self._type_to_atomtype[type_] = atomtype
            logger.debug(f"Created atom type ({type_}, *)")
            return atomtype
        
        # 情况3: 任意 type，具体 class ("*", class)
        if type_ == "*" and canonical_class != "*":
            if canonical_class in self._class_to_atomtype:
                return self._class_to_atomtype[canonical_class]
            atomstyle = self._ff.def_atomstyle("full")
            atomtype = atomstyle.def_type(name=name, **kwargs)
            ov = kwargs.get("overrides")
            if ov and type_ != "*":
                self._overrides[type_] = ov
            self._class_to_atomtype[canonical_class] = atomtype
            logger.debug(f"Created atom type (*, {class_})")
            return atomtype
        
        # 情况4: 具体 type 和 class (type, class)
        # 优先查找 type mapping
        if type_ in self._type_to_atomtype:
            return self._type_to_atomtype[type_]
        # 创建新的 AtomType
        atomstyle = self._ff.def_atomstyle("full")
        atomtype = atomstyle.def_type(name=name, **kwargs)
        ov = kwargs.get("overrides")
        if ov and type_ != "*":
            self._overrides[type_] = ov
        # 只存入 type_to_atomtype，不存入 class_to_atomtype
        # 因为 class_to_atomtype 只用于存储 (*, class) 形式的通用 AtomType
        self._type_to_atomtype[type_] = atomtype
        logger.debug(f"Created atom type ({type_}, {canonical_class})")
        return atomtype
        self._type_to_atomtype[type_] = atomtype
        self._class_to_atomtype[class_] = atomtype
        logger.debug(f"Created atom type ({type_}, {class_})")
        return atomtype


    def _parse_atomtypes(self, element: ET.Element) -> None:
        """
        Parse AtomTypes section.

        Args:
            element: AtomTypes XML element
        """
        atomstyle = self._ff.def_atomstyle("full")
        count = 0

        for type_elem in element:
            if type_elem.tag != "Type":
                continue

            # Extract attributes
            # In XML: "name" is the type, "class" is the class
            type_name = type_elem.get("name")  # 可能是 None
            class_name = type_elem.get("class")  # 可能是 None
            element_sym = type_elem.get("element")
            mass_str = type_elem.get("mass")
            def_str = type_elem.get("def")
            desc_str = type_elem.get("desc")
            doi_str = type_elem.get("doi")
            overrides_str = type_elem.get("overrides")

            # Parse mass
            mass = float(mass_str) if mass_str else 0.0

            # Create atom type using _get_or_create_atomtype
            # 使用 _normalize_to_wildcard 将 None 转换为 "*"
            atomtype = self._get_or_create_atomtype(
                type_=type_name or "",  # 空字符串会被规范化为 "*"
                class_=class_name or "",
                element=element_sym,
                mass=mass,
                def_=def_str,
                desc=desc_str,
                doi=doi_str,
                overrides=overrides_str,
            )

            count += 1

        logger.info(f"Parsed {count} atom types")

    def _parse_bonds(self, element: ET.Element) -> None:
        """
        Parse HarmonicBondForce section.

        Args:
            element: HarmonicBondForce XML element
        """
        bondstyle = self._ff.def_bondstyle("harmonic")
        count = 0

        for bond_elem in element:
            if bond_elem.tag != "Bond":
                continue

            # Get atom type identifiers
            class1 = bond_elem.get("class1", "*")
            class2 = bond_elem.get("class2", "*")
            type1 = bond_elem.get("type1", "*")
            type2 = bond_elem.get("type2", "*")

            # Get bond parameters
            length_str = bond_elem.get("length")
            k_str = bond_elem.get("k")

            # Check if types are present (None means missing, "" or "*" are valid wildcards)
            if type1 is None or type2 is None:
                logger.warning("Skipping bond without type information")
                continue

            # Get or create atom types (handles wildcards)
            at1 = self._get_or_create_atomtype(type1, class1)
            at2 = self._get_or_create_atomtype(type2, class2)

            # Parse parameters
            r0 = float(length_str) if length_str else 0.0
            k = float(k_str) if k_str else 0.0

            # Define bond type
            bondstyle.def_type(at1, at2, r0=r0, k=k)
            count += 1

        logger.info(f"Parsed {count} bond types")

    def _parse_angles(self, element: ET.Element) -> None:
        """
        Parse HarmonicAngleForce section.

        Args:
            element: HarmonicAngleForce XML element
        """
        anglestyle = self._ff.def_anglestyle("harmonic")
        count = 0

        for angle_elem in element:
            if angle_elem.tag != "Angle":
                continue

            # Get atom type identifiers
            class1 = angle_elem.get("class1", "*")
            class2 = angle_elem.get("class2", "*")
            class3 = angle_elem.get("class3", "*")
            type1 = angle_elem.get("type1", "*")
            type2 = angle_elem.get("type2", "*")
            type3 = angle_elem.get("type3", "*")

            # Get angle parameters
            angle_str = angle_elem.get("angle")
            k_str = angle_elem.get("k")

            # Check if types are present (None means missing, "" or "*" are valid wildcards)
            if type1 is None or type2 is None or type3 is None:
                logger.warning("Skipping angle without type information")
                continue

            # Get or create atom types (handles wildcards)
            at1 = self._get_or_create_atomtype(type1, class1)
            at2 = self._get_or_create_atomtype(type2, class2)
            at3 = self._get_or_create_atomtype(type3, class3)

            # Parse parameters
            theta0 = float(angle_str) if angle_str else 0.0
            k = float(k_str) if k_str else 0.0

            # Define angle type
            anglestyle.def_type(at1, at2, at3, theta0=theta0, k=k)
            count += 1

        logger.info(f"Parsed {count} angle types")

    def _parse_dihedrals(self, element: ET.Element) -> None:
        """
        Parse RBTorsionForce section (Ryckaert-Bellemans dihedrals).

        Args:
            element: RBTorsionForce XML element
        """
        dihedralstyle = self._ff.def_dihedralstyle("opls")
        count = 0

        for dihedral_elem in element:
            if dihedral_elem.tag != "Proper":
                continue

            # Get atom type identifiers
            class1 = dihedral_elem.get("class1", "*")
            class2 = dihedral_elem.get("class2", "*")
            class3 = dihedral_elem.get("class3", "*")
            class4 = dihedral_elem.get("class4", "*")
            type1 = dihedral_elem.get("type1", "*")
            type2 = dihedral_elem.get("type2", "*")
            type3 = dihedral_elem.get("type3", "*")
            type4 = dihedral_elem.get("type4", "*")

            # Check if types are present (None means missing, "" or "*" are valid wildcards)
            if type1 is None or type2 is None or type3 is None or type4 is None:
                logger.warning("Skipping dihedral without type information")
                continue

            # Get or create atom types (handles wildcards)
            at1 = self._get_or_create_atomtype(type1, class1)
            at2 = self._get_or_create_atomtype(type2, class2)
            at3 = self._get_or_create_atomtype(type3, class3)
            at4 = self._get_or_create_atomtype(type4, class4)

            # Parse RB coefficients (c0-c5)
            params = {}
            for i in range(6):
                c_str = dihedral_elem.get(f"c{i}")
                if c_str:
                    params[f"c{i}"] = float(c_str)

            # Define dihedral type - create and add manually
            # Use more specific name: if type is "*", use class name instead
            name1 = type1 if type1 != "*" else class1
            name2 = type2 if type2 != "*" else class2
            name3 = type3 if type3 != "*" else class3
            name4 = type4 if type4 != "*" else class4
            dihedral_name = f"{name1}-{name2}-{name3}-{name4}"
            
            # If still not unique (same name already exists), append counter
            # Note: we check by trying to find existing type with same name
            base_name = dihedral_name
            counter = 1
            existing_names = {dt.name for dt in dihedralstyle.types.bucket(DihedralType)}
            while dihedral_name in existing_names:
                dihedral_name = f"{base_name}#{counter}"
                counter += 1
            
            dihedral_type = DihedralType(dihedral_name, at1, at2, at3, at4, **params)
            dihedralstyle.types.add(dihedral_type)
            count += 1

        logger.info(f"Parsed {count} dihedral types")

    def _parse_nonbonded(self, element: ET.Element) -> None:
        """
        Parse NonbondedForce section (LJ and Coulomb parameters).

        Args:
            element: NonbondedForce XML element
        """
        # Get scaling factors
        coulomb14scale = element.get("coulomb14scale", "0.5")
        lj14scale = element.get("lj14scale", "0.5")

        pairstyle = self._ff.def_pairstyle(
            "lj/cut/coul/cut",
            coulomb14scale=float(coulomb14scale),
            lj14scale=float(lj14scale)
        )
        count = 0

        for atom_elem in element:
            if atom_elem.tag != "Atom":
                continue

            # Get atom type identifier
            type_name = atom_elem.get("type")
            if not type_name:
                logger.warning("Skipping nonbonded atom without type")
                continue

            # Get or create atom type
            # 在 NonbondedForce 中，type 通常对应实际的 type，没有 class
            atomtype = self._get_or_create_atomtype(type_=type_name, class_="*")

            # Parse parameters
            charge_str = atom_elem.get("charge")
            sigma_str = atom_elem.get("sigma")
            epsilon_str = atom_elem.get("epsilon")

            charge = float(charge_str) if charge_str else 0.0
            sigma = float(sigma_str) if sigma_str else 0.0
            epsilon = float(epsilon_str) if epsilon_str else 0.0

            # Define pair parameters (self-interaction)
            pairstyle.def_type(
                atomtype,
                atomtype,
                charge=charge,
                sigma=sigma,
                epsilon=epsilon
            )
            count += 1

        logger.info(f"Parsed {count} nonbonded parameters")

    def _is_wildcard(self, type_name: str) -> bool:
        """检查类型名是否为通配符 ("*")"""
        return type_name == "*" or type_name == ""


def read_xml_forcefield(
    filepath: str | Path, forcefield: AtomisticForcefield | None = None
) -> AtomisticForcefield:
    """
    Convenience function to read an XML force field file.

    Args:
        filepath: Path to the XML force field file, or filename for built-in files
                 (e.g., "oplsaa.xml" will load from molpy/data/forcefield/)
        forcefield: Optional existing force field to populate

    Returns:
        Populated AtomisticForcefield object

    Example:
        >>> # Load built-in OPLS-AA force field
        >>> ff = read_xml_forcefield("oplsaa.xml")
        >>> 
        >>> # Load custom force field from path
        >>> from pathlib import Path
        >>> ff = read_xml_forcefield(Path("/path/to/custom.xml"))
    """
    reader = XMLForceFieldReader(filepath)
    return reader.read(forcefield)
