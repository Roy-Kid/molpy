import csv
from collections import defaultdict
from io import StringIO
from pathlib import Path
from typing import Any, TypeVar, Generic, Iterator, cast

from .frame import Block

# --- TypeBucket 通用实现 ---

T = TypeVar("T")

def get_nearest_type(item: T) -> type[T]:
    """获取对象的具体类型"""
    return type(item)


class TypeBucket(Generic[T]):
    """
    通用的类型桶（TypeBucket）实现，按对象的具体类型分组存储。
    
    支持两种存储模式：
    - 集合模式（set）：用于存储唯一的类型定义，自动去重
    - 列表模式（list）：用于存储实体对象，保持顺序
    
    使用示例：
        # 集合模式（forcefield 使用）
        bucket = TypeBucket(container_type=set)
        
        # 列表模式（entity 使用）
        bucket = TypeBucket(container_type=list)
    """
    
    def __init__(self, container_type: type = set) -> None:
        """
        初始化 TypeBucket
        
        Args:
            container_type: 容器类型，set（去重）或 list（保持顺序）
        """
        self.container_type = container_type
        if container_type == set:
            self._items: dict[type[T], Any] = defaultdict(set)
        else:
            self._items: dict[type[T], Any] = defaultdict(list)
    
    def add(self, item: T) -> None:
        """添加一个对象到其对应类型的桶中"""
        cls = get_nearest_type(item)
        if self.container_type == set:
            self._items[cls].add(item)
        else:
            self._items[cls].append(item)
    
    def add_many(self, items: Any) -> None:
        """添加多个对象"""
        for item in items:
            self.add(item)
    
    def remove(self, item: T) -> bool:
        """
        从桶中移除对象
        
        Returns:
            如果成功移除返回 True，否则返回 False
        """
        cls = get_nearest_type(item)
        bucket = self._items.get(cls)
        if not bucket:
            return False
        
        if self.container_type == set:
            if item in bucket:
                bucket.discard(item)
                if not bucket:
                    self._items.pop(cls, None)
                return True
            return False
        else:
            # 列表模式：使用身份比较（is）而非相等比较（==）
            for i, obj in enumerate(bucket):
                if obj is item:
                    bucket.pop(i)
                    if not bucket:
                        self._items.pop(cls, None)
                    return True
            return False
    
    def bucket(self, cls: type[T]) -> list[T]:
        """
        获取指定类型及其所有子类的对象
        
        Args:
            cls: 目标类型
            
        Returns:
            包含所有匹配对象的列表
        """
        result: list[T] = []
        for k, items in self._items.items():
            if isinstance(k, type) and issubclass(k, cls):
                result.extend(items)
        return result
    
    def exact_bucket(self, cls: type[T]) -> list[T]:
        """
        获取精确类型的对象（不包括子类）
        
        Args:
            cls: 目标类型
            
        Returns:
            包含精确匹配对象的列表
        """
        bucket = self._items.get(cls)
        return list(bucket) if bucket else []
    
    def classes(self) -> Iterator[type[T]]:
        """返回当前存储的所有具体类型"""
        return iter(self._items.keys())
    
    def all(self) -> list[T]:
        """返回所有桶中的所有对象"""
        result: list[T] = []
        for bucket in self._items.values():
            result.extend(bucket)
        return result
    
    def __len__(self) -> int:
        """返回所有桶中对象的总数"""
        return sum(len(b) for b in self._items.values())
    
    def __getitem__(self, cls: type[T]) -> list[T]:
        """获取精确类型的桶（不包括子类）"""
        return self.exact_bucket(cls)
    
    def __setitem__(self, cls: type[T], items: Any) -> None:
        """设置指定类型的桶"""
        if self.container_type == set:
            self._items[cls] = set(items)
        else:
            self._items[cls] = list(items)


def to_dict_of_list(list_of_dict: list[dict[str, Any]]) -> dict[str, list[Any]]:
    result = defaultdict(list)
    for item in list_of_dict:
        for key, value in item.items():
            result[key].append(value)
    return dict(result)


def to_list_of_dict(list_of_dict: dict[str, list[Any]]) -> list[dict[str, Any]]:
    keys = list(list_of_dict.keys())
    length = len(next(iter(list_of_dict.values())))
    return [{key: list_of_dict[key][i] for key in keys} for i in range(length)]


def read_csv(file: Path | StringIO, delimiter: str = ",") -> Block:
    """
    Read a CSV file or StringIO object and return a Block object.

    Args:
        file: Path to the CSV file or a StringIO object containing CSV data.
        delimiter: Delimiter used in the CSV file (default is comma).

    Returns:
        Block: A Block object containing the data from the CSV.
    """
    if isinstance(file, StringIO):
        file.seek(0)  # Ensure we read from the start
        reader = csv.DictReader(file, delimiter=delimiter)
        data = [row for row in reader]
    else:
        file = Path(file)
        if not file.exists():
            raise FileNotFoundError(f"File {file} does not exist.")

        with open(file, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            data = [row for row in reader]

    # Convert list of dicts to dict of lists for Block
    dict_of_lists = to_dict_of_list(data)

    # Handle empty data case - create empty arrays for each column
    if not dict_of_lists and isinstance(file, StringIO):
        file.seek(0)
        reader = csv.DictReader(file, delimiter=delimiter)
        fieldnames = reader.fieldnames or []
        dict_of_lists = {field: [] for field in fieldnames}
    elif not dict_of_lists and not isinstance(file, StringIO):
        with open(file, "r", newline="") as csvfile:
            reader = csv.DictReader(csvfile, delimiter=delimiter)
            fieldnames = reader.fieldnames or []
            dict_of_lists = {field: [] for field in fieldnames}

    return Block(dict_of_lists)
