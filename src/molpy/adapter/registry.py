from __future__ import annotations
from typing import Any, Callable, TypeVar, overload

# ---- 基础类型 ----
_To = TypeVar("_To")
_ConverterFn = Callable[[Any], Any]
_Pred = Callable[[Any], bool]

# ---- 注册表 ----
class ConverterRegistry:
    def __init__(self) -> None:
        self._by_pair: dict[tuple[type, type], _ConverterFn] = {}
        self._preds: list[tuple[_Pred, type, _ConverterFn]] = []  # (pred, dst_type, fn)

    def register(self, src: type, dst: type, fn: _ConverterFn) -> None:
        self._by_pair[(src, dst)] = fn

    def register_predicate(self, pred: _Pred, dst: type, fn: _ConverterFn) -> None:
        self._preds.append((pred, dst, fn))

    def resolve(self, obj_or_src: Any | type, dst: type) -> _ConverterFn | None:
        # 允许传实例或源类型
        if isinstance(obj_or_src, type):
            return self._by_pair.get((obj_or_src, dst))
        obj = obj_or_src
        # 1) 按 MRO 精确/上溯匹配
        for src in type(obj).mro():
            fn = self._by_pair.get((src, dst))
            if fn:
                return fn
        # 2) 谓词匹配（按注册顺序）
        for pred, d, fn in self._preds:
            if d is dst:
                try:
                    if pred(obj):
                        return fn
                except Exception:
                    pass
        return None

    def get_converter(self, src: type, dst: type) -> _ConverterFn | None:
        return self._by_pair.get((src, dst))


REG = ConverterRegistry()

# ---- 统一入口：既可调用也可下标 ----
class _Convert:
    def __call__(self, obj: Any, to: type[_To]) -> _To:
        fn = REG.resolve(obj, to)
        if fn is None:
            raise TypeError(f"No converter registered for {type(obj)} -> {to}, alternatives: {list(REG._by_pair.keys())}")
        return fn(obj)  # type: ignore[return-value]

    def __getitem__(self, dst: type[_To]) -> Callable[[Any], _To]:
        def _f(obj: Any) -> _To:
            return self(obj, dst)
        return _f

convert = _Convert()  # ★ 统一的单例

# ---- 装饰器糖 ----
def register(src: type, dst: type):
    def _wrap(fn: _ConverterFn) -> _ConverterFn:
        REG.register(src, dst, fn)
        return fn
    return _wrap

def register_if(pred: _Pred, dst: type):
    def _wrap(fn: _ConverterFn) -> _ConverterFn:
        REG.register_predicate(pred, dst, fn)
        return fn
    return _wrap
