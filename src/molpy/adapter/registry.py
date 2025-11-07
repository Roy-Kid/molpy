from typing import Any, Callable, TypeVar

T = TypeVar("T")
U = TypeVar("U")

_Pred = Callable[[T], bool]
_ConverterFn = Callable[[T], U]

class ConverterRegistry:
    """Registry for converter functions between different molecular representations."""

    def __init__(self) -> None:
        self._by_pair: dict[tuple[type, type], _ConverterFn] = {}
        self._preds: list[tuple[_Pred, type, _ConverterFn]] = []  # (pred, dst_type, fn)

    def register(self, src: type, dst: type, fn: _ConverterFn) -> None:
        self._by_pair[(src, dst)] = fn

    def register_predicate(self, pred: _Pred, dst: type, fn: _ConverterFn) -> None:
        self._preds.append((pred, dst, fn))

    def resolve(self, obj: Any, dst: type) -> _ConverterFn | None:
        for src in type(obj).mro():
            fn = self._by_pair.get((src, dst))
            if fn: return fn
        for pred, d, fn in self._preds:
            if d is dst:
                try:
                    if pred(obj): return fn
                except Exception:
                    pass
        return None
    
    def get_converter(self, src: type, dst: type) -> _ConverterFn | None:
        return self._by_pair.get((src, dst))

REG = ConverterRegistry()