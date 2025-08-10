from abc import ABC, abstractmethod
from dataclasses import dataclass
import numpy as np
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .frame import Block

__all__ = [
    "MaskPredicate",
    "AtomType",
    "AtomIndex",
]

class MaskPredicate(ABC):
    """Boolean mask producer combinable with &, |, ~."""

    @abstractmethod
    def mask(self, block: "Block") -> np.ndarray: ...

    def __call__(self, block: "Block") -> "Block":
        return block[self.mask(block)]

    # compositional logic
    def __and__(self, other: "MaskPredicate") -> "MaskPredicate":
        return _And(self, other)
    def __or__(self, other: "MaskPredicate") -> "MaskPredicate":
        return _Or(self, other)
    def __invert__(self) -> "MaskPredicate":
        return _Not(self)

    __rand__ = __and__
    __ror__ = __or__


@dataclass(frozen=True)
class AtomType(MaskPredicate):
    atom_type: int | str
    field: str = "type"
    def mask(self, block: "Block") -> np.ndarray:  # type: ignore[override]
        return (block[self.field] == self.atom_type)


@dataclass(frozen=True)
class AtomIndex(MaskPredicate):
    indices: list[int]
    id_field: str = "id"
    def mask(self, block: "Block") -> np.ndarray:  # type: ignore[override]
        return np.isin(block[self.id_field], self.indices)


# ------------------------------------------------------------------ combinators
@dataclass(frozen=True)
class _And(MaskPredicate):
    a: MaskPredicate
    b: MaskPredicate
    def mask(self, block: "Block") -> np.ndarray:  # type: ignore[override]
        return self.a.mask(block) & self.b.mask(block)

@dataclass(frozen=True)
class _Or(MaskPredicate):
    a: MaskPredicate
    b: MaskPredicate
    def mask(self, block: "Block") -> np.ndarray:  # type: ignore[override]
        return self.a.mask(block) | self.b.mask(block)

@dataclass(frozen=True)
class _Not(MaskPredicate):
    a: MaskPredicate
    def mask(self, block: "Block") -> np.ndarray:  # type: ignore[override]
        return ~self.a.mask(block)

# Backward compatible alias
Selection = MaskPredicate