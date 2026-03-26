"""Dihedral potentials."""

from .opls import DihedralOPLSStyle, DihedralOPLSType
from .periodic import DihedralPeriodicStyle, DihedralPeriodicType

__all__ = [
    "DihedralOPLSStyle",
    "DihedralOPLSType",
    "DihedralPeriodicStyle",
    "DihedralPeriodicType",
]
