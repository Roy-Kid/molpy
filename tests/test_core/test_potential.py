import numpy as np
import pytest

import molpy as mp
from molpy.core.forcefield import ForceField
from molpy.potential.base import KernelMeta, Potential


class TestRegisterPotential:

    def test_metaclass_register(self):

        class MetaclassRegisterPotential(metaclass=KernelMeta):
            name = "PotentialA"

        assert (
            ForceField._kernel_registry["root"]["PotentialA"]
            is MetaclassRegisterPotential
        )

    def test_base_potential_register(self):

        class BasePotential(Potential):
            name = "BasePotential"

        assert ForceField._kernel_registry["root"]["BasePotential"] is BasePotential

    def test_type_potential_register(self):

        class TypedPotential(Potential):
            name = "harmonic"
            type = "angle"

        assert ForceField._kernel_registry["angle"]["harmonic"] is TypedPotential
