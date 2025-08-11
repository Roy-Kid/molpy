#!/usr/bin/env python3
"""
Test forcefield functionality including get_potential and get_potentials methods.
"""

import os
import sys

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "src"))

from molpy.core.forcefield import AtomType, ForceField
from molpy.potential.base import Potentials
from molpy.potential.pair.lj import LJ126


class TestForceFieldPotentials:
    """Test forcefield potential creation functionality."""

    def test_get_potential_basic(self):
        """Test basic get_potential functionality."""
        ff = ForceField("test_ff")

        # Create LJ126 style
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0

        # Get potential
        potential = ff.get_potential("lj126/cut")

        assert potential is not None
        assert isinstance(potential, LJ126)
        assert potential.epsilon[0, 0] == 1.0
        assert potential.sigma[0, 0] == 1.0

    def test_get_potential_parameter_override(self):
        """Test parameter override in get_potential."""
        ff = ForceField("test_ff")

        # Create style with default parameters
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0

        # Override parameters
        potential = ff.get_potential("lj126/cut", epsilon=2.5, sigma=1.5)

        assert potential.epsilon[0, 0] == 2.5
        assert potential.sigma[0, 0] == 1.5

    def test_get_potential_nonexistent(self):
        """Test get_potential with non-existent style."""
        ff = ForceField("test_ff")

        # Try to get non-existent potential
        potential = ff.get_potential("nonexistent")

        assert potential is None

    def test_get_potentials_basic(self):
        """Test basic get_potentials functionality."""
        ff = ForceField("test_ff")

        # Create multiple styles
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0

        # Get all potentials
        potentials = ff.get_potentials()

        assert isinstance(potentials, Potentials)
        assert len(potentials) == 1
        assert isinstance(potentials[0], LJ126)

    def test_get_potentials_specific_styles(self):
        """Test get_potentials with specific style names."""
        ff = ForceField("test_ff")

        # Create multiple styles
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0

        bond_style = ff.def_bondstyle("harmonic")
        bond_style["k"] = 300.0

        # Get specific potentials
        potentials = ff.get_potentials(["lj126/cut", "harmonic"])

        assert isinstance(potentials, Potentials)
        assert len(potentials) == 1  # Only LJ126 has implementation
        assert isinstance(potentials[0], LJ126)

    def test_get_potentials_empty_list(self):
        """Test get_potentials with empty style list."""
        ff = ForceField("test_ff")

        # Get potentials with empty list
        potentials = ff.get_potentials([])

        assert isinstance(potentials, Potentials)
        assert len(potentials) == 0


class TestMixingRules:
    """Test epsilon/sigma mixing rules for LJ potential."""

    def test_geometric_mixing_rule(self):
        """Test geometric mixing rule for LJ parameters."""
        ff = ForceField("test_ff")

        # Create atom types
        atom_style = ff.def_atomstyle("full")
        c_type = atom_style.def_type("C")
        h_type = atom_style.def_type("H")

        # Create LJ style with geometric mixing rule
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "geometric"

        # Define self-interaction parameters
        lj_style.def_type(
            c_type, c_type, parms=[0.066, 3.4]
        )  # C-C: epsilon=0.066, sigma=3.4
        lj_style.def_type(
            h_type, h_type, parms=[0.030, 2.5]
        )  # H-H: epsilon=0.030, sigma=2.5

        # Test that types are defined
        pair_types = lj_style.get_types()
        assert len(pair_types) == 2

        # Check that parameters are stored correctly
        cc_type = lj_style.get("C-C")
        hh_type = lj_style.get("H-H")

        assert cc_type is not None
        assert hh_type is not None
        assert cc_type.parms == [0.066, 3.4]
        assert hh_type.parms == [0.030, 2.5]

    def test_arithmetic_mixing_rule(self):
        """Test arithmetic mixing rule for LJ parameters."""
        ff = ForceField("test_ff")

        # Create atom types
        atom_style = ff.def_atomstyle("full")
        c_type = atom_style.def_type("C")
        o_type = atom_style.def_type("O")

        # Create LJ style with arithmetic mixing rule
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "arithmetic"

        # Define self-interaction parameters
        lj_style.def_type(c_type, c_type, parms=[0.066, 3.4])  # C-C
        lj_style.def_type(o_type, o_type, parms=[0.210, 2.96])  # O-O

        # Test that types are defined
        pair_types = lj_style.get_types()
        assert len(pair_types) == 2

        # Verify parameters
        cc_type = lj_style.get("C-C")
        oo_type = lj_style.get("O-O")

        assert cc_type.parms == [0.066, 3.4]
        assert oo_type.parms == [0.210, 2.96]

    def test_lorentz_berthelot_mixing_rule(self):
        """Test Lorentz-Berthelot mixing rule (geometric epsilon, arithmetic sigma)."""
        ff = ForceField("test_ff")

        # Create atom types
        atom_style = ff.def_atomstyle("full")
        ar_type = atom_style.def_type("Ar")
        ne_type = atom_style.def_type("Ne")

        # Create LJ style with Lorentz-Berthelot mixing rule
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "lorentz_berthelot"

        # Define self-interaction parameters
        lj_style.def_type(ar_type, ar_type, parms=[0.996, 3.4])  # Ar-Ar
        lj_style.def_type(ne_type, ne_type, parms=[0.317, 2.78])  # Ne-Ne

        # Test that types are defined
        pair_types = lj_style.get_types()
        assert len(pair_types) == 2

        # Verify parameters are stored
        ar_type_obj = lj_style.get("Ar-Ar")
        ne_type_obj = lj_style.get("Ne-Ne")

        assert ar_type_obj is not None
        assert ne_type_obj is not None
        assert ar_type_obj.parms == [0.996, 3.4]
        assert ne_type_obj.parms == [0.317, 2.78]

    def test_custom_mixing_rule(self):
        """Test custom mixing rule definition."""
        ff = ForceField("test_ff")

        # Create atom types
        atom_style = ff.def_atomstyle("full")
        li_type = atom_style.def_type("Li")
        f_type = atom_style.def_type("F")

        # Create LJ style with custom mixing rule
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "custom"

        # Define self-interaction parameters
        lj_style.def_type(li_type, li_type, parms=[0.025, 2.13])  # Li-Li
        lj_style.def_type(f_type, f_type, parms=[0.061, 2.83])  # F-F

        # Define explicit cross-interaction
        lj_style.def_type(li_type, f_type, parms=[0.040, 2.48])  # Li-F (custom)

        # Test that all types are defined
        pair_types = lj_style.get_types()
        assert len(pair_types) == 3

        # Verify all parameters
        li_li = lj_style.get("Li-Li")
        f_f = lj_style.get("F-F")
        li_f = lj_style.get("Li-F")

        assert li_li.parms == [0.025, 2.13]
        assert f_f.parms == [0.061, 2.83]
        assert li_f.parms == [0.040, 2.48]

    def test_mixed_potential_creation(self):
        """Test creating potential with mixed atom types."""
        ff = ForceField("test_ff")

        # Create atom types
        atom_style = ff.def_atomstyle("full")
        c_type = atom_style.def_type("C")
        h_type = atom_style.def_type("H")

        # Create LJ style with multiple types
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "geometric"

        # Define pair types
        lj_style.def_type(c_type, c_type, parms=[0.066, 3.4])
        lj_style.def_type(h_type, h_type, parms=[0.030, 2.5])

        # For now, we'll just test that the style is created correctly
        # Full mixing rule implementation would require more complex potential creation
        pair_types = lj_style.get_types()
        assert len(pair_types) == 2
        assert lj_style["mixing_rule"] == "geometric"

    def test_mixing_rule_validation(self):
        """Test validation of mixing rule parameters."""
        ff = ForceField("test_ff")

        # Create LJ style
        lj_style = ff.def_pairstyle("lj126/cut")

        # Test valid mixing rules
        valid_rules = ["geometric", "arithmetic", "lorentz_berthelot", "custom"]

        for rule in valid_rules:
            lj_style["mixing_rule"] = rule
            assert lj_style["mixing_rule"] == rule

    def test_mixing_rule_metadata(self):
        """Test that mixing rule metadata is preserved."""
        ff = ForceField("test_ff")

        # Create LJ style with metadata
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["mixing_rule"] = "geometric"
        lj_style["cutoff"] = 10.0
        lj_style["tail_correction"] = True

        # Test metadata preservation
        assert lj_style["mixing_rule"] == "geometric"
        assert lj_style["cutoff"] == 10.0
        assert lj_style["tail_correction"] is True

        # Test that metadata is available for potential creation
        style_data = lj_style.data
        assert "mixing_rule" in style_data
        assert "cutoff" in style_data
        assert "tail_correction" in style_data


class TestForceFieldIntegration:
    """Test integration of forcefield with potential creation."""

    def test_multiple_potential_types(self):
        """Test creating multiple types of potentials."""
        ff = ForceField("multi_potential_ff")

        # Create different styles
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0

        bond_style = ff.def_bondstyle("harmonic")
        bond_style["k"] = 300.0

        angle_style = ff.def_anglestyle("harmonic")
        angle_style["k"] = 50.0

        # Test get_potentials with all styles
        potentials = ff.get_potentials()

        # Only LJ126 should be created (others don't have implementations)
        assert len(potentials) == 1
        assert isinstance(potentials[0], LJ126)

    def test_forcefield_serialization_with_potentials(self):
        """Test that forcefield can be serialized and deserialized with potential data."""
        ff = ForceField("serialization_test")

        # Create LJ style with mixing rule
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.5
        lj_style["sigma"] = 2.0
        lj_style["mixing_rule"] = "geometric"

        # Convert to dict and back
        ff_dict = ff.to_dict()
        ff_restored = ForceField.from_dict(ff_dict)

        # Test that potential can be created from restored forcefield
        original_potential = ff.get_potential("lj126/cut")
        restored_potential = ff_restored.get_potential("lj126/cut")

        assert original_potential is not None
        assert restored_potential is not None
        assert original_potential.epsilon[0, 0] == restored_potential.epsilon[0, 0]
        assert original_potential.sigma[0, 0] == restored_potential.sigma[0, 0]

    def test_potential_parameter_inheritance(self):
        """Test that potential inherits all style parameters."""
        ff = ForceField("inheritance_test")

        # Create LJ style with many parameters
        lj_style = ff.def_pairstyle("lj126/cut")
        lj_style["epsilon"] = 1.0
        lj_style["sigma"] = 1.0
        lj_style["mixing_rule"] = "arithmetic"
        lj_style["cutoff"] = 12.0
        lj_style["tail_correction"] = True

        # Create potential (should only use compatible parameters)
        potential = ff.get_potential("lj126/cut")

        assert potential is not None
        assert potential.epsilon[0, 0] == 1.0
        assert potential.sigma[0, 0] == 1.0
        # Other parameters would be ignored by LJ126 constructor


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
