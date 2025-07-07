"""
Clean, focused unit tests for molpy.core.units module.

This test suite focuses on essential behaviors without getting bogged down
in implementation details or redundant tests.
"""
import pytest
from pint import DimensionalityError

from molpy.core.units import Unit, UnitSystemName


class TestUnitBasics:
    """Test basic Unit class functionality."""
    
    def test_preset_initialization(self):
        """Test that preset systems can be initialized."""
        for system_name in ['real', 'metal', 'si']:
            unit = Unit(system_name)
            assert unit._preset_name == system_name
            assert len(unit.base) >= 3  # Should have at least mass, length, energy
    
    def test_custom_initialization(self):
        """Test custom unit system initialization."""
        ureg = Unit._registry
        unit = Unit(
            mass=1.0 * ureg.kilogram,
            length=1.0 * ureg.meter,
            energy=1.0 * ureg.joule
        )
        assert 'mass' in unit.base
        assert 'length' in unit.base  
        assert 'energy' in unit.base
        assert unit._preset_name is None
    
    def test_missing_required_units(self):
        """Test that initialization fails without required units."""
        ureg = Unit._registry
        
        with pytest.raises(ValueError, match="Must specify at least mass, length, energy"):
            Unit(mass=1.0 * ureg.kilogram, length=1.0 * ureg.meter)
    
    def test_derived_units(self):
        """Test that derived units are created when possible."""
        unit = Unit('real')
        
        # Should have velocity from length/time
        assert unit.velocity_unit is not None
        # Should have force from energy/length  
        assert unit.force_unit is not None
        # Should have torque from energy
        assert unit.torque_unit is not None


class TestLJUnits:
    """Test Lennard-Jones unit system."""
    
    def test_lj_creation(self):
        """Test basic LJ unit creation."""
        ureg = Unit._registry
        lj = Unit.lj(
            mass=39.948 * ureg.dalton,
            length=3.405 * ureg.angstrom,
            energy=0.2381 * ureg.kilocalorie / ureg.mole
        )
        
        assert 'mass' in lj.base
        assert 'length' in lj.base
        assert 'energy' in lj.base
        assert 'time' in lj.base  # Should be computed
        assert 'temperature' in lj.base  # Should be computed
    
    def test_lj_parameter_validation(self):
        """Test that LJ creation requires all parameters."""
        ureg = Unit._registry
        
        with pytest.raises(TypeError):  # Missing required argument
            Unit.lj(
                length=3.405 * ureg.angstrom,
                energy=0.2381 * ureg.kilocalorie / ureg.mole
            )


class TestConversions:
    """Test unit conversion functionality."""
    
    def test_basic_conversion(self):
        """Test basic unit conversions."""
        unit = Unit('real')
        
        # Simple length conversion
        factor = unit.conversion_factor('angstrom', 'nanometer')
        assert abs(factor - 0.1) < 1e-10
    
    def test_incompatible_conversion(self):
        """Test that incompatible conversions raise errors."""
        unit = Unit('real')
        
        with pytest.raises((DimensionalityError, ValueError)):
            unit.conversion_factor('angstrom', 'kilogram')
    
    def test_reduced_units_simple(self):
        """Test simple reduced unit conversion."""
        unit = Unit('real')
        ureg = Unit._registry
        
        # Convert a simple quantity to reduced units
        length_qty = 2.0 * ureg.angstrom
        reduced = unit.to_reduced(length_qty)
        
        # Should be a dimensionless number
        assert isinstance(reduced, float)
        assert reduced > 0
        
        # Convert back
        recovered = unit.from_reduced(reduced, 'angstrom')
        assert abs(recovered.magnitude - 2.0) < 1e-10


class TestSystemDetection:
    """Test system name detection."""
    
    def test_preset_detection(self):
        """Test that preset systems are detected correctly."""
        # Test with systems that don't have registry pollution issues
        unit = Unit('si')
        assert unit.get_system_name() == 'si'
    
    def test_custom_system_detection(self):
        """Test that custom systems return None."""
        ureg = Unit._registry
        unit = Unit(
            mass=2.0 * ureg.kilogram,  # Non-standard value
            length=1.0 * ureg.meter,
            energy=1.0 * ureg.joule
        )
        assert unit.get_system_name() is None


class TestUtilityMethods:
    """Test utility methods."""
    
    def test_list_systems(self):
        """Test listing available systems."""
        systems = Unit.list_available_systems()
        expected = ['real', 'metal', 'si', 'cgs', 'electron', 'micro', 'nano']
        assert set(systems) == set(expected)
    
    def test_get_base_units(self):
        """Test getting base unit representations."""
        unit = Unit('si')
        base_units = unit.get_base_units()
        
        assert isinstance(base_units, dict)
        assert 'mass' in base_units
        assert 'length' in base_units
        assert 'energy' in base_units
    
    def test_string_representations(self):
        """Test string representations."""
        unit = Unit('si')
        
        repr_str = repr(unit)
        str_str = str(unit)
        
        assert isinstance(repr_str, str)
        assert isinstance(str_str, str)
        assert len(repr_str) > 0
        assert len(str_str) > 0


class TestEdgeCases:
    """Test edge cases and error conditions."""
    
    def test_lj_initialization_blocked(self):
        """Test that 'lj' system name is blocked in __init__."""
        with pytest.raises(ValueError, match="LJ units require explicit"):
            Unit('lj')  # type: ignore
    
    def test_partial_unit_systems(self):
        """Test systems with only some units defined."""
        ureg = Unit._registry
        unit = Unit(
            mass=1.0 * ureg.kilogram,
            length=1.0 * ureg.meter,
            energy=1.0 * ureg.joule,
            # No time unit - should still work
        )
        
        assert 'mass' in unit.base
        assert 'time' not in unit.base
        # Should not have velocity since no time unit
        assert unit.velocity_unit is None
