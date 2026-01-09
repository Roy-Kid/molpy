import os
import tempfile
import math

import pytest

import molpy as mp
from molpy import AngleStyle, AtomStyle, BondStyle, DihedralStyle, PairStyle
from molpy.io.forcefield import LAMMPSForceFieldWriter, ParameterFormatterRegistry
from molpy.potential.bond import BondHarmonicStyle


# ===================================================================
#               Basic Writer Tests
# ===================================================================


class TestLAMMPSForceFieldWriter:
    """Tests for basic LAMMPSForceFieldWriter functionality."""

    def test_full_forcefield_write(self):
        """Test writing a complete force field with all style types."""
        ff = mp.ForceField("testff")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_H = atomstyle.def_type("H", mass=1.008)

        bondstyle = ff.def_style(BondStyle("harmonic"))
        bondstyle.def_type(atomtype_C, atomtype_H, k=100.0, r0=1.09)

        anglestyle = ff.def_style(AngleStyle("harmonic"))
        anglestyle.def_type(atomtype_C, atomtype_H, atomtype_C, k=50.0, theta0=109.5)

        dihedralstyle = ff.def_style(DihedralStyle("opls"))
        dihedralstyle.def_type(
            atomtype_C,
            atomtype_H,
            atomtype_C,
            atomtype_H,
            k1=0.5,
            k2=1.0,
            k3=0.0,
            k4=0.0,
        )

        pairstyle = ff.def_style(PairStyle("lj/cut"))
        pairstyle.def_type(atomtype_C, atomtype_C, epsilon=0.2, sigma=3.4)
        pairstyle.def_type(atomtype_H, atomtype_H, epsilon=0.05, sigma=2.5)
        pairstyle.def_type(atomtype_C, atomtype_H, epsilon=0.1, sigma=3.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "lammps.ff")
            writer = LAMMPSForceFieldWriter(outpath)
            writer.write(ff)
            with open(outpath) as f:
                content = f.read()
                print(content)
            assert "bond_style harmonic" in content
            assert "angle_style harmonic" in content
            assert "dihedral_style opls" in content
            assert "pair_style lj/cut" in content
            assert "bond_coeff" in content
            assert "angle_coeff" in content
            assert "dihedral_coeff" in content
            assert "pair_coeff" in content
            assert "C" in content and "H" in content
            assert "bond_coeff" in content
            assert "50" in content and "109.5" in content
            assert "0.5" in content and "1" in content
            assert "pair_coeff C" in content and "3.4" in content
            assert "pair_coeff H" in content and "2.5" in content
            assert "pair_coeff C H" in content or "pair_coeff H C" in content

    def test_angle_radians_to_degrees_conversion(self):
        """Test that angle theta0 is correctly converted from radians to degrees.

        Force field stores angles internally in radians, but LAMMPS requires degrees.
        This test verifies the conversion is performed correctly.
        """
        from molpy.potential.angle import AngleHarmonicStyle
        import re

        ff = mp.ForceField("test_angle_units")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_O = atomstyle.def_type("O", mass=16.00)

        # Create angle style with theta0 in RADIANS (as stored internally)
        theta0_radians = 1.9111355  # ~109.5 degrees in radians
        expected_degrees = math.degrees(theta0_radians)

        anglestyle = ff.def_style(AngleHarmonicStyle())
        anglestyle.def_type(
            atomtype_C, atomtype_O, atomtype_C, k=50.0, theta0=theta0_radians
        )

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "angle_test.ff")
            writer = LAMMPSForceFieldWriter(outpath)
            writer.write(ff)

            with open(outpath) as f:
                content = f.read()

            print(f"Written content:\n{content}")

            # Extract the theta0 value from the angle_coeff line
            match = re.search(r"angle_coeff\s+\S+\s+([\d.]+)\s+([\d.]+)", content)
            assert match, f"Could not find angle_coeff in output: {content}"

            written_k = float(match.group(1))
            written_theta0 = float(match.group(2))

            print(f"Input theta0: {theta0_radians} rad = {expected_degrees}°")
            print(f"Written theta0: {written_theta0}°")

            # Verify the theta0 was converted to degrees
            assert abs(written_theta0 - expected_degrees) < 0.1, (
                f"theta0 not converted correctly: expected ~{expected_degrees:.1f}° "
                f"but got {written_theta0:.1f}°. "
                f"Input was {theta0_radians} rad."
            )

            # Also verify k is preserved
            assert (
                abs(written_k - 50.0) < 0.01
            ), f"k value changed: expected 50.0, got {written_k}"


# ===================================================================
#               Formatter Registry Tests
# ===================================================================


class CustomBondStyle(BondStyle):
    """Custom bond style for testing formatter registration."""

    def __init__(self):
        super().__init__("custom")


class MockBondType:
    """Mock bond type for testing."""

    def __init__(self, k: float, r0: float, custom: float = 0.0):
        self.params = type(
            "Params", (), {"kwargs": {"k": k, "r0": r0, "custom": custom}}
        )()
        self.itom = type("AtomType", (), {"name": "C"})()
        self.jtom = type("AtomType", (), {"name": "H"})()
        self.name = "C-H"


class TestLAMMPSForceFieldWriterFormatters:
    """Tests for LAMMPSForceFieldWriter formatter registry integration."""

    def test_has_class_level_registry(self):
        """Test that LAMMPSForceFieldWriter has a class-level formatter registry."""
        assert hasattr(LAMMPSForceFieldWriter, "_formatters")
        assert isinstance(
            LAMMPSForceFieldWriter._formatters, ParameterFormatterRegistry
        )

    def test_default_formatters_registered(self):
        """Test that default formatters are automatically registered."""
        from molpy.potential.angle import AngleHarmonicStyle
        from molpy.potential.dihedral import DihedralOPLSStyle
        from molpy.potential.pair import PairLJ126CoulCutStyle

        formatters = LAMMPSForceFieldWriter._formatters.list_formatters()

        # Check specialized formatters
        assert BondHarmonicStyle in formatters
        assert AngleHarmonicStyle in formatters
        assert DihedralOPLSStyle in formatters
        assert PairLJ126CoulCutStyle in formatters

        # Check generic fallback formatters
        assert BondStyle in formatters
        assert AngleStyle in formatters
        assert DihedralStyle in formatters
        assert PairStyle in formatters

    def test_register_custom_formatter(self):
        """Test registering a custom formatter for a custom style."""

        def custom_formatter(typ):
            return [
                typ.params.kwargs["k"],
                typ.params.kwargs["r0"],
                typ.params.kwargs["custom"],
            ]

        # Register custom formatter
        LAMMPSForceFieldWriter.register_formatter(CustomBondStyle, custom_formatter)

        # Verify it's registered
        formatter = LAMMPSForceFieldWriter._formatters.get_formatter(CustomBondStyle)
        assert formatter is custom_formatter

        # Test formatting
        mock_type = MockBondType(k=300.0, r0=1.5, custom=999.0)
        params = formatter(mock_type)
        assert params == [300.0, 1.5, 999.0]

    def test_instance_uses_class_registry(self):
        """Test that writer instances use the class-level registry by default."""
        writer = LAMMPSForceFieldWriter("/tmp/test.lmp")

        # Should use class-level registry
        assert writer.formatters is LAMMPSForceFieldWriter._formatters

    def test_instance_level_registry_override(self):
        """Test that writer instances can override with custom registry."""
        # Create custom registry
        custom_registry = ParameterFormatterRegistry()

        def custom_formatter(typ):
            return [999.0, 888.0]

        custom_registry.register(CustomBondStyle, custom_formatter)

        # Create writer with custom registry
        writer = LAMMPSForceFieldWriter(
            "/tmp/test.lmp", formatter_registry=custom_registry
        )

        # Should use instance registry
        assert writer.formatters is custom_registry
        assert writer.formatters is not LAMMPSForceFieldWriter._formatters

        # Verify custom formatter is used
        formatter = writer.formatters.get_formatter(CustomBondStyle)
        assert formatter is custom_formatter

    def test_formats_parameters_correctly(self):
        """Test that writer correctly formats parameters using the registry."""
        ff = mp.ForceField("testff")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_H = atomstyle.def_type("H", mass=1.008)

        bondstyle = ff.def_style(BondHarmonicStyle())
        bond_type = bondstyle.def_type(atomtype_C, atomtype_H, k=300.0, r0=1.5)

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "test.ff")
            writer = LAMMPSForceFieldWriter(outpath)
            writer.write(ff)

            with open(outpath) as f:
                content = f.read()

            # Verify bond parameters are formatted correctly
            assert "bond_coeff" in content
            assert "300" in content
            assert "1.5" in content

    def test_missing_formatter_error(self):
        """Test that helpful error is raised when formatter is missing."""

        class UnregisteredStyle(BondStyle):
            def __init__(self):
                super().__init__("unregistered")

        ff = mp.ForceField("testff")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_H = atomstyle.def_type("H", mass=1.008)

        # Create style without registered formatter
        unregistered_style = ff.def_style(UnregisteredStyle())
        unregistered_style.def_type(atomtype_C, atomtype_H, k=100.0, r0=1.0)

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "test.ff")
            writer = LAMMPSForceFieldWriter(outpath)

            # Should raise helpful error
            with pytest.raises(ValueError, match="Failed to format parameters"):
                writer.write(ff)

    def test_custom_formatter_integration(self):
        """Test end-to-end integration with custom formatter."""

        class SpecialBondStyle(BondStyle):
            def __init__(self):
                super().__init__("special")

        def special_formatter(typ):
            # Custom format: multiply k by 2
            return [typ.params.kwargs["k"] * 2, typ.params.kwargs["r0"]]

        # Register custom formatter
        LAMMPSForceFieldWriter.register_formatter(SpecialBondStyle, special_formatter)

        ff = mp.ForceField("testff")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_H = atomstyle.def_type("H", mass=1.008)

        special_style = ff.def_style(SpecialBondStyle())
        special_style.def_type(atomtype_C, atomtype_H, k=100.0, r0=1.5)

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "test.ff")
            writer = LAMMPSForceFieldWriter(outpath)
            writer.write(ff)

            with open(outpath) as f:
                content = f.read()

            # Verify custom formatter was used (k should be 200.0 = 100.0 * 2)
            assert "bond_coeff" in content
            assert "200" in content  # k * 2
            assert "1.5" in content  # r0 unchanged

    def test_registry_isolation(self):
        """Test that different writer classes have isolated registries."""
        # This test verifies the concept even though we only have LAMMPSForceFieldWriter
        # If we had GROMACSForceFieldWriter, they should have separate registries

        original_formatters = LAMMPSForceFieldWriter._formatters.list_formatters()

        # Register a custom formatter
        LAMMPSForceFieldWriter.register_formatter(CustomBondStyle, lambda typ: [1.0])

        # Verify it's in the registry
        assert CustomBondStyle in LAMMPSForceFieldWriter._formatters.list_formatters()

        # Original formatters should still be there
        for fmt in original_formatters:
            assert fmt in LAMMPSForceFieldWriter._formatters.list_formatters()

    def test_backward_compatibility(self):
        """Test that existing code works without changes (backward compatibility)."""
        ff = mp.ForceField("testff")
        atomstyle = ff.def_style(AtomStyle("full"))
        atomtype_C = atomstyle.def_type("C", mass=12.01)
        atomtype_H = atomstyle.def_type("H", mass=1.008)

        bondstyle = ff.def_style(BondStyle("harmonic"))
        bondstyle.def_type(atomtype_C, atomtype_H, k=100.0, r0=1.09)

        with tempfile.TemporaryDirectory() as tmpdir:
            outpath = os.path.join(tmpdir, "test.ff")

            # Old-style usage: no formatter_registry parameter
            writer = LAMMPSForceFieldWriter(outpath)
            writer.write(ff)

            with open(outpath) as f:
                content = f.read()

            # Should work exactly as before
            assert "bond_style harmonic" in content
            assert "bond_coeff" in content
            assert "100" in content
            assert "1.09" in content
