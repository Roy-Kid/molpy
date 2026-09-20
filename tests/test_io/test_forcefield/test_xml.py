#!/usr/bin/env python3
"""Unit tests for XML force field reader.

This module contains comprehensive tests for:
- XML force field file reading functionality
- Force field parameter extraction and validation
- Atom types, bond types, angle types, dihedral types, and pair types
- Error handling and edge cases

Uses pytest framework with modern Python 3.10+ type hints and Google-style docstrings.
"""

import math
from pathlib import Path

import pytest

from molpy import ForceField, AtomType, BondType
from molpy.data import get_forcefield_path
from molpy.io.forcefield.xml import XMLForceFieldReader, read_xml_forcefield


class TestXMLForceFieldReader:
    """Test suite for XML force field reader."""

    def test_file_not_found_error(self) -> None:
        """Test that FileNotFoundError is raised for nonexistent files."""
        with pytest.raises(FileNotFoundError):
            read_xml_forcefield("nonexistent_forcefield.xml")

    def test_read_xml_forcefield_convenience_function(
        self, TEST_DATA_DIR: Path
    ) -> None:
        ff = read_xml_forcefield(TEST_DATA_DIR / "xml" / "pf6.xml")
        assert isinstance(ff, ForceField)
        assert {t.name for t in ff.get_types(AtomType)} == {"P", "F1", "F2", "F3"}

    def test_read_xml_forcefield_with_existing_forcefield(
        self, TEST_DATA_DIR: Path
    ) -> None:
        existing_ff = ForceField(name="test", units="real")
        ff = read_xml_forcefield(
            TEST_DATA_DIR / "xml" / "pf6.xml", forcefield=existing_ff
        )
        assert ff is existing_ff
        assert len(ff.get_types(AtomType)) == 4

    def test_parse_bond_with_class_only_creates_wildcard_atomtypes(
        self, tmp_path: Path
    ) -> None:
        """Class-only bond endpoints create wildcard AtomTypes (type="*", class=…).

        A three-element fixture is enough — loading full ``oplsaa.xml`` and
        walking every AtomType.params was multi-second noise.
        """
        xml = tmp_path / "class_bond.xml"
        xml.write_text(
            """<?xml version="1.0"?>
<ForceField>
  <AtomTypes>
    <Type name="opls_c" class="C" element="C" mass="12.01"/>
    <Type name="opls_o" class="O_3" element="O" mass="16.00"/>
  </AtomTypes>
  <HarmonicBondForce>
    <Bond class1="C" class2="O_3" length="0.14" k="1000.0"/>
  </HarmonicBondForce>
</ForceField>
"""
        )
        ff = read_xml_forcefield(xml)

        wildcards = {
            (
                at.params.kwargs.get("class_", ""),
                at.params.kwargs.get("type_", ""),
            ): at
            for at in ff.get_types(AtomType)
        }
        assert ("O_3", "*") in wildcards
        assert ("C", "*") in wildcards
        assert wildcards[("O_3", "*")].name == "O_3"
        assert wildcards[("C", "*")].name == "C"

        # molrs rebuilds endpoint AtomTypes; class-only wildcards use name==class.
        bond_found = False
        for bt in ff.get_types(BondType):
            if {bt.itom.name, bt.jtom.name} == {"C", "O_3"}:
                bond_found = True
                break
        assert bond_found, "C - O_3 bond type should exist and use wildcard AtomTypes"

    def test_class_based_bond_typing_works(self, tmp_path: Path) -> None:
        """Bonds type via class match once wildcard AtomTypes exist."""
        from molpy import Atomistic
        from molpy.typifier import ForceFieldParams

        xml = tmp_path / "class_bond.xml"
        xml.write_text(
            """<?xml version="1.0"?>
<ForceField>
  <AtomTypes>
    <Type name="opls_267" class="C" element="C" mass="12.01"/>
    <Type name="opls_269" class="O_3" element="O" mass="16.00"/>
  </AtomTypes>
  <HarmonicBondForce>
    <Bond class1="C" class2="O_3" length="0.14" k="250000.0"/>
  </HarmonicBondForce>
  <NonbondedForce coulomb14scale="0.5" lj14scale="0.5">
    <Atom type="opls_267" charge="0.0" sigma="0.35" epsilon="0.3"/>
    <Atom type="opls_269" charge="-0.5" sigma="0.3" epsilon="0.7"/>
  </NonbondedForce>
</ForceField>
"""
        )
        ff = read_xml_forcefield(xml)

        asm = Atomistic()
        atom1 = asm.def_atom(symbol="O", type="opls_269")  # class="O_3"
        atom2 = asm.def_atom(symbol="C", type="opls_267")  # class="C"
        bond = asm.def_bond(atom1, atom2)

        typed = ForceFieldParams(ff, strict=False).assign(asm)
        typed_bond = next(iter(typed.bonds))

        assert typed_bond.get("type") is not None, "Bond should have a type assigned"
        assert "k" in typed_bond.data or "r0" in typed_bond.data, (
            "Bond should have parameters"
        )
        assert bond.data == {}, "assign must not mutate its input"


class TestAngleUnitOption:
    """Input angle unit is configurable; internal storage is always radians."""

    def _theta0_rad(self, ff):
        from molpy import AngleStyle

        for style in ff.get_styles(AngleStyle):
            for typ in style.types:
                v = typ.params.kwargs.get("theta0")
                if v:
                    return v
        return None

    def test_radian_input_is_kept_as_the_internal_unit(self):
        """The default (radian) XML input needs no conversion — radians are internal.

        This asserted degrees, which is what let the reader ship a 104.52 that
        molrs's LAMMPS writer then multiplied by 180/π into 5988.55.
        """
        ff = XMLForceFieldReader(
            get_forcefield_path("oplsaa.xml"), angle_unit="radian"
        ).read()
        theta0 = self._theta0_rad(ff)
        assert theta0 is not None
        assert 1.4 < theta0 < math.pi  # radians — the molrs internal unit

    def test_degree_input_is_converted_and_radian_input_is_not(self):
        """A degree file converts at the boundary; a radian file passes through."""
        from molpy.io.forcefield.xml import _angle_to_internal

        assert abs(_angle_to_internal(109.5, "degree") - math.radians(109.5)) < 1e-12
        assert _angle_to_internal(math.radians(109.5), "radian") == math.radians(109.5)

    def test_writer_inverts_to_output_unit(self):
        """Reader(unit) -> internal radians -> writer(unit) round-trips."""
        from molpy.io.forcefield.xml import _angle_from_internal, _angle_to_internal

        for unit in ("radian", "degree"):
            internal = _angle_to_internal(1.91 if unit == "radian" else 109.5, unit)
            back = _angle_from_internal(internal, unit)
            assert abs(back - (1.91 if unit == "radian" else 109.5)) < 1e-9

    def test_invalid_unit_raises(self):
        with pytest.raises(ValueError, match="angle_unit"):
            XMLForceFieldReader(get_forcefield_path("oplsaa.xml"), angle_unit="grad")


class TestAngleUnitDetection:
    """Warnings/errors that flag a likely angle-unit mismatch."""

    def test_degrees_read_as_radians_warns(self):
        """The classic bug: a degree value (104.52) declared radian -> warn."""
        from molpy.io.forcefield.xml import AngleUnitWarning, _normalize_angle

        with pytest.warns(AngleUnitWarning):
            _normalize_angle(104.52, "radian", kind="equilibrium", label="theta0")

    def test_plausible_value_does_not_warn(self):
        """A genuine radian equilibrium angle is silent."""
        import warnings

        from molpy.io.forcefield.xml import _normalize_angle

        with warnings.catch_warnings():
            warnings.simplefilter("error")  # any warning becomes a failure
            rad = _normalize_angle(
                math.radians(109.5), "radian", kind="equilibrium", label="theta0"
            )
        assert abs(rad - math.radians(109.5)) < 1e-9

    def test_phase_out_of_range_warns(self):
        from molpy.io.forcefield.xml import AngleUnitWarning, _normalize_angle

        with pytest.warns(AngleUnitWarning):
            _normalize_angle(400.0, "degree", kind="phase", label="phase1")

    def test_phase_radian_converts_without_warning(self):
        import warnings

        from molpy.io.forcefield.xml import _normalize_angle

        with warnings.catch_warnings():
            warnings.simplefilter("error")
            rad = _normalize_angle(
                math.radians(180.0), "radian", kind="phase", label="phase1"
            )
        assert abs(rad - math.pi) < 1e-9


class TestAbsentChargeAttribute:
    """A `charge` the file never states must not become an explicit ``0.0``.

    TIP3P takes charge from the residue (`UseAttributeFromResidue`), so its
    `<Atom>` entries under `NonbondedForce` carry sigma and epsilon and nothing
    else. Recording a fabricated `charge=0.0` there made
    `ForceFieldParams.assign` overwrite the charges already on the graph.
    """

    def test_nonbonded_type_has_no_charge_when_the_file_states_none(self):
        ff = read_xml_forcefield(get_forcefield_path("tip3p.xml"))
        pairstyle = next(iter(ff.get_styles("pair")))
        for typ in pairstyle.types:
            assert "charge" not in typ.params.kwargs, (
                f"{typ.name} invented a charge the force-field file never gave"
            )

    def test_assign_leaves_the_graphs_own_charges_alone(self):
        import molpy as mp
        from molpy.typifier import ForceFieldParams

        ff = read_xml_forcefield(get_forcefield_path("tip3p.xml"))
        water = mp.Atomistic()
        o = water.def_atom(
            element="O", type="tip3p-O", x=0.0, y=0.0, z=0.0, charge=-0.834
        )
        h1 = water.def_atom(
            element="H", type="tip3p-H", x=0.9572, y=0.0, z=0.0, charge=0.417
        )
        h2 = water.def_atom(
            element="H", type="tip3p-H", x=-0.24, y=0.927, z=0.0, charge=0.417
        )
        water.def_bond(o, h1)
        water.def_bond(o, h2)

        typed = ForceFieldParams(ff).assign(water.get_topo(gen_angle=True))

        assert [a["charge"] for a in typed.atoms] == [-0.834, 0.417, 0.417]
        # File states epsilon in kJ/mol (OpenMM); molrs OPLS path converts to
        # kcal/mol for the real-unit force field surface (÷ 4.184).
        assert typed.atoms[0]["epsilon"] == pytest.approx(0.635968 / 4.184)
