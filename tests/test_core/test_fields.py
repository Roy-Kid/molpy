"""Field vocabulary: molrs keys re-exported by name, and the two formatters."""

import numpy as np
import pytest

import molrs
import molrs.keys

import molpy as mp
from molpy.core import fields
from molpy.core.fields import CHARGE, FieldFormatter, ForceFieldFormatter
from molpy.core.forcefield import BondHarmonicStyle, BondStyle


def test_every_molrs_key_is_a_named_string_constant():
    for name in (n for n in dir(molrs.keys) if n.isupper()):
        key = getattr(molrs.keys, name)
        expected = key.key if hasattr(key, "key") else key
        assert getattr(fields, name) == expected


def test_site_is_a_molpy_owned_column():
    assert fields.SITE == "site"
    assert "SITE" in fields.__all__


class _AcLike(FieldFormatter):
    _field_formatters = {"q": CHARGE}


class TestFieldFormatter:
    def test_canonicalize_and_localize_rename_in_place_and_invert(self):
        block = molrs.Block()
        block["q"] = np.array([0.1, -0.1])
        block["x"] = np.array([0.0, 1.0])
        _AcLike().canonicalize(block)
        assert sorted(block.keys()) == ["charge", "x"]
        _AcLike().localize(block)
        assert sorted(block.keys()) == ["q", "x"]

    def test_canonicalize_frame_walks_every_block(self):
        frame = molrs.Frame()
        frame["atoms"] = {"q": np.array([0.5]), "x": np.array([0.0])}
        assert _AcLike().canonicalize_frame(frame) is frame
        assert "charge" in frame["atoms"] and "q" not in frame["atoms"]

    def test_register_field_extends_the_mapping_at_runtime(self):
        class _Fmt(FieldFormatter):
            _field_formatters = {}

        _Fmt.register_field("qq", CHARGE)
        block = molrs.Block()
        block["qq"] = np.array([1.0])
        _Fmt().canonicalize(block)
        assert "charge" in block


class TestForceFieldFormatter:
    def _harmonic(self):
        ff = mp.ForceField(name="t", units="real")
        style = ff.def_bondstyle("harmonic")
        typ = style.def_type("CT", "CT", k=300.0, r0=1.5)
        return typ, style, ff

    def test_specialised_formatter_beats_the_category_fallback(self):
        class _Fmt(ForceFieldFormatter):
            _param_formatters = {
                BondStyle: lambda t: ["generic"],
                BondHarmonicStyle: lambda t: [t.get("k"), t.get("r0")],
            }

        typ, style, _ = self._harmonic()
        assert _Fmt().format_params(typ, style) == [300.0, 1.5]

    def test_category_fallback_when_no_specialised_match(self):
        class _Fmt(ForceFieldFormatter):
            _param_formatters = {BondStyle: lambda t: ["generic"]}

        typ, style, _ = self._harmonic()
        assert _Fmt().format_params(typ, style) == ["generic"]

    def test_unregistered_style_raises(self):
        class _Fmt(ForceFieldFormatter):
            _param_formatters = {}

        typ, style, _ = self._harmonic()
        with pytest.raises(ValueError, match="No param formatter registered"):
            _Fmt().format_params(typ, style)

    def test_subclasses_get_their_own_registry(self):
        class _Base(ForceFieldFormatter):
            _param_formatters = {}

        class _Child(_Base):
            pass

        _Child.register_param_formatter(BondStyle, lambda t: [1.0])
        assert BondStyle in _Child._param_formatters
        assert BondStyle not in _Base._param_formatters
