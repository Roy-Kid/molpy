"""
Test suite for molpy.core.struct module.

This test suite provides tests for the base Struct class and MolecularStructure.
Atomic structure tests have been moved to test_atoms.py.
"""

import numpy as np
import pytest

from molpy.core._atomistic import Angle, Atom, Atomistic, Bond, Dihedral


class TestEntity:
    def test_dict_behavior(self):
        e = Atomistic(name="foo", bar=123)
        assert e["name"] == "foo"
        e["baz"] = 456
        assert e["baz"] == 456
        d = e.to_dict()
        assert d["bar"] == 123

    def test_copy(self):
        e = Atomistic(name="foo", bar=123)
        e2 = e.copy(bar=999)
        assert e2["bar"] == 999
        assert e["bar"] == 123
        assert e2 is not e

    def test_call(self):
        e = Atomistic(name="foo", bar=123)
        e2 = e(bar=888)
        assert e2["bar"] == 888
        assert e["bar"] == 123


class TestStruct:
    def test_init_basic(self):
        struct = Atomistic(name="test_struct")
        assert struct["name"] == "test_struct"
        assert "test_struct" in repr(struct)
        unnamed = Atomistic()
        assert "name" not in repr(unnamed)

    def test_copy(self):
        struct = Atomistic(name="foo")
        struct2 = struct.copy(name="bar")
        assert struct2["name"] == "bar"
        assert struct["name"] == "foo"
        assert struct2 is not struct
