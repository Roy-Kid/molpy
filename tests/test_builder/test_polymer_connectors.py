"""Tests for the unified Connector class."""

import pytest

from molpy import Atomistic
from molpy.builder.polymer.connectors import Connector, ConnectorContext
from molpy.builder.polymer.errors import AmbiguousPortsError
from molpy.builder.polymer.port_utils import PortInfo
from molpy.core.atomistic import Atom
from molpy.reacter import Reacter
from molpy.reacter.transformers import form_single_bond


def _make_reacter():
    return Reacter(
        name="test",
        anchor_selector_left=lambda a, port_atom: port_atom,
        anchor_selector_right=lambda a, port_atom: port_atom,
        leaving_selector_left=lambda a, anchor: [],
        leaving_selector_right=lambda a, anchor: [],
        bond_former=form_single_bond,
    )


def _make_monomer_with_ports(port_names, roles=None):
    asm = Atomistic()
    atoms = []
    for i, pname in enumerate(port_names):
        a = Atom(symbol="C")
        a["port"] = pname
        if roles and i < len(roles):
            a["port_role"] = roles[i]
        asm.add_entity(a)
        atoms.append(a)
    return asm, atoms


@pytest.fixture
def reacter():
    return _make_reacter()


@pytest.fixture
def ctx():
    return ConnectorContext(step=0, left_label="A", right_label="B")


class TestConnectorInit:
    def test_default_empty_port_map(self, reacter):
        c = Connector(reacter=reacter)
        assert c.port_map == {}
        assert c.overrides == {}

    def test_with_port_map(self, reacter):
        pm = {("A", "B"): (">", "<")}
        c = Connector(reacter=reacter, port_map=pm)
        assert c.port_map == pm

    def test_with_overrides(self, reacter):
        other = _make_reacter()
        c = Connector(reacter=reacter, overrides={("X", "Y"): other})
        assert c.get_reacter("X", "Y") is other
        assert c.get_reacter("A", "B") is reacter


class TestSelectPortsExplicit:
    def test_explicit_port_map(self, reacter, ctx):
        c = Connector(reacter=reacter, port_map={("A", "B"): (">", "<")})
        left_ports = {">": [PortInfo(">", Atom(symbol="C"))]}
        right_ports = {"<": [PortInfo("<", Atom(symbol="C"))]}
        result = c.select_ports(Atomistic(), Atomistic(), left_ports, right_ports, ctx)
        assert result[0] == ">"
        assert result[2] == "<"

    def test_explicit_port_not_found_raises(self, reacter, ctx):
        c = Connector(reacter=reacter, port_map={("A", "B"): (">", "<")})
        left_ports = {"x": [PortInfo("x", Atom(symbol="C"))]}
        right_ports = {"<": [PortInfo("<", Atom(symbol="C"))]}
        with pytest.raises(AmbiguousPortsError):
            c.select_ports(Atomistic(), Atomistic(), left_ports, right_ports, ctx)


class TestSelectPortsRoleBased:
    def test_role_right_left(self, reacter, ctx):
        c = Connector(reacter=reacter)
        lp = PortInfo(">", Atom(symbol="C"), role="right")
        rp = PortInfo("<", Atom(symbol="C"), role="left")
        result = c.select_ports(Atomistic(), Atomistic(), {">": [lp]}, {"<": [rp]}, ctx)
        assert result[0] == ">"
        assert result[2] == "<"


class TestSelectPortsSinglePort:
    def test_single_port_each_side(self, reacter, ctx):
        c = Connector(reacter=reacter)
        left_ports = {"$": [PortInfo("$", Atom(symbol="C"))]}
        right_ports = {"$": [PortInfo("$", Atom(symbol="C"))]}
        result = c.select_ports(Atomistic(), Atomistic(), left_ports, right_ports, ctx)
        assert result[0] == "$"
        assert result[2] == "$"


class TestSelectPortsAmbiguous:
    def test_ambiguous_raises(self, reacter, ctx):
        c = Connector(reacter=reacter)
        left_ports = {
            "a": [PortInfo("a", Atom(symbol="C"))],
            "b": [PortInfo("b", Atom(symbol="C"))],
        }
        right_ports = {
            "x": [PortInfo("x", Atom(symbol="C"))],
            "y": [PortInfo("y", Atom(symbol="C"))],
        }
        with pytest.raises(AmbiguousPortsError):
            c.select_ports(Atomistic(), Atomistic(), left_ports, right_ports, ctx)
