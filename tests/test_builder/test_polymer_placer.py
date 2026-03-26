"""Unit tests for placer module: PortInfo, port utilities, and Placer components."""

import numpy as np
import pytest

from molpy.builder.polymer.placer import (
    CovalentSeparator,
    LinearOrienter,
    Placer,
    VdWSeparator,
)
from molpy.builder.polymer.port_utils import (
    PortInfo,
    get_all_ports,
    get_port_atom,
    get_port_info,
    get_ports,
    get_ports_on_node,
)
from molpy.core.atomistic import Atomistic


# ---- PortInfo Tests ----


class TestPortInfo:
    def test_immutable(self):
        struct = Atomistic()
        atom = struct.def_atom(symbol="C")
        port = PortInfo(name=">", target=atom)
        with pytest.raises(AttributeError):
            port.name = "<"  # type: ignore[misc]

    def test_connects_to_right_left(self):
        struct = Atomistic()
        a = struct.def_atom(symbol="C")
        b = struct.def_atom(symbol="C")
        assert PortInfo(">", a).connects_to(PortInfo("<", b))
        assert PortInfo("<", b).connects_to(PortInfo(">", a))

    def test_connects_to_same_name(self):
        struct = Atomistic()
        a = struct.def_atom(symbol="C")
        b = struct.def_atom(symbol="C")
        assert PortInfo("$", a).connects_to(PortInfo("$", b))
        assert PortInfo("$1", a).connects_to(PortInfo("$1", b))

    def test_does_not_connect_incompatible(self):
        struct = Atomistic()
        a = struct.def_atom(symbol="C")
        b = struct.def_atom(symbol="C")
        assert not PortInfo(">", a).connects_to(PortInfo(">", b))
        assert not PortInfo("<", a).connects_to(PortInfo("<", b))
        assert not PortInfo("$1", a).connects_to(PortInfo("$2", b))


# ---- Port Utility Tests ----


class TestPortUtilities:
    def _make_struct(self) -> Atomistic:
        struct = Atomistic()
        struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0, port=">")
        struct.def_atom(symbol="C", x=1.5, y=0.0, z=0.0)
        struct.def_atom(symbol="O", x=3.0, y=0.0, z=0.0, port="<")
        return struct

    def test_get_ports(self):
        struct = self._make_struct()
        ports = get_ports(struct)
        assert ">" in ports
        assert "<" in ports
        assert len(ports[">"]) == 1
        assert len(ports["<"]) == 1

    def test_get_port_atom(self):
        struct = self._make_struct()
        atom = get_port_atom(struct, ">")
        assert atom is not None
        assert atom.get("symbol") == "C"

    def test_get_port_atom_not_found(self):
        struct = self._make_struct()
        assert get_port_atom(struct, "nonexistent") is None

    def test_get_port_info(self):
        struct = self._make_struct()
        info = get_port_info(struct, ">")
        assert info is not None
        assert info.name == ">"
        assert info.target.get("symbol") == "C"

    def test_get_port_info_not_found(self):
        struct = self._make_struct()
        assert get_port_info(struct, "nope") is None

    def test_get_all_ports(self):
        struct = self._make_struct()
        all_ports = get_all_ports(struct)
        assert len(all_ports) == 2
        for name, plist in all_ports.items():
            for p in plist:
                assert isinstance(p, PortInfo)
                assert p.name == name

    def test_get_ports_on_node(self):
        struct = Atomistic()
        a1 = struct.def_atom(symbol="C", port=">", monomer_node_id=0)
        a2 = struct.def_atom(symbol="C", port="<", monomer_node_id=1)

        node0_ports = get_ports_on_node(struct, 0)
        assert ">" in node0_ports
        assert "<" not in node0_ports

        node1_ports = get_ports_on_node(struct, 1)
        assert "<" in node1_ports
        assert ">" not in node1_ports

    def test_get_ports_on_node_empty(self):
        struct = Atomistic()
        struct.def_atom(symbol="C", monomer_node_id=0)
        assert get_ports_on_node(struct, 0) == {}

    def test_multiple_ports_same_name(self):
        struct = Atomistic()
        struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0, port="$")
        struct.def_atom(symbol="C", x=1.0, y=0.0, z=0.0, port="$")
        ports = get_all_ports(struct)
        assert len(ports["$"]) == 2


# ---- Separator Tests ----


class TestVdWSeparator:
    def test_basic_separation(self):
        sep = VdWSeparator(buffer=0.0)
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="C", x=1.0, y=0.0, z=0.0)
        left_port = PortInfo(">", a)
        right_port = PortInfo("<", b)

        dist = sep.get_separation(struct, struct, left_port, right_port)
        assert dist > 0

    def test_buffer(self):
        sep0 = VdWSeparator(buffer=0.0)
        sep1 = VdWSeparator(buffer=1.0)
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="C", x=1.0, y=0.0, z=0.0)
        lp = PortInfo(">", a)
        rp = PortInfo("<", b)

        d0 = sep0.get_separation(struct, struct, lp, rp)
        d1 = sep1.get_separation(struct, struct, lp, rp)
        assert d1 == pytest.approx(d0 + 1.0)


class TestCovalentSeparator:
    def test_cc_bond_length(self):
        sep = CovalentSeparator(buffer=0.0)
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="C", x=1.0, y=0.0, z=0.0)
        lp = PortInfo(">", a)
        rp = PortInfo("<", b)

        dist = sep.get_separation(struct, struct, lp, rp)
        assert dist == pytest.approx(1.54, abs=0.01)

    def test_co_bond_length(self):
        sep = CovalentSeparator(buffer=0.0)
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="O", x=1.0, y=0.0, z=0.0)
        lp = PortInfo(">", a)
        rp = PortInfo("<", b)

        dist = sep.get_separation(struct, struct, lp, rp)
        assert dist == pytest.approx(1.43, abs=0.01)

    def test_buffer(self):
        sep = CovalentSeparator(buffer=-0.1)
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="C", x=1.0, y=0.0, z=0.0)
        lp = PortInfo(">", a)
        rp = PortInfo("<", b)

        dist = sep.get_separation(struct, struct, lp, rp)
        assert dist == pytest.approx(1.54 - 0.1, abs=0.01)

    def test_unknown_elements_fallback(self):
        sep = CovalentSeparator(buffer=0.0)
        struct = Atomistic()
        a = struct.def_atom(symbol="Xe", x=0.0, y=0.0, z=0.0)
        b = struct.def_atom(symbol="Kr", x=1.0, y=0.0, z=0.0)
        lp = PortInfo(">", a)
        rp = PortInfo("<", b)

        dist = sep.get_separation(struct, struct, lp, rp)
        assert dist == pytest.approx(1.54, abs=0.01)  # Default C-C


# ---- Orienter Tests ----


class TestLinearOrienter:
    def _make_struct_with_port(self, x, port_name):
        struct = Atomistic()
        a = struct.def_atom(symbol="C", x=x, y=0.0, z=0.0, port=port_name)
        b = struct.def_atom(symbol="C", x=x - 1.0, y=0.0, z=0.0)
        struct.def_bond(a, b, order=1)
        return struct, PortInfo(port_name, a)

    def test_orientation_returns_correct_shapes(self):
        orienter = LinearOrienter()
        left, lp = self._make_struct_with_port(0.0, ">")
        right, rp = self._make_struct_with_port(5.0, "<")

        translation, rotation = orienter.get_orientation(left, right, lp, rp, 1.5)
        assert translation.shape == (3,)
        assert rotation.shape == (3, 3)

    def test_rotation_is_orthogonal(self):
        orienter = LinearOrienter()
        left, lp = self._make_struct_with_port(0.0, ">")
        right, rp = self._make_struct_with_port(5.0, "<")

        _, rotation = orienter.get_orientation(left, right, lp, rp, 1.5)
        # R @ R^T should be identity
        identity = rotation @ rotation.T
        np.testing.assert_allclose(identity, np.eye(3), atol=1e-10)


# ---- Placer Integration Tests ----


class TestPlacer:
    def test_place_monomer_modifies_right(self):
        left = Atomistic()
        la = left.def_atom(symbol="C", x=0.0, y=0.0, z=0.0, port=">")
        left.def_atom(symbol="C", x=-1.0, y=0.0, z=0.0)
        left.def_bond(list(left.atoms)[0], list(left.atoms)[1], order=1)

        right = Atomistic()
        ra = right.def_atom(symbol="C", x=10.0, y=0.0, z=0.0, port="<")
        right.def_atom(symbol="C", x=11.0, y=0.0, z=0.0)
        right.def_bond(list(right.atoms)[0], list(right.atoms)[1], order=1)

        lp = PortInfo(">", la)
        rp = PortInfo("<", ra)

        placer = Placer(
            separator=CovalentSeparator(buffer=0.0),
            orienter=LinearOrienter(),
        )

        # Record original position
        orig_x = ra["x"]

        placer.place_monomer(left, right, lp, rp)

        # Position should have changed
        assert ra["x"] != orig_x
