"""
Phase 1 Tests: Topology-only linear polymer assembly.

These tests focus on connector logic and port selection,
not on 3D geometry or coordinates.
"""

import pytest
from molpy.builder.polymer import (
    PolymerBuilder,
    AutoConnector,
    TableConnector,
    CallbackConnector,
    ChainConnector,
)
from molpy.builder.errors import (
    SequenceError,
    AmbiguousPortsError,
    BondKindConflictError,
    PortReuseError,
)
from molpy.core.wrappers.monomer import Monomer
from molpy.core.atomistic import Atomistic, Atom


# ========== Test Fixtures ==========

@pytest.fixture
def simple_monomer_with_roles():
    """
    Create a simple monomer with 'left' and 'right' role ports.
    Simulates BigSMILES [<]CC[>] pattern.
    """
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    atomistic.add_bond(c1, c2)
    
    monomer = Monomer(atomistic)
    monomer.set_port("left", c1, role="left")
    monomer.set_port("right", c2, role="right")
    
    return monomer


@pytest.fixture
def monomer_A_roles():
    """Monomer A with left/right roles."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    o = atomistic.add_atom(symbol="O", xyz=[3, 0, 0])
    atomistic.add_bond(c1, c2)
    atomistic.add_bond(c2, o)
    
    monomer = Monomer(atomistic)
    monomer.set_port("left", c1, role="left")
    monomer.set_port("right", o, role="right")
    
    return monomer


@pytest.fixture
def monomer_B_roles():
    """Monomer B with left/right roles."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    c3 = atomistic.add_atom(symbol="C", xyz=[1.5, 1.5, 0])
    o = atomistic.add_atom(symbol="O", xyz=[3, 0, 0])
    atomistic.add_bond(c1, c2)
    atomistic.add_bond(c2, c3)
    atomistic.add_bond(c2, o)
    
    monomer = Monomer(atomistic)
    monomer.set_port("left", c1, role="left")
    monomer.set_port("right", o, role="right")
    
    return monomer


@pytest.fixture
def terminus_monomer():
    """Terminus with single port."""
    atomistic = Atomistic()
    h = atomistic.add_atom(symbol="H", xyz=[0, 0, 0])
    o = atomistic.add_atom(symbol="O", xyz=[1, 0, 0])
    atomistic.add_bond(h, o)
    
    monomer = Monomer(atomistic)
    monomer.set_port("t", o, role="left")  # Connects to 'right' of chain
    
    return monomer


@pytest.fixture
def numeric_port_A():
    """Monomer with numeric port '1' that can be used twice."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    o = atomistic.add_atom(symbol="O", xyz=[3, 0, 0])
    atomistic.add_bond(c1, c2)
    atomistic.add_bond(c2, o)
    
    monomer = Monomer(atomistic)
    monomer.set_port("1", o, multiplicity=2)  # Can be used twice for middle monomers
    
    return monomer


@pytest.fixture
def numeric_port_B_two_ports():
    """Monomer B with two numeric ports '2' and '3'."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    c3 = atomistic.add_atom(symbol="C", xyz=[1.5, 1.5, 0])
    o = atomistic.add_atom(symbol="O", xyz=[3, 0, 0])
    atomistic.add_bond(c1, c2)
    atomistic.add_bond(c2, c3)
    atomistic.add_bond(c2, o)
    
    monomer = Monomer(atomistic)
    monomer.set_port("2", c3)
    monomer.set_port("3", o)
    
    return monomer


# ========== Test 1: Auto Connector with BigSMILES Roles ==========

def test_linear_auto_bigsmiles_roles(monomer_A_roles, monomer_B_roles, terminus_monomer):
    """
    Test automatic connection using BigSMILES-style roles.
    Monomers with 'left' and 'right' roles should connect automatically.
    """
    library = {
        "A": monomer_A_roles,
        "B": monomer_B_roles,
        "T": terminus_monomer,
    }
    
    poly = PolymerBuilder.linear(
        sequence="TABBA",
        library=library,
        connector=AutoConnector(),
    )
    
    # Should have successfully connected
    assert poly is not None
    assert len(poly.unwrap().atoms) > 0
    
    # Check that we have bonds (5 labels → 4 connections + internal bonds)
    assert len(poly.unwrap().bonds) >= 4


def test_linear_auto_simple_sequence(simple_monomer_with_roles):
    """Test simple AA sequence with auto connector."""
    library = {"A": simple_monomer_with_roles}
    
    poly = PolymerBuilder.linear(
        sequence="AAA",
        library=library,
        connector=AutoConnector(),
    )
    
    # 3 monomers with 2 atoms each = 6 atoms
    # (assuming deep copy works correctly)
    assert len(poly.unwrap().atoms) == 6
    
    # 3 monomers × 1 internal bond + 2 inter-monomer bonds = 5 bonds
    assert len(poly.unwrap().bonds) == 5


# ========== Test 2: Table Connector with Numeric Ports ==========

def test_linear_table_numeric_ports(numeric_port_A, numeric_port_B_two_ports, terminus_monomer):
    """
    Test table connector with explicit port mappings.
    B has two ports, so AutoConnector would be ambiguous.
    """
    library = {
        "A": numeric_port_A,
        "B": numeric_port_B_two_ports,
        "T": terminus_monomer,
    }
    
    rules = {
        ("T", "A"): ("t", "1"),
        ("A", "B"): ("1", "2"),
        ("B", "B"): ("3", "2"),  # B-B connection
        ("B", "A"): ("3", "1"),
    }
    
    poly = PolymerBuilder.linear(
        sequence="TABBA",
        library=library,
        connector=TableConnector(rules),
    )
    
    assert poly is not None
    assert len(poly.unwrap().atoms) > 0


# ========== Test 3: Ambiguous Ports Without Rules ==========

def test_ambiguous_without_rules(numeric_port_B_two_ports):
    """
    Test that ambiguous ports raise error with AutoConnector.
    B has two ports and no roles → AutoConnector should fail.
    """
    library = {"B": numeric_port_B_two_ports}
    
    with pytest.raises(AmbiguousPortsError, match="Cannot auto-select ports"):
        PolymerBuilder.linear(
            sequence="BB",
            library=library,
            connector=AutoConnector(),
        )


# ========== Test 4: Bond Kind Resolution Priority ==========

def test_bond_kind_resolution_global_override():
    """Test that global bond_kind overrides everything."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    atomistic.add_bond(c1, c2)
    
    monomer = Monomer(atomistic)
    monomer.set_port("left", c1, role="left", bond_kind="-")
    monomer.set_port("right", c2, role="right", bond_kind="-")
    
    library = {"A": monomer}
    details = {}
    
    poly = PolymerBuilder.linear(
        sequence="AAA",
        library=library,
        connector=AutoConnector(),
        bond_kind="=",  # Global override
        details=details,
    )
    
    # All connections should use '=' bond
    for conn in details["connections"]:
        assert conn["bond_kind"] == "="


def test_bond_kind_conflict_raises():
    """Test that conflicting port bond kinds raise error."""
    atomistic_a = Atomistic()
    c1 = atomistic_a.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic_a.add_atom(symbol="C", xyz=[1.5, 0, 0])
    atomistic_a.add_bond(c1, c2)
    
    atomistic_b = Atomistic()
    c3 = atomistic_b.add_atom(symbol="C", xyz=[0, 0, 0])
    c4 = atomistic_b.add_atom(symbol="C", xyz=[1.5, 0, 0])
    atomistic_b.add_bond(c3, c4)
    
    mono_a = Monomer(atomistic_a)
    mono_a.set_port("left", c1, role="left")
    mono_a.set_port("right", c2, role="right", bond_kind="-")
    
    mono_b = Monomer(atomistic_b)
    mono_b.set_port("left", c3, role="left", bond_kind="=")  # Conflict!
    mono_b.set_port("right", c4, role="right")
    
    library = {"A": mono_a, "B": mono_b}
    
    with pytest.raises(BondKindConflictError, match="Conflicting bond kinds"):
        PolymerBuilder.linear(
            sequence="AB",
            library=library,
            connector=AutoConnector(),
        )


# ========== Test 5: Port Consumption and Reuse Prevention ==========

def test_consumption_forbids_reuse():
    """Test that consumed ports (multiplicity=0) cannot be reused."""
    atomistic = Atomistic()
    c = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    
    monomer = Monomer(atomistic)
    monomer.set_port("single", c, multiplicity=1)  # Can only be used once
    
    library = {"A": monomer}
    
    # This should fail: trying to connect AAA but each A only has one port
    # After first connection, left A's port is consumed
    with pytest.raises(Exception):  # Could be NoCompatiblePortsError or PortReuseError
        PolymerBuilder.linear(
            sequence="AAA",
            library=library,
            connector=AutoConnector(),
        )


# ========== Test 6: Deep Copy - No Reference to Source ==========

def test_deep_copy_no_reference(simple_monomer_with_roles):
    """Test that modifying source monomer after build doesn't affect polymer."""
    library = {"A": simple_monomer_with_roles}
    
    poly = PolymerBuilder.linear(
        sequence="AA",
        library=library,
        connector=AutoConnector(),
    )
    
    original_atom_count = len(poly.unwrap().atoms)
    
    # Modify the source monomer
    simple_monomer_with_roles.unwrap().add_atom(symbol="N", xyz=[5, 5, 5])
    
    # Polymer should be unchanged
    assert len(poly.unwrap().atoms) == original_atom_count


# ========== Test 7: Sequence Validation ==========

def test_sequence_too_short():
    """Test that sequence with < 2 labels raises SequenceError."""
    library = {"A": simple_monomer_with_roles}
    
    with pytest.raises(SequenceError, match="at least 2 labels"):
        PolymerBuilder.linear(
            sequence="A",  # Too short
            library=library,
            connector=AutoConnector(),
        )


def test_sequence_unknown_label(simple_monomer_with_roles):
    """Test that unknown label in sequence raises SequenceError."""
    library = {"A": simple_monomer_with_roles}
    
    with pytest.raises(SequenceError, match="not found in library"):
        PolymerBuilder.linear(
            sequence="AXA",  # X not in library
            library=library,
            connector=AutoConnector(),
        )


# ========== Test 8: Chain Connector ==========

def test_chain_connector_fallback(numeric_port_A, numeric_port_B_two_ports):
    """Test ChainConnector trying table first, then fallback."""
    library = {
        "A": numeric_port_A,
        "B": numeric_port_B_two_ports,
    }
    
    # Only define one rule, let AutoConnector handle the rest
    rules = {("A", "B"): ("1", "2")}
    
    connector = ChainConnector([
        TableConnector(rules),
        AutoConnector(),  # Will fail for B-B but that's ok for this test
    ])
    
    # This should work: A-B uses table rule
    poly = PolymerBuilder.linear(
        sequence="AB",
        library=library,
        connector=connector,
    )
    
    assert poly is not None


# ========== Test 9: Callback Connector ==========

def test_callback_connector():
    """Test custom callback for port selection."""
    atomistic = Atomistic()
    c1 = atomistic.add_atom(symbol="C", xyz=[0, 0, 0])
    c2 = atomistic.add_atom(symbol="C", xyz=[1.5, 0, 0])
    atomistic.add_bond(c1, c2)
    
    monomer = Monomer(atomistic)
    monomer.set_port("p1", c1)
    monomer.set_port("p2", c2)
    
    def my_selector(left, right, left_ports, right_ports, ctx):
        # Always connect p2 to p1
        return ("p2", "p1")
    
    library = {"A": monomer}
    
    poly = PolymerBuilder.linear(
        sequence="AAA",
        library=library,
        connector=CallbackConnector(my_selector),
    )
    
    assert poly is not None


# ========== Test 10: Details Audit Trail ==========

def test_details_audit_trail(simple_monomer_with_roles):
    """Test that details dict is filled with connection audit."""
    library = {"A": simple_monomer_with_roles}
    details = {}
    
    poly = PolymerBuilder.linear(
        sequence="AAA",
        library=library,
        connector=AutoConnector(),
        details=details,
    )
    
    # Should have recorded 2 connections (AAA → 2 steps)
    assert "connections" in details
    assert len(details["connections"]) == 2
    
    # Each connection should have required fields
    for conn in details["connections"]:
        assert "step" in conn
        assert "left_label" in conn
        assert "right_label" in conn
        assert "left_port" in conn
        assert "right_port" in conn
        assert "bond_kind" in conn
    
    # Should have leftover ports info
    assert "leftover_ports" in details


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
