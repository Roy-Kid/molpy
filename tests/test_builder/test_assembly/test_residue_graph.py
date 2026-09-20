"""Unit tests for :mod:`molpy.builder.assembly._residue_graph`."""

from molpy.builder.assembly import (
    linear_topology,
    ring_topology,
    star_topology,
)


class TestLinearTopology:
    def test_builds_path_nodes_and_edges(self):
        g = linear_topology(["A", "B", "C"])
        assert [n.label for n in g.nodes] == ["A", "B", "C"]
        assert len(g.bonds) == 2


class TestRingTopology:
    def test_closes_cycle(self):
        g = ring_topology("EO", 4)
        assert len(g.nodes) == 4
        assert len(g.bonds) == 4


class TestStarTopology:
    def test_core_and_arms(self):
        g = star_topology("X", "A", n_arms=3, arm_length=2)
        assert g.nodes[0].label == "X"
        assert len(g.nodes) == 1 + 3 * 2
