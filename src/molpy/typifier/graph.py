"""Module for SMARTSGraph and SMARTS matching logic."""

import itertools
from collections import OrderedDict, defaultdict

from igraph import Graph, plot

from molpy.core.element import Element
from molpy.parser.smarts import SmartsParser, SmartsIR, AtomExpressionIR, AtomPrimitiveIR


class SMARTSGraph(Graph):
    """A graph representation of a SMARTS pattern.

    Attributes
    ----------
    smarts_string : str
        The SMARTS string outlined in the force field
    parser : SmartsParser
        The parser that converts SMARTS string into IR
    name : str
    overrides : set
        Rules or SMARTSGraph over which this SMARTSGraph takes precedence
    ir : SmartsIR
        The intermediate representation of the SMARTS pattern

    Attributes
    ----------
    graph_matcher : smarts_graph.SMARTSMatcher
        implementation of VF2 that handles subgraph matching

    Notes
    -----
    SMARTSGraph inherits from igraph.Graph
    """

    def __init__(
        self,
        smarts_string: str,
        parser: SmartsParser | None = None,
        name=None,
        overrides=None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        self.smarts_string = smarts_string
        self.name = name
        self.overrides = overrides

        if parser is None:
            self.ir = SmartsParser().parse_smarts(smarts_string)
        else:
            self.ir = parser.parse_smarts(smarts_string)

        self._atom_indices = OrderedDict()
        self._add_nodes()
        self._add_edges()
        self._graph_matcher = None

    def __repr__(self):
        return f"<SmartsGraph({self.smarts_string})>"

    def plot(self, *args, **kwargs):
        """Plot the SMARTS graph."""
        graph = Graph(edges=self.get_edgelist())
        graph.vs["label"] = [v.index for v in self.vs]
        return plot(graph, *args, **kwargs)

    def override(self, overrides):
        """Set the priority of this SMART"""
        self.overrides = overrides
        # self.priority = max([override.priority for override in overrides]) + 1

    @property
    def priority(self):
        if self.overrides is None:
            return 0
        return max([override.priority for override in self.overrides]) + 1

    def _add_nodes(self):
        """Add all atoms in the SMARTS IR as nodes in the graph."""
        atoms = self.ir.atoms
        self.add_vertices(len(atoms), {"atom": atoms})
        for i, atom in enumerate(atoms):
            self._atom_indices[id(atom)] = i

    def _add_edges(self):
        """Add all bonds in the SMARTS IR as edges in the graph."""
        atom_indices = self._atom_indices
        for bond in self.ir.bonds:
            start_idx = atom_indices[id(bond.start)]
            end_idx = atom_indices[id(bond.end)]
            self.add_edge(start_idx, end_idx, bond_type=bond.bond_type)

    def _node_match_fn(self, g1, g2, v1, v2):
        """Determine if two graph nodes are equal."""
        host = g1.vs[v1]
        pattern = g2.vs[v2]
        atom_ir = pattern["atom"]  # This is a SmartsAtomIR
        neighbors = g1.neighbors(v1)
        result = self._atom_expr_matches(atom_ir.expression, host, neighbors, g1)
        return result

    def _atom_expr_matches(self, atom_expr: AtomExpressionIR | AtomPrimitiveIR, atom, bond_partners, graph):
        """Evaluate SMARTS IR expressions."""
        # Handle AtomPrimitiveIR directly
        if isinstance(atom_expr, AtomPrimitiveIR):
            return self._atom_primitive_matches(atom_expr, atom, bond_partners, graph)
        
        # Handle AtomExpressionIR
        if atom_expr.op == "not":
            return not self._atom_expr_matches(
                atom_expr.children[0], atom, bond_partners, graph
            )
        elif atom_expr.op in ("and", "weak_and"):
            result = True
            for child in atom_expr.children:
                result = result and self._atom_expr_matches(
                    child, atom, bond_partners, graph
                )
            return result
        elif atom_expr.op == "or":
            for child in atom_expr.children:
                if self._atom_expr_matches(child, atom, bond_partners, graph):
                    return True
            return False
        elif atom_expr.op == "primitive":
            if atom_expr.children:
                return self._atom_expr_matches(
                    atom_expr.children[0], atom, bond_partners, graph
                )
            return True
        else:
            raise TypeError(
                f"Unexpected atom expression op: {atom_expr.op}"
            )

    @staticmethod
    def _atom_primitive_matches(atom_primitive: AtomPrimitiveIR, atom, bond_partners, graph):
        """Compare atomic primitives against atom properties."""
        atomic_num = atom.attributes().get("number", None)
        atom_name = atom.attributes().get("name", None)
        atom_idx = atom.index
        assert atomic_num or atom_name, f"Atom {atom_idx} has no atomic number or name."

        if atom_primitive.type == "atomic_num":
            assert isinstance(atom_primitive.value, int)
            return atomic_num == atom_primitive.value
        elif atom_primitive.type == "symbol":
            symbol_val = str(atom_primitive.value)
            if symbol_val == "*":
                return True
            elif symbol_val.startswith("_"):
                # Store non-element elements in .name
                return atom_name == symbol_val
            else:
                # Handle lowercase (aromatic) symbols
                if symbol_val.islower():
                    # Check both element and aromaticity
                    # For now, just check element
                    return atomic_num == Element(symbol_val.upper()).number
                return atomic_num == Element(symbol_val).number
        elif atom_primitive.type == "wildcard":
            return True
        elif atom_primitive.type == "has_label":
            label = str(atom_primitive.value)[1:]  # Strip the % sign
            return label in graph.vs[atom_idx].get("type", [])
        elif atom_primitive.type == "neighbor_count":
            assert isinstance(atom_primitive.value, int)
            return len(bond_partners) == atom_primitive.value
        elif atom_primitive.type == "ring_size":
            assert isinstance(atom_primitive.value, int)
            cycle_len = atom_primitive.value
            for cycle in graph.vs[atom_idx].get("cycles", []):
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_primitive.type == "ring_count":
            assert isinstance(atom_primitive.value, int)
            n_cycles = len(graph.vs[atom_idx].get("cycles", []))
            return n_cycles == atom_primitive.value
        elif atom_primitive.type == "matches_smarts":
            raise NotImplementedError("Recursive SMARTS (matches_smarts) is not yet implemented")
        else:
            raise ValueError(f"Unknown atom primitive type: {atom_primitive.type}")

    def find_matches(self, graph):
        """Return sets of atoms that match this SMARTS pattern in a topology.

        Parameters
        ----------
        structure : TopologyGraph
            The topology that we are trying to atomtype.
        typemap : dict
            The target typemap being used/edited

        Notes
        -----
        When this function gets used in atomtyper.py, we actively modify the
        white- and blacklists of the atoms in `topology` after finding a match.
        This means that between every successive call of
        `subgraph_isomorphisms_iter()`, the topology against which we are
        matching may have actually changed. Currently, we take advantage of this
        behavior in some edges cases (e.g. see `test_hexa_coordinated` in
        `test_smarts.py`).

        """

        self.calc_signature(graph)

        self._graph_matcher = SMARTSMatcher(
            graph, self, node_match_fn=self._node_match_fn
        )

        matches = self._graph_matcher.subgraph_isomorphisms()
        match_index = set([match[0] for match in matches])
        return match_index

    def calc_signature(self, graph):
        """Calculate graph signatures for pattern matching."""
        # Check if any atoms have ring-related properties
        def check_expr_for_rings(expr):
            """Recursively check expression for ring-related primitives."""
            if isinstance(expr, AtomPrimitiveIR):
                return expr.type in ("ring_size", "ring_count")
            if isinstance(expr, AtomExpressionIR):
                return any(check_expr_for_rings(child) for child in expr.children)
            return False
        
        has_ring_rules = any(
            check_expr_for_rings(atom.expression) 
            for atom in self.ir.atoms
        )

        if has_ring_rules:
            graph.vs["cycles"] = [set() for _ in graph.vs]
            all_cycles = _find_chordless_cycles(graph, max_cycle_size=6)
            for i, cycles in enumerate(all_cycles):
                for cycle in cycles:
                    graph.vs[i]["cycles"].add(tuple(cycle))


class SMARTSMatcher:
    """Inherits and implements VF2 for a SMARTSGraph."""

    def __init__(self, G1: Graph, G2: Graph, node_match_fn):
        self.G1 = G1
        self.G2 = G2
        self.node_match_fn = node_match_fn

    @property
    def is_isomorphic(self):
        """Return True if the two graphs are isomorphic."""
        return self.G1.isomorphic(self.G2)

    def subgraph_isomorphisms(self):
        """Iterate over all subgraph isomorphisms between G1 and G2."""
        matches = self.G1.get_subisomorphisms_vf2(
            self.G2, node_compat_fn=self.node_match_fn
        )
        results = []
        for sgi in matches:
            sg = self.G1.subgraph(sgi)
            if sg.get_isomorphisms_vf2(self.G2):
                results.append(sgi)
        return results

    def candidate_pairs_iter(self):
        """Iterate over candidate pairs of nodes in G1 and G2."""
        # All computations are done using the current state!
        G2_nodes = self.G2_nodes

        # First we compute the inout-terminal sets.
        T1_inout = set(self.inout_1.keys()) - set(self.core_1.keys())
        T2_inout = set(self.inout_2.keys()) - set(self.core_2.keys())

        # If T1_inout and T2_inout are both nonempty.
        # P(s) = T1_inout x {min T2_inout}
        if T1_inout and T2_inout:
            for node in T1_inout:
                yield node, min(T2_inout)
        else:
            # First we determine the candidate node for G2
            other_node = min(G2_nodes - set(self.core_2))
            host_nodes = self.valid_nodes if other_node == 0 else self.G1.nodes()
            for node in host_nodes:
                if node not in self.core_1:
                    yield node, other_node

        # For all other cases, we don't have any candidate pairs.


def _find_chordless_cycles(graph, max_cycle_size):
    """Find all chordless cycles (i.e. rings) in the bond graph.

    Traverses the bond graph to determine all cycles (i.e. rings) each
    atom is contained within. Algorithm has been adapted from:
    https://stackoverflow.com/questions/4022662/find-all-chordless-cycles-in-an-undirected-graph/4028855#4028855
    """
    cycles = [[] for _ in graph.vs]

    """
    For all nodes we need to find the cycles that they are included within.
    """
    for i, node in enumerate(graph.vs):
        node_idx = node.index
        neighbors = list(graph.neighbors(node_idx))
        pairs = list(itertools.combinations(neighbors, 2))
        """
        Loop over all pairs of neighbors of the node. We will see if a ring
        exists that includes these branches.
        """
        for pair in pairs:
            """
            We need to store all node sequences that could be rings. We will
            update this as we traverse the graph.
            """
            connected = False
            possible_rings = []

            last_node = pair[0]
            ring = [last_node, node_idx, pair[1]]
            possible_rings.append(ring)

            if graph.are_adjacent(last_node, pair[1]):
                cycles[i].append(ring)
                connected = True

            while not connected:
                """
                Branch and create a new list of possible rings
                """
                new_possible_rings = []
                for possible_ring in possible_rings:
                    next_neighbors = graph.neighbors(possible_ring[-1])
                    for next_neighbor in next_neighbors:
                        if next_neighbor != possible_ring[-2]:
                            new_possible_rings.append(possible_ring + [next_neighbor])
                possible_rings = new_possible_rings

                for possible_ring in possible_rings:
                    if graph.are_adjacent(possible_ring[-1], last_node):
                        if any(
                            [
                                graph.are_adjacent(possible_ring[-1], internal_node)
                                for internal_node in possible_ring[1:-2]
                            ]
                        ):
                            pass
                        else:
                            cycles[i].append(possible_ring)
                            connected = True

                if not possible_rings or len(possible_rings[0]) == max_cycle_size:
                    break

    return cycles
