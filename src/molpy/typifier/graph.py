"""Module for SMARTSGraph and SMARTS matching logic."""

import itertools
import re
from collections import OrderedDict, defaultdict
from typing import Any, Callable

from igraph import Graph, plot

from molpy.core.element import Element
from molpy.parser.smarts import SmartsParser, SmartsIR, AtomExpressionIR, AtomPrimitiveIR


class SMARTSGraph(Graph):
    """A graph representation of a SMARTS pattern.

    Constructed from SMARTS strings and parsed into an intermediate representation (IR)
    for pattern matching against molecular graphs.

    Attributes
    ----------
    atomtype_name : str
        The atom type this pattern assigns
    priority : int
        Priority for conflict resolution (higher wins)
    target_vertices : list[int]
        Which pattern vertices should receive the atom type (default: [0] for root only)
    source : str
        Source identifier for debugging
    smarts_string : str
        The SMARTS string
    ir : SmartsIR
        The intermediate representation from SMARTS parser
    dependencies : set[str]
        Set of atom type names this pattern depends on (from %opls_XXX references)
    level : int | None
        Topological level for dependency-aware typing (0=no deps, higher=has deps)

    Notes
    -----
    SMARTSGraph inherits from igraph.Graph
    
    Vertex attributes:
        - atom: AtomExpressionIR - SMARTS IR node for matching
    
    Edge attributes:
        - bond_type: BondTypeIR - SMARTS bond type specification
    """

    def __init__(
        self,
        smarts_string: str,
        parser: SmartsParser | None = None,
        name: str | None = None,
        atomtype_name: str | None = None,
        priority: int = 0,
        target_vertices: list[int] | None = None,
        source: str = "",
        overrides: set | None = None,
        *args,
        **kwargs,
    ):
        super().__init__(*args, **kwargs)

        # Metadata
        self.atomtype_name = atomtype_name or name or ""
        self.priority = priority
        self.target_vertices = target_vertices if target_vertices is not None else [0]
        self.source = source
        self.overrides = overrides
        
        # Dependency tracking
        self.dependencies: set[str] = set()
        self.level: int | None = None
        
        # Parse SMARTS string
        self.smarts_string = smarts_string
        if parser is None:
            self.ir = SmartsParser().parse_smarts(smarts_string)
        else:
            self.ir = parser.parse_smarts(smarts_string)

        self._atom_indices = OrderedDict()
        self._add_nodes()
        self._add_edges()
        
        # Extract dependencies from SMARTS IR
        self.dependencies = self.extract_dependencies()


    def __repr__(self):
        return f"<SmartsGraph({self.smarts_string})>"
    
    def plot(self, *args, **kwargs):
        """Plot the SMARTS graph."""
        graph = Graph(edges=self.get_edgelist())
        graph.vs["label"] = [v.index for v in self.vs]
        return plot(graph, *args, **kwargs)
    
    def extract_dependencies(self) -> set[str]:
        """Extract type references from SMARTS IR.
        
        Finds all has_label primitives that reference atom types (e.g., %opls_154).
        These are parsed by Lark as AtomPrimitiveIR(type="has_label", value="%opls_154").
        
        Returns:
            Set of referenced atom type names (e.g., {'opls_154', 'opls_135'})
        """
        if not self.ir or not self.ir.atoms:
            return set()
        
        dependencies = set()
        
        def extract_from_expr(expr):
            """Recursively extract dependencies from expression."""
            if isinstance(expr, AtomPrimitiveIR):
                if expr.type == "has_label" and isinstance(expr.value, str):
                    # has_label value is like "%opls_154"
                    label = expr.value
                    if label.startswith('%opls_'):
                        # Strip the % to get "opls_154"
                        dependencies.add(label[1:])
            elif isinstance(expr, AtomExpressionIR):
                for child in expr.children:
                    extract_from_expr(child)
        
        # Extract from all atoms
        for atom in self.ir.atoms:
            extract_from_expr(atom.expression)
        
        return dependencies

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
        """Determine if pattern node matches host graph node.
        
        Evaluates the SMARTS IR expression for the pattern node against
        the host node's attributes.
        
        Args:
            g1: Host graph
            g2: Pattern graph (self)
            v1: Host vertex index
            v2: Pattern vertex index
            
        Returns:
            True if pattern node matches host node
        """
        host = g1.vs[v1]
        pattern = g2.vs[v2]
        
        # Evaluate SMARTS IR
        if "atom" in pattern.attributes():
            atom_ir = pattern["atom"]
            neighbors = g1.neighbors(v1)
            return self._atom_expr_matches(atom_ir.expression, host, neighbors, g1)
        
        # No constraints - match anything
        return True
    
    def _edge_match_fn(self, g1, g2, e1, e2):
        """Determine if pattern edge matches host graph edge.
        
        Currently uses simple bond type matching. Full bond type logic
        can be implemented as needed.
        
        Args:
            g1: Host graph
            g2: Pattern graph (self)
            e1: Host edge index
            e2: Pattern edge index
            
        Returns:
            True if pattern edge matches host edge
        """
        # Simple bond type matching - always match for now
        return True

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
            # Type reference (e.g., %opls_154)
            label = str(atom_primitive.value)
            if label.startswith('%'):
                # Strip % to get type name (e.g., "opls_154")
                required_type = label[1:]
                # Check atomtype attribute (set by LayeredTypingEngine)
                assigned_type = atom.attributes().get("atomtype")
                if assigned_type is None:
                    # Type not yet assigned - no match
                    return False
                return assigned_type == required_type
            else:
                # Check if label is in type attribute (list)
                atom_vertex = graph.vs[atom_idx]
                type_attr = atom_vertex["type"] if "type" in atom_vertex.attributes() else []
                return label in type_attr
        elif atom_primitive.type == "neighbor_count":
            assert isinstance(atom_primitive.value, int)
            return len(bond_partners) == atom_primitive.value
        elif atom_primitive.type == "ring_size":
            assert isinstance(atom_primitive.value, int)
            cycle_len = atom_primitive.value
            atom_vertex = graph.vs[atom_idx]
            cycles = atom_vertex["cycles"] if "cycles" in atom_vertex.attributes() else []
            for cycle in cycles:
                if len(cycle) == cycle_len:
                    return True
            return False
        elif atom_primitive.type == "ring_count":
            assert isinstance(atom_primitive.value, int)
            atom_vertex = graph.vs[atom_idx]
            cycles = atom_vertex["cycles"] if "cycles" in atom_vertex.attributes() else []
            n_cycles = len(cycles)
            return n_cycles == atom_primitive.value
        elif atom_primitive.type == "matches_smarts":
            raise NotImplementedError("Recursive SMARTS (matches_smarts) is not yet implemented")
        else:
            raise ValueError(f"Unknown atom primitive type: {atom_primitive.type}")
