"""Pattern builder for constructing SmartsGraph programmatically.

This module provides a fluent API for building SMARTS patterns without
requiring a SMARTS string parser.
"""

from typing import Any
from igraph import Graph
from .predicates import VertexPredicate, EdgePredicate, is_element


class PatternBuilder:
    """Builder for creating SmartsGraph patterns programmatically.
    
    Example:
        >>> pb = PatternBuilder("t_CO_single", priority=0)
        >>> c = pb.add_vertex(preds=[is_element("C")])
        >>> o = pb.add_vertex(preds=[is_element("O")])
        >>> pb.add_edge(c, o, preds=[bond_order(1)])
        >>> graph = pb.build()
    """
    
    def __init__(self, atomtype_name: str, priority: int = 0, 
                 source: str = ""):
        """Initialize pattern builder.
        
        Args:
            atomtype_name: Name of the atom type this pattern assigns
            priority: Priority for conflict resolution (higher wins)
            source: Source identifier (e.g., "File:Line" or "test:pattern_1")
        """
        self.atomtype_name = atomtype_name
        self.priority = priority
        self.source = source or f"pattern:{atomtype_name}"
        
        # Internal graph for building
        self._vertices: list[int] = []
        self._vertex_preds: dict[int, list[VertexPredicate]] = {}
        self._edges: list[tuple[int, int]] = []
        self._edge_preds: dict[tuple[int, int], list[EdgePredicate]] = {}
        
        # Target vertices (which vertices should be typed)
        self._target_vertices: list[int] = []
        
        # Vertex counter
        self._next_vertex_id = 0
    
    def add_vertex(self, preds: list[VertexPredicate] | None = None) -> int:
        """Add a vertex with predicates.
        
        Args:
            preds: List of vertex predicates (default: wildcard)
        
        Returns:
            Vertex ID (integer index)
        """
        if preds is None:
            preds = [is_element("*")]
        
        vertex_id = self._next_vertex_id
        self._next_vertex_id += 1
        
        self._vertices.append(vertex_id)
        self._vertex_preds[vertex_id] = preds
        
        return vertex_id
    
    def add_edge(self, u: int, v: int, 
                 preds: list[EdgePredicate] | None = None) -> None:
        """Add an edge between two vertices.
        
        Args:
            u: Source vertex ID
            v: Target vertex ID
            preds: List of edge predicates (default: any bond)
        """
        if preds is None:
            preds = []
        
        if u not in self._vertices:
            raise ValueError(f"Vertex {u} not in pattern")
        if v not in self._vertices:
            raise ValueError(f"Vertex {v} not in pattern")
        
        # Store as undirected edge (canonicalized)
        edge = (min(u, v), max(u, v))
        self._edges.append(edge)
        self._edge_preds[edge] = preds
    
    def set_target_vertices(self, vertices: list[int]) -> None:
        """Set which vertices should receive the atom type.
        
        Args:
            vertices: List of vertex IDs that should be typed
        
        If not called, all matched vertices receive the type.
        """
        for v in vertices:
            if v not in self._vertices:
                raise ValueError(f"Vertex {v} not in pattern")
        self._target_vertices = vertices
    
    def build(self) -> "SMARTSGraph":
        """Build the SMARTSGraph.
        
        Returns:
            SMARTSGraph instance with predicates and metadata
        """
        from .graph import SMARTSGraph
        
        # Create igraph Graph
        g = Graph(n=len(self._vertices), directed=False)
        
        # Add edges
        edge_list = [(self._vertices.index(u), self._vertices.index(v))
                     for u, v in self._edges]
        if edge_list:
            g.add_edges(edge_list)
        
        # Convert vertex IDs to indices
        vertex_id_to_idx = {vid: idx for idx, vid in enumerate(self._vertices)}
        
        # Attach predicates to vertices
        for idx, vertex_id in enumerate(self._vertices):
            g.vs[idx]["preds"] = self._vertex_preds.get(vertex_id, [])
        
        # Attach predicates to edges
        for eid in range(g.ecount()):
            edge = g.es[eid]
            source_id = self._vertices[edge.source]
            target_id = self._vertices[edge.target]
            edge_key = (min(source_id, target_id), max(source_id, target_id))
            edge["preds"] = self._edge_preds.get(edge_key, [])
        
        # Convert target vertices to indices
        target_indices = [vertex_id_to_idx[v] for v in self._target_vertices]
        
        # Create SMARTSGraph wrapper
        smarts_graph = SMARTSGraph.from_igraph(
            g,
            atomtype_name=self.atomtype_name,
            priority=self.priority,
            target_vertices=target_indices,
            source=self.source
        )
        
        return smarts_graph
    
    def __repr__(self) -> str:
        return (f"PatternBuilder(name={self.atomtype_name!r}, "
                f"vertices={len(self._vertices)}, edges={len(self._edges)})")


def quick_pattern(atomtype_name: str, element: str, 
                   priority: int = 0, **predicates) -> "SMARTSGraph":
    """Quick helper to create single-atom patterns.
    
    Args:
        atomtype_name: Atom type name
        element: Element symbol
        priority: Priority
        **predicates: Additional predicates as kwargs
            - is_aromatic: bool
            - charge: int
            - degree: int
            - hyb: int
            - in_ring: bool
    
    Example:
        >>> pattern = quick_pattern("t_O_neg", "O", charge=-1)
    """
    from .predicates import (is_element, is_aromatic, charge, degree, 
                              hyb, in_ring)
    
    pb = PatternBuilder(atomtype_name, priority=priority)
    
    preds = [is_element(element)]
    
    if predicates.get("is_aromatic") is not None:
        preds.append(is_aromatic(predicates["is_aromatic"]))
    if predicates.get("charge") is not None:
        preds.append(charge(predicates["charge"]))
    if predicates.get("degree") is not None:
        preds.append(degree(predicates["degree"]))
    if predicates.get("hyb") is not None:
        preds.append(hyb(predicates["hyb"]))
    if predicates.get("in_ring") is not None:
        preds.append(in_ring(predicates["in_ring"]))
    
    pb.add_vertex(preds=preds)
    
    return pb.build()
