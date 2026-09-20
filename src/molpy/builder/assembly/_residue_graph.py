"""Residue topology constructors for polymer assembly.

Helpers such as :func:`linear_topology` build a :class:`CGSmilesGraphIR`
directly. There is no CGSmiles string round-trip and no Lark parser.
"""

from __future__ import annotations

from collections.abc import Sequence

from molpy.builder.assembly._cgsmiles_ir import (
    CGSmilesBondIR,
    CGSmilesGraphIR,
    CGSmilesNodeIR,
)


def linear_topology(labels: Sequence[str]) -> CGSmilesGraphIR:
    """Path topology: one node per label, edges between consecutive residues.

    Raises:
        ValueError: if ``labels`` is empty.
    """
    if not labels:
        raise ValueError("linear topology needs at least one residue label")
    labels = [str(lab) for lab in labels]
    nodes = [CGSmilesNodeIR(label=lab) for lab in labels]
    bonds = [
        CGSmilesBondIR(node_i=nodes[i], node_j=nodes[i + 1])
        for i in range(len(nodes) - 1)
    ]
    return CGSmilesGraphIR(nodes=nodes, bonds=bonds)


def ring_topology(label: str, n: int) -> CGSmilesGraphIR:
    """Cycle of ``n`` identical residues (``n >= 3``)."""
    if n < 3:
        raise ValueError(f"a residue ring needs n >= 3, got {n}")
    graph = linear_topology([str(label)] * n)
    graph.bonds.append(CGSmilesBondIR(node_i=graph.nodes[-1], node_j=graph.nodes[0]))
    return graph


def star_topology(
    core: str,
    arm: str,
    *,
    n_arms: int,
    arm_length: int,
    cap: str | None = None,
) -> CGSmilesGraphIR:
    """Star: one core node bonded to ``n_arms`` linear arms of ``arm_length``.

    Raises:
        ValueError: if ``n_arms < 2`` or ``arm_length < 1``.
    """
    if n_arms < 2:
        raise ValueError(f"a star needs n_arms >= 2, got {n_arms}")
    if arm_length < 1:
        raise ValueError(f"arm_length must be >= 1, got {arm_length}")

    core_node = CGSmilesNodeIR(label=str(core))
    nodes: list[CGSmilesNodeIR] = [core_node]
    bonds: list[CGSmilesBondIR] = []

    for _ in range(n_arms):
        prev = core_node
        for _j in range(arm_length):
            node = CGSmilesNodeIR(label=str(arm))
            nodes.append(node)
            bonds.append(CGSmilesBondIR(node_i=prev, node_j=node))
            prev = node
        if cap is not None:
            cap_node = CGSmilesNodeIR(label=str(cap))
            nodes.append(cap_node)
            bonds.append(CGSmilesBondIR(node_i=prev, node_j=cap_node))

    return CGSmilesGraphIR(nodes=nodes, bonds=bonds)


# ---------------------------------------------------------------------------
