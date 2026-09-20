"""``molpy.compute`` package facade — the Compute contract is molrs-owned."""

from __future__ import annotations

from molrs.compute import Compute as MolrsCompute

from molpy.compute import Compute


def test_compute_protocol_is_the_molrs_one():
    assert Compute is MolrsCompute
