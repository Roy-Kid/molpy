"""PMFTXY — reads its query-point axes from the frame's ``orientations`` block.

Regression guard: the op reads its per-particle orientation axis from the
frame's core ``orientations`` topology block — one ``(head, tail)`` atom pair
per row, the same on-disk schema as ``bonds`` (endpoint columns ``atomi`` /
``atomj``). The molpy wrapper therefore forwards ``(frames, nlists)`` ONLY; no
separate angle / orientation array is passed. A prior signature passed such an
external array — these tests pin the no-external-array contract (the axis is
the internal expansion ``normalize(pos[head] - pos[tail])``).
"""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
import pytest

import molrs

from molpy.compute import PMFTXY, NeighborList


def _pmft_frame(
    random_frame: Callable[..., molrs.Frame],
    attach_orientations: Callable[..., molrs.Frame],
    n: int = 20,
    box_len: float = 12.0,
    seed: int = 1,
) -> molrs.Frame:
    frame = random_frame(n=n, box_len=box_len, seed=seed)
    idx = np.arange(n, dtype=np.uint32)
    # One (head, tail) row per particle (query-point index order).
    attach_orientations(frame, heads=(idx + 1) % n, tails=idx)
    return frame


def test_pmftxy_reads_orientations_from_frame(
    random_periodic_frame, orientations_frame
):
    frame = _pmft_frame(random_periodic_frame, orientations_frame)
    nlist = NeighborList(cutoff=3.0).compute(frame)
    out = PMFTXY(x_max=5.0, y_max=5.0, n_x=20, n_y=20).compute(frame, nlist)
    assert isinstance(out, list) and len(out) == 1
    counts, _density, _pmf = out[0]
    assert np.asarray(counts).shape == (20, 20)


def test_pmftxy_lab_frame_without_block(random_periodic_frame):
    # No orientations block => lab frame (the old `orientations=None` path).
    frame = random_periodic_frame(n=20, box_len=12.0, seed=2)
    nlist = NeighborList(cutoff=3.0).compute(frame)
    out = PMFTXY(x_max=5.0, y_max=5.0, n_x=8, n_y=8).compute(frame, nlist)
    assert isinstance(out, list) and len(out) == 1


def test_pmftxy_rejects_external_orientations(
    random_periodic_frame, orientations_frame
):
    frame = _pmft_frame(random_periodic_frame, orientations_frame)
    nlist = NeighborList(cutoff=3.0).compute(frame)
    with pytest.raises(TypeError):
        PMFTXY(5.0, 5.0, 20, 20).compute(frame, nlist, [[0.0] * 20])
