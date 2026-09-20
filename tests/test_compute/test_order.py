"""Nematic order parameter — reads its axes from the frame's ``orientations`` block.

Regression guard: the op reads its per-particle orientation axis from the
frame's core ``orientations`` topology block — one ``(head, tail)`` atom pair
per row, the same on-disk schema as ``bonds`` (endpoint columns ``atomi`` /
``atomj``). The molpy wrapper therefore forwards ``(frames)`` ONLY; no separate
director array is passed. A prior signature passed such an external array —
these tests pin the no-external-array contract (the director/axis is the
internal expansion ``normalize(pos[head] - pos[tail])``).
"""

from __future__ import annotations

import numpy as np
import pytest

from molpy.compute import Nematic


def test_frame_carries_orientations_block(axis_frame):
    frame = axis_frame()
    assert "orientations" in frame.keys()


def test_nematic_reads_orientations_from_frame(axis_frame):
    # No director array — directors are the unit head-tail vectors of the
    # `orientations` block. All axes point +z, so order ~ 1, director ~ z.
    frame = axis_frame()
    order, eigenvalues, director, q_tensor = Nematic().compute(frame)
    assert np.asarray(eigenvalues).shape == (3,)
    assert np.asarray(q_tensor).shape == (3, 3)
    assert order > 0.9
    assert abs(np.asarray(director)[2]) > 0.9


def test_nematic_rejects_external_directors(axis_frame):
    frame = axis_frame()
    directors = np.zeros((8, 3))
    with pytest.raises(TypeError):
        Nematic().compute(frame, directors)
