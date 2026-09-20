"""SpatialDistribution — reads its target axes from the frame's ``orientations`` block.

Regression guard: the op reads its per-particle orientation axis from the
frame's core ``orientations`` topology block — one ``(head, tail)`` atom pair
per row, the same on-disk schema as ``bonds`` (endpoint columns ``atomi`` /
``atomj``). The molpy wrapper therefore forwards ``(frames)`` ONLY; no separate
orientation-pair array is passed. A prior signature passed such an external
array — these tests pin the no-external-array contract (the axis is the
internal expansion ``normalize(pos[head] - pos[tail])``).
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pytest

from molpy.compute import SpatialDistribution


_TEMPLATE = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])


def _sdf_kwargs(target: Sequence[int]) -> dict[str, object]:
    return dict(
        reference=[0, 1, 2],
        template=_TEMPLATE,
        target=target,
        n=(8, 8, 8),
        extent=(6.0, 6.0, 6.0),
    )


def test_sdf_reads_orientations_from_frame(random_periodic_frame, orientations_frame):
    n = 12
    frame = random_periodic_frame(n=n, box_len=10.0, seed=3)
    target = list(range(3, n))
    # One (head, tail) row per target atom, in target order.
    orientations_frame(frame, heads=target, tails=[(t + 1) % n for t in target])
    res = SpatialDistribution(**_sdf_kwargs(target)).compute([frame])
    assert np.asarray(res.density).size > 0
    assert res.orientation is not None  # per-voxel mean-orientation field present


def test_sdf_without_orientations_block_is_isotropic(random_periodic_frame):
    n = 12
    frame = random_periodic_frame(n=n, box_len=10.0, seed=3)
    res = SpatialDistribution(**_sdf_kwargs(list(range(3, n)))).compute([frame])
    assert res.orientation is None


def test_sdf_rejects_orientation_pairs_kwarg():
    n = 12
    with pytest.raises(TypeError):
        SpatialDistribution(
            orientation_pairs=np.zeros((n - 3, 2), dtype=np.int64),
            **_sdf_kwargs(list(range(3, n))),
        )
