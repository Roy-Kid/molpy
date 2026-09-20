"""Mean squared displacement of particle positions."""

from __future__ import annotations

from molrs.compute.msd import MSD as _MolrsMSD
from molrs.compute.msd import MSDTimeSeries as _MolrsMSDTimeSeries


class MSD:
    """Mean squared displacement.

    Two estimators, chosen by ``method`` — they are not interchangeable:

    - ``"direct"`` (default): ``<|r(t) - r(0)|^2>``, frame 0 the one time origin.
    - ``"window"``: ``<|r(tau+t) - r(tau)|^2>`` averaged over **every** time
      origin. Better statistics at long lag, which is what a diffusion
      coefficient needs, and O(T log T) rather than O(T^2).

    Examples
    --------
    >>> series = MSD(method="window").compute(trajectory_frames)
    >>> series.mean.shape    # (n_frames,)
    """

    def __init__(self, method: str = "direct") -> None:
        self.method = method
        self._impl = _MolrsMSD(method=method)

    def compute(self, frames) -> _MolrsMSDTimeSeries:
        return self._impl.compute(frames)
