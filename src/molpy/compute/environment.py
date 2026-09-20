"""Local-environment observables — native-backed.

``BondOrder`` histograms neighbor bond directions on a (theta, phi) grid.
Thin shell over ``BondOrder``; takes
``(frames, nlists)`` like ``RDF``.

References
----------
- V. Ramasubramani et al., *Comput. Phys. Commun.* **254**, 107275 (2020) — the
  freud library, whose ``environment.BondOrder`` this mirrors.
"""

from __future__ import annotations

import molrs


class BondOrder:
    """Bond-orientational order diagram on a spherical (theta, phi) grid.

    Parameters
    ----------
    n_theta : int
        Number of polar-angle bins.
    n_phi : int
        Number of azimuthal-angle bins.
    """

    def __init__(self, n_theta: int, n_phi: int):
        self._inner = molrs.compute.environment.BondOrder(n_theta, n_phi)

    def compute(self, frames, nlists):
        return self._inner.compute(frames, nlists)
