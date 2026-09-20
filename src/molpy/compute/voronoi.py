"""Radical (Laguerre) Voronoi tessellation, domains, voids & integration — native-backed.

The **radical** (power / Laguerre) Voronoi tessellation partitions space by
radius-weighted planes, so atoms of different size get cells proportional to
their radii — the physically correct partition for polydisperse systems. On top
of the tessellation:

- :func:`voronoi_domains` merges cells sharing a label into connected domains
  (e.g. polar vs. apolar nanostructuring in ionic liquids).
- :func:`voronoi_voids` aggregates the cells flagged as empty into void volumes.
- :class:`VoronoiIntegration` integrates an electron density over the cells to
  yield per-molecule charges and dipoles (Voronoi/atomic-charge partitioning),
  the basis for predicting infrared spectra from *ab initio* MD.

Every name here is the native analysis-parity object itself, re-exported
unchanged: :class:`RadicalVoronoi` tessellates through ``build(positions,
radii, box)`` and :class:`VoronoiIntegration` reduces a density through
``integrate(positions, radii, atomic_numbers, atom_to_mol, n_mol, grid, box)``.

References
----------
- B. J. Gellatly, J. L. Finney, *J. Non-Cryst. Solids* **50**, 313 (1982) — radical
  (power) Voronoi tessellation.
- M. Thomas, M. Brehm, B. Kirchner, *Phys. Chem. Chem. Phys.* **17**, 3207 (2015)
  — Voronoi integration of the electron density for molecular dipoles.
- M. Brehm, M. Thomas, S. Gehrke, B. Kirchner, *J. Chem. Phys.* **152**, 164105
  (2020) — reference implementation; domain and void analysis.
"""

from __future__ import annotations

import molrs

# Re-export the tessellation/integration operators, the per-cell result type and
# the domain/void reductions.
RadicalVoronoi = molrs.compute.voronoi.RadicalVoronoi
VoronoiIntegration = molrs.compute.voronoi.VoronoiIntegration
VoronoiCells = molrs.compute.voronoi.VoronoiCells
voronoi_domains = molrs.compute.voronoi.voronoi_domains
voronoi_voids = molrs.compute.voronoi.voronoi_voids

__all__ = [
    "RadicalVoronoi",
    "VoronoiCells",
    "VoronoiIntegration",
    "voronoi_domains",
    "voronoi_voids",
]
