"""Dimensionality reduction + clustering ML primitives — native-backed.

- ``Pca`` (alias for the native core ``Pca2``): 2-component PCA over a list of
  ``DescriptorRow`` objects.
- ``KMeans``: k-means clustering over a ``PcaResult``.

``DescriptorRow`` is re-exported as the input wrapper: each row is a
1-D float ndarray passed through ``DescriptorRow(row)``.
"""

from __future__ import annotations

from molrs.compute.ml import (
    DescriptorRow as DescriptorRow,
    KMeans as _MolrsKMeans,
    Pca2 as _MolrsPca2,
)


class Pca:
    """Two-component PCA. Input: list of ``DescriptorRow``."""

    def __init__(self) -> None:
        self._impl = _MolrsPca2()

    def compute(self, rows):
        return self._impl.compute(rows)


class KMeans:
    """k-means over a ``PcaResult``.

    Parameters
    ----------
    k : int
        Number of clusters.
    max_iter : int, default 100
        Lloyd-iteration cap.
    seed : int, default 0
        RNG seed for centroid initialization.
    """

    def __init__(self, k: int, max_iter: int = 100, seed: int = 0) -> None:
        self._impl = _MolrsKMeans(k, max_iter, seed)

    def compute(self, pca_result):
        return self._impl.compute(pca_result)


__all__ = ["DescriptorRow", "Pca", "KMeans"]
