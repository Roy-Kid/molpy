import freud

from ..base import Compute


class RDF(Compute):

    def __init__(self, bins, r_max, r_min: int = 0, normalization_mode="exact"):
        super().__init__()
        self._kernel = freud.density.RDF(
            bins=bins,
            r_max=r_max,
            r_min=r_min,
            normalization_mode=normalization_mode,
        )

    def compute(self, points, query_points, box, nblist=None, reset: bool = True):

        self._kernel.compute(
            (box, points), query_points=query_points, neighbors=nblist, reset=reset
        )

        return self._kernel.rdf
