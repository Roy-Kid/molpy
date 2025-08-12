from abc import ABC, abstractmethod
from typing import Any, Literal

import freud
import numpy as np

import molpy as mp


class _NeighborQuery(ABC):

    @abstractmethod
    def query(self, points: np.ndarray, **query_kwargs) -> Any:
        pass


class AABBQuery(_NeighborQuery):

    def __init__(self, box: mp.Box, points: np.ndarray):
        self.box = box
        self.points = points
        self.kernel = freud.locality.AABBQuery(box, points)

    def query(self, points: np.ndarray, **query_kwargs) -> Any:
        return self.kernel.query(points, query_kwargs)


class LinkedCell(_NeighborQuery):

    def __init__(self, box: mp.Box, points: np.ndarray):
        self.box = box
        self.points = points
        self.kernel = freud.locality.LinkCell(box, points)

    def query(self, points: np.ndarray, **query_kwargs):
        return self.kernel.query(points, query_kwargs)


class Neighborlist:

    def __init__(
        self,
        method: Literal["aabb", "linkedcell"] = "aabb",
        points: np.ndarray | None = None,
        box: mp.Box | None = None,
    ):
        self._method = method
        self._query = None
        if points is not None and box is not None:
            self._query = self.build(points, box)

    def build(self, points: np.ndarray, box: mp.Box) -> _NeighborQuery:
        """
        Build a neighbor query for the given box and points.

        Parameters
        ----------
        box : mp.Box
            The simulation box.
        points : np.ndarray
            The points for which neighbors are to be queried.

        Returns
        -------
        _NeighborQuery
            An instance of _NeighborQuery containing the neighbor information.
        """
        if self._method == "aabb":
            kernel = AABBQuery(box, points)
        elif self._method == "linked_cell":
            kernel = LinkedCell(box, points)
        else:
            raise ValueError(f"Unknown method: {self._method}")
        return kernel

    def query_self(self, points: np.ndarray, box: mp.Box, **query_kwargs):
        """
        Compute the self-neighbors for the given box and points.

        Parameters
        ----------
        box : mp.Box
            The simulation box.
        points : np.ndarray
            The points for which self-neighbors are to be computed.
        """
        nblist = self.build(points, box)
        for pair in nblist.query(points, **query_kwargs):
            yield pair

    def query_neighbors(
        self, query_points: np.ndarray, nblist: _NeighborQuery, **query_kwargs
    ):
        """
        Compute the neighbors for the given box and points.

        Parameters
        ----------
        box : mp.Box
            The simulation box.
        points : np.ndarray
            The points for which neighbors are to be computed.
        """
        if self._query is None:
            raise ValueError("Neighbor list has not been built. Call build() first.")

        for pair in self._query.query(query_points, **query_kwargs):
            yield pair

    def add_points(self, points: np.ndarray):
        """
        Add points to the neighbor list.

        Parameters
        ----------
        points : np.ndarray
            The points to be added.
        """
        raise NotImplementedError(
            "Adding points is not supported for this neighbor list."
        )
