import freud

import molpy as mp

from ..base import Compute


class ChainFinder(Compute):
    """
    Find connected components (chains) in a molecular structure.

    This class analyzes the topology of a molecular frame to identify
    connected components, which represent separate molecular chains
    or disconnected fragments.
    """

    def compute(self, topo: mp.Topology):
        """
        Compute connected components from the frame's topology.

        Args:
            frame: Molecular frame containing atoms and bonds

        Returns:
            ChainResult containing the connected components

        Raises:
            ValueError: If frame has no topology information
        """
        components = topo.connected_components()
        membership = components.membership
        return membership


class GyrationTensor(Compute):
    """
    Compute the gyration tensor of a molecular structure.
    """

    def compute(self, xyz, box, cluster_id):
        """
        Compute the gyration tensor of the frame.
        """
        cl = freud.cluster.ClusterProperties()
        cl.compute((box.matrix, xyz), cluster_id)
        return cl.gyrations


class RadiiOfGyration(Compute):
    """
    Compute the radii of gyration of a molecular structure.
    """

    def compute(self, xyz, box, cluster_id):
        """
        Compute the radii of gyration of the frame.
        """
        cl = freud.cluster.ClusterProperties()
        cl.compute((box, xyz), cluster_id)
        return cl.radii_of_gyration
