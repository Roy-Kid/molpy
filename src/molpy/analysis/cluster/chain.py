import numpy as np

import freud

from molpy.core.frame import Frame
from ..base import Compute, ComputeContext

class ChainFinder(Compute):
    """
    Find connected components (chains) in a molecular structure.
    
    This class analyzes the topology of a molecular frame to identify
    connected components, which represent separate molecular chains
    or disconnected fragments.
    """

    def compute(self, context: ComputeContext) -> ComputeContext:
        """
        Compute connected components from the frame's topology.
        
        Args:
            frame: Molecular frame containing atoms and bonds
            
        Returns:
            ChainResult containing the connected components
            
        Raises:
            ValueError: If frame has no topology information
        """
        frame = context.frame
        topo = frame.get_topology()
        components = topo.connected_components()
        membership = components.membership
        context.result[f"{self.name}_chain_idx"] = membership
        return context


class GyrationTensor(Compute):
    """
    Compute the gyration tensor of a molecular structure.
    """

    def __init__(self, name: str, cluster_field: str = "chain_idx"):
        super().__init__(name)
        self.cluster_field = cluster_field

    def compute(self, context: ComputeContext) -> ComputeContext:
        """
        Compute the gyration tensor of the frame.
        """
        xyz = context.frame["atoms"]["xyz"]
        box = context.frame.box.matrix
        cl = freud.cluster.ClusterProperties()
        cl.compute((box, xyz), context.result[self.cluster_field])
        context.result[f"{self.name}_gyrations"] = cl.gyrations
        return context

class RadiiOfGyration(Compute):
    """
    Compute the radii of gyration of a molecular structure.
    """

    def __init__(self, name: str, cluster_field: str = "chain_idx"):
        super().__init__(name)
        self.cluster_field = cluster_field
        
    def compute(self, context: ComputeContext) -> ComputeContext:
        """
        Compute the radii of gyration of the frame.
        """
        xyz = context.frame["atoms"][["xu", "yu", "zu"]]
        box = context.frame.box.matrix
        cl = freud.cluster.ClusterProperties()
        print(xyz.shape)
        print(context.result[self.cluster_field].shape)
        cl.compute((box, xyz), context.result[self.cluster_field])
        context.result[f"{self.name}_radii_of_gyration"] = cl.radii_of_gyration
        return context