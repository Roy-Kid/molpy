# author: Roy Kid
# contact: lijichen365@126.com
# date: 2022-06-10
# version: 0.0.1

from typing import Dict, Optional, Union
import numpy as np
import numpy.typing as npt
from molpy.topo import Topo


class Graph:
    
    def __init__(self, withTopo:bool=False):
        
        self._n_nodes = 0
        self._topo: Optional[Topo] = Topo() if withTopo else None
        self._nodes: Dict[str, npt.NDArray] = {}  # 
        self.withTopo = withTopo
        
    def set_node(self, key:str, value:npt.ArrayLike, ref:npt.ArrayLike=None):
        
        v = np.array(value)

        if ref is not None:
            args = np.argsort(ref)
            v = v[args]

        self._nodes[key] = v
        
    def get_node(self, key:str)->npt.NDArray:
        
        return self._nodes[key]
        
    def get_subgraph(self, index: Union[slice, int, list, npt.NDArray]) -> 'Graph':
        
        graph = Graph.from_graph(self)
        nodes = {key: value[index] for key, value in self._nodes.items()}
        graph._nodes = nodes
        return graph
        
    def __getitem__(self, o: Union[str, slice, int, npt.NDArray]) -> Union['Graph', npt.NDArray]:
        
        if isinstance(o, str):
            return self.get_node(o)
        
        return self.get_subgraph(o)
    
    def __setitem__(self, key, value):
        
        self.set_node(key, value)
    
    @classmethod
    def from_graph(cls, graph: 'Graph') -> 'Graph':
        
        ins = cls()
        ins._nodes = graph._nodes
        return ins
    
    @property
    def n_nodes(self):
        
        n = 0
        for value in self._nodes.values():
            n = max(n, len(value))
        return n
    
    @property
    def n_edges(self):
        if self._topo is None:
            return 0
        return self._topo.n_edges
    
    def __repr__(self)->str:
        
        if self._topo is None:
            return f'<Graph: {self._n_nodes} nodes>'
        else:
            return f'<Graph: {self._n_nodes} nodes, {self._topo.n_edges} edges>'
        
    def append(self, graph: 'Graph'):
        
        for node in graph._nodes:
            self.set_node(node, np.concatenate((self._nodes[node], graph._nodes[node])))
        
        if graph.withTopo:
            self.topo.append(graph._topo)
        
    def replace_nodes(self, nodes: Dict[str, npt.ArrayLike]):

        for key, value in nodes.items():
            self.set_node(key, value)