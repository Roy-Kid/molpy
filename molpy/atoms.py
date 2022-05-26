# author: Roy Kid
# contact: lijichen365@126.com
# date: 2022-01-08
# version: 0.0.2

from molpy.bond import Bond
from molpy.angle import Angle
from molpy.dihedral import Dihedral
from .model import Graph
from .atom import Atom
import numpy as np
from typing import List
from freud import cluster


class Atoms(Graph):
    def __init__(self, name=None):
        super().__init__(name)
        
    @staticmethod
    def fromDict(nodes, edges, globals, topo, name=None):
        atoms = Atoms(name)
        atoms.nodes = nodes
        atoms.edges = edges
        atoms.globals = globals
        atoms.topo.setTopo(topo)
        return atoms
    
    @classmethod
    def fromAtoms(cls, atoms):
        ins = cls.fromCopy(atoms)
        return ins

    def __iter__(self):
        return iter(self.atoms)

    @property
    def atoms(self):
        return self.getAtoms()
    
    def getPositions(self, *fields):
        
        if 'positions' in self.nodes:
            return self['positions']
        
        elif 'x' in self.nodes:
            return np.array([self['x'], self['y'], self['z']]).T
        
        elif fields:
            return np.array([self[field] for field in fields]).T

    def getAtoms(self):
        atomList = []
        for i in range(self.natoms):
            atomInfo = {key: value[i] for key, value in self.nodes.items()}
            atomList.append(Atom(atomInfo))
        return atomList

    @property
    def natoms(self):
        return self.nnodes

    def getAtomByIdx(self, idx):
        fields = {key: value[idx] for key, value in self.nodes.items()}
        return Atom.fromAtom(fields)

    def setTopo(self, connection):
        self.topo.setTopo(connection)

    def getBonds(self) -> List[Bond]:
        bondIdx = self.topo.bonds
        if bondIdx is None:
            return []
        atoms = np.array(self.atoms)
        itoms = atoms[bondIdx[:, 0]]
        jtoms = atoms[bondIdx[:, 1]]
        bondTypes = getattr(self.topo, '_bondTypes', [None]*len(bondIdx))
        bonds = []
        for i in range(len(bondIdx)):
            bonds.append(Bond(itoms[i], jtoms[i], bondTypes[i]))
        return bonds

    @property
    def nbonds(self) -> int:
        return self.topo.nbonds

    def getBondIdx(self) -> List[List]:
        return self.topo.bonds

    def getAngles(self) -> List[Angle]:
        angleIdx = self.topo.angles
        angleTypes = getattr(self.topo, '_angleTypes', [None]*len(angleIdx))
        if angleIdx is None:
            return []
        atoms = np.array(self.atoms)
        itoms = atoms[angleIdx[:, 0]]
        jtoms = atoms[angleIdx[:, 1]]
        ktoms = atoms[angleIdx[:, 2]]
        
        angles = []
        for id, i in enumerate(1, range(len(angleIdx))):
            angles.append(Angle(id, itoms[i], jtoms[i], ktoms[i], angleTypes[i]))
        
        return angles

    @property
    def nangles(self) -> int:
        return self.topo.nangles

    def getAnglesIdx(self) -> List[List]:
        return self.topo.angles

    def getDihedrals(self) -> List[Dihedral]:

        dihedralIdx = self.topo.dihedrals
        dihedralTypes = getattr(self.topo, '_dihedralTypes', [None]*len(dihedralIdx))
        if dihedralIdx is None:
            return []
        atoms = np.array(self.atoms)
        itoms = atoms[dihedralIdx[:, 0]]
        jtoms = atoms[dihedralIdx[:, 1]]
        ktoms = atoms[dihedralIdx[:, 2]]
        ltoms = atoms[dihedralIdx[:, 3]]
        dihedrals = []
        for id, i in enumerate(1, range(len(dihedralIdx))):
            dihedrals.append(Dihedral(id, itoms[i], jtoms[i], ktoms[i], ltoms[i], dihedralTypes[i]))

        return dihedrals

    @property
    def ndihedrals(self) -> int:
        return self.topo.ndihedrals

    def getDihedralIdx(self) -> List[List]:
        return self.topo.dihedrals     
    
    def append(self, atoms):
        pass
    
    def splitToCluster(self, box, settings:dict):
        
        if box is None:
            raise ValueError('box is required for cluster method')
        
        clusterKernel = cluster.Cluster()
        positions = self.getPositions()
        clusterKernel.compute((box, positions), neighbors=settings)
        
        atoms_list = self.split(clusterKernel.cluster_keys)
        
        
        return atoms_list
    
    def split(self, keys:List[List]):
        
        atoms_list = []
        for key in keys:
            
            # nodes
            nodes = {field: value[key] for field, value in self.nodes.items()}
            
            # edges is not supported yet
            atoms_list.append(Atoms.fromDict(nodes, None, None, None, f'atoms with keys {keys}'))
            
        return atoms_list
        
        
            
class Cluster(Atoms):
    
    def __init__(self, name=None):
        super().__init__(name)
        self.clusterPropertiesKernel = cluster.ClusterProperties()
        self.clusterPropertiesKernel.compute()
    
    def splitToCluster(self, box, settings: dict):
        atoms_list = super().splitToCluster(box, settings)
        return list(map(Cluster.fromAtoms, atoms_list))
    
        

class AtomVec:
    
    def __init__(self):
        
        self.atoms = Atoms()
        
    def append(self):
        
        self.atoms.append()