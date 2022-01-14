# author: Roy Kid
# contact: lijichen365@126.com
# date: 2022-01-13
# version: 0.0.2

from typing import Literal

import numpy as np

class TypeBase:
    
    def __init__(self, id, name, attr):
        self.id = id
        self.name = name
        self.update(attr)
        
    def update(self, attr:dict):
        for k, v in attr.items():
            setattr(self, k, v)
            
    def __repr__(self):
        return f'< {self.__class__.__name__} {self.id}: {self.name} >'

class AtomType(TypeBase):
    
    def __init__(self, id, name, atomClass, attr:dict):
        super().__init__(id, name, attr)
        self.atomClass = atomClass

        
class BondType(TypeBase):
    
    def __init__(self, id, name, bondClass, attr:dict):
        super().__init__(id, name, attr)
        self.bondClass = bondClass


# below

class TypeManagement:
    
    pass

class AtomTypeManagement(TypeManagement):
    
    def __init__(self):
        self._id2AtomType = {}  # id -> atomType
        self._name2AtomType = {}  # name -> atomType
        self._class2AtomTypes = {}  # class -> [atomType]
        self._nAtomType = 0
        
    def defAtomType(self, name, atomClass=None, attr=None):
        
        atomTypeId = self._nAtomType + 1  # start from 1
        at = AtomType(atomTypeId, name, atomClass, attr)
        
        if name in self._name2AtomType:
            raise KeyError(f'AtomTypeName {name} is already defined!')
        self._name2AtomType[name] = at
        
        if atomTypeId in self._id2AtomType:
            raise KeyError(f'AtomTypeId {atomTypeId} is already defined! This error is a bug of the atomTypeId management, please contact to developers.')
        self._id2AtomType[atomTypeId] = at
        
        if atomClass in self._id2AtomType:
            self._class2AtomTypes[atomClass].add(at) 
        else:
            self._class2AtomTypes[atomClass]= set()
            
        self._nAtomType += 1
        return at
    
    def getAtomTypeByName(self, name):
        return self._name2AtomType[name]
    
    def getAtomTypeById(self, id):
        return self._id2AtomType[id]
    
    def getAtomTypeByClass(self, class_):
        return list(self._class2AtomTypes[class_])
    
    def matchAtomType(self, ):
        pass
    
class BondTypeManagement(TypeManagement):
    
    def __init__(self):
        self._id2BondType = {}
        self._name2BondType = {}
        self._class2BondType = {}
        
        self._atomTypes2BondType = {}  # {atomType1: {atomType2: bondType}}
        
        self._nBondType = 0
        
    def defBondType(self, name, bondClass=None, attr=None):
        
        bondTypeId = self._nBondType + 1
        bt = BondType(bondTypeId, name, bondClass, attr)
        
        if name in self._name2BondType:
            raise KeyError(f'BondTypeName {name} is already defined!')
        
        if bondTypeId in self._id2BondType:
            raise KeyError(f'BondTypeId {bondTypeId} is already defined! This error is a bug of the atomTypeId management, please contact to developers.')
        
        if bondClass in self._class2BondType:
            self._class2BondType[bondClass] = set()
        else:
            self._class2BondType[bondClass].add(bondClass)
        
        self._nBondType += 1
        return bt
    
    def getBondTypeByName(self, name):
        return self._name2BondType[name]
    
    def getBondTypeById(self, id):
        return self._id2BondType[id]
    
    def getBondTypeByClass(self, class_):
        return list(self._class2BondTypes[class_])
    
    def matchBondType(self, ):
        pass
        
class AngleTypeManagement(TypeManagement):
    pass

class DihedralTypeManagement(TypeManagement):
    pass

class ImproperTypeManagement(TypeManagement):
    pass
    
class ForceField(AtomTypeManagement, BondTypeManagement):
    
    def __init__(self, unit='SI'):
        
        AtomTypeManagement.__init__(self)
        BondTypeManagement.__init__(self)
        
        self._unit = unit
    
    def matchAtomTypeOfAtoms(self, atoms, field='type', ref:Literal['id', 'name', 'class']='name'):
        
        types = atoms.fields[field]
        # ats = map(self._atomTypeManagement.getAtomTypeByName, types)
        ats = np.vectorize(getattr(self, f'getAtomTypeBy{ref.capitalize()}'))(types)
        atoms.appendFields({'atomType': ats})
        return atoms
    
    def matchBondTypeOfAtoms(self, atoms, ):
        bonds = atoms.getBondIdx()
        return bonds