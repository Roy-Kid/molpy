# author: Roy Kid
# contact: lijichen365@126.com
# date: 2022-01-07
# version: 0.0.1

import numpy as np
     
class Angle:
    
    def __init__(self, itom, jtom, ktom):
        self.itom = itom
        self.jtom = jtom
        self.ktom = ktom
        
    def getAngle(self, ):
        
        va = self.jtom.position - self.itom.position
        vb = self.jtom.position - self.ktom.position
        
        return np.arccos(np.dot(va, vb) / (np.linalg.norm(va) * np.linalg.norm(vb)))
    
    def setAngleType(self, angleTypes):
        pass

class Angles:
    
    def __init__(self, angles, atoms=None):
        
        self.angleIdx = Angles.unique(angles)
        self.atoms = atoms
        
    @staticmethod
    def unique(angles):
        angles = np.array(angles)
        angles = np.where((angles[:,0]>angles[:,2]).reshape((-1, 1)), angles[:, ::-1], angles)
        return np.unique(angles, axis=0)
    
    @property
    def nangles(self):
        return len(self.angleIdx)
    
    def __len__(self):
        return len(self.angleIdx)
    
    def getAngleInstances(self):
        
        if self.atoms is None:
            raise ValueError(f'need atoms to generate Angle instances')
    
        itoms = self.atoms[self.angleIdx[:, 0]]
        jtoms = self.atoms[self.angleIdx[:, 1]]
        ktoms = self.atoms[self.angleIdx[:, 2]]
        
        angles = [Angle(itom, jtom, ktom) for itom, jtom, ktom in zip(itoms, jtoms, ktoms)]      
        return angles  