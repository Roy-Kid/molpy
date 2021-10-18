# author: Roy Kid
# contact: lijichen365@126.com
# date: 2021-10-17
# version: 0.0.1

import pytest
import molpy as mp
import numpy as np

class TestGroup:
    
    @pytest.fixture(scope='class')
    def CH4(self):
        CH4 = mp.Group('CH4')
        C = mp.Atom('C')
        Hs = [mp.Atom(f'H{i}') for i in range(4)]
        CH4.add(C)
        for H in Hs:
            CH4.add(H)
        yield CH4
        
    @pytest.fixture(scope='class')
    def C6(self):
        C6 = mp.Group('C6')
        for C in [mp.Atom(f'C{i}') for i in range(6)]:
            C6.add(C)
        covalentMap = np.array([[0, 1, 2, 3, 2, 1],
                                [1, 0, 1, 2, 3, 2],
                                [2, 1, 0, 1, 2, 3],
                                [3, 2, 1, 0, 1, 2],
                                [2, 3, 2, 1, 0, 1],
                                [1, 2, 3, 2, 1, 0]], dtype=int)
        C6.setTopoByCovalentMap(covalentMap)
        C6.reference_covalentMap = covalentMap
        yield C6
            
    def test_setTopoByCovalentMap(self, CH4):
        covalentMap = np.zeros((CH4.natoms, CH4.natoms), dtype=int)
        covalentMap[0, 1:] = covalentMap[1:, 0] = 1
        CH4.setTopoByCovalentMap(covalentMap)
        assert len(CH4['C'].bondedAtoms) == 4
        assert CH4['H0'].bondedAtoms[0] == CH4['C']
        
    def test_getCovalentMap(self, CH4):
        co = CH4.getCovalentMap()
        print(co)
        assert co[0, 0] == 0
        assert all(co[0, 1:] == 1)
        assert all(co[1:, 0] == 1)
        
    def test_getRingCovalentMap(self, C6):
        assert(C6.reference_covalentMap == C6.getCovalentMap()).all()
        
    def test_setRingTopoByCovalentMap(self, C6):
        atoms = C6.getAtoms()
        assert len(atoms[0].bondedAtoms) == 2
        assert len(atoms[1].bondedAtoms) == 2
        assert len(atoms[2].bondedAtoms) == 2
        assert len(atoms[3].bondedAtoms) == 2