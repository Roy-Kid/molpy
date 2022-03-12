# author: Roy Kid
# contact: lijichen365@126.com
# date: 2022-03-06
# version: 0.0.1

from molpy.modeller.randomWalk import RandomWalkOnFcc
from molpy.modeller.utils import toAtoms
from molpy.system import System
import pytest

class TestSystem:
    
    @pytest.fixture(scope='class', name='system')
    def test_init(self):
        
        system = System()
        
        rw = RandomWalkOnFcc(10, 10, 10)
        positions, bonds = rw.linear(10)
        atoms = toAtoms(positions, bonds)
        system.atomManager.append(atoms)
        
        return system
    
    def test_properties(self, system):
        
        assert system.natoms == 10
        assert system.nbonds == 9
        assert system.nangles == 8
        assert system.ndihedrals == 7
        
    def test_add_atoms(self, system):
        
        rw = RandomWalkOnFcc(10, 10, 10)
        positions, bonds = rw.linear(10)
        atoms = toAtoms(positions, bonds)
        system.atomManager.append(atoms)
        
        assert system.natoms == 20
        assert system.nbonds == 18
        assert system.nangles == 16
        assert system.ndihedrals == 14