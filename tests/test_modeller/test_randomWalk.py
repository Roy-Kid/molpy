import pytest
import numpy as np
from molpy.modeller.randomWalk import SimpleRW

class TestSimpleRandomWalk:
    
    def test_simpleRW(self):
        
        rw = SimpleRW()
        rw.walk(3, 2)
        assert len(rw.positions) == 6
        assert len(rw.links) == 4
        
        positions = rw.positions
        del rw
        assert len(positions) == 6
        
    
    def test_walkOnce(self):
        
        rw = SimpleRW()
        rw.walkOnce(start=[1,2,3], lchain=3)
        rw.walkOnce(np.array([1,2,3]), 3)
    
    