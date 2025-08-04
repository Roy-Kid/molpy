"""
Tests for chain analysis tools.
"""

import random
import numpy as np
import pytest
import igraph as ig

import molpy as mp
from molpy.analysis.cluster.chain import ChainFinder, ChainResult


@pytest.fixture
def mock_frame():
    """Create a simple topology with two disconnected chains."""
    random.seed(0)
    n_nodes = 50
    g = ig.Graph.GRG(n_nodes, 0.15)
    bonds = np.array(g.get_edgelist())
    frame = mp.Frame()
    frame["atoms"] = mp.Block()
    frame["atoms"]["id"] = np.arange(n_nodes)
    frame["bonds"] = mp.Block()
    frame["bonds"]["i"] = bonds[:, 0]
    frame["bonds"]["j"] = bonds[:, 1]
    
    return frame

@pytest.fixture
def chain_finder():
    """Create a ChainFinder instance."""
    return ChainFinder()

class TestChainResult:
    """Test the ChainResult class."""
    
    def test_init_with_components(self):
        """Test ChainResult initialization with components."""
        chain_idx = [np.array([0, 1, 2]), np.array([3, 4])]
        result = ChainResult(chain_idx=chain_idx)
        
        assert len(result.chain_idx) == 2
        assert np.array_equal(result.chain_idx[0], np.array([0, 1, 2]))
        assert np.array_equal(result.chain_idx[1], np.array([3, 4]))
    
    def test_init_empty(self):
        """Test ChainResult initialization with empty components."""
        result = ChainResult(chain_idx=[])
        assert len(result.chain_idx) == 0


class TestChainFinder:
    """Test the ChainFinder class."""
    
    def test_init(self, chain_finder):
        """Test ChainFinder initialization."""
        assert isinstance(chain_finder, ChainFinder)
    
    def test_compute(self, chain_finder, mock_frame):
        """Test ChainFinder computation."""
        result = chain_finder.compute(mock_frame)
        assert isinstance(result, ChainResult)
        assert len(result.chain_idx) == 8