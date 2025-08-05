"""
Tests for chain analysis tools.
"""

import random
import numpy as np
import pytest
import igraph as ig

import molpy as mp
from molpy.analysis.cluster.chain import ChainFinder


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
    return ChainFinder(name="test")

class TestChainFinder:
    """Test the ChainFinder class."""
    
    def test_init(self, chain_finder):
        """Test ChainFinder initialization."""
        assert isinstance(chain_finder, ChainFinder)
    
    def test_compute(self, chain_finder, mock_frame):
        """Test ChainFinder computation."""
        context = chain_finder(mock_frame)
        assert "test_chain_idx" in context.result
        result = context.result["test_chain_idx"]
        assert isinstance(result, list)
        assert np.unique(result).size == 8