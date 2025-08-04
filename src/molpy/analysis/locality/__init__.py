"""
Locality analysis tools for molpy.

This module provides tools for analyzing spatial relationships between atoms,
including neighbor list generation and spatial clustering.
"""

from dataclasses import dataclass
from typing import Optional, Tuple, List
import numpy as np
import molpy as mp
import freud

@dataclass
class NeighborListResult:
    """
    Result from neighbor list computation.
    
    Attributes:
        pairs: List of (i, j) atom index pairs within cutoff
        distances: Optional array of distances for each pair
    """
    pairs: List[Tuple[int, int]]
    distances: np.array | None = None


class Neighborlist:
    """
    Stateless neighbor list calculator for molpy.Frame objects.
    
    This class computes atom pairs within a specified cutoff distance,
    supporting periodic boundary conditions. It uses freud-analysis as the
    backend for efficient neighbor finding with O(N log N) performance.
    The class is stateless and does not cache any per-frame data.
    
    Attributes:
        cutoff: Distance cutoff for neighbor pairs (in same units as positions)
    """
    
    def __init__(self, cutoff: float):
        """
        Initialize the neighbor list calculator.
        
        Args:
            cutoff: Distance cutoff for neighbor pairs
        """
        self.cutoff = cutoff
    
    def compute(self, frame: mp.Frame) -> NeighborListResult:
        """
        Compute neighbor pairs for a given frame.
        
        Args:
            frame: molpy.Frame object containing atomic positions
            
        Returns:
            NeighborListResult containing atom pairs and optional distances
            
        Raises:
            ValueError: If frame doesn't contain position data
            ValueError: If cutoff is negative
        """
        if self.cutoff < 0:
            raise ValueError("Cutoff must be non-negative")
        
        # Extract positions from frame
        positions = self._get_positions(frame)
        if positions is None:
            raise ValueError("Frame must contain position data ('xyz' or 'positions')")
        
        # Extract box information
        box = frame.box
        pbc = box.pbc if box else np.array([False, False, False])
        cell = box.matrix if box else None
        
        # Compute neighbor pairs using freud backend
        pairs, distances = self._compute_pairs_freud(positions, cell, pbc)
        
        return NeighborListResult(pairs=pairs, distances=distances)
    
    def _get_positions(self, frame: mp.Frame) -> np.array:
        """
        Extract positions from frame, checking common variable names.
        
        Args:
            frame: molpy.Frame object
            
        Returns:
            Position array if found, None otherwise
        """
        return frame["atoms"]["xyz"]
    
    def _compute_pairs_freud(
        self, 
        positions: np.ndarray, 
        cell: np.array | None, 
        pbc: np.ndarray
    ) -> Tuple[List[Tuple[int, int]], np.array | None]:
        """
        Compute atom pairs within cutoff distance using freud backend.
        
        Args:
            positions: Atomic positions (N, 3)
            cell: Cell matrix (3, 3) or None for non-periodic
            pbc: Periodic boundary conditions (3,) boolean array
            
        Returns:
            Tuple of (pairs, distances)
        """
        # Create freud box
        if cell is not None and np.any(pbc):
            # Periodic system
            freud_box = freud.box.Box.from_matrix(cell)
        else:
            # Non-periodic system
            freud_box = freud.box.Box.cube(1e6)  # Large box for non-periodic
        
        # Create freud neighbor query
        query = freud.locality.AABBQuery(freud_box, positions)
        
        # Query for neighbors within cutoff
        query_args = dict(mode='ball', r_max=self.cutoff, exclude_ii=True)
        result = query.query(positions, query_args)
        
        # Convert to our format
        nlist = result.toNeighborList()

        return NeighborListResult(
            pairs=nlist.pairs,
            distances=nlist.distances,
        )