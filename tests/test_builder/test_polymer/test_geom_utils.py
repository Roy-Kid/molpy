"""
Tests for geometry utility functions (Rodrigues, VDW, orientation inference).
"""

import pytest
import math
from molpy.builder.polymer.geom_utils import rodrigues, get_vdw_radius
from molpy.core.ops.geometry import _dot, _cross, _norm, _unit
from molpy.core.entity import Entity


class TestRodrigues:
    """Test Rodrigues rotation matrix computation."""
    
    def test_identity_rotation(self):
        """Test that aligning a vector to itself gives identity."""
        v = [1.0, 0.0, 0.0]
        R = rodrigues(v, v)
        
        # Should be identity matrix
        assert R[0] == pytest.approx([1.0, 0.0, 0.0], abs=1e-10)
        assert R[1] == pytest.approx([0.0, 1.0, 0.0], abs=1e-10)
        assert R[2] == pytest.approx([0.0, 0.0, 1.0], abs=1e-10)
    
    def test_90_degree_rotation(self):
        """Test 90° rotation from x to y axis."""
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        R = rodrigues(a, b)
        
        # Apply rotation to a
        result = [
            R[0][0] * a[0] + R[0][1] * a[1] + R[0][2] * a[2],
            R[1][0] * a[0] + R[1][1] * a[1] + R[1][2] * a[2],
            R[2][0] * a[0] + R[2][1] * a[1] + R[2][2] * a[2],
        ]
        
        # Result should be b
        assert result == pytest.approx(b, abs=1e-10)
    
    def test_180_degree_rotation(self):
        """Test 180° rotation (antiparallel vectors)."""
        a = [1.0, 0.0, 0.0]
        b = [-1.0, 0.0, 0.0]
        R = rodrigues(a, b)
        
        # Apply rotation to a
        result = [
            R[0][0] * a[0] + R[0][1] * a[1] + R[0][2] * a[2],
            R[1][0] * a[0] + R[1][1] * a[1] + R[1][2] * a[2],
            R[2][0] * a[0] + R[2][1] * a[1] + R[2][2] * a[2],
        ]
        
        # Result should be b
        assert result == pytest.approx(b, abs=1e-10)
    
    def test_arbitrary_rotation(self):
        """Test rotation between arbitrary unit vectors."""
        a = _unit([1.0, 1.0, 1.0])
        b = _unit([1.0, -1.0, 0.0])
        R = rodrigues(a, b)
        
        # Apply rotation
        result = [
            R[0][0] * a[0] + R[0][1] * a[1] + R[0][2] * a[2],
            R[1][0] * a[0] + R[1][1] * a[1] + R[1][2] * a[2],
            R[2][0] * a[0] + R[2][1] * a[1] + R[2][2] * a[2],
        ]
        
        # Result should be b
        assert result == pytest.approx(b, abs=1e-10)
    
    def test_rotation_determinant_is_one(self):
        """Test that rotation matrix has determinant +1."""
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        R = rodrigues(a, b)
        
        # Compute determinant
        det = (
            R[0][0] * (R[1][1] * R[2][2] - R[1][2] * R[2][1])
            - R[0][1] * (R[1][0] * R[2][2] - R[1][2] * R[2][0])
            + R[0][2] * (R[1][0] * R[2][1] - R[1][1] * R[2][0])
        )
        
        assert det == pytest.approx(1.0, abs=1e-10)


class TestVectorOps:
    """Test basic vector operations."""
    
    def test__unit(self):
        """Test vector normalization."""
        v = [3.0, 4.0, 0.0]
        u = _unit(v)
        assert u == pytest.approx([0.6, 0.8, 0.0], abs=1e-10)
        assert _norm(u) == pytest.approx(1.0, abs=1e-10)
    
    def test_dot_product(self):
        """Test dot product."""
        a = [1.0, 2.0, 3.0]
        b = [4.0, 5.0, 6.0]
        assert _dot(a, b) == 32.0
    
    def test_cross_product(self):
        """Test cross product."""
        a = [1.0, 0.0, 0.0]
        b = [0.0, 1.0, 0.0]
        c = _cross(a, b)
        assert c == pytest.approx([0.0, 0.0, 1.0], abs=1e-10)


class TestVDWRadius:
    """Test VDW radius lookup."""
    
    def test_explicit_vdw_in_entity(self):
        """Test that explicit vdw in entity.data is used."""
        entity = Entity()
        entity.data['vdw'] = 2.5
        assert get_vdw_radius(entity) == 2.5
    
    def test_element_lookup_carbon(self):
        """Test Element class lookup for carbon."""
        entity = Entity()
        entity.data['symbol'] = 'C'
        radius = get_vdw_radius(entity)
        assert radius == pytest.approx(1.70, abs=0.01)
    
    def test_element_lookup_hydrogen(self):
        """Test Element class lookup for hydrogen."""
        entity = Entity()
        entity.data['symbol'] = 'H'
        radius = get_vdw_radius(entity)
        assert radius == pytest.approx(1.20, abs=0.01)
    
    def test_element_lookup_oxygen(self):
        """Test Element class lookup for oxygen."""
        entity = Entity()
        entity.data['symbol'] = 'O'
        radius = get_vdw_radius(entity)
        assert radius == pytest.approx(1.52, abs=0.01)
    
    def test_fallback_to_carbon(self):
        """Test fallback to carbon for unknown elements."""
        entity = Entity()
        # No symbol - should default to carbon
        radius = get_vdw_radius(entity)
        assert radius == pytest.approx(1.70, abs=0.01)
