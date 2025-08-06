"""
Test the new flexible map interface for Compute class.
"""

import pytest
import numpy as np
from molpy.analysis.base import Compute, ComputeContext, get_selected_entities, get_entities_from_frame
from molpy.analysis.selection.expression import SelectResult, ExpressionSelection
from molpy.core.frame import Frame


class SimpleCompute(Compute):
    """Simple compute that just returns the entity_id squared."""
    
    def compute(self, context: ComputeContext) -> ComputeContext:
        entity_id = context.result.get('entity_id', 0)
        context.result['squared'] = entity_id ** 2
        context.result['doubled'] = entity_id * 2
        return context


class MockFrame(Frame):
    """Mock frame for testing."""
    
    def __init__(self, n_entities: int = 5):
        super().__init__()
        self["molecules"] = {
            "id": np.arange(n_entities),
            "name": [f"MOL_{i}" for i in range(n_entities)]
        }
        self["atoms"] = {
            "id": np.arange(n_entities * 3),  # 3 atoms per molecule
            "mol_id": np.repeat(np.arange(n_entities), 3)
        }


def test_map_with_explicit_entity_ids():
    """Test map with explicitly provided entity IDs."""
    frame = MockFrame(5)
    context = ComputeContext.attach_frame(frame)
    compute = SimpleCompute("test")
    
    # Map over specific entities
    entity_ids = np.array([1, 3, 4])
    result = compute.map(context, entity_ids=entity_ids, selection_field="molecules")
    
    # Check results
    assert "molecules_results" in result.result
    assert "mapped_molecules" in result.result
    
    molecules_results = result.result["molecules_results"]
    mapped_molecules = result.result["mapped_molecules"]
    
    assert len(molecules_results) == 3
    assert np.array_equal(mapped_molecules, entity_ids)
    
    # Check individual results
    assert molecules_results[1]['squared'] == 1
    assert molecules_results[3]['squared'] == 9
    assert molecules_results[4]['squared'] == 16
    
    assert molecules_results[1]['doubled'] == 2
    assert molecules_results[3]['doubled'] == 6
    assert molecules_results[4]['doubled'] == 8


def test_map_all_entities():
    """Test map over all entities when no entity_ids provided."""
    frame = MockFrame(3)
    context = ComputeContext.attach_frame(frame)
    compute = SimpleCompute("test")
    
    # Map over all molecules
    result = compute.map(context, selection_field="molecules")
    
    molecules_results = result.result["molecules_results"]
    mapped_molecules = result.result["mapped_molecules"]
    
    assert len(molecules_results) == 3
    assert np.array_equal(mapped_molecules, [0, 1, 2])
    
    for i in range(3):
        assert molecules_results[i]['squared'] == i ** 2
        assert molecules_results[i]['doubled'] == i * 2


def test_map_with_selection_result():
    """Test map using entity IDs from SelectResult."""
    frame = MockFrame(5)
    context = ComputeContext.attach_frame(frame)
    
    # Create selection
    def select_even(ctx: ComputeContext) -> np.ndarray:
        n_mols = len(ctx.frame["molecules"]["id"])
        return np.arange(n_mols) % 2 == 0
    
    selection = ExpressionSelection("even_molecules", select_even, select_field="molecules")
    selected_context = selection.compute(context)
    
    # Get selected IDs using utility function
    selected_ids = get_selected_entities(selected_context)
    assert selected_ids is not None
    assert np.array_equal(selected_ids, [0, 2, 4])
    
    # Map using selected IDs
    compute = SimpleCompute("test")
    result = compute.map(context, entity_ids=selected_ids, selection_field="molecules")
    
    molecules_results = result.result["molecules_results"]
    assert len(molecules_results) == 3
    assert 0 in molecules_results
    assert 2 in molecules_results
    assert 4 in molecules_results
    assert 1 not in molecules_results
    assert 3 not in molecules_results


def test_utility_functions():
    """Test utility functions for getting entity IDs."""
    frame = MockFrame(4)
    context = ComputeContext.attach_frame(frame)
    
    # Test get_entities_from_frame
    mol_ids = get_entities_from_frame(context, "molecules")
    assert mol_ids is not None
    assert np.array_equal(mol_ids, [0, 1, 2, 3])
    
    atom_ids = get_entities_from_frame(context, "atoms")
    assert atom_ids is not None
    assert np.array_equal(atom_ids, np.arange(12))  # 4 molecules * 3 atoms
    
    # Test non-existent field
    none_ids = get_entities_from_frame(context, "nonexistent")
    assert none_ids is None
    
    # Test get_selected_entities with no selection
    selected = get_selected_entities(context)
    assert selected is None


def test_map_different_fields():
    """Test mapping over different entity fields."""
    frame = MockFrame(3)
    context = ComputeContext.attach_frame(frame)
    compute = SimpleCompute("test")
    
    # Map over atoms
    atom_result = compute.map(context, selection_field="atoms")
    atom_results = atom_result.result["atoms_results"]
    
    assert len(atom_results) == 9  # 3 molecules * 3 atoms
    assert atom_results[0]['squared'] == 0
    assert atom_results[5]['squared'] == 25
    assert atom_results[8]['squared'] == 64


def test_map_error_cases():
    """Test error cases for map method."""
    frame = MockFrame(3)
    context = ComputeContext.attach_frame(frame)
    compute = SimpleCompute("test")
    
    # Test with empty entity_ids
    with pytest.raises(ValueError, match="No entities found for mapping"):
        compute.map(context, entity_ids=np.array([]), selection_field="molecules")
    
    # Test with non-existent field
    with pytest.raises(ValueError, match="No entities found for mapping"):
        compute.map(context, selection_field="nonexistent")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
