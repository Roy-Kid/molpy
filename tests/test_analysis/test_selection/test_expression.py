"""
Unit tests for ComputeContext and SelectResult classes.
"""

import pytest
import numpy as np
from molpy.analysis.base import ComputeContext
from molpy.analysis.selection.expression import SelectResult, ExpressionSelection
from molpy.core.frame import Frame


class MockFrame(Frame):
    """Mock frame for testing."""
    
    def __init__(self, n_atoms: int = 100):
        super().__init__()
        # Create atoms block with positions
        self["atoms"] = {
            "xyz": np.random.random((n_atoms, 3)),
            "id": np.arange(n_atoms)
        }
        self.n_atoms = n_atoms
        

class TestComputeContext:
    """Test cases for ComputeContext class."""
    
    def test_compute_context_initialization(self):
        """Test ComputeContext initialization."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        assert context.frame is frame
        assert isinstance(context.result, dict)
        assert len(context.result) == 0
    
    def test_compute_context_result_storage(self):
        """Test storing results in ComputeContext."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        # Store some results
        context.result["test_key"] = "test_value"
        context.result["atoms"] = np.array([1, 2, 3])
        
        assert context.result["test_key"] == "test_value"
        assert np.array_equal(context.result["atoms"], np.array([1, 2, 3]))
    
    def test_compute_context_frame_property(self):
        """Test frame property accessor."""
        frame = MockFrame(30)
        context = ComputeContext.attach_frame(frame)
        
        assert context.frame == frame
        assert frame.n_atoms == 30
    
    def test_compute_context_chain_operations(self):
        """Test chain operations in ComputeContext."""
        frame = MockFrame(50)
        root_context = ComputeContext.attach_frame(frame)
        
        # Create a chain of contexts
        context1 = ComputeContext(root_context)
        context2 = ComputeContext(context1)
        context3 = ComputeContext(context2)
        
        # Test chain depth
        assert root_context.get_depth() == 0
        assert context1.get_depth() == 1
        assert context2.get_depth() == 2
        assert context3.get_depth() == 3
        
        # Test pop operations
        assert context3.pop() is context2
        assert context2.pop() is context1
        assert context1.pop() is root_context
        
        # Test get_stack
        chain = context3.get_stack()
        assert len(chain) == 4
        assert chain[0] is context3
        assert chain[1] is context2
        assert chain[2] is context1
        assert chain[3] is root_context
        
        # Test frame access through chain
        assert context3.frame is frame


class TestSelectResult:
    """Test cases for SelectResult class."""
    
    def test_select_result_is_compute_context(self):
        """Test that SelectResult is a ComputeContext."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        result = SelectResult(context)
        
        # SelectResult should be an instance of ComputeContext
        assert isinstance(result, ComputeContext)
        assert isinstance(result, SelectResult)
    
    def test_select_result_initialization(self):
        """Test SelectResult initialization."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        result = SelectResult(context)
        
        # Should access frame through context chain
        assert result.frame is frame
        assert result.selected is None
        
        # Should store previous context
        assert result.unwrap() is context
    
    def test_select_result_without_selection(self):
        """Test SelectResult initialization without selection."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        result = SelectResult(context)
        
        # Should not have any selections by default
        assert "atoms" not in result.result
    
    def test_select_result_pop(self):
        """Test popping selection to get previous context."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        result = SelectResult(context)
        popped = result.pop()
        
        assert popped is context
    
    def test_select_result_chain_depth_single(self):
        """Test chain depth with single SelectResult."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        result = SelectResult(context)
        
        assert result.get_depth() == 1
    
    def test_select_result_chain_depth_multiple(self):
        """Test chain depth with multiple SelectResults."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        # Chain multiple selections
        result1 = SelectResult(context)
        result1.selected = np.array([1, 3, 5, 7, 9])
        
        result2 = SelectResult(result1)
        result2.selected = np.array([1, 5, 9])
        
        result3 = SelectResult(result2)
        result3.selected = np.array([1, 9])
        
        assert result1.get_depth() == 1
        assert result2.get_depth() == 2
        assert result3.get_depth() == 3
    
    def test_select_result_select_stack(self):
        """Test getting selection stack from context chain."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        # Create mixed chain with regular context and SelectResults
        regular_context = ComputeContext(context)
        result1 = SelectResult(regular_context)
        result1.selected = np.array([1, 3, 5])
        
        result2 = SelectResult(result1)
        result2.selected = np.array([1, 5])
        
        # Get selection stack
        select_stack = result2.select_stack()
        
        # Should only contain SelectResult instances
        assert len(select_stack) == 2
        assert all(isinstance(ctx, SelectResult) for ctx in select_stack)
        assert select_stack[0] is result2  # Most recent first
        assert select_stack[1] is result1


class TestExpressionSelection:
    """Test cases for ExpressionSelection class."""
    
    def test_expression_selection_basic(self):
        """Test basic expression selection."""
        frame = MockFrame(100)
        context = ComputeContext.attach_frame(frame)
        
        # Expression to select even indices
        def select_even(ctx: ComputeContext) -> np.ndarray:
            n_atoms = len(ctx.frame["atoms"]["id"])
            return np.arange(n_atoms) % 2 == 0
        
        selection = ExpressionSelection("even_atoms", select_even)
        result = selection.compute(context)
        
        # Should return SelectResult
        assert isinstance(result, SelectResult)
        assert isinstance(result, ComputeContext)
        
        # Should have selected even indices
        expected = np.arange(0, 100, 2)  # [0, 2, 4, 6, ...]
        assert result.selected is not None
        assert np.array_equal(result.selected, expected)
    
    def test_expression_selection_chained(self):
        """Test chained expression selections."""
        frame = MockFrame(100)
        context = ComputeContext.attach_frame(frame)
        
        # First selection: indices < 50
        def select_first_half(ctx: ComputeContext) -> np.ndarray:
            n_atoms = len(ctx.frame["atoms"]["id"])
            return np.arange(n_atoms) < 50
        
        # Second selection: even indices from frame (independent of first selection)
        def select_even(ctx: ComputeContext) -> np.ndarray:
            n_atoms = len(ctx.frame["atoms"]["id"])
            return np.arange(n_atoms) % 2 == 0
        
        selection1 = ExpressionSelection("first_half", select_first_half)
        selection2 = ExpressionSelection("even_indices", select_even)
        
        # Chain selections
        result1 = selection1.compute(context)
        result2 = selection2.compute(result1)
        
        # Verify chain structure
        assert isinstance(result2, SelectResult)
        assert result2.get_depth() == 2
        
        # Get selection stack to access historical selections
        select_stack = result2.select_stack()
        assert len(select_stack) == 2
        
        # First selection should be [0, 1, 2, ..., 49]
        first_result = select_stack[1]  # Older selection
        assert isinstance(first_result, SelectResult)
        expected_first = np.arange(50)
        assert first_result.selected is not None
        assert np.array_equal(first_result.selected, expected_first)
        
        # Second selection should be all even numbers from frame [0, 2, 4, ..., 98]
        second_result = select_stack[0]  # Current selection
        assert isinstance(second_result, SelectResult)
        expected_second = np.arange(0, 100, 2)  # [0, 2, 4, ..., 98]
        assert second_result.selected is not None
        assert np.array_equal(second_result.selected, expected_second)
    
    def test_expression_selection_custom_field(self):
        """Test expression selection with custom field name."""
        frame = MockFrame(50)
        context = ComputeContext.attach_frame(frame)
        
        def select_small(ctx: ComputeContext) -> np.ndarray:
            n_atoms = len(ctx.frame["atoms"]["id"])
            return np.arange(n_atoms) < 10
        
        selection = ExpressionSelection("small_indices", select_small, select_field="atoms")
        result = selection.compute(context)
        
        # Should store selected IDs in selected attribute
        assert isinstance(result, SelectResult)
        expected = np.arange(10)
        assert result.selected is not None
        assert np.array_equal(result.selected, expected)
    
    def test_compute_context_error_handling(self):
        """Test error handling in ComputeContext operations."""
        frame = MockFrame(10)
        context = ComputeContext.attach_frame(frame)
        
        # Test pop on root context
        with pytest.raises(ValueError, match="No previous context found"):
            context.pop()
    
    def test_select_result_nested_chains(self):
        """Test deeply nested selection chains."""
        frame = MockFrame(20)
        context = ComputeContext.attach_frame(frame)
        
        # Create a deep chain of selections
        current = context
        results = []
        for i in range(5):
            result = SelectResult(current)
            result.selected = np.arange(i, 20-i)  # Different selections
            results.append(result)
            current = result
        
        # Test chain properties
        assert current.get_depth() == 5
        
        # Test select_stack filtering
        select_stack = current.select_stack()
        assert len(select_stack) == 5
        assert all(isinstance(r, SelectResult) for r in select_stack)
        
        # Check order: most recent first
        for i, result in enumerate(select_stack):
            expected_selected = np.arange(4-i, 20-(4-i))  # Reverse order
            assert result.selected is not None
            assert np.array_equal(result.selected, expected_selected)
    
    def test_mixed_context_chain(self):
        """Test chains with mixed ComputeContext and SelectResult."""
        frame = MockFrame(30)
        root = ComputeContext.attach_frame(frame)
        
        # Mixed chain: Context -> SelectResult -> Context -> SelectResult
        ctx1 = ComputeContext(root)
        sel1 = SelectResult(ctx1)
        sel1.selected = np.array([1, 3, 5])
        
        ctx2 = ComputeContext(sel1)
        sel2 = SelectResult(ctx2)
        sel2.selected = np.array([7, 9, 11])
        
        # Test total chain depth
        assert sel2.get_depth() == 4
        
        # Test full stack
        full_stack = sel2.get_stack()
        assert len(full_stack) == 5  # sel2, ctx2, sel1, ctx1, root
        
        # Test select_stack filtering
        select_stack = sel2.select_stack()
        assert len(select_stack) == 2  # Only SelectResults
        assert select_stack[0] is sel2
        assert select_stack[1] is sel1


if __name__ == "__main__":
    pytest.main([__file__])
