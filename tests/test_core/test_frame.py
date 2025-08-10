import numpy as np
import pytest
import tempfile
import os
from io import StringIO
from molpy.core.frame import Frame, Block

@pytest.fixture
def simple_frame() -> Frame:
    f = Frame()
    f["atoms", "xyz"] = np.arange(9).reshape(3, 3)
    f["atoms", "charge"] = np.array([-1.0, 0.5, 0.5])
    f["bonds", "i"] = np.arange(3)
    return f

class TestBlock:
    def test_block_set_and_get(self):
        blk = Block()
        blk["foo"] = np.arange(5)
        assert np.array_equal(blk["foo"], np.arange(5))

    def test_block_len_and_iter(self):
        blk = Block({"a": np.arange(2), "b": np.ones(2)})
        assert set(blk) == {"a", "b"}
        assert len(blk) == 2

    def test_block_to_from_dict(self):
        blk = Block({"a": np.arange(2), "b": np.ones(2)})
        dct = blk.to_dict()
        restored = Block.from_dict(dct)
        for k in ("a", "b"):
            assert np.array_equal(blk[k], restored[k])

    def test_block_from_csv_file(self):
        """Test Block.from_csv with file path."""
        csv_content = """x,y,z,atom_type
0.0,0.0,0.0,C
1.0,1.0,1.0,O
2.0,2.0,2.0,N"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            temp_file = f.name
        
        try:
            block = Block.from_csv(temp_file)
            
            # Test basic structure
            assert set(block.keys()) == {"x", "y", "z", "atom_type"}
            assert block.nrows == 3
            assert block.shape == (3, 4)
            
            # Test numeric columns
            assert np.array_equal(block["x"], np.array([0.0, 1.0, 2.0]))
            assert np.array_equal(block["y"], np.array([0.0, 1.0, 2.0]))
            assert np.array_equal(block["z"], np.array([0.0, 1.0, 2.0]))
            
            # Test string column
            assert np.array_equal(block["atom_type"], np.array(["C", "O", "N"]))
            
        finally:
            os.unlink(temp_file)

    def test_block_from_csv_stringio(self):
        """Test Block.from_csv with StringIO."""
        csv_content = """name,age,height
Alice,25,1.65
Bob,30,1.80
Charlie,35,1.75"""
        
        csv_io = StringIO(csv_content)
        block = Block.from_csv(csv_io)
        
        # Test basic structure
        assert set(block.keys()) == {"name", "age", "height"}
        assert block.nrows == 3
        assert block.shape == (3, 3)
        
        # Test mixed data types
        assert np.array_equal(block["name"], np.array(["Alice", "Bob", "Charlie"]))
        assert np.array_equal(block["age"], np.array([25.0, 30.0, 35.0]))
        assert np.array_equal(block["height"], np.array([1.65, 1.80, 1.75]))

    def test_block_from_csv_delimiter(self):
        """Test Block.from_csv with custom delimiter."""
        csv_content = """x;y;z;atom_type
0.0;0.0;0.0;C
1.0;1.0;1.0;O"""
        
        csv_io = StringIO(csv_content)
        block = Block.from_csv(csv_io, delimiter=";")
        
        assert set(block.keys()) == {"x", "y", "z", "atom_type"}
        assert block.nrows == 2
        assert np.array_equal(block["x"], np.array([0.0, 1.0]))

    def test_block_from_csv_empty_file(self):
        """Test Block.from_csv with empty file."""
        csv_io = StringIO("")
        with pytest.raises(ValueError, match="CSV file is empty"):
            Block.from_csv(csv_io)

    def test_block_from_csv_no_header(self):
        """Test Block.from_csv with no header CSV."""
        csv_content = """0.0,0.0,0.0,C
1.0,1.0,1.0,O
2.0,2.0,2.0,N"""
        
        csv_io = StringIO(csv_content)
        block = Block.from_csv(csv_io, header=["x", "y", "z", "atom_type"])
        
        # Test basic structure
        assert set(block.keys()) == {"x", "y", "z", "atom_type"}
        assert block.nrows == 3
        assert block.shape == (3, 4)
        
        # Test data
        assert np.array_equal(block["x"], np.array([0.0, 1.0, 2.0]))
        assert np.array_equal(block["y"], np.array([0.0, 1.0, 2.0]))
        assert np.array_equal(block["z"], np.array([0.0, 1.0, 2.0]))
        assert np.array_equal(block["atom_type"], np.array(["C", "O", "N"]))

    def test_block_from_csv_no_header_file(self):
        """Test Block.from_csv with no header CSV file."""
        csv_content = """0.0,0.0,0.0
1.0,1.0,1.0
2.0,2.0,2.0"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.csv', delete=False) as f:
            f.write(csv_content)
            temp_file = f.name
        
        try:
            block = Block.from_csv(temp_file, header=["x", "y", "z"])
            
            assert set(block.keys()) == {"x", "y", "z"}
            assert block.nrows == 3
            assert np.array_equal(block["x"], np.array([0.0, 1.0, 2.0]))
            assert np.array_equal(block["y"], np.array([0.0, 1.0, 2.0]))
            assert np.array_equal(block["z"], np.array([0.0, 1.0, 2.0]))
            
        finally:
            os.unlink(temp_file)

    def test_block_from_csv_header_mismatch(self):
        """Test Block.from_csv with header length mismatch."""
        csv_content = """0.0,0.0,0.0
1.0,1.0,1.0"""
        
        csv_io = StringIO(csv_content)
        # Header has more columns than data
        with pytest.raises(IndexError):
            Block.from_csv(csv_io, header=["x", "y", "z", "extra"])
        
        csv_io.seek(0)
        # Header has fewer columns than data
        block = Block.from_csv(csv_io, header=["x", "y"])
        assert set(block.keys()) == {"x", "y"}
        assert block.nrows == 2

    def test_block_from_csv_type_inference(self):
        """Test Block.from_csv with automatic type inference."""
        csv_content = """id,name,age,height,active
1,Alice,25,1.65,True
2,Bob,30,1.80,False
3,Charlie,35,1.75,True"""
        
        csv_io = StringIO(csv_content)
        block = Block.from_csv(csv_io)
        
        # Test type inference
        assert block["id"].dtype == np.dtype('int64')  # Should be int
        assert block["name"].dtype == np.dtype('<U7')  # Should be string
        assert block["age"].dtype == np.dtype('int64')  # Should be int (25, 30, 35 are integers)
        assert block["height"].dtype == np.dtype('float64')  # Should be float
        assert block["active"].dtype == np.dtype('<U5')  # Should be string
        
        # Test values
        assert np.array_equal(block["id"], np.array([1, 2, 3]))
        assert np.array_equal(block["name"], np.array(["Alice", "Bob", "Charlie"]))
        assert np.array_equal(block["age"], np.array([25, 30, 35]))
        assert np.array_equal(block["height"], np.array([1.65, 1.80, 1.75]))
        assert np.array_equal(block["active"], np.array(["True", "False", "True"]))

    def test_block_from_csv_mixed_types_no_header(self):
        """Test Block.from_csv with mixed types and no header."""
        csv_content = """1,Alice,25.5,True
2,Bob,30.0,False
3,Charlie,35.7,True"""
        
        csv_io = StringIO(csv_content)
        block = Block.from_csv(csv_io, header=["id", "name", "score", "active"])
        
        # Test type inference
        assert block["id"].dtype == np.dtype('int64')  # Should be int
        assert block["name"].dtype == np.dtype('<U7')  # Should be string
        assert block["score"].dtype == np.dtype('float64')  # Should be float
        assert block["active"].dtype == np.dtype('<U5')  # Should be string
        
        # Test values
        assert np.array_equal(block["id"], np.array([1, 2, 3]))
        assert np.array_equal(block["name"], np.array(["Alice", "Bob", "Charlie"]))
        assert np.array_equal(block["score"], np.array([25.5, 30.0, 35.7]))
        assert np.array_equal(block["active"], np.array(["True", "False", "True"]))

    def test_block_sort_basic(self):
        """Test basic sorting functionality."""
        blk = Block({
            "x": np.array([3, 1, 2]),
            "y": np.array([30, 10, 20]),
            "z": np.array([300, 100, 200])
        })
        
        sorted_blk = blk.sort("x")
        
        # Check that x is sorted
        assert np.array_equal(sorted_blk["x"], np.array([1, 2, 3]))
        # Check that other variables are sorted accordingly
        assert np.array_equal(sorted_blk["y"], np.array([10, 20, 30]))
        assert np.array_equal(sorted_blk["z"], np.array([100, 200, 300]))

    def test_block_sort_reverse(self):
        """Test reverse sorting."""
        blk = Block({
            "x": np.array([1, 2, 3]),
            "y": np.array([10, 20, 30])
        })
        
        sorted_blk = blk.sort("x", reverse=True)
        
        assert np.array_equal(sorted_blk["x"], np.array([3, 2, 1]))
        assert np.array_equal(sorted_blk["y"], np.array([30, 20, 10]))

    def test_block_sort_strings(self):
        """Test sorting with string data."""
        blk = Block({
            "name": np.array(["Charlie", "Alice", "Bob"]),
            "age": np.array([35, 25, 30])
        })
        
        sorted_blk = blk.sort("name")
        
        assert np.array_equal(sorted_blk["name"], np.array(["Alice", "Bob", "Charlie"]))
        assert np.array_equal(sorted_blk["age"], np.array([25, 30, 35]))

    def test_block_sort_float(self):
        """Test sorting with float data."""
        blk = Block({
            "score": np.array([3.14, 1.41, 2.71]),
            "name": np.array(["pi", "sqrt2", "e"])
        })
        
        sorted_blk = blk.sort("score")
        
        assert np.array_equal(sorted_blk["score"], np.array([1.41, 2.71, 3.14]))
        assert np.array_equal(sorted_blk["name"], np.array(["sqrt2", "e", "pi"]))

    def test_block_sort_empty(self):
        """Test sorting empty block."""
        blk = Block()
        sorted_blk = blk.sort("x")
        assert len(sorted_blk) == 0

    def test_block_sort_single_variable(self):
        """Test sorting block with single variable."""
        blk = Block({"x": np.array([3, 1, 2])})
        sorted_blk = blk.sort("x")
        assert np.array_equal(sorted_blk["x"], np.array([1, 2, 3]))

    def test_block_sort_key_error(self):
        """Test sorting with non-existent key."""
        blk = Block({"x": np.array([1, 2, 3])})
        with pytest.raises(KeyError, match="Variable 'y' not found in block"):
            blk.sort("y")

    def test_block_sort_length_mismatch(self):
        """Test sorting with variables of different lengths."""
        blk = Block({
            "x": np.array([1, 2, 3]),
            "y": np.array([10, 20])  # Different length
        })
        
        with pytest.raises(ValueError, match="Variable 'y' has different length than 'x'"):
            blk.sort("x")

    def test_block_sort_immutable(self):
        """Test that sorting doesn't modify the original block."""
        blk = Block({
            "x": np.array([3, 1, 2]),
            "y": np.array([30, 10, 20])
        })
        
        original_x = blk["x"].copy()
        original_y = blk["y"].copy()
        
        sorted_blk = blk.sort("x")
        
        # Original should remain unchanged
        assert np.array_equal(blk["x"], original_x)
        assert np.array_equal(blk["y"], original_y)
        
        # Sorted should be different
        assert not np.array_equal(sorted_blk["x"], original_x)

    def test_block_boolean_mask_indexing(self):
        """Test Block boolean mask indexing."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3]),
            "mass": np.array([12.0, 14.0, 12.0, 16.0, 14.0])
        })
        
        # Simple boolean mask
        mask = np.array([True, False, True, False, True])
        masked_blk = blk[mask]
        
        assert masked_blk.nrows == 3
        assert np.array_equal(masked_blk["id"], np.array([1, 3, 5]))
        assert np.array_equal(masked_blk["type"], np.array(['A', 'A', 'B']))
        assert np.array_equal(masked_blk["mol"], np.array([1, 2, 3]))

    def test_block_condition_mask_indexing(self):
        """Test Block condition-based mask indexing."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3]),
            "mass": np.array([12.0, 14.0, 12.0, 16.0, 14.0])
        })
        
        # Condition-based mask
        mol_mask = blk["mol"] < 3
        masked_blk = blk[mol_mask]
        
        assert masked_blk.nrows == 4
        assert np.array_equal(masked_blk["id"], np.array([1, 2, 3, 4]))
        assert np.array_equal(masked_blk["mol"], np.array([1, 1, 2, 2]))

    def test_block_complex_mask_indexing(self):
        """Test Block complex condition mask indexing."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3]),
            "mass": np.array([12.0, 14.0, 12.0, 16.0, 14.0])
        })
        
        # Complex condition mask
        complex_mask = (blk["mol"] < 3) & (blk["mass"] > 12.0)
        masked_blk = blk[complex_mask]
        
        assert masked_blk.nrows == 2
        assert np.array_equal(masked_blk["id"], np.array([2, 4]))
        assert np.array_equal(masked_blk["type"], np.array(['B', 'C']))

    def test_block_integer_array_indexing(self):
        """Test Block integer array indexing."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3]),
            "mass": np.array([12.0, 14.0, 12.0, 16.0, 14.0])
        })
        
        # Integer array indexing
        indices = np.array([0, 2, 4])
        indexed_blk = blk[indices]
        
        assert indexed_blk.nrows == 3
        assert np.array_equal(indexed_blk["id"], np.array([1, 3, 5]))
        assert np.array_equal(indexed_blk["type"], np.array(['A', 'A', 'B']))

    def test_block_reverse_integer_indexing(self):
        """Test Block reverse integer array indexing."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3]),
            "mass": np.array([12.0, 14.0, 12.0, 16.0, 14.0])
        })
        
        # Reverse order indexing
        reverse_indices = np.array([4, 3, 2, 1, 0])
        reverse_blk = blk[reverse_indices]
        
        assert reverse_blk.nrows == 5
        assert np.array_equal(reverse_blk["id"], np.array([5, 4, 3, 2, 1]))
        assert np.array_equal(reverse_blk["type"], np.array(['B', 'C', 'A', 'B', 'A']))

    def test_block_mask_wrong_length(self):
        """Test Block mask with wrong length raises ValueError."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B'])
        })
        
        # Wrong length boolean mask
        wrong_mask = np.array([True, False, True])
        
        with pytest.raises(IndexError):
            blk[wrong_mask]

    def test_block_mask_empty_result(self):
        """Test Block mask that results in empty selection."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3])
        })
        
        # Mask that selects nothing
        empty_mask = np.array([False, False, False, False, False])
        empty_blk = blk[empty_mask]
        
        assert empty_blk.nrows == 0
        assert len(empty_blk) == 3  # Still has all variables
        assert "id" in empty_blk
        assert "type" in empty_blk
        assert "mol" in empty_blk

    def test_block_mask_all_true(self):
        """Test Block mask that selects all rows."""
        blk = Block({
            "id": np.array([1, 2, 3, 4, 5]),
            "type": np.array(['A', 'B', 'A', 'C', 'B']),
            "mol": np.array([1, 1, 2, 2, 3])
        })
        
        # Mask that selects everything
        all_true_mask = np.array([True, True, True, True, True])
        all_blk = blk[all_true_mask]
        
        assert all_blk.nrows == 5
        assert np.array_equal(all_blk["id"], blk["id"])
        assert np.array_equal(all_blk["type"], blk["type"])
        assert np.array_equal(all_blk["mol"], blk["mol"])

    def test_block_mask_with_xyz_coordinates(self):
        """Test Block mask with xyz coordinates."""
        blk = Block({
            "id": np.array([1, 2, 3]),
            "xyz": np.array([
                [0.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                [0.0, 1.0, 0.0]
            ]),
            "type": np.array(['C', 'H', 'H'])
        })
        
        # Select atoms with x > 0.5
        x_mask = blk["xyz"][:, 0] > 0.5
        masked_blk = blk[x_mask]
        
        assert masked_blk.nrows == 1
        assert np.array_equal(masked_blk["id"], np.array([2]))
        assert np.array_equal(masked_blk["type"], np.array(['H']))
        assert np.array_equal(masked_blk["xyz"], np.array([[1.0, 0.0, 0.0]]))

    def test_block_mask_real_world_example(self):
        """Test Block mask with real-world molecular example."""
        # Create a molecular system with multiple molecules
        blk = Block({
            "id": np.arange(1, 21),  # 20 atoms
            "mol": np.array([1, 1, 1, 1, 1,  # molecule 1 (5 atoms)
                             2, 2, 2, 2, 2,  # molecule 2 (5 atoms)
                             3, 3, 3, 3, 3,  # molecule 3 (5 atoms)
                             4, 4, 4, 4, 4]), # molecule 4 (5 atoms)
            "type": np.array(['C', 'H', 'H', 'H', 'H'] * 4),  # 4 methane molecules
            "mass": np.array([12.0, 1.0, 1.0, 1.0, 1.0] * 4)
        })
        
        # Select atoms from molecules 1 and 2 (like in the notebook example)
        molid = blk["mol"]
        poly_atoms = blk[molid < 3]
        
        assert poly_atoms.nrows == 10
        assert np.array_equal(np.unique(poly_atoms["mol"]), np.array([1, 2]))
        assert np.array_equal(poly_atoms["id"], np.arange(1, 11))
        
        # Select heavy atoms from molecules 1 and 2
        heavy_atoms = blk[(molid < 3) & (blk["mass"] > 5.0)]
        
        assert heavy_atoms.nrows == 2
        assert np.array_equal(heavy_atoms["type"], np.array(['C', 'C']))
        assert np.array_equal(heavy_atoms["mol"], np.array([1, 2]))

class TestFrame:
    def test_set_and_get_variable(self, simple_frame):
        assert np.isclose(simple_frame["atoms", "charge"][0], -1.0)
        assert np.array_equal(simple_frame["bonds", "i"], np.arange(3))

    def test_setitem_creates_block(self):
        f = Frame()
        f["foo"] = {"bar": np.ones(4)}
        assert "foo" in list(f.blocks())
        assert np.array_equal(f["foo", "bar"], np.ones(4))

    def test_get_block(self, simple_frame):
        blk = simple_frame["atoms"]
        assert isinstance(blk, Block)
        assert set(blk) == {"xyz", "charge"}

    def test_variables(self, simple_frame):
        assert set(simple_frame.variables("atoms")) == {"xyz", "charge"}
        assert set(simple_frame.variables("bonds")) == {"i"}

    def test_blocks_iter_and_len(self, simple_frame):
        blocks = set(simple_frame.blocks())
        assert blocks == {"atoms", "bonds"}
        assert len(list(simple_frame.blocks())) == 2

    def test_delete_variable(self, simple_frame):
        del simple_frame["atoms"]["charge"]
        assert "charge" not in simple_frame.variables("atoms")

    def test_delete_block(self, simple_frame):
        del simple_frame._blocks["bonds"]
        assert "bonds" not in set(simple_frame.blocks())

    def test_to_from_dict_roundtrip(self, simple_frame):
        dct = simple_frame.to_dict()
        restored = Frame.from_dict(dct)
        for g in restored.blocks():
            for v in restored.variables(g):
                assert np.array_equal(restored[g, v], simple_frame[g, v])

    def test_forbid_assigning_non_block(self):
        f = Frame()
        with pytest.raises(ValueError):
            f["invalid"] = np.array([1, 2, 3])  # type: ignore[assignment]