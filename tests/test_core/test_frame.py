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