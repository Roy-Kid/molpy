import csv
import tempfile
from io import StringIO
from pathlib import Path

import pytest

from molpy.core.frame import Block
from molpy.core.utils import read_csv, to_dict_of_list, to_list_of_dict


class TestReadCSV:
    """Test the read_csv function."""

    def test_read_csv_from_file(self):
        """Test reading CSV from a file."""
        # Create a temporary CSV file
        csv_content = """name,age,city
John,30,New York
Jane,25,London
Bob,35,Paris"""

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False) as f:
            f.write(csv_content)
            temp_file = Path(f.name)

        try:
            # Test reading the CSV file
            result = read_csv(temp_file)

            # Verify the result is a Block
            assert isinstance(result, Block)

            # Verify the data content
            assert len(result) == 3  # 3 columns
            assert "name" in result
            assert "age" in result
            assert "city" in result

            # Check specific values
            assert result["name"].tolist() == ["John", "Jane", "Bob"]
            assert result["age"].tolist() == ["30", "25", "35"]
            assert result["city"].tolist() == ["New York", "London", "Paris"]

        finally:
            # Clean up
            temp_file.unlink()

    def test_read_csv_from_stringio(self):
        """Test reading CSV from StringIO."""
        csv_content = """id,value,status
1,100,active
2,200,inactive
3,300,active"""

        string_io = StringIO(csv_content)
        result = read_csv(string_io)

        # Verify the result is a Block
        assert isinstance(result, Block)

        # Verify the data content
        assert len(result) == 3  # 3 columns
        assert "id" in result
        assert "value" in result
        assert "status" in result

        # Check specific values
        assert result["id"].tolist() == ["1", "2", "3"]
        assert result["value"].tolist() == ["100", "200", "300"]
        assert result["status"].tolist() == ["active", "inactive", "active"]

    def test_read_csv_with_custom_delimiter(self):
        """Test reading CSV with custom delimiter."""
        csv_content = """name;score;grade
Alice;95;A
Bob;87;B
Charlie;92;A"""

        string_io = StringIO(csv_content)
        result = read_csv(string_io, delimiter=";")

        # Verify the result is a Block
        assert isinstance(result, Block)

        # Verify the data content
        assert len(result) == 3  # 3 columns
        assert "name" in result
        assert "score" in result
        assert "grade" in result

        # Check specific values
        assert result["name"].tolist() == ["Alice", "Bob", "Charlie"]
        assert result["score"].tolist() == ["95", "87", "92"]
        assert result["grade"].tolist() == ["A", "B", "A"]

    def test_read_csv_empty_file(self):
        """Test reading an empty CSV file."""
        csv_content = """name,age
"""

        string_io = StringIO(csv_content)
        result = read_csv(string_io)

        # Verify the result is a Block
        assert isinstance(result, Block)

        # Should have columns but no data
        assert len(result) == 2  # 2 columns
        assert "name" in result
        assert "age" in result

        # Check that arrays are empty
        assert len(result["name"]) == 0
        assert len(result["age"]) == 0

    def test_read_csv_file_not_found(self):
        """Test reading a non-existent file."""
        non_existent_file = Path("/tmp/non_existent_file.csv")

        with pytest.raises(FileNotFoundError, match="File .* does not exist"):
            read_csv(non_existent_file)

    def test_read_csv_with_spaces_in_values(self):
        """Test reading CSV with spaces in values."""
        csv_content = """name,description,value
Product A,High quality product,100
Product B,Low cost option,50
Product C,Premium choice,200"""

        string_io = StringIO(csv_content)
        result = read_csv(string_io)

        # Verify the result is a Block
        assert isinstance(result, Block)

        # Check specific values with spaces
        assert result["name"].tolist() == ["Product A", "Product B", "Product C"]
        assert result["description"].tolist() == [
            "High quality product",
            "Low cost option",
            "Premium choice",
        ]
        assert result["value"].tolist() == ["100", "50", "200"]

    def test_read_csv_with_quoted_values(self):
        """Test reading CSV with quoted values containing commas."""
        csv_content = """name,address,phone
"John Doe","123 Main St, City",555-1234
"Jane Smith","456 Oak Ave, Town",555-5678"""

        string_io = StringIO(csv_content)
        result = read_csv(string_io)

        # Verify the result is a Block
        assert isinstance(result, Block)

        # Check that quoted values are properly parsed
        assert result["name"].tolist() == ["John Doe", "Jane Smith"]
        assert result["address"].tolist() == ["123 Main St, City", "456 Oak Ave, Town"]
        assert result["phone"].tolist() == ["555-1234", "555-5678"]


class TestUtilityFunctions:
    """Test other utility functions."""

    def test_to_dict_of_list(self):
        """Test converting list of dicts to dict of lists."""
        list_of_dict = [
            {"name": "Alice", "age": 30, "city": "NYC"},
            {"name": "Bob", "age": 25, "city": "LA"},
            {"name": "Charlie", "age": 35, "city": "Chicago"},
        ]

        result = to_dict_of_list(list_of_dict)

        expected = {
            "name": ["Alice", "Bob", "Charlie"],
            "age": [30, 25, 35],
            "city": ["NYC", "LA", "Chicago"],
        }

        assert result == expected

    def test_to_list_of_dict(self):
        """Test converting dict of lists to list of dicts."""
        dict_of_list = {
            "name": ["Alice", "Bob", "Charlie"],
            "age": [30, 25, 35],
            "city": ["NYC", "LA", "Chicago"],
        }

        result = to_list_of_dict(dict_of_list)

        expected = [
            {"name": "Alice", "age": 30, "city": "NYC"},
            {"name": "Bob", "age": 25, "city": "LA"},
            {"name": "Charlie", "age": 35, "city": "Chicago"},
        ]

        assert result == expected

    def test_roundtrip_conversion(self):
        """Test that converting back and forth preserves data."""
        original = [{"x": 1, "y": 2}, {"x": 3, "y": 4}, {"x": 5, "y": 6}]

        # Convert to dict of lists and back
        dict_version = to_dict_of_list(original)
        back_to_list = to_list_of_dict(dict_version)

        assert back_to_list == original
