"""
Tests for DataReader base class, focusing on BytesIO and file-like object support.
"""

from io import BytesIO, StringIO
from pathlib import Path

import pytest

from molpy.core.frame import Frame
from molpy.io.data.base import DataReader


# Create a concrete implementation of DataReader for testing
class SimpleTestReader(DataReader):
    """Simple test reader that reads a single integer from the first line."""

    def read(self, frame: Frame | None = None) -> Frame:
        """Read a simple format: first line is number of atoms, rest is ignored."""
        frame = frame or Frame()
        lines = self.read_lines()
        if lines:
            # Store the first line as metadata
            frame.metadata["first_line"] = lines[0].strip()
            if len(lines) > 1:
                frame.metadata["second_line"] = lines[1].strip()
        return frame


class TestDataReaderBytesIO:
    """Test DataReader base class with BytesIO support."""

    def test_bytesio_basic_reading(self):
        """Test basic reading from BytesIO."""
        data = b"Line 1\nLine 2\nLine 3\n"
        bytesio = BytesIO(data)

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        assert frame.metadata["first_line"] == "Line 1"
        assert frame.metadata["second_line"] == "Line 2"

    def test_bytesio_with_utf8_encoding(self):
        """Test BytesIO with UTF-8 encoded text."""
        data = "第一行\n第二行\n".encode("utf-8")
        bytesio = BytesIO(data)

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        assert frame.metadata["first_line"] == "第一行"
        assert frame.metadata["second_line"] == "第二行"

    def test_bytesio_empty(self):
        """Test reading from empty BytesIO."""
        bytesio = BytesIO(b"")

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        # Should not crash, just return empty frame
        assert "first_line" not in frame.metadata

    def test_bytesio_single_line(self):
        """Test reading from BytesIO with single line."""
        bytesio = BytesIO(b"Only one line")

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        assert frame.metadata["first_line"] == "Only one line"
        assert "second_line" not in frame.metadata

    def test_bytesio_rewind_and_reread(self):
        """Test that BytesIO can be rewound and read multiple times."""
        data = b"Test data\nSecond line\n"
        bytesio = BytesIO(data)

        # First read
        reader1 = SimpleTestReader(bytesio)
        frame1 = reader1.read()

        # Rewind
        bytesio.seek(0)

        # Second read
        reader2 = SimpleTestReader(bytesio)
        frame2 = reader2.read()

        # Both reads should produce the same result
        assert frame1.metadata["first_line"] == frame2.metadata["first_line"]
        assert frame1.metadata["second_line"] == frame2.metadata["second_line"]

    def test_stringio_support(self):
        """Test that StringIO is also supported."""
        stringio = StringIO("String line 1\nString line 2\n")

        reader = SimpleTestReader(stringio)
        frame = reader.read()

        assert frame.metadata["first_line"] == "String line 1"
        assert frame.metadata["second_line"] == "String line 2"

    def test_file_path_still_works(self, tmp_path):
        """Test that traditional file path reading still works."""
        # Create a temporary file
        test_file = tmp_path / "test.txt"
        test_file.write_text("Path line 1\nPath line 2\n")

        reader = SimpleTestReader(test_file)
        frame = reader.read()

        assert frame.metadata["first_line"] == "Path line 1"
        assert frame.metadata["second_line"] == "Path line 2"

    def test_file_path_as_string(self, tmp_path):
        """Test reading from file path as string."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("String path 1\nString path 2\n")

        reader = SimpleTestReader(str(test_file))
        frame = reader.read()

        assert frame.metadata["first_line"] == "String path 1"
        assert frame.metadata["second_line"] == "String path 2"

    def test_bytesio_with_context_manager(self):
        """Test BytesIO with context manager usage."""
        bytesio = BytesIO(b"Context line 1\nContext line 2\n")

        with SimpleTestReader(bytesio) as reader:
            frame = reader.read()

        assert frame.metadata["first_line"] == "Context line 1"
        assert frame.metadata["second_line"] == "Context line 2"

    def test_file_path_with_context_manager(self, tmp_path):
        """Test file path with context manager usage."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("CM line 1\nCM line 2\n")

        with SimpleTestReader(test_file) as reader:
            frame = reader.read()

        assert frame.metadata["first_line"] == "CM line 1"
        assert frame.metadata["second_line"] == "CM line 2"

    def test_bytesio_read_lines_method(self):
        """Test the read_lines() method with BytesIO."""
        bytesio = BytesIO(b"Line A\nLine B\nLine C\n")

        reader = SimpleTestReader(bytesio)
        lines = reader.read_lines()

        assert len(lines) == 3
        assert lines[0].strip() == "Line A"
        assert lines[1].strip() == "Line B"
        assert lines[2].strip() == "Line C"

    def test_bytesio_iteration(self):
        """Test iterating over non-blank lines with BytesIO."""
        bytesio = BytesIO(b"Line 1\n\nLine 3\n  \nLine 5\n")

        reader = SimpleTestReader(bytesio)
        lines = list(reader)

        # Should skip blank lines
        assert len(lines) == 3
        assert lines[0] == "Line 1"
        assert lines[1] == "Line 3"
        assert lines[2] == "Line 5"

    def test_bytesio_with_special_characters(self):
        """Test BytesIO with special characters and whitespace."""
        data = "  Spaces  \n\tTabs\t\nMixed \t \n".encode("utf-8")
        bytesio = BytesIO(data)

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        # read_lines() preserves whitespace, but we strip in our reader
        assert frame.metadata["first_line"] == "Spaces"
        assert frame.metadata["second_line"] == "Tabs"

    def test_bytesio_multiline_content(self):
        """Test BytesIO with realistic multi-line content."""
        content = """# Comment line
Data line 1
Data line 2
Data line 3
# Another comment
Data line 4
"""
        bytesio = BytesIO(content.encode("utf-8"))

        reader = SimpleTestReader(bytesio)
        lines = reader.read_lines()

        assert len(lines) == 6
        assert lines[0].strip() == "# Comment line"
        assert lines[1].strip() == "Data line 1"

    def test_bytesio_preserves_underlying_object(self):
        """Test that BytesIO object is not closed by the reader."""
        bytesio = BytesIO(b"Test content\n")

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        # BytesIO should still be usable after reading
        bytesio.seek(0)
        content = bytesio.read()
        assert content == b"Test content\n"

    def test_type_detection_with_read_attribute(self):
        """Test that objects with 'read' attribute are detected as file-like."""
        # BytesIO has 'read' attribute
        bytesio = BytesIO(b"Test\n")
        reader = SimpleTestReader(bytesio)
        assert reader._is_file_object is True

        # StringIO has 'read' attribute
        stringio = StringIO("Test\n")
        reader2 = SimpleTestReader(stringio)
        assert reader2._is_file_object is True

    def test_type_detection_with_path(self, tmp_path):
        """Test that paths are detected correctly."""
        test_file = tmp_path / "test.txt"
        test_file.write_text("Test\n")

        # Path object
        reader = SimpleTestReader(test_file)
        assert reader._is_file_object is False

        # String path
        reader2 = SimpleTestReader(str(test_file))
        assert reader2._is_file_object is False

    def test_bytesio_with_binary_content(self):
        """Test BytesIO with content that needs UTF-8 decoding."""
        # Create content with various UTF-8 characters
        content = "Hello 世界 🌍\nTest ñ é ü\n".encode("utf-8")
        bytesio = BytesIO(content)

        reader = SimpleTestReader(bytesio)
        frame = reader.read()

        assert "世界" in frame.metadata["first_line"]
        assert "🌍" in frame.metadata["first_line"]
        assert "ñ" in frame.metadata["second_line"]


class TestDataReaderEdgeCases:
    """Test edge cases and error conditions."""

    def test_bytesio_with_windows_line_endings(self):
        """Test BytesIO with Windows-style line endings."""
        bytesio = BytesIO(b"Line 1\r\nLine 2\r\nLine 3\r\n")

        reader = SimpleTestReader(bytesio)
        lines = reader.read_lines()

        # Should handle \r\n correctly
        assert len(lines) == 3

    def test_bytesio_with_mac_line_endings(self):
        """Test BytesIO with old Mac-style line endings."""
        bytesio = BytesIO(b"Line 1\rLine 2\rLine 3\r")

        reader = SimpleTestReader(bytesio)
        lines = reader.read_lines()

        # Should handle \r correctly
        assert len(lines) >= 1

    def test_bytesio_with_mixed_line_endings(self):
        """Test BytesIO with mixed line endings."""
        bytesio = BytesIO(b"Line 1\nLine 2\r\nLine 3\rLine 4\n")

        reader = SimpleTestReader(bytesio)
        lines = reader.read_lines()

        # Should handle mixed line endings
        assert len(lines) >= 3

    def test_large_bytesio_content(self):
        """Test BytesIO with large content."""
        # Create a large content
        lines = [f"Line {i}\n" for i in range(10000)]
        content = "".join(lines).encode("utf-8")
        bytesio = BytesIO(content)

        reader = SimpleTestReader(bytesio)
        read_lines = reader.read_lines()

        assert len(read_lines) == 10000
        assert read_lines[0].strip() == "Line 0"
        assert read_lines[9999].strip() == "Line 9999"
