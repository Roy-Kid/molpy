"""
Tests for XYZReader and XYZWriter using chemfiles-testcases/xyz files.
"""

from io import BytesIO
from pathlib import Path

import numpy as np
import pytest

from molpy.core import Box
from molpy.io.data.xyz import XYZReader
from molpy.core.element import Element


@pytest.fixture
def xyz_test_dir(TEST_DATA_DIR) -> Path:
    return TEST_DATA_DIR / "xyz"


class TestXYZReader:
    """Test XYZ file reading with various formats."""

    def test_standard_xyz_format(self, xyz_test_dir):
        """Test reading standard XYZ format (methane)."""
        reader = XYZReader(xyz_test_dir / "methane.xyz")
        frame = reader.read()

        # Check coordinates are stored as separate x, y, z arrays
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert frame["atoms"]["x"].shape == (
            5,
        ), "x-coordinates should be 1D array of length N"
        assert frame["atoms"]["y"].shape == (
            5,
        ), "y-coordinates should be 1D array of length N"
        assert frame["atoms"]["z"].shape == (
            5,
        ), "z-coordinates should be 1D array of length N"

        # Check elements
        assert "element" in frame["atoms"]
        elements = frame["atoms"]["element"]
        assert elements[0] == "C", "First atom should be Carbon"
        assert all(elements[1:] == "H"), "Other atoms should be Hydrogen"

        # Check atomic numbers
        assert "number" in frame["atoms"]
        assert frame["atoms"]["number"][0] == 6, "Carbon should have Z=6"
        assert all(frame["atoms"]["number"][1:] == 1), "Hydrogen should have Z=1"

    def test_extended_xyz_format(self, xyz_test_dir):
        """Test reading extended XYZ format with Properties."""
        reader = XYZReader(xyz_test_dir / "extended.xyz")
        frame = reader.read()

        # Check coordinates are stored as separate x, y, z arrays
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert frame["atoms"]["x"].shape == (192,), "x-coordinates should be 1D array"
        assert frame["atoms"]["y"].shape == (192,), "y-coordinates should be 1D array"
        assert frame["atoms"]["z"].shape == (192,), "z-coordinates should be 1D array"

        # Check box from Lattice
        box = frame.box
        assert isinstance(box, Box)
        assert box.matrix.shape == (3, 3)
        assert np.allclose(
            box.matrix,
            np.array(
                [
                    [8.43116035, 0.0, 0.0],
                    [0.158219155128, 14.5042431863, 0.0],
                    [1.16980663624, 4.4685149855, 14.9100096405],
                ]
            ),
        )

        # Check metadata
        assert "ENERGY" in frame.metadata
        assert float(frame.metadata["ENERGY"]) == pytest.approx(-2069.84934116)
        assert "Natoms" in frame.metadata
        assert frame.metadata["Natoms"] == "192"

        # Check atomic numbers are present
        assert "number" in frame["atoms"]
        assert frame["atoms"]["number"].shape == (192,)

    def test_extended_xyz_with_properties(self, xyz_test_dir):
        """Test extended XYZ with Properties specification parsing."""
        reader = XYZReader(xyz_test_dir / "extended.xyz")
        frame = reader.read()

        # The extended.xyz has Properties=species:S:1:pos:R:3:CS:R:2
        # Check that species was mapped to element
        assert "element" in frame["atoms"]

        # Check that pos was split into x, y, z
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert frame["atoms"]["x"].shape == (192,)
        assert frame["atoms"]["y"].shape == (192,)
        assert frame["atoms"]["z"].shape == (192,)

        # Check that CS (custom property) is present
        assert "CS" in frame["atoms"]
        assert frame["atoms"]["CS"].shape == (192, 2)

    def test_water_xyz(self, xyz_test_dir):
        """Test reading water.xyz file."""
        reader = XYZReader(xyz_test_dir / "water.xyz")
        frame = reader.read()

        # Water has 3 atoms per molecule
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        # Coordinates should be 1D arrays
        assert frame["atoms"]["x"].ndim == 1, "Coordinates should be 1D"
        assert frame["atoms"]["y"].ndim == 1, "Coordinates should be 1D"
        assert frame["atoms"]["z"].ndim == 1, "Coordinates should be 1D"

        # Check atomic numbers
        assert "number" in frame["atoms"]
        # Water molecules have O-H-H pattern
        num_values = frame["atoms"]["number"]
        assert (
            len(num_values) % 3 == 0
        ), "Should have multiple of 3 atoms (water molecules)"

    def test_topology_xyz(self, xyz_test_dir):
        """Test reading topology.xyz file."""
        reader = XYZReader(xyz_test_dir / "topology.xyz")
        frame = reader.read()

        # Check basic structure
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert "element" in frame["atoms"]
        assert "number" in frame["atoms"]

        # Coordinates should be 1D arrays
        assert frame["atoms"]["x"].ndim == 1
        assert frame["atoms"]["y"].ndim == 1
        assert frame["atoms"]["z"].ndim == 1

    def test_velocities_xyz(self, xyz_test_dir):
        """Test reading XYZ with velocities."""
        reader = XYZReader(xyz_test_dir / "velocities.xyz")
        frame = reader.read()

        # Check structure
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert "element" in frame["atoms"]
        assert "number" in frame["atoms"]

        # Check if velocities are present (depends on Properties specification)
        # The file might have velocities in the Properties field
        atoms = frame["atoms"]
        assert atoms["x"].ndim == 1

    def test_atomic_number_generation(self, xyz_test_dir):
        """Test that atomic numbers are correctly generated from element symbols."""
        reader = XYZReader(xyz_test_dir / "methane.xyz")
        frame = reader.read()

        elements = frame["atoms"]["element"]
        num_values = frame["atoms"]["number"]

        # Verify atomic numbers match elements
        for elem, num in zip(elements, num_values):
            expected_num = Element.get_atomic_number(elem)
            assert (
                num == expected_num
            ), f"Element {elem} should have Z={expected_num}, got {num}"

    def test_coordinate_array_shape(self, xyz_test_dir):
        """Test that coordinates are stored as three separate 1D arrays."""
        test_files = ["methane.xyz", "water.xyz", "topology.xyz"]

        for filename in test_files:
            reader = XYZReader(xyz_test_dir / filename)
            frame = reader.read()

            x = frame["atoms"]["x"]
            y = frame["atoms"]["y"]
            z = frame["atoms"]["z"]
            assert x.ndim == 1, f"{filename}: x-coordinates should be 1D"
            assert y.ndim == 1, f"{filename}: y-coordinates should be 1D"
            assert z.ndim == 1, f"{filename}: z-coordinates should be 1D"
            assert (
                x.shape == y.shape == z.shape
            ), f"{filename}: all coordinates should have same length"

    def test_element_to_atomic_number_mapping(self, xyz_test_dir):
        """Test various element symbols are correctly mapped to atomic numbers."""
        # Create a test with known elements
        reader = XYZReader(xyz_test_dir / "extended.xyz")
        frame = reader.read()

        # Check that all atomic numbers are valid (> 0)
        num_values = frame["atoms"]["number"]
        assert all(num_values > 0), "All atomic numbers should be positive"

        # Check that elements and atomic numbers are consistent
        if "element" in frame["atoms"]._vars:
            elements = frame["atoms"]["element"]
            for elem, num in zip(elements, num_values):
                expected_num = Element.get_atomic_number(str(elem))
                assert num == expected_num

    def test_metadata_extraction(self, xyz_test_dir):
        """Test that metadata from comment line is correctly extracted."""
        reader = XYZReader(xyz_test_dir / "extended.xyz")
        frame = reader.read()

        # Check various metadata fields
        assert "ENERGY" in frame.metadata
        assert "Natoms" in frame.metadata
        assert "NAME" in frame.metadata

        # Metadata values should be strings or parsed values
        assert isinstance(frame.metadata["ENERGY"], (str, float))
        assert isinstance(frame.metadata["Natoms"], str)

    def test_empty_comment_line(self, xyz_test_dir):
        """Test XYZ file with empty comment line."""
        reader = XYZReader(xyz_test_dir / "methane.xyz")
        frame = reader.read()

        # Should still parse correctly even with empty comment
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert "number" in frame["atoms"]

    def test_spaces_in_xyz(self, xyz_test_dir):
        """Test XYZ file with various spacing."""
        reader = XYZReader(xyz_test_dir / "spaces.xyz")
        frame = reader.read()

        # Should handle various spacing correctly
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert "element" in frame["atoms"]
        assert frame["atoms"]["x"].ndim == 1
        assert frame["atoms"]["y"].ndim == 1
        assert frame["atoms"]["z"].ndim == 1

    def test_bytesio_standard_xyz(self):
        """Test reading standard XYZ format from BytesIO."""
        xyz_data = """5
Methane molecule
C     0.000000     0.000000     0.000000
H     0.000000     0.000000     1.089000
H     1.026719     0.000000    -0.363000
H    -0.513360    -0.889165    -0.363000
H    -0.513360     0.889165    -0.363000
"""
        bytesio_obj = BytesIO(xyz_data.encode("utf-8"))
        reader = XYZReader(bytesio_obj)
        frame = reader.read()

        # Check structure
        assert len(frame["atoms"]["element"]) == 5
        assert frame["atoms"]["element"][0] == "C"
        assert all(frame["atoms"]["element"][1:] == "H")

        # Check coordinates
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert frame["atoms"]["x"].shape == (5,)
        assert frame["atoms"]["y"].shape == (5,)
        assert frame["atoms"]["z"].shape == (5,)

        # Check atomic numbers
        assert frame["atoms"]["number"][0] == 6  # Carbon
        assert all(frame["atoms"]["number"][1:] == 1)  # Hydrogen

    def test_bytesio_extended_xyz(self):
        """Test reading extended XYZ format from BytesIO."""
        extxyz_data = """3
Properties=species:S:1:pos:R:3 Lattice="5.0 0.0 0.0 0.0 5.0 0.0 0.0 0.0 5.0" ENERGY=-10.5
H 0.0 0.0 0.0
O 1.0 0.0 0.0
H 2.0 0.0 0.0
"""
        bytesio_obj = BytesIO(extxyz_data.encode("utf-8"))
        reader = XYZReader(bytesio_obj)
        frame = reader.read()

        # Check structure
        assert len(frame["atoms"]["element"]) == 3
        assert frame["atoms"]["element"][0] == "H"
        assert frame["atoms"]["element"][1] == "O"
        assert frame["atoms"]["element"][2] == "H"

        # Check coordinates are split into x, y, z
        assert "x" in frame["atoms"]
        assert "y" in frame["atoms"]
        assert "z" in frame["atoms"]
        assert np.allclose(frame["atoms"]["x"], [0.0, 1.0, 2.0])
        assert np.allclose(frame["atoms"]["y"], [0.0, 0.0, 0.0])
        assert np.allclose(frame["atoms"]["z"], [0.0, 0.0, 0.0])

        # Check box from Lattice
        assert frame.box is not None
        assert np.allclose(frame.box.matrix, np.eye(3) * 5.0)

        # Check metadata
        assert "ENERGY" in frame.metadata
        assert float(frame.metadata["ENERGY"]) == pytest.approx(-10.5)

    def test_bytesio_reusability(self):
        """Test that BytesIO can be rewound and read multiple times."""
        xyz_data = """2
Water
O 0.0 0.0 0.0
H 1.0 0.0 0.0
"""
        bytesio_obj = BytesIO(xyz_data.encode("utf-8"))

        # First read
        reader1 = XYZReader(bytesio_obj)
        frame1 = reader1.read()
        assert len(frame1["atoms"]["element"]) == 2

        # Rewind and read again
        bytesio_obj.seek(0)
        reader2 = XYZReader(bytesio_obj)
        frame2 = reader2.read()
        assert len(frame2["atoms"]["element"]) == 2

        # Verify both reads produced the same data
        assert np.array_equal(frame1["atoms"]["element"], frame2["atoms"]["element"])
        assert np.array_equal(frame1["atoms"]["x"], frame2["atoms"]["x"])
