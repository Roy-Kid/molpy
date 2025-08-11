import tempfile
from pathlib import Path

import numpy as np
import pytest

import molpy as mp
from molpy.io.trajectory.lammps import LammpsTrajectoryReader, LammpsTrajectoryWriter


class TestReadLammpsTrajectory:

    def test_read_frame_simple(self, TEST_DATA_DIR):
        reader = LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        )
        frame = reader.read_frame(0)

        assert frame.metadata["timestep"] is not None
        assert frame.metadata["box"] is not None
        assert frame.metadata["box"].matrix.shape == (3, 3)
        assert frame["atoms"].nrows > 0

    def test_read_frame_with_properties(self, TEST_DATA_DIR):
        reader = LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/properties.lammpstrj"
        )
        frame = reader.read_frame(0)

        assert frame.metadata["timestep"] is not None
        assert frame.metadata["box"] is not None
        assert frame.metadata["box"].matrix.shape == (3, 3)
        assert frame["atoms"].nrows > 0

        # Check that atoms block exists and has content
        atoms = frame["atoms"]
        assert len(atoms) > 0  # Should have variable names

    def test_trajectory_iteration(self, TEST_DATA_DIR):
        """Test iterating through trajectory frames."""
        reader = LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        )

        # Test iterator protocol
        frames = []
        for frame in reader:
            frames.append(frame)
            if len(frames) >= 3:  # Limit to avoid long test
                break

        assert len(frames) > 0
        for frame in frames:
            assert isinstance(frame, mp.Frame)
            assert frame.metadata["timestep"] is not None
            assert "atoms" in frame
            assert frame.metadata["box"] is not None

    def test_frame_properties(self, TEST_DATA_DIR):
        """Test that frames have correct properties."""
        reader = LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        )
        frame = reader.read_frame(0)

        # Check timestep
        assert frame.metadata["timestep"] is not None
        assert isinstance(frame.metadata["timestep"], (int, np.integer))

        # Check box
        assert frame.metadata["box"] is not None
        assert hasattr(frame.metadata["box"], "matrix")
        assert frame.metadata["box"].matrix.shape == (3, 3)

    def test_context_manager(self, TEST_DATA_DIR):
        """Test using trajectory reader as context manager."""
        with LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        ) as reader:
            frame = reader.read_frame(0)
            assert frame.metadata["timestep"] is not None
            assert "atoms" in frame

    def test_multiple_files_support(self, TEST_DATA_DIR):
        """Test reading from multiple trajectory files."""
        file1 = TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        file2 = TEST_DATA_DIR / "trajectory/lammps/properties.lammpstrj"

        # Test with multiple files if they exist
        if file1.exists() and file2.exists():
            reader = LammpsTrajectoryReader([file1, file2])

            # Should have frames from both files
            assert reader.n_frames > 0

            # Read some frames
            frame1 = reader.read_frame(0)
            assert frame1.metadata["timestep"] is not None
            assert "atoms" in frame1

            # Test reading range
            if reader.n_frames > 1:
                frames = reader.read_range(0, min(3, reader.n_frames))
                assert len(frames) >= 1
                for frame in frames:
                    assert isinstance(frame, mp.Frame)


class TestWriteLammpsTrajectory:

    def test_write_simple_trajectory(self):
        """Test writing a simple trajectory."""
        # Create test frames
        frames = []
        for i in range(3):
            frame = mp.Frame()

            # Create atoms data using Block structure
            atoms_data = {
                "id": [0, 1, 2],
                "type": [1, 1, 2],
                "x": [0.0 + i * 0.1, 1.0 + i * 0.1, 0.5 + i * 0.1],
                "y": [0.0, 0.0, 1.0],
                "z": [0.0, 0.0, 0.0],
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = i * 100
            frame.metadata["box"] = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)

        # Write trajectory
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            for frame in frames:
                timestep = frame.metadata["timestep"]
                writer.write_frame(frame, timestep=timestep)
            writer.close()

            # Read back and verify
            reader = LammpsTrajectoryReader(tmp.name)

            # Check that we can read the frames back
            for i, frame_read in enumerate(reader):
                if i >= len(frames):
                    break
                assert frame_read.metadata["timestep"] == frames[i].metadata["timestep"]
                assert "atoms" in frame_read
                # Check that positions changed over time
                if i > 0:
                    # x coordinates should be different between frames
                    x_vals = frame_read["atoms"]["x"]
                    x_vals_prev = frames[0]["atoms"]["x"]
                    assert not np.allclose(x_vals, x_vals_prev)

    def test_write_with_context_manager(self):
        """Test writing trajectory using context manager."""
        frame = mp.Frame()
        atoms_data = {
            "id": [0, 1],
            "type": [1, 1],
            "x": [0.0, 1.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 0
        frame.metadata["box"] = mp.Box(np.eye(3) * 5.0)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            with LammpsTrajectoryWriter(tmp.name) as writer:
                writer.write_frame(frame)

    def test_trajectory_roundtrip(self):
        """Test writing and reading back maintains data integrity."""
        # Create a more complex frame
        frame = mp.Frame()

        atoms_data = {
            "id": [0, 1, 2, 3],
            "type": [1, 1, 2, 2],
            "x": [0.0, 1.0, 0.5, 1.5],
            "y": [0.0, 0.0, 1.0, 1.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "vx": [0.1, -0.1, 0.2, -0.2],
            "vy": [0.0, 0.0, 0.1, -0.1],
            "vz": [0.0, 0.0, 0.0, 0.0],
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 1000
        frame.metadata["box"] = mp.Box(np.diag([5.0, 5.0, 5.0]))

        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            # Write
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame)
            writer.close()

            # Read back
            reader = LammpsTrajectoryReader(tmp.name)
            frame_read = reader.read_frame(0)

            # Verify timestep
            assert frame_read.metadata["timestep"] == 1000

            # Verify atoms data exists
            assert "atoms" in frame_read
            atoms = frame_read["atoms"]
            assert atoms.nrows == 4

            # Verify box
            assert frame_read.metadata["box"] is not None
            assert np.allclose(
                frame_read.metadata["box"].matrix.diagonal(), [5.0, 5.0, 5.0]
            )


class TestErrorHandling:

    def test_read_nonexistent_file(self):
        """Test reading non-existent trajectory file."""
        with pytest.raises((FileNotFoundError, IOError)):
            reader = LammpsTrajectoryReader("nonexistent.dump")

    def test_read_empty_file(self):
        """Test reading empty trajectory file."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            tmp.write("")  # Empty file

            # Should handle empty file gracefully
            with pytest.raises((IndexError, ValueError, EOFError)):
                reader = LammpsTrajectoryReader(tmp.name)

    def test_malformed_trajectory(self):
        """Test handling malformed trajectory files."""
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            tmp.write("ITEM: TIMESTEP\n")
            tmp.write("0\n")
            tmp.write("ITEM: NUMBER OF ATOMS\n")
            tmp.write("2\n")
            tmp.write("ITEM: BOX BOUNDS pp pp pp\n")
            tmp.write("0 10\n")
            tmp.write("0 10\n")
            tmp.write("0 10\n")
            tmp.write("ITEM: ATOMS id type x y z\n")
            tmp.write("1 1 0.0 0.0 0.0\n")
            # Missing second atom - malformed
            tmp.flush()  # Ensure data is written to disk

            reader = LammpsTrajectoryReader(tmp.name)
            # Should handle malformed file gracefully
            try:
                frame = reader.read_frame(0)
                # If it doesn't raise an error, check that the frame is reasonable
                assert isinstance(frame, mp.Frame)
            except (ValueError, IndexError, EOFError):
                # These exceptions are acceptable for malformed files
                pass

    def test_invalid_frame_index(self, TEST_DATA_DIR):
        """Test reading invalid frame indices."""
        reader = LammpsTrajectoryReader(
            TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj"
        )

        # Test too large index
        with pytest.raises((IndexError, ValueError)):
            reader.read_frame(999999)


class TestMultipleFilesTrajectory:
    """Test multiple file trajectory reading capabilities."""

    def test_create_multi_file_trajectory(self):
        """Test creating and reading a multi-file trajectory."""
        # Create multiple trajectory files
        filenames = []

        for file_idx in range(3):
            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".dump", delete=False
            ) as tmp:
                filenames.append(tmp.name)

                writer = LammpsTrajectoryWriter(tmp.name)

                # Write 2 frames per file
                for frame_idx in range(2):
                    frame = mp.Frame()

                    atoms_data = {
                        "id": [0, 1],
                        "type": [1, 2],
                        "x": [
                            0.0 + file_idx * 0.5 + frame_idx * 0.1,
                            1.0 + file_idx * 0.5 + frame_idx * 0.1,
                        ],
                        "y": [0.0, 0.0],
                        "z": [0.0, 0.0],
                    }
                    frame["atoms"] = atoms_data
                    frame.metadata["timestep"] = file_idx * 100 + frame_idx * 10
                    frame.metadata["box"] = mp.Box(np.eye(3) * 10.0)

                    writer.write_frame(frame)

                writer.close()

        # Test reading from multiple files
        reader = LammpsTrajectoryReader(filenames)

        # Should have 6 frames total (3 files × 2 frames each)
        assert reader.n_frames == 6

        # Test reading all frames
        all_frames = reader.read_all()
        assert len(all_frames) == 6

        # Check that timesteps are in expected order
        expected_timesteps = [0, 10, 100, 110, 200, 210]
        for i, frame in enumerate(all_frames):
            assert frame.metadata["timestep"] == expected_timesteps[i]

        # Test reading specific frames
        frame_2 = reader.read_frame(2)  # First frame of second file
        assert frame_2.metadata["timestep"] == 100

        frame_5 = reader.read_frame(5)  # Last frame
        assert frame_5.metadata["timestep"] == 210

        # Test reading range
        frames_1_to_3 = reader.read_range(1, 4)
        assert len(frames_1_to_3) == 3
        assert frames_1_to_3[0].metadata["timestep"] == 10
        assert frames_1_to_3[2].metadata["timestep"] == 110

        # Cleanup
        for filename in filenames:
            Path(filename).unlink()

    def test_multi_file_with_different_formats(self):
        """Test reading multiple files with potentially different atom formats."""
        filenames = []

        # File 1: basic format
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dump", delete=False
        ) as tmp1:
            filenames.append(tmp1.name)

            writer = LammpsTrajectoryWriter(tmp1.name)
            frame = mp.Frame()

            atoms_data = {
                "id": [0, 1],
                "type": [1, 1],
                "x": [0.0, 1.0],
                "y": [0.0, 0.0],
                "z": [0.0, 0.0],
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = 0
            frame.metadata["box"] = mp.Box(np.eye(3) * 10.0)

            writer.write_frame(frame)
            writer.close()

        # File 2: format with velocities
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".dump", delete=False
        ) as tmp2:
            filenames.append(tmp2.name)

            writer = LammpsTrajectoryWriter(tmp2.name)
            frame = mp.Frame()

            atoms_data = {
                "id": [0, 1, 2],
                "type": [1, 1, 2],
                "x": [0.0, 1.0, 0.5],
                "y": [0.0, 0.0, 1.0],
                "z": [0.0, 0.0, 0.0],
                "vx": [0.1, -0.1, 0.0],
                "vy": [0.0, 0.0, 0.1],
                "vz": [0.0, 0.0, 0.0],
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = 100
            frame.metadata["box"] = mp.Box(np.eye(3) * 10.0)

            writer.write_frame(frame)
            writer.close()

        # Read both files together
        reader = LammpsTrajectoryReader(filenames)

        assert reader.n_frames == 2

        # Read first frame (from file 1)
        frame1 = reader.read_frame(0)
        assert frame1.metadata["timestep"] == 0
        assert frame1["atoms"].nrows == 2

        # Read second frame (from file 2)
        frame2 = reader.read_frame(1)
        assert frame2.metadata["timestep"] == 100
        assert frame2["atoms"].nrows == 3

        # Cleanup
        for filename in filenames:
            Path(filename).unlink()


class TestTrajectoryIntegration:

    def test_data_to_trajectory_conversion(self):
        """Test converting data format to trajectory format."""
        # Create a frame in data format
        frame = mp.Frame()

        atoms_data = {
            "id": [0, 1, 2],
            "type": [1, 1, 2],  # Use numeric types for LAMMPS
            "x": [0.0, 0.816, -0.816],
            "y": [0.0, 0.577, 0.577],
            "z": [0.0, 0.0, 0.0],
            "q": [-0.8476, 0.4238, 0.4238],
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 0
        frame.metadata["box"] = mp.Box(np.diag([10.0, 10.0, 10.0]))

        # Write as trajectory
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame)
            writer.close()

            # Read back as trajectory
            reader = LammpsTrajectoryReader(tmp.name)
            frame_read = reader.read_frame(0)

            assert frame_read.metadata["timestep"] is not None
            assert "atoms" in frame_read
            assert frame_read.metadata["box"] is not None

    def test_multiple_formats_consistency(self):
        """Test that data and trajectory formats are consistent."""
        # This test ensures that a frame written in one format
        # can be meaningfully compared with the other format
        frame_original = mp.Frame()

        atoms_data = {
            "id": [0, 1],
            "type": [1, 2],
            "x": [0.0, 1.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        }
        frame_original["atoms"] = atoms_data
        frame_original.metadata["timestep"] = 100
        frame_original.metadata["box"] = mp.Box(np.eye(3) * 5.0)

        # Write as trajectory and read back
        with tempfile.NamedTemporaryFile(mode="w", suffix=".dump", delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame_original)
            writer.close()

            reader = LammpsTrajectoryReader(tmp.name)
            frame_traj = reader.read_frame(0)

            # Both should have same basic structure
            assert frame_traj.metadata["timestep"] is not None
            assert "atoms" in frame_traj
            assert frame_traj.metadata["box"] is not None
            assert frame_original.metadata["box"] is not None

            # Box dimensions should be similar
            assert np.allclose(
                frame_traj.metadata["box"].matrix.diagonal(),
                frame_original.metadata["box"].matrix.diagonal(),
            )
