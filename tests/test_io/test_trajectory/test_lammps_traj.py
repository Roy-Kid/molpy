import pytest
import tempfile
import numpy as np
import molpy as mp
from molpy.io.trajectory.lammps import LammpsTrajectoryReader, LammpsTrajectoryWriter
from pathlib import Path


class TestReadLammpsTrajectory:

    def test_read_frame_simple(self, TEST_DATA_DIR):
        reader = LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj")
        frame = reader.read_frame(0)

        assert frame.metadata["timestep"] is not None
        assert frame.box.matrix.shape == (3, 3)  # type: ignore
        assert frame["atoms"].nrows > 0
        assert frame.box is not None

    def test_read_frame_with_properties(self, TEST_DATA_DIR):
        reader = LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/properties.lammpstrj")
        frame = reader.read_frame(0)

        assert frame.metadata["timestep"] is not None
        assert frame.box.matrix.shape == (3, 3)  # type: ignore
        assert frame["atoms"].nrows > 0
        
        # Check that atoms block exists and has content
        atoms = frame["atoms"]
        assert len(atoms) > 0  # Should have variable names

    def test_read_multi_traj(self, TEST_DATA_DIR):
        reader = LammpsTrajectoryReader([TEST_DATA_DIR/"trajectory/lammps/unwrapped.lammpstrj", TEST_DATA_DIR/"trajectory/lammps/wrapped.lammpstrj"])

        assert reader.n_frames >= 2
        frame0 = reader.read_frame(0)
        assert frame0.metadata["timestep"] is not None
        assert len(frame0["atoms"]) > 0  # Check atoms exist
        frame1 = reader.read_frame(1)
        assert frame1.metadata["timestep"] is not None
        assert len(frame1["atoms"]) > 0  # Check atoms exist

    def test_trajectory_iteration(self, TEST_DATA_DIR):
        """Test iterating through trajectory frames."""
        reader = LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj")
        
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
            assert frame.box is not None

    def test_frame_properties(self, TEST_DATA_DIR):
        """Test that frames have correct properties."""
        reader = LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj")
        frame = reader.read_frame(0)
        
        # Check timestep
        assert frame.metadata["timestep"] is not None
        assert isinstance(frame.metadata["timestep"], (int, np.integer))
        
        # Check box
        assert frame.box is not None
        assert hasattr(frame.box, 'matrix')
        assert frame.box.matrix.shape == (3, 3)

    def test_context_manager(self, TEST_DATA_DIR):
        """Test using trajectory reader as context manager."""
        with LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj") as reader:
            frame = reader.read_frame(0)
            assert frame.metadata["timestep"] is not None
            assert "atoms" in frame


class TestWriteLammpsTrajectory:

    def test_write_simple_trajectory(self):
        """Test writing a simple trajectory."""
        # Create test frames
        frames = []
        for i in range(3):
            frame = mp.Frame()
            
            # Create atoms data using Block structure
            atoms_data = {
                'id': [0, 1, 2],
                'type': [1, 1, 2],
                'x': [0.0 + i*0.1, 1.0 + i*0.1, 0.5 + i*0.1],
                'y': [0.0, 0.0, 1.0],
                'z': [0.0, 0.0, 0.0]
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = i * 100
            frame.box = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)
        
        # Write trajectory
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            for frame in frames:
                writer.write_frame(frame)
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
            'id': [0, 1],
            'type': [1, 1],
            'x': [0.0, 1.0],
            'y': [0.0, 0.0],
            'z': [0.0, 0.0]
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 0
        frame.box = mp.Box(np.eye(3) * 5.0)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            with LammpsTrajectoryWriter(tmp.name) as writer:
                writer.write_frame(frame)

    def test_trajectory_roundtrip(self):
        """Test writing and reading back maintains data integrity."""
        # Create a more complex frame
        frame = mp.Frame()
        
        atoms_data = {
            'id': [0, 1, 2, 3],
            'type': [1, 1, 2, 2],
            'x': [0.0, 1.0, 0.5, 1.5],
            'y': [0.0, 0.0, 1.0, 1.0],
            'z': [0.0, 0.0, 0.0, 0.0],
            'vx': [0.1, -0.1, 0.2, -0.2],
            'vy': [0.0, 0.0, 0.1, -0.1],
            'vz': [0.0, 0.0, 0.0, 0.0]
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 1000
        frame.box = mp.Box(np.diag([5.0, 5.0, 5.0]))
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
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
            assert frame_read.box is not None
            assert np.allclose(frame_read.box.matrix.diagonal(), [5.0, 5.0, 5.0])


class TestErrorHandling:

    def test_read_nonexistent_file(self):
        """Test reading non-existent trajectory file."""
        with pytest.raises((FileNotFoundError, IOError)):
            reader = LammpsTrajectoryReader("nonexistent.dump")
            reader.read_frame(0)

    def test_read_empty_file(self):
        """Test reading empty trajectory file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            tmp.write("")  # Empty file
            
            reader = LammpsTrajectoryReader(tmp.name)
            # Should handle empty file gracefully
            with pytest.raises((IndexError, ValueError, EOFError)):
                reader.read_frame(0)

    def test_malformed_trajectory(self):
        """Test handling malformed trajectory files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
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
        reader = LammpsTrajectoryReader(TEST_DATA_DIR / "trajectory/lammps/unwrapped.lammpstrj")
        
        # Test negative index
        with pytest.raises((IndexError, ValueError)):
            reader.read_frame(-1)
        
        # Test too large index
        with pytest.raises((IndexError, ValueError)):
            reader.read_frame(999999)


class TestTrajectoryIntegration:

    def test_data_to_trajectory_conversion(self):
        """Test converting data format to trajectory format."""
        # Create a frame in data format
        frame = mp.Frame()
        
        atoms_data = {
            'id': [0, 1, 2],
            'molid': [1, 1, 2],
            'type': ['O', 'H', 'H'],
            'q': [-0.8476, 0.4238, 0.4238],
            'xyz': [[0.0, 0.0, 0.0], [0.816, 0.577, 0.0], [-0.816, 0.577, 0.0]]
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 0
        frame.box = mp.Box(np.diag([10.0, 10.0, 10.0]))
        
        # Write as trajectory
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame)
            writer.close()
            
            # Read back as trajectory
            reader = LammpsTrajectoryReader(tmp.name)
            frame_read = reader.read_frame(0)
            
            assert frame_read.metadata["timestep"] is not None
            assert "atoms" in frame_read
            assert frame_read.box is not None

    def test_multiple_formats_consistency(self):
        """Test that data and trajectory formats are consistent."""
        # This test ensures that a frame written in one format
        # can be meaningfully compared with the other format
        frame_original = mp.Frame()
        
        atoms_data = {
            'id': [0, 1],
            'type': [1, 2],
            'xyz': [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]
        }
        frame_original["atoms"] = atoms_data
        frame_original.metadata["timestep"] = 100
        frame_original.box = mp.Box(np.eye(3) * 5.0)
        
        # Write as trajectory and read back
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame_original)
            writer.close()
            
            reader = LammpsTrajectoryReader(tmp.name)
            frame_traj = reader.read_frame(0)
            
            # Both should have same basic structure
            assert frame_traj.metadata["timestep"] is not None
            assert "atoms" in frame_traj
            assert frame_traj.box is not None
            assert frame_original.box is not None
            
            # Box dimensions should be similar
            assert np.allclose(frame_traj.box.matrix.diagonal(), 
                             frame_original.box.matrix.diagonal())

# ===== Merged tests from test_lammps_trajectory_fixes.py =====

class TestLammpsTrajectoryFixes:
    """Test LAMMPS trajectory functionality with fixes."""

    def test_basic_trajectory_write_read(self):
        """Test basic trajectory writing and reading."""
        # Create frames
        frames = []
        for timestep in [0, 100, 200]:
            frame = mp.Frame()
            
            atoms_data = {
                'id': [0, 1],
                'type': [1, 2],  # Use numeric types for trajectory
                'x': [0.0 + timestep*0.01, 1.0 + timestep*0.01],
                'y': [0.0, 0.0],
                'z': [0.0, 0.0],
                'q': [0.5, -0.5]
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = timestep
            frame.box = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            # Write trajectory
            writer = LammpsTrajectoryWriter(tmp.name)
            for frame in frames:
                timestep = frame.metadata["timestep"]
                writer.write_frame(frame, timestep=timestep)
            writer.close()
            
            # Verify file exists and has content
            assert Path(tmp.name).exists()
            assert Path(tmp.name).stat().st_size > 0
            
            # Read back and verify
            reader = LammpsTrajectoryReader(tmp.name)
            
            # Test reading first frame
            frame_read = reader.read_frame(0)
            assert frame_read.metadata["timestep"] == 0
            assert frame_read["atoms"].nrows == 2
            
            # Verify coordinates are preserved
            x_values = frame_read["atoms"]["x"]
            assert np.isclose(x_values[0], 0.0)
            assert np.isclose(x_values[1], 1.0)

    def test_trajectory_timestep_handling(self):
        """Test that timesteps are correctly handled in trajectory."""
        frames = []
        timesteps = [0, 25, 50, 75, 100]
        
        for timestep in timesteps:
            frame = mp.Frame()
            
            atoms_data = {
                'id': [0, 1, 2],
                'type': ['A', 'B', 'C'],
                'x': [0.0 + timestep*0.001, 1.0 + timestep*0.001, 2.0 + timestep*0.001],
                'y': [0.0, 1.0, 2.0],
                'z': [0.0, 0.0, 0.0]
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = timestep
            frame.box = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            
            for frame in frames:
                timestep = frame.metadata["timestep"]
                writer.write_frame(frame, timestep=timestep)
            
            writer.close()
            
            # Read back and verify timesteps
            reader = LammpsTrajectoryReader(tmp.name)
            
            for i, expected_timestep in enumerate(timesteps):
                frame_read = reader.read_frame(i)
                assert frame_read.metadata["timestep"] == expected_timestep

    def test_multiple_frame_trajectory(self):
        """Test writing and reading multiple frames in a trajectory."""
        n_frames = 5
        frames = []
        
        for i in range(n_frames):
            frame = mp.Frame()
            
            atoms_data = {
                'id': [0, 1, 2],
                'type': [1, 1, 2],
                'x': [0.0 + i*0.1, 1.0 + i*0.1, 0.5 + i*0.1],
                'y': [0.0, 0.0, 1.0],
                'z': [0.0, 0.0, 0.0]
            }
            frame["atoms"] = atoms_data
            frame.metadata["timestep"] = i * 10
            frame.box = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            
            for frame in frames:
                timestep = frame.metadata["timestep"]
                writer.write_frame(frame, timestep=timestep)
            
            writer.close()
            
            # Verify file content
            with open(tmp.name, 'r') as f:
                content = f.read()
            
            # Should contain all timesteps
            assert "ITEM: TIMESTEP" in content
            assert "0\n" in content  # First timestep
            assert "40\n" in content  # Last timestep
            
            # Should contain atom data for all frames
            assert content.count("ITEM: ATOMS") == 5

    def test_trajectory_box_handling(self):
        """Test that box information is correctly written to trajectory."""
        frame = mp.Frame()
        
        atoms_data = {
            'id': [0],
            'type': [1],
            'x': [0.0],
            'y': [0.0],
            'z': [0.0]
        }
        frame["atoms"] = atoms_data
        frame.metadata["timestep"] = 0
        
        # Test with custom box dimensions
        frame.box = mp.Box(np.diag([5.0, 7.5, 10.0]))
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.dump', delete=False) as tmp:
            writer = LammpsTrajectoryWriter(tmp.name)
            writer.write_frame(frame, timestep=0)
            writer.close()
            
            with open(tmp.name, 'r') as f:
                content = f.read()
            
            # Should contain box bounds
            assert "ITEM: BOX BOUNDS" in content
            assert "5.000000" in content  # x dimension
            assert "7.500000" in content  # y dimension
            assert "10.000000" in content  # z dimension
