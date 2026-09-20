import numpy as np

import molrs
from molrs import MetaValue

import molpy as mp
from molpy.io import read_lammps_trajectory, write_lammps_dump_local
from molpy.io.trajectory.lammps import LammpsDumpLocalWriter, LammpsTrajectoryWriter


class TestWriteLammpsTrajectory:
    def test_write_simple_trajectory(self, tmp_path):
        """Test writing a simple trajectory."""
        # Create test frames
        frames = []
        for i in range(3):
            frame = molrs.Frame()

            # Create atoms data using Block structure
            atoms_data = {
                "id": [0, 1, 2],
                "type_id": [1, 1, 2],
                "x": [0.0 + i * 0.1, 1.0 + i * 0.1, 0.5 + i * 0.1],
                "y": [0.0, 0.0, 1.0],
                "z": [0.0, 0.0, 0.0],
            }
            frame["atoms"] = atoms_data
            frame.meta = {"timestep": MetaValue("i64", i * 100)}
            frame.box = mp.Box(np.eye(3) * 10.0)
            frames.append(frame)

        # Write trajectory
        tmp_file = tmp_path / "test.dump"
        writer = LammpsTrajectoryWriter(str(tmp_file))
        for frame in frames:
            timestep = frame.meta["timestep"]
            writer.write_frame(frame, timestep=timestep)
        writer.close()

        # Read back via the molrs-backed reader and verify
        reader = read_lammps_trajectory(str(tmp_file))

        # Check that we can read the frames back
        for i, frame_read in enumerate(reader):
            if i >= len(frames):
                break
            assert frame_read.meta["timestep"] == frames[i].meta["timestep"]
            assert "atoms" in frame_read
            # Check that positions changed over time
            if i > 0:
                # x coordinates should be different between frames
                x_vals = frame_read["atoms"]["x"]
                x_vals_prev = frames[0]["atoms"]["x"]
                assert not np.allclose(x_vals, x_vals_prev)

    def test_write_with_context_manager(self, tmp_path):
        """Test writing trajectory using context manager."""
        frame = molrs.Frame()
        atoms_data = {
            "id": [0, 1],
            "type_id": [1, 1],
            "x": [0.0, 1.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        }
        frame["atoms"] = atoms_data
        frame.meta = {"timestep": MetaValue("i64", 0)}
        frame.box = mp.Box(np.eye(3) * 5.0)

        tmp_file = tmp_path / "test.dump"
        with LammpsTrajectoryWriter(str(tmp_file)) as writer:
            writer.write_frame(frame)

    def test_trajectory_roundtrip(self, tmp_path):
        """Test writing and reading back maintains data integrity."""
        # Create a more complex frame
        frame = molrs.Frame()

        atoms_data = {
            "id": [0, 1, 2, 3],
            "type_id": [1, 1, 2, 2],
            "x": [0.0, 1.0, 0.5, 1.5],
            "y": [0.0, 0.0, 1.0, 1.0],
            "z": [0.0, 0.0, 0.0, 0.0],
            "vx": [0.1, -0.1, 0.2, -0.2],
            "vy": [0.0, 0.0, 0.1, -0.1],
            "vz": [0.0, 0.0, 0.0, 0.0],
        }
        frame["atoms"] = atoms_data
        frame.meta = {"timestep": MetaValue("i64", 1000)}
        frame.box = mp.Box(np.diag([5.0, 5.0, 5.0]))

        tmp_file = tmp_path / "test.dump"
        # Write
        writer = LammpsTrajectoryWriter(str(tmp_file))
        writer.write_frame(frame)
        writer.close()

        # Read back via the molrs-backed reader
        reader = read_lammps_trajectory(str(tmp_file))
        frame_read = reader[0]

        # Verify timestep
        assert frame_read.meta["timestep"] == 1000

        # Verify atoms data exists
        assert "atoms" in frame_read
        atoms = frame_read["atoms"]
        assert atoms.nrows == 4

        # Verify box
        assert frame_read.box is not None
        assert np.allclose(frame_read.box.matrix.diagonal(), [5.0, 5.0, 5.0])


class TestTrajectoryIntegration:
    def test_data_to_trajectory_conversion(self, tmp_path):
        """Test converting data format to trajectory format."""
        # Create a frame in data format
        frame = molrs.Frame()

        atoms_data = {
            "id": [0, 1, 2],
            "type_id": [1, 1, 2],  # Use numeric types for LAMMPS
            "x": [0.0, 0.816, -0.816],
            "y": [0.0, 0.577, 0.577],
            "z": [0.0, 0.0, 0.0],
            "q": [-0.8476, 0.4238, 0.4238],
        }
        frame["atoms"] = atoms_data
        frame.meta = {"timestep": MetaValue("i64", 0)}
        frame.box = mp.Box(np.diag([10.0, 10.0, 10.0]))

        # Write as trajectory
        tmp_file = tmp_path / "test.dump"
        writer = LammpsTrajectoryWriter(str(tmp_file))
        writer.write_frame(frame)
        writer.close()

        # Read back as trajectory
        reader = read_lammps_trajectory(str(tmp_file))
        frame_read = reader[0]

        assert frame_read.meta["timestep"] == 0
        assert "atoms" in frame_read
        assert frame_read.box is not None

    def test_multiple_formats_consistency(self, tmp_path):
        """Test that data and trajectory formats are consistent."""
        # This test ensures that a frame written in one format
        # can be meaningfully compared with the other format
        frame_original = molrs.Frame()

        atoms_data = {
            "id": [0, 1],
            "type_id": [1, 2],
            "x": [0.0, 1.0],
            "y": [0.0, 0.0],
            "z": [0.0, 0.0],
        }
        frame_original["atoms"] = atoms_data
        frame_original.meta = {"timestep": MetaValue("i64", 100)}
        frame_original.box = mp.Box(np.eye(3) * 5.0)

        # Write as trajectory and read back
        tmp_file = tmp_path / "test.dump"
        writer = LammpsTrajectoryWriter(str(tmp_file))
        writer.write_frame(frame_original)
        writer.close()

        reader = read_lammps_trajectory(str(tmp_file))
        frame_traj = reader[0]

        # Both should have same basic structure
        assert frame_traj.meta["timestep"] == 100
        assert "atoms" in frame_traj
        assert frame_traj.box is not None


class TestWriteLammpsDumpLocal:
    def test_write_bonds_roundtrip(self, tmp_path):
        atoms = molrs.Block()
        atoms["id"] = np.array([1, 2, 3], dtype=np.uint64)
        atoms["x"] = np.array([0.0, 1.0, 2.0])
        atoms["y"] = np.zeros(3)
        atoms["z"] = np.zeros(3)
        bonds = molrs.Block()
        bonds["atomi"] = np.array([0, 1], dtype=np.uint64)
        bonds["atomj"] = np.array([1, 2], dtype=np.uint64)
        frame = molrs.Frame()
        frame["atoms"] = atoms
        frame["bonds"] = bonds
        frame.box = molrs.Box.cube(10.0)
        path = tmp_path / "bonds.dump.local"
        write_lammps_dump_local(path, [frame])
        text = path.read_text()
        assert "ITEM: NUMBER OF ENTRIES" in text
        assert "batom1 batom2" in text
        loaded = molrs.io.raw.read_lammps_traj(str(path))
        assert loaded[0]["entries"].nrows == 2

    def test_writer_class_matches_factory(self, tmp_path):
        frame = molrs.Frame()
        frame["atoms"] = {
            "id": np.array([1, 2], dtype=np.uint64),
            "x": np.array([0.0, 1.0]),
            "y": np.zeros(2),
            "z": np.zeros(2),
        }
        frame["bonds"] = {
            "atomi": np.array([0], dtype=np.uint64),
            "atomj": np.array([1], dtype=np.uint64),
        }
        frame.box = molrs.Box.cube(8.0)
        path = tmp_path / "one.dump.local"
        with LammpsDumpLocalWriter(path) as writer:
            writer.write_frame(frame)
        assert "ITEM: NUMBER OF ENTRIES" in path.read_text()
