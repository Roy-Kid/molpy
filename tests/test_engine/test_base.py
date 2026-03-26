"""Unit tests for engine base classes."""

import tempfile
from pathlib import Path

import pytest

from molpy import Script
from molpy.engine import CP2KEngine, LAMMPSEngine


class TestEngineInit:
    """Test engine initialization."""

    def test_init_with_defaults(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert engine.executable == "lmp"
        assert engine.work_dir is None
        assert engine.env_vars == {}
        assert engine.env is None
        assert engine.env_manager is None

    def test_init_with_workdir(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )
            assert engine.work_dir == Path(tmpdir)

    def test_init_with_env_vars(self):
        engine = LAMMPSEngine(
            executable="lmp",
            env_vars={"OMP_NUM_THREADS": "4"},
            check_executable=False,
        )
        assert engine.env_vars == {"OMP_NUM_THREADS": "4"}

    def test_init_env_validation(self):
        # Both None OK
        LAMMPSEngine(executable="lmp", check_executable=False)

        # Both set OK
        LAMMPSEngine(
            executable="lmp", env="myenv", env_manager="conda", check_executable=False
        )

        # Only env set -> raises
        with pytest.raises(ValueError, match="environment configuration is incomplete"):
            LAMMPSEngine(executable="lmp", env="myenv", check_executable=False)

        # Only env_manager set -> raises
        with pytest.raises(ValueError, match="environment configuration is incomplete"):
            LAMMPSEngine(executable="lmp", env_manager="conda", check_executable=False)

    def test_check_executable_missing(self):
        with pytest.raises(FileNotFoundError):
            LAMMPSEngine(executable="nonexistent_lammps_binary_xyz123")

    def test_repr(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert "lmp" in repr(engine)
        assert "LAMMPSEngine" in repr(engine)


class TestEngineRun:
    """Test engine run method."""

    def test_run_no_scripts_raises(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        with pytest.raises(ValueError, match="At least one script is required"):
            engine.run()

    def test_run_empty_list_raises(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        with pytest.raises(ValueError, match="At least one script is required"):
            engine.run([])

    def test_run_with_script_saves_files(self):
        script = Script.from_text("input", "units real\natom_style full\n")

        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            try:
                engine.run(script, capture_output=True, check=False)
            except FileNotFoundError:
                pass  # LAMMPS not installed

            # Script should have been saved
            assert (Path(tmpdir) / "input.lmp").exists()

    def test_run_with_string(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            try:
                engine.run("units real\n", capture_output=True, check=False)
            except FileNotFoundError:
                pass

            assert (Path(tmpdir) / "input.lmp").exists()


class TestCP2KEngine:
    """Test CP2K engine specifics."""

    def test_name(self):
        engine = CP2KEngine(executable="cp2k", check_executable=False)
        assert engine.name == "CP2K"

    def test_extension(self):
        engine = CP2KEngine(executable="cp2k", check_executable=False)
        assert engine._get_default_extension() == ".inp"


class TestLAMMPSEngine:
    """Test LAMMPS engine specifics."""

    def test_name(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert engine.name == "LAMMPS"

    def test_extension(self):
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert engine._get_default_extension() == ".lmp"
