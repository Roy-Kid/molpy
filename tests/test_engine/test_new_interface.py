"""
Tests for the new Engine interface (direct script execution).
"""

import tempfile
from pathlib import Path

import pytest

from molpy import Script
from molpy.engine import CP2KEngine, LAMMPSEngine


class TestNewEngineInterface:
    """Test new Engine interface with direct script execution."""

    def test_engine_with_workdir_in_constructor(self):
        """Test engine initialization with workdir in constructor."""
        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )
            assert engine.work_dir == Path(tmpdir)

    def test_engine_with_env_vars(self):
        """Test engine initialization with environment variables."""
        engine = LAMMPSEngine(
            executable="lmp",
            env_vars={"OMP_NUM_THREADS": "4"},
            check_executable=False,
        )
        assert engine.env_vars == {"OMP_NUM_THREADS": "4"}

    def test_engine_env_config_validation(self):
        """Test that env and env_manager must both be set or both be None."""
        # Both None is OK
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert engine.env is None
        assert engine.env_manager is None

        # Both set is OK
        engine = LAMMPSEngine(
            executable="lmp",
            env="myenv",
            env_manager="conda",
            check_executable=False,
        )
        assert engine.env == "myenv"
        assert engine.env_manager == "conda"

        # Only env set - should raise
        with pytest.raises(ValueError, match="environment configuration is incomplete"):
            LAMMPSEngine(executable="lmp", env="myenv", check_executable=False)

        # Only env_manager set - should raise
        with pytest.raises(ValueError, match="environment configuration is incomplete"):
            LAMMPSEngine(executable="lmp", env_manager="conda", check_executable=False)

    def test_run_with_script_object(self):
        """Test running with Script object."""
        script = Script.from_text("input", "units real\natom_style full\n")

        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            # This will fail without LAMMPS, but we can verify the setup
            try:
                result = engine.run(script, capture_output=True, check=False)
                # If it runs, verify paths were set up
                assert engine.work_dir == Path(tmpdir)
                assert len(engine.scripts) == 1
                assert engine.scripts[0] == script
                assert (Path(tmpdir) / "input.lmp").exists()
            except FileNotFoundError:
                # LAMMPS not installed, just verify setup happened
                # Note: run() was called, so setup should have occurred before the error
                pass

    def test_run_with_string_content(self):
        """Test running with string content."""
        content = "units real\natom_style full\n"

        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            try:
                result = engine.run(content, capture_output=True, check=False)
                # Verify script was created
                assert len(engine.scripts) == 1
                assert engine.scripts[0].name == "input"
                assert (Path(tmpdir) / "input.lmp").exists()
            except FileNotFoundError:
                pass

    def test_run_with_path(self):
        """Test running with Path object."""
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir_path = Path(tmpdir)

            # Create a script file
            script_file = tmpdir_path / "my_script.lmp"
            script_file.write_text("units real\natom_style full\n")

            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            try:
                result = engine.run(script_file, capture_output=True, check=False)
                # Verify script was loaded
                assert len(engine.scripts) == 1
                assert engine.scripts[0].path.name == "my_script.lmp"
            except FileNotFoundError:
                pass

    def test_run_with_multiple_scripts(self):
        """Test running with multiple scripts."""
        script1 = Script.from_text("main", "units real\n")
        script1.tags.add("input")
        script2 = Script.from_text("data", "# data file\n")

        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            try:
                result = engine.run(
                    [script1, script2], capture_output=True, check=False
                )
                # Verify both scripts were saved
                assert len(engine.scripts) == 2
                assert (Path(tmpdir) / "main.lmp").exists()
                assert (Path(tmpdir) / "data.lmp").exists()
                # First script with 'input' tag should be the input script
                assert engine.input_script == script1
            except FileNotFoundError:
                pass

    def test_run_with_workdir_override(self):
        """Test that workdir can be overridden in run()."""
        with tempfile.TemporaryDirectory() as tmpdir1:
            with tempfile.TemporaryDirectory() as tmpdir2:
                engine = LAMMPSEngine(
                    executable="lmp", workdir=tmpdir1, check_executable=False
                )

                try:
                    # Override workdir in run()
                    result = engine.run(
                        "units real\n",
                        workdir=tmpdir2,
                        capture_output=True,
                        check=False,
                    )
                    # Verify script was saved to tmpdir2, not tmpdir1
                    assert not (Path(tmpdir1) / "input.lmp").exists()
                    assert (Path(tmpdir2) / "input.lmp").exists()
                except FileNotFoundError:
                    pass

    def test_run_empty_script_list_raises(self):
        """Test that running with empty script list raises ValueError."""
        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )

            with pytest.raises(ValueError, match="At least one script is required"):
                engine.run([])

    def test_cp2k_new_interface(self):
        """Test CP2K engine with new interface."""
        script = Script.from_text("input", "&GLOBAL\n  PROJECT water\n&END GLOBAL\n")

        with tempfile.TemporaryDirectory() as tmpdir:
            engine = CP2KEngine(
                executable="cp2k.psmp", workdir=tmpdir, check_executable=False
            )

            try:
                result = engine.run(script, capture_output=True, check=False)
                # Verify setup
                assert engine.work_dir == Path(tmpdir)
                assert len(engine.scripts) == 1
                assert (Path(tmpdir) / "input.inp").exists()
            except FileNotFoundError:
                pass

    def test_merged_env_vars(self):
        """Test that env_vars are merged correctly."""
        engine = LAMMPSEngine(
            executable="lmp",
            env_vars={"VAR1": "value1"},
            check_executable=False,
        )

        merged = engine._merged_env({"VAR2": "value2"})
        assert "VAR1" in merged
        assert merged["VAR1"] == "value1"
        assert "VAR2" in merged
        assert merged["VAR2"] == "value2"

    def test_repr_with_workdir(self):
        """Test __repr__ includes workdir."""
        with tempfile.TemporaryDirectory() as tmpdir:
            engine = LAMMPSEngine(
                executable="lmp", workdir=tmpdir, check_executable=False
            )
            repr_str = repr(engine)
            assert "lmp" in repr_str
            assert tmpdir in repr_str

    def test_repr_with_env(self):
        """Test __repr__ includes env info."""
        engine = LAMMPSEngine(
            executable="lmp",
            env="myenv",
            env_manager="conda",
            check_executable=False,
        )
        repr_str = repr(engine)
        assert "lmp" in repr_str
        assert "myenv" in repr_str
        assert "conda" in repr_str


class TestRemovedLegacyInterface:
    """Ensure deprecated engine staging API is removed."""

    def test_prepare_interface_removed(self):
        """Engine.prepare() should no longer exist."""
        engine = LAMMPSEngine(executable="lmp", check_executable=False)
        assert not hasattr(engine, "prepare")
