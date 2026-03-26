"""
Engine base classes for molecular simulation engines.

This module provides abstract base classes for running external computational
chemistry programs like LAMMPS, CP2K, etc. It integrates with the core Script
class for input file management.
"""

import shutil
import subprocess
from abc import ABC, abstractmethod
from collections.abc import Sequence
from pathlib import Path
from typing import TYPE_CHECKING, Any, Self

if TYPE_CHECKING:
    pass

from molpy.core.script import Script


class Engine(ABC):
    """
    Abstract base class for computational chemistry engines.

    Provides a common interface for running external programs like LAMMPS, CP2K, etc.
    Each engine handles setup, execution, and output processing for its specific program.

    The Engine class integrates with the core Script class for input file management.
    Scripts can be created from text, loaded from files, or loaded from URLs.

    Attributes:
        executable: Path or command to the executable
        work_dir: Working directory for calculations
        scripts: List of Script objects for input files
        input_script: Primary input script (first script or script with 'input' tag)
        output_file: Primary output file from the calculation

    Example:
        >>> from molpy.core.script import Script
        >>> from molpy.engine import LAMMPSEngine
        >>>
        >>> # Create input script
        >>> script = Script.from_text(
        ...     name="input",
        ...     text="units real\\natom_style full\\n",
        ...     language="other"
        ... )
        >>>
        >>> # Create engine and prepare
        >>> engine = LAMMPSEngine(executable="lmp")
        >>> engine.run(script, workdir="./calc", check=False)
        >>>
        >>> # Run calculation
        >>> result = engine.run()
        >>> print(result.returncode)
        0
    """

    def __init__(
        self,
        executable: str,
        *,
        workdir: str | Path | None = None,
        env_vars: dict[str, str] | None = None,
        env: str | None = None,
        env_manager: str | None = None,
        check_executable: bool = True,
    ):
        """
        Initialize the engine.

        Args:
            executable: Path or command to the executable.
            workdir: Default working directory for calculations.
            env_vars: Environment variables to set during execution.
            env: Conda/virtual environment name.
            env_manager: Environment manager type (e.g. "conda").
            check_executable: Whether to check if executable exists in PATH.

        Raises:
            FileNotFoundError: If check_executable is True and executable not found.
            ValueError: If env and env_manager are not both set or both None.
        """
        self.executable = executable
        self.work_dir = Path(workdir) if workdir is not None else None
        self.env_vars = env_vars or {}
        self.env = env
        self.env_manager = env_manager

        if (env is None) != (env_manager is None):
            raise ValueError(
                "Both 'env' and 'env_manager' must be set together, or both None. "
                "Your environment configuration is incomplete."
            )

        if check_executable:
            self.check_executable()

    @property
    @abstractmethod
    def name(self) -> str:
        """Return the name of the engine."""
        pass

    def check_executable(self) -> None:
        """
        Check if the executable is available in the system PATH.

        Raises:
            FileNotFoundError: If the executable is not found
        """
        if not shutil.which(self.executable):
            raise FileNotFoundError(
                f"Executable '{self.executable}' not found in PATH. "
                f"Please install the engine or set the correct executable path."
            )

    def __repr__(self) -> str:
        parts = [f"executable='{self.executable}'"]
        if self.work_dir is not None:
            parts.append(f"workdir='{self.work_dir}'")
        if self.env is not None:
            parts.append(f"env='{self.env}'")
        if self.env_manager is not None:
            parts.append(f"env_manager='{self.env_manager}'")
        return f"<{self.__class__.__name__}({', '.join(parts)})>"

    def _merged_env(self, extra: dict[str, str] | None = None) -> dict[str, str]:
        """Merge base env_vars with optional extra vars.

        Args:
            extra: Additional environment variables to merge.

        Returns:
            Merged dict (base vars overridden by extra).
        """
        import os

        merged = dict(os.environ)
        merged.update(self.env_vars)
        if extra:
            merged.update(extra)
        return merged

    def _find_input_script(self) -> Script | None:
        """
        Find the primary input script.

        Returns:
            Script object with 'input' tag, or first script if no tag found
        """
        # Look for script with 'input' tag
        for script in self.scripts:
            if "input" in script.tags:
                return script

        # Return first script if no tag found
        return self.scripts[0] if self.scripts else None

    @abstractmethod
    def _get_default_extension(self) -> str:
        """
        Get the default file extension for input files.

        Returns:
            Default file extension (e.g., '.inp', '.lmp')
        """
        pass

    def run(
        self,
        scripts: "Script | str | Path | Sequence[Script] | None" = None,
        *,
        workdir: str | Path | None = None,
        capture_output: bool = False,
        check: bool = True,
        **kwargs: Any,
    ) -> subprocess.CompletedProcess:
        """Execute the engine calculation.

        Accepts scripts as Script objects, string content, Path to a file,
        or a list of Scripts. If workdir is given, it overrides the default.

        Args:
            scripts: Input script(s) — Script, str, Path, or list[Script].
            workdir: Override working directory for this run.
            capture_output: Capture stdout/stderr.
            check: Raise on non-zero exit code.
            **kwargs: Additional arguments for subclass run.

        Returns:
            CompletedProcess with execution results.

        Raises:
            ValueError: If no scripts provided and none previously prepared.
        """
        # Resolve workdir
        run_dir = Path(workdir) if workdir is not None else self.work_dir
        if run_dir is None:
            import tempfile

            run_dir = Path(tempfile.mkdtemp())
        run_dir.mkdir(parents=True, exist_ok=True)
        self.work_dir = run_dir

        # Normalize scripts
        if scripts is not None:
            if isinstance(scripts, str):
                scripts = [Script.from_text("input", scripts)]
            elif isinstance(scripts, Path):
                scripts = [Script.from_path(scripts)]
            elif isinstance(scripts, Script):
                scripts = [scripts]
            else:
                scripts = list(scripts)

            if not scripts:
                raise ValueError("At least one script is required")

            self.scripts = scripts
        elif not hasattr(self, "scripts") or not self.scripts:
            raise ValueError("At least one script is required")

        # Save scripts to workdir
        for script in self.scripts:
            if script.path is not None:
                script_path = run_dir / script.path.name
            else:
                ext = self._get_default_extension()
                script_path = run_dir / f"{script.name}{ext}"
            script.save(script_path)

        self.input_script = self._find_input_script()

        return self._execute(capture_output=capture_output, check=check, **kwargs)

    @abstractmethod
    def _execute(
        self,
        capture_output: bool = False,
        check: bool = True,
        **kwargs: Any,
    ) -> subprocess.CompletedProcess:
        """Run the engine subprocess. Subclasses implement this."""
        pass
