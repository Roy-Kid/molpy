"""
LAMMPS molecular dynamics engine.

This module provides the LAMMPSEngine class for running LAMMPS calculations.
"""

import subprocess
from pathlib import Path
from typing import Any

from .base import Engine


class LAMMPSEngine(Engine):
    """
    LAMMPS molecular dynamics engine.

    This engine runs LAMMPS calculations with input scripts.

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

    @property
    def name(self) -> str:
        """Return engine name."""
        return "LAMMPS"

    def _get_default_extension(self) -> str:
        """Get default file extension for LAMMPS input files."""
        return ".lmp"

    def _execute(
        self,
        capture_output: bool = False,
        check: bool = True,
        **kwargs: Any,
    ) -> subprocess.CompletedProcess:
        """Execute LAMMPS calculation."""
        if self.input_script is None or self.input_script.path is None:
            raise RuntimeError("No input script found.")

        input_file = self.input_script.path.name
        command = [self.executable, "-in", input_file, "-log", "log.lammps"]

        return subprocess.run(
            command,
            cwd=self.work_dir,
            capture_output=capture_output,
            text=True,
            check=check,
            env=self._merged_env(),
        )
