"""
CP2K quantum chemistry engine.

This module provides the CP2KEngine class for running CP2K calculations.
"""

import subprocess
from pathlib import Path
from typing import Any

from .base import Engine


class CP2KEngine(Engine):
    """
    CP2K quantum chemistry engine.

    This engine runs CP2K calculations with input scripts.

    Example:
        >>> from molpy.core.script import Script
        >>> from molpy.engine import CP2KEngine
        >>>
        >>> # Create input script
        >>> script = Script.from_text(
        ...     name="input",
        ...     text="&GLOBAL\\n  PROJECT water\\n&END GLOBAL\\n",
        ...     language="other"
        ... )
        >>>
        >>> # Create engine and prepare
        >>> engine = CP2KEngine(executable="cp2k.psmp")
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
        return "CP2K"

    def _get_default_extension(self) -> str:
        """Get default file extension for CP2K input files."""
        return ".inp"

    def _execute(
        self,
        capture_output: bool = False,
        check: bool = True,
        **kwargs: Any,
    ) -> subprocess.CompletedProcess:
        """Execute CP2K calculation."""
        if self.input_script is None or self.input_script.path is None:
            raise RuntimeError("No input script found.")

        input_file = self.input_script.path.name
        command = [self.executable, "-i", input_file, "-o", "cp2k.out"]

        return subprocess.run(
            command,
            cwd=self.work_dir,
            capture_output=capture_output,
            text=True,
            check=check,
            env=self._merged_env(),
        )
