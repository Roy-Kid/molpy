"""XSF (XCrySDen Structure File) I/O — native-backed.

Read/write go through the native ``read_xsf`` / the native ``write_xsf``.
Atoms carry ``atomic_number``, ``element``, and ``x``/``y``/``z``; crystal
structures attach a periodic box, molecules a free box.
"""

from __future__ import annotations

from pathlib import Path

from molrs import Frame

from .base import DataReader, DataWriter


class XsfReader(DataReader):
    """Read an XSF file into a :class:`~molpy.Frame` natively."""

    def __init__(self, file: str | Path) -> None:
        super().__init__(Path(file))
        self._file = Path(file)

    def read(self, frame: Frame | None = None) -> Frame:
        """Read the XSF path.

        Args:
            frame: Accepted for API parity; ignored (the native core always returns a
                new Frame).

        Returns:
            Frame with atoms and box.
        """
        del frame
        import molrs.io

        if not self._file.exists():
            raise FileNotFoundError(f"XSF file not found: {self._file}")
        try:
            return molrs.io.read_xsf(str(self._file))
        except OSError as exc:
            msg = str(exc)
            lower = msg.lower()
            if "no such file" in lower or "not found" in lower:
                raise FileNotFoundError(msg) from exc
            raise ValueError(msg) from exc


class XsfWriter(DataWriter):
    """Write a Frame to an XSF file natively."""

    def __init__(self, file: str | Path) -> None:
        super().__init__(Path(file))
        self._file = Path(file)

    def write(self, frame: Frame) -> None:
        """Write *frame* as XSF (``atomic_number`` + coordinates required)."""
        import molrs.io

        try:
            molrs.io.write_xsf(str(self._file), frame)
        except OSError as exc:
            raise ValueError(str(exc)) from exc
