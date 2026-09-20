import contextlib
from collections.abc import Callable, Iterator
from pathlib import Path

import mollog
import pytest


class _RecordingHandler(mollog.Handler):
    """mollog handler that collects records for assertions in tests."""

    def __init__(self) -> None:
        super().__init__()
        self.records: list[mollog.LogRecord] = []

    def emit(self, record: mollog.LogRecord) -> None:
        self.records.append(record)


@pytest.fixture
def mollog_capture() -> Callable[[str], "contextlib.AbstractContextManager"]:
    """Factory fixture to capture records on a named mollog logger.

    MolPy logs through :mod:`mollog`, not the stdlib ``logging`` module, so
    pytest's ``caplog`` does not see its records. Usage::

        with mollog_capture("molpy.io.data.lammps_bond_react") as records:
            do_something()
        assert any(r.level >= mollog.Level.WARNING for r in records)
    """

    @contextlib.contextmanager
    def _capture(name: str) -> Iterator[list[mollog.LogRecord]]:
        logger = mollog.get_logger(name)
        handler = _RecordingHandler()
        logger.add_handler(handler)
        try:
            yield handler.records
        finally:
            logger.remove_handler(handler)

    return _capture


@pytest.fixture(scope="session", name="TEST_DATA_DIR")
def test_data_dir() -> Path:
    """Fixture files committed under tests/tests-data."""
    return Path(__file__).resolve().parent / "tests-data"
