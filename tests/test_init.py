"""Unit test for the ``molpy`` facade module (``src/molpy/__init__.py``)."""

from __future__ import annotations

import pytest

import molpy


def test_unregistered_name_raises_attribute_error() -> None:
    """``__getattr__`` resolves lazy submodules only; anything else fails.

    A neutral name that was never registered must raise ``AttributeError``
    rather than falling through to an import attempt.
    """
    with pytest.raises(AttributeError, match="not_a_submodule"):
        molpy.not_a_submodule
