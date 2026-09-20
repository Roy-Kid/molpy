"""Antechamber .ac reader seam: error mapping and the ``q`` -> ``charge`` rename."""

import numpy as np
import pytest

import molrs

from molpy.io.data.ac import AcFieldFormatter, AcReader


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        AcReader(tmp_path / "nope.ac").read()


def test_formatter_canonicalises_the_charge_column():
    block = molrs.Block()
    block["q"] = np.array([-0.1, 0.1])
    AcFieldFormatter().canonicalize(block)
    assert "charge" in block and "q" not in block
