"""AMBER frcmod I/O seam: missing files and the section round trip."""

import pytest

from molpy.io.forcefield.frcmod import read_frcmod, write_frcmod


def test_missing_file_raises(tmp_path):
    with pytest.raises(FileNotFoundError):
        read_frcmod(tmp_path / "nope.frcmod")


def test_sections_round_trip(tmp_path):
    path = tmp_path / "deep" / "mol.frcmod"
    write_frcmod(path, {"remark": "hand-written test", "mass": "CT 12.011 0.878"})
    sections = read_frcmod(path)
    assert "hand-written test" in sections["raw_text"]
    assert "CT" in sections["raw_text"]
    assert set(sections) >= {
        "remark",
        "mass",
        "bond",
        "angle",
        "dihe",
        "improper",
        "nonbon",
    }
