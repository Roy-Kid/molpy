"""Tests for GROMACS topology (.top) data file reader and writer."""

from pathlib import Path

import numpy as np
import pytest

import molrs
from molrs import MetaValue

import molpy as mp
from molpy.io.data.top import TopReader, TopWriter


@pytest.fixture
def top_dir(TEST_DATA_DIR: Path) -> Path:
    return TEST_DATA_DIR / "top"


class TestTopReader:
    """TopReader parses a GROMACS topology into per-section blocks.

    Fixtures: benzene.top (12 atoms, 12 bonds, an #include the reader must
    skip) and chain.top (four atoms with bonds, pairs, angles, dihedrals).
    """

    def test_atoms_section_columns_and_values(self, top_dir: Path) -> None:
        atoms = TopReader(top_dir / "benzene.top").read()["atoms"]
        assert atoms.nrows == 12
        for key in ("id", "type", "charge", "mass", "name"):
            assert key in atoms
        first = atoms[0]
        assert int(first["id"]) == 1
        assert str(first["type"]) == "opls_145"
        assert str(first["name"]) == "C"
        assert float(first["charge"]) == pytest.approx(-0.115)
        assert float(first["mass"]) == pytest.approx(12.011)
        assert list(atoms["type"][6:]) == ["opls_146"] * 6

    def test_bond_indices_stay_one_based(self, top_dir: Path) -> None:
        bonds = TopReader(top_dir / "benzene.top").read()["bonds"]
        assert bonds.nrows == 12
        assert int(bonds[0]["atomi"]) == 1
        assert int(bonds[0]["atomj"]) == 2
        assert bonds["atomi"].min() == 1

    def test_every_bonded_section_is_read(self, top_dir: Path) -> None:
        frame = TopReader(top_dir / "chain.top").read()
        assert frame["atoms"].nrows == 4
        assert frame["bonds"].nrows == 3
        assert frame["pairs"].nrows == 1
        assert frame["angles"].nrows == 2
        assert frame["dihedrals"].nrows == 1
        dihedral = frame["dihedrals"][0]
        assert [int(dihedral[k]) for k in ("atomi", "atomj", "atomk", "atoml")] == [
            1,
            2,
            3,
            4,
        ]

    def test_section_headers_without_spaces_are_accepted(self, tmp_path: Path) -> None:
        top_file = tmp_path / "nospaces.top"
        top_file.write_text(
            "[moleculetype]\nMOL  3\n\n[atoms]\n1  CT  1  MOL  C  1  -0.1  12.011\n"
        )
        assert TopReader(top_file).read()["atoms"].nrows == 1

    def test_empty_frame_when_no_sections(self, tmp_path: Path) -> None:
        top_file = tmp_path / "empty.top"
        top_file.write_text("; just a comment\n")
        assert "atoms" not in TopReader(top_file).read()


class TestTopWriter:
    """Tests for TopWriter producing valid GROMACS topology files."""

    def _make_minimal_frame(self) -> molrs.Frame:
        """Create a minimal two-atom frame with one bond."""
        frame = molrs.Frame()
        frame.meta = {"name": MetaValue("string", "MOL")}
        frame["atoms"] = {
            "id": np.array([1, 2]),
            "type": np.array(["CT", "HC"]),
            "resnr": np.array([1, 1]),
            "residu": np.array(["MOL", "MOL"]),
            "name": np.array(["C", "H"]),
            "cgnr": np.array([1, 2]),
            "charge": np.array([-0.1, 0.1]),
            "mass": np.array([12.011, 1.008]),
        }
        frame["bonds"] = {
            "atomi": np.array([1]),
            "atomj": np.array([2]),
            "type_id": np.array([1]),
        }
        return frame

    def test_write_creates_file(self, tmp_path: Path) -> None:
        """TopWriter.write() creates a file at the given path."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "out.top"
        writer = TopWriter(out_file)
        writer.write(frame)
        assert out_file.exists()

    def test_write_contains_sections(self, tmp_path: Path) -> None:
        """Written file contains expected GROMACS section headers."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ moleculetype ]" in content
        assert "[ atoms ]" in content
        assert "[ bonds ]" in content
        assert "[ system ]" in content
        assert "[ molecules ]" in content

    def test_write_molecule_name(self, tmp_path: Path) -> None:
        """Written file uses frame.meta['name'] as molecule name."""
        frame = self._make_minimal_frame()
        frame.meta = {**frame.meta, "name": MetaValue("string", "BENZENE")}
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "BENZENE" in content

    def test_roundtrip_atoms(self, tmp_path: Path) -> None:
        """Atoms written by TopWriter can be read back by TopReader."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "roundtrip.top"
        TopWriter(out_file).write(frame)

        frame2 = TopReader(out_file).read()
        assert "atoms" in frame2
        assert frame2["atoms"].nrows == 2

        # Check first atom
        a0 = frame2["atoms"][0]
        assert int(a0["id"]) == 1
        assert str(a0["type"]) == "CT"
        assert pytest.approx(float(a0["charge"]), abs=1e-4) == -0.1
        assert pytest.approx(float(a0["mass"]), abs=1e-3) == 12.011

    def test_roundtrip_bonds(self, tmp_path: Path) -> None:
        """Bonds written by TopWriter can be read back by TopReader."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "roundtrip.top"
        TopWriter(out_file).write(frame)

        frame2 = TopReader(out_file).read()
        assert "bonds" in frame2
        assert frame2["bonds"].nrows == 1
        bond = frame2["bonds"][0]
        assert int(bond["atomi"]) == 1
        assert int(bond["atomj"]) == 2

    def test_write_pairs_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ pairs ] section when present in frame."""
        frame = self._make_minimal_frame()
        frame["pairs"] = {
            "atomi": np.array([1]),
            "atomj": np.array([2]),
            "type_id": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ pairs ]" in content

    def test_write_angles_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ angles ] section when present in frame."""
        frame = self._make_minimal_frame()
        frame["angles"] = {
            "atomi": np.array([1]),
            "atomj": np.array([2]),
            "atomk": np.array([3]),
            "type_id": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ angles ]" in content

    def test_write_dihedrals_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ dihedrals ] section when present in frame."""
        frame = self._make_minimal_frame()
        frame["dihedrals"] = {
            "atomi": np.array([1]),
            "atomj": np.array([2]),
            "atomk": np.array([3]),
            "atoml": np.array([4]),
            "type_id": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ dihedrals ]" in content

    def test_write_via_factory(self, tmp_path: Path) -> None:
        """write_top factory function writes topology correctly."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "factory.top"
        mp.io.write_top(str(out_file), frame)
        assert out_file.exists()
        content = out_file.read_text()
        assert "[ atoms ]" in content
