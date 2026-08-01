"""Tests for GROMACS topology (.top) data file reader and writer."""

from pathlib import Path

import numpy as np
import pytest

import molrs
from molrs import MetaValue

import molpy as mp
from molpy.io.data.top import TopFormatError, TopReader, TopWriter


class TestTopReader:
    """Tests for TopReader parsing GROMACS topology files."""

    def test_read_benzene_atoms(self, TEST_DATA_DIR: Path) -> None:
        """TopReader should parse [atoms] section correctly."""
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        reader = TopReader(top_file)
        frame = reader.read()

        assert "atoms" in frame
        atoms = frame["atoms"]
        assert atoms.nrows == 12  # 6 C + 6 H in benzene

    def test_read_benzene_bonds(self, TEST_DATA_DIR: Path) -> None:
        """TopReader should parse [bonds] section correctly."""
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        reader = TopReader(top_file)
        frame = reader.read()

        assert "bonds" in frame
        bonds = frame["bonds"]
        assert bonds.nrows == 12

    def test_atom_fields(self, TEST_DATA_DIR: Path) -> None:
        """Atom block should contain expected fields."""
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        frame = TopReader(top_file).read()
        atoms = frame["atoms"]

        assert "id" in atoms
        assert "type" in atoms
        assert "charge" in atoms
        assert "mass" in atoms
        assert "name" in atoms

    def test_first_atom_values(self, TEST_DATA_DIR: Path) -> None:
        """First atom should have correct values from benzene.top."""
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        frame = TopReader(top_file).read()
        atom = frame["atoms"][0]

        assert int(atom["id"]) == 1
        assert str(atom["type"]) == "opls_145"
        assert pytest.approx(float(atom["charge"]), abs=1e-4) == -0.115
        assert pytest.approx(float(atom["mass"]), abs=1e-3) == 12.011

    def test_bond_indices_are_zero_based_row_indices(self, TEST_DATA_DIR: Path) -> None:
        """Endpoints are 0-based row indices, not the file's 1-based serials.

        This test previously asserted the opposite (``atomi == 1`` for the
        record ``1 2 1``), pinning the off-by-one as if it were the contract.
        """
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        frame = TopReader(top_file).read()
        first_bond = frame["bonds"][0]

        # First bond record in benzene.top is "1 2 1" (serials 1 and 2).
        assert int(first_bond["atomi"]) == 0
        assert int(first_bond["atomj"]) == 1

    def test_all_endpoints_index_within_atoms_block(self, TEST_DATA_DIR: Path) -> None:
        """Every endpoint in every relation block addresses a real atom row.

        Under the old 1-based storage the *last* atom's index equalled
        ``n_atoms`` and was silently out of range. Benzene's ring bonds do
        reach atom 12 of 12, so this fails loudly on the unfixed reader.
        """
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        frame = TopReader(top_file).read()
        n_atoms = frame["atoms"].nrows

        for key in ("bonds", "pairs", "angles", "dihedrals"):
            if key not in frame:
                continue
            block = frame[key]
            for col in ("atomi", "atomj", "atomk", "atoml"):
                if col not in block:
                    continue
                values = np.asarray(block[col])
                assert values.min() >= 0, f"{key}.{col} has a negative index"
                assert values.max() < n_atoms, (
                    f"{key}.{col} indexes atom {values.max()} but the frame "
                    f"only has {n_atoms} atoms"
                )

    def test_endpoint_corpus_touches_first_and_last_atom(
        self, TEST_DATA_DIR: Path
    ) -> None:
        """Meta-test: the fixture must exercise both ends of the index range.

        A 1-based index set is ``[1..N]``, so only the value ``N`` is out of
        range. Without a record touching the final atom, the test above would
        pass on the broken reader and prove nothing.
        """
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        frame = TopReader(top_file).read()
        n_atoms = frame["atoms"].nrows
        endpoints = np.concatenate(
            [
                np.asarray(frame["bonds"][c]).ravel()
                for c in ("atomi", "atomj")
                if c in frame["bonds"]
            ]
        )
        assert endpoints.min() == 0, "corpus never references the first atom"
        assert endpoints.max() == n_atoms - 1, (
            "corpus never references the last atom, so it cannot detect an "
            "off-by-one in endpoint indexing"
        )

    def test_zero_serial_is_rejected_not_dropped(self, tmp_path: Path) -> None:
        """A 0 serial is invalid GROMACS and must raise, not vanish.

        ``0 - 1 == -1`` would wrap to the last atom; returning ``None`` would
        silently drop the bond. Both are worse than failing.
        """
        top_file = tmp_path / "zero_serial.top"
        top_file.write_text(
            "[moleculetype]\nMOL  3\n\n"
            "[atoms]\n1  CT  1  MOL  C  1  -0.1  12.011\n"
            "2  HC  1  MOL  H  2  0.1  1.008\n\n"
            "[bonds]\n0  2  1\n"
        )

        with pytest.raises(TopFormatError, match="1-based"):
            TopReader(top_file).read()

    def test_read_bromobutane_all_sections(self, TEST_DATA_DIR: Path) -> None:
        """1-bromobutane has atoms, bonds, pairs, angles, and dihedrals."""
        top_file = TEST_DATA_DIR / "top/1-bromobutane.top"
        if not top_file.exists():
            pytest.skip("1-bromobutane.top test data not available")

        frame = TopReader(top_file).read()

        assert "atoms" in frame
        assert "bonds" in frame
        assert "pairs" in frame
        assert "angles" in frame
        assert "dihedrals" in frame

    def test_bromobutane_atom_count(self, TEST_DATA_DIR: Path) -> None:
        """1-bromobutane has 14 atoms."""
        top_file = TEST_DATA_DIR / "top/1-bromobutane.top"
        if not top_file.exists():
            pytest.skip("1-bromobutane.top test data not available")

        frame = TopReader(top_file).read()
        assert frame["atoms"].nrows == 14

    def test_bromobutane_bond_count(self, TEST_DATA_DIR: Path) -> None:
        """1-bromobutane has 13 bonds."""
        top_file = TEST_DATA_DIR / "top/1-bromobutane.top"
        if not top_file.exists():
            pytest.skip("1-bromobutane.top test data not available")

        frame = TopReader(top_file).read()
        assert frame["bonds"].nrows == 13

    def test_section_normalization_spaces(self, tmp_path: Path) -> None:
        """Reader handles both `[ atoms ]` and `[atoms]` section styles."""
        # Write file with spacing-less headers
        top_content = """; test
[moleculetype]
MOL  3

[atoms]
1  CT  1  MOL  C  1  -0.1  12.011

[bonds]
"""
        top_file = tmp_path / "test_nospaces.top"
        top_file.write_text(top_content)

        frame = TopReader(top_file).read()
        assert "atoms" in frame
        assert frame["atoms"].nrows == 1

    def test_empty_frame_when_no_sections(self, tmp_path: Path) -> None:
        """Reader returns an empty frame for a file with no known sections."""
        top_file = tmp_path / "empty.top"
        top_file.write_text("; just a comment\n")
        frame = TopReader(top_file).read()
        assert "atoms" not in frame

    def test_via_factory_function(self, TEST_DATA_DIR: Path) -> None:
        """mp.io.read_top with a frame argument reads topology data."""
        top_file = TEST_DATA_DIR / "top/benzene.top"
        if not top_file.exists():
            pytest.skip("benzene.top test data not available")

        # read_top reads into a ForceField (forcefield reader)
        # For data, use TopReader directly
        frame = TopReader(top_file).read()
        assert "atoms" in frame


class TestTopWriter:
    """Tests for TopWriter producing valid GROMACS topology files."""

    def _make_minimal_frame(self, n_atoms: int = 2) -> molrs.Frame:
        """Create a small frame with one bond between the first two atoms.

        Endpoints are 0-based row indices. The previous fixture used ``atomi=1,
        atomj=2`` on a two-atom frame — index 2 does not exist, so the fixture
        itself carried the off-by-one it was meant to test.
        """
        frame = molrs.Frame()
        frame.meta = {"name": MetaValue("string", "MOL")}
        idx = np.arange(n_atoms)
        frame["atoms"] = {
            "id": idx + 1,
            "type": np.array(["CT", "HC"] * n_atoms)[:n_atoms],
            "resnr": np.ones(n_atoms, dtype=int),
            "residu": np.array(["MOL"] * n_atoms),
            "name": np.array(["C", "H"] * n_atoms)[:n_atoms],
            "cgnr": idx + 1,
            "charge": np.linspace(-0.1, 0.1, n_atoms),
            "mass": np.array([12.011, 1.008] * n_atoms)[:n_atoms],
        }
        frame["bonds"] = {
            "atomi": np.array([0]),
            "atomj": np.array([1]),
            "type": np.array([1]),
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
        # Written as GROMACS serials 1/2, read back as row indices 0/1.
        assert int(bond["atomi"]) == 0
        assert int(bond["atomj"]) == 1

    def test_writer_emits_one_based_serials(self, tmp_path: Path) -> None:
        """The [ bonds ] record on disk must carry GROMACS 1-based serials."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        section = out_file.read_text().split("[ bonds ]")[1]
        record = [
            ln.strip()
            for ln in section.splitlines()
            if ln.strip() and not ln.strip().startswith(";")
        ][0]
        # Frame endpoints are rows 0 and 1 → serials 1 and 2 on disk.
        assert record.split()[:2] == ["1", "2"]

    def test_write_pairs_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ pairs ] section when present in frame."""
        frame = self._make_minimal_frame()
        frame["pairs"] = {
            "atomi": np.array([0]),
            "atomj": np.array([1]),
            "type": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ pairs ]" in content

    def test_write_angles_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ angles ] section when present in frame."""
        frame = self._make_minimal_frame(n_atoms=3)
        frame["angles"] = {
            "atomi": np.array([0]),
            "atomj": np.array([1]),
            "atomk": np.array([2]),
            "type": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ angles ]" in content

    def test_write_dihedrals_section(self, tmp_path: Path) -> None:
        """TopWriter writes [ dihedrals ] section when present in frame."""
        frame = self._make_minimal_frame(n_atoms=4)
        frame["dihedrals"] = {
            "atomi": np.array([0]),
            "atomj": np.array([1]),
            "atomk": np.array([2]),
            "atoml": np.array([3]),
            "type": np.array([1]),
        }
        out_file = tmp_path / "out.top"
        TopWriter(out_file).write(frame)

        content = out_file.read_text()
        assert "[ dihedrals ]" in content

    @pytest.mark.parametrize("missing", ["type", "charge", "mass"])
    def test_write_rejects_missing_forcefield_column(
        self, tmp_path: Path, missing: str
    ) -> None:
        """Force-field data is never fabricated: the write fails instead.

        The writer used to substitute ``type="X"``, ``charge=0.0`` and
        ``mass=0.0``, producing a topology that loads cleanly and is physically
        meaningless.
        """
        frame = self._make_minimal_frame()
        atoms = {k: np.asarray(frame["atoms"][k]) for k in frame["atoms"].keys()}
        del atoms[missing]
        frame["atoms"] = atoms

        with pytest.raises(ValueError, match=missing):
            TopWriter(tmp_path / "out.top").write(frame)

    def test_write_via_factory(self, tmp_path: Path) -> None:
        """write_top factory function writes topology correctly."""
        frame = self._make_minimal_frame()
        out_file = tmp_path / "factory.top"
        mp.io.write_top(str(out_file), frame)
        assert out_file.exists()
        content = out_file.read_text()
        assert "[ atoms ]" in content
