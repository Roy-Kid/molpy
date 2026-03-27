"""Tests for AMBER prmtop file reading."""

import numpy as np
import pytest

from molpy.core.frame import Frame
from molpy.io import read_amber
from molpy.io.forcefield.amber import AmberPrmtopReader, CHARGE_CONVERSION_FACTOR


@pytest.fixture
def litfsi_prmtop(TEST_DATA_DIR):
    """Path to LiTFSI test prmtop file."""
    return TEST_DATA_DIR / "prmtop" / "LiTFSI.prmtop"


@pytest.fixture
def litfsi_inpcrd(TEST_DATA_DIR):
    """Path to LiTFSI test inpcrd file."""
    return TEST_DATA_DIR / "inpcrd" / "LiTFSI.inpcrd"


def test_prmtop_file_exists(litfsi_prmtop):
    """Test that the test prmtop file exists."""
    assert litfsi_prmtop.exists(), f"Test file not found: {litfsi_prmtop}"


def test_prmtop_reader_initialization(litfsi_prmtop):
    """Test AmberPrmtopReader initialization."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    assert reader.file == litfsi_prmtop
    assert reader.raw_data == {}
    assert reader.meta == {}


def test_prmtop_sanitizer():
    """Test the sanitizer static method."""
    assert AmberPrmtopReader.sanitizer("  test  \n") == "test"
    assert AmberPrmtopReader.sanitizer("\t  data  ") == "data"
    assert AmberPrmtopReader.sanitizer("nowhitespace") == "nowhitespace"


def test_prmtop_read_section_int():
    """Test read_section with integer conversion."""
    lines = ["1 2 3", "4 5 6"]
    result = AmberPrmtopReader.read_section(lines, int)
    assert result == [1, 2, 3, 4, 5, 6]


def test_prmtop_read_section_float():
    """Test read_section with float conversion."""
    lines = ["1.0 2.5", "3.7 4.2"]
    result = AmberPrmtopReader.read_section(lines, float)
    assert result == [1.0, 2.5, 3.7, 4.2]


def test_prmtop_read_section_str():
    """Test read_section with string conversion."""
    lines = ["ABC DEF", "GHI"]
    result = AmberPrmtopReader.read_section(lines, str)
    assert result == ["ABC", "DEF", "GHI"]


def test_prmtop_read_basic(litfsi_prmtop):
    """Test basic reading of prmtop file."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    # Check that frame and forcefield were returned
    assert frame is not None
    assert ff is not None

    # Check that basic blocks exist
    assert "atoms" in frame
    assert "bonds" in frame
    assert "angles" in frame
    assert "dihedrals" in frame


def test_prmtop_read_pointers(litfsi_prmtop):
    """Test reading POINTERS section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    # LiTFSI has 16 atoms (15 TFSI + 1 Li)
    assert frame.metadata["n_atoms"] == 16
    assert frame["atoms"].nrows == 16

    # Check that meta contains expected fields
    assert "n_bonds" in frame.metadata
    assert "n_angles" in frame.metadata
    assert "n_dihedrals" in frame.metadata
    assert "n_atomtypes" in frame.metadata


def test_prmtop_read_atom_names(litfsi_prmtop):
    """Test reading ATOM_NAME section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "name" in atoms

    # Check some expected atom names
    names = atoms["name"]
    assert len(names) == 16
    # First few atoms should be F, C, F1, F2, S, ...
    assert names[0] == "F"
    assert names[1] == "C"
    # Last atom should be LI
    assert names[-1] == "LI"


def test_prmtop_read_charges(litfsi_prmtop):
    """Test reading and converting CHARGE section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "charge" in atoms

    charges = atoms["charge"]
    assert len(charges) == 16

    # Charges should be converted by dividing by CHARGE_CONVERSION_FACTOR
    # Last atom (Li+) should have charge close to +1
    assert np.isclose(charges[-1], 1.0, atol=0.01)

    # Check that charges sum to approximately 0 (neutral system)
    # Note: LiTFSI is Li+ with TFSI-, so should sum to 0
    total_charge = np.sum(charges)
    assert np.isclose(total_charge, 0.0, atol=0.01)


def test_prmtop_read_atomic_numbers(litfsi_prmtop):
    """Test reading ATOMIC_NUMBER section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "atomic_number" in atoms

    atomic_numbers = atoms["atomic_number"]
    assert len(atomic_numbers) == 16

    # Li has atomic number 3
    assert atomic_numbers[-1] == 3
    # F has atomic number 9
    assert atomic_numbers[0] == 9


def test_prmtop_read_masses(litfsi_prmtop):
    """Test reading MASS section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "mass" in atoms

    masses = atoms["mass"]
    assert len(masses) == 16

    # Li mass ~6.94
    assert np.isclose(masses[-1], 6.94, atol=0.1)
    # F mass ~19.0
    assert np.isclose(masses[0], 19.0, atol=0.5)


def test_prmtop_read_atom_types(litfsi_prmtop):
    """Test reading AMBER_ATOM_TYPE section."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "type" in atoms

    types = atoms["type"]
    assert len(types) == 16
    # All should be strings
    assert all(isinstance(t, str) for t in types)


def test_prmtop_read_bonds(litfsi_prmtop):
    """Test reading bond information."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    bonds = frame["bonds"]
    assert "atomi" in bonds
    assert "atomj" in bonds
    assert "type" in bonds
    assert "type_id" in bonds
    assert "id" in bonds

    n_bonds = frame.metadata["n_bonds"]
    assert len(bonds["atomi"]) == n_bonds
    assert len(bonds["atomj"]) == n_bonds

    # Bond indices should be valid (0-indexed, less than n_atoms)
    assert all(0 <= i < 16 for i in bonds["atomi"])
    assert all(0 <= j < 16 for j in bonds["atomj"])


def test_prmtop_read_angles(litfsi_prmtop):
    """Test reading angle information."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    angles = frame["angles"]
    assert "atomi" in angles
    assert "atomj" in angles
    assert "atomk" in angles
    assert "type" in angles
    assert "type_id" in angles
    assert "id" in angles

    n_angles = frame.metadata["n_angles"]
    assert len(angles["atomi"]) == n_angles
    assert len(angles["atomj"]) == n_angles
    assert len(angles["atomk"]) == n_angles

    # Angle indices should be valid
    assert all(0 <= i < 16 for i in angles["atomi"])
    assert all(0 <= j < 16 for j in angles["atomj"])
    assert all(0 <= k < 16 for k in angles["atomk"])


def test_prmtop_read_dihedrals(litfsi_prmtop):
    """Test reading dihedral information."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    dihedrals = frame["dihedrals"]
    assert "atomi" in dihedrals
    assert "atomj" in dihedrals
    assert "atomk" in dihedrals
    assert "atoml" in dihedrals
    assert "type" in dihedrals
    assert "type_id" in dihedrals
    assert "id" in dihedrals

    n_dihedrals = frame.metadata["n_dihedrals"]
    assert len(dihedrals["atomi"]) == n_dihedrals
    assert len(dihedrals["atomj"]) == n_dihedrals
    assert len(dihedrals["atomk"]) == n_dihedrals
    assert len(dihedrals["atoml"]) == n_dihedrals

    # Dihedral indices should be valid
    assert all(0 <= i < 16 for i in dihedrals["atomi"])
    assert all(0 <= j < 16 for j in dihedrals["atomj"])
    assert all(0 <= k < 16 for k in dihedrals["atomk"])
    assert all(0 <= l < 16 for l in dihedrals["atoml"])


def test_prmtop_read_residues(litfsi_prmtop):
    """Test reading residue information."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    assert "residue" in atoms

    residues = atoms["residue"]
    assert len(residues) == 16

    # Should have residue assignments for all atoms
    assert all(isinstance(r, (int, np.integer)) for r in residues)


def test_prmtop_forcefield_structure(litfsi_prmtop):
    """Test that forcefield has correct structure."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    # Check forcefield units
    assert ff.units == "real"

    # Check that forcefield has styles attribute
    assert hasattr(ff, "styles")
    assert ff.styles is not None

    # Check that we can get atom types
    atomtypes = ff.get_atomtypes()
    assert len(atomtypes) > 0

    # Check that we can get bond types
    bondtypes = ff.get_bondtypes()
    assert len(bondtypes) > 0

    # Check that we can get angle types
    angletypes = ff.get_angletypes()
    assert len(angletypes) > 0


def test_prmtop_charge_conversion_constant():
    """Test that CHARGE_CONVERSION_FACTOR is defined correctly."""
    assert CHARGE_CONVERSION_FACTOR == 18.2223


def test_prmtop_nonexistent_file():
    """Test error handling for nonexistent file."""
    reader = AmberPrmtopReader("/nonexistent/file.prmtop")
    frame = Frame()

    with pytest.raises(FileNotFoundError):
        reader.read(frame)


def test_prmtop_get_bond_with_H(litfsi_prmtop):
    """Test get_bond_with_H method."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    reader.read(frame)

    # Should return list of tuples
    bonds_with_h = reader.get_bond_with_H()
    assert isinstance(bonds_with_h, list)

    # Each tuple should have 5 elements: (type, i, j, force_constant, equil_length)
    for bond in bonds_with_h:
        assert len(bond) == 5
        assert isinstance(bond[0], int)  # type
        assert isinstance(bond[1], int)  # i
        assert isinstance(bond[2], int)  # j
        assert isinstance(bond[3], float)  # force constant
        assert isinstance(bond[4], float)  # equilibrium length


def test_prmtop_get_bond_without_H(litfsi_prmtop):
    """Test get_bond_without_H method."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    reader.read(frame)

    # Should return list of tuples
    bonds_without_h = reader.get_bond_without_H()
    assert isinstance(bonds_without_h, list)

    # Each tuple should have 5 elements
    for bond in bonds_without_h:
        assert len(bond) == 5
        assert isinstance(bond[0], int)  # type
        assert isinstance(bond[1], int)  # i
        assert isinstance(bond[2], int)  # j
        assert isinstance(bond[3], float)  # force constant
        assert isinstance(bond[4], float)  # equilibrium length


def test_prmtop_parse_angle_params(litfsi_prmtop):
    """Test parse_angle_params method."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    reader.read(frame)

    # Should return list of tuples
    angles = reader.parse_angle_params()
    assert isinstance(angles, list)

    # Each tuple should have 6 elements: (type, i, j, k, force_constant, equil_angle)
    for angle in angles:
        assert len(angle) == 6
        assert isinstance(angle[0], int)  # type
        assert isinstance(angle[1], int)  # i
        assert isinstance(angle[2], int)  # j
        assert isinstance(angle[3], int)  # k
        assert isinstance(angle[4], float)  # force constant
        assert isinstance(angle[5], float)  # equilibrium angle (in degrees)
        # Angle should be in reasonable range (0-180 degrees)
        assert 0 <= angle[5] <= 180


def test_prmtop_parse_dihedral_params(litfsi_prmtop):
    """Test parse_dihedral_params method."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    reader.read(frame)

    # Should return list of tuples
    dihedrals = reader.parse_dihedral_params()
    assert isinstance(dihedrals, list)

    # Each tuple should have 8 elements: (type, i, j, k, l, force_constant, phase, periodicity)
    for dihedral in dihedrals:
        assert len(dihedral) == 8
        assert isinstance(dihedral[0], int)  # type
        assert isinstance(dihedral[1], int)  # i
        assert isinstance(dihedral[2], int)  # j
        assert isinstance(dihedral[3], int)  # k
        assert isinstance(dihedral[4], int)  # l
        assert isinstance(dihedral[5], float)  # force constant
        assert isinstance(dihedral[6], float)  # phase
        assert isinstance(dihedral[7], int)  # periodicity


def test_prmtop_parse_nonbond_params(litfsi_prmtop):
    """Test parse_nonbond_params method."""
    reader = AmberPrmtopReader(litfsi_prmtop)
    frame = Frame()
    frame, ff = reader.read(frame)

    atoms = frame["atoms"]
    nonbond_params = reader.parse_nonbond_params(atoms)

    assert isinstance(nonbond_params, list)
    assert len(nonbond_params) == 16  # One per atom

    # Each tuple should have 3 elements: (atom_index, sigma, epsilon)
    for param in nonbond_params:
        assert len(param) == 3
        assert isinstance(param[0], int)  # atom index (1-based)
        assert isinstance(param[1], float)  # sigma
        assert isinstance(param[2], float)  # epsilon
        # Sigma and epsilon should be non-negative
        assert param[1] >= 0
        assert param[2] >= 0


def test_read_amber_helper_reads_prmtop_and_inpcrd(litfsi_prmtop, litfsi_inpcrd):
    """read_amber helper should return frame+ff with coordinates loaded."""
    frame, ff = read_amber(litfsi_prmtop, litfsi_inpcrd)

    assert frame["atoms"].nrows == 16
    assert ff is not None
    assert "x" in frame["atoms"]
    assert "y" in frame["atoms"]
    assert "z" in frame["atoms"]
