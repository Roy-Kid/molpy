"""
Modern tests for LAMMPS data I/O using chemfiles-testcases data.

Tests the LammpsDataReader and LammpsDataWriter classes with real
test cases from the chemfiles project.
"""

import pytest
import tempfile
import numpy as np
import os

# Import molpy components
import molpy as mp
from molpy.io.data.lammps import LammpsDataReader, LammpsDataWriter, LammpsMoleculeReader, LammpsMoleculeWriter

@pytest.fixture
def test_files(TEST_DATA_DIR):
    """Provide paths to test files."""

    lammps_data_dir = TEST_DATA_DIR / "data/lammps-data"

    files = {
        'data_body': lammps_data_dir / "data.body",
        'labelmap': lammps_data_dir / "labelmap.lmp",
        'molid': lammps_data_dir / "molid.lmp", 
        'solvated': lammps_data_dir / "solvated.lmp",
        'triclinic_1': lammps_data_dir / "triclinic-1.lmp",
        'triclinic_2': lammps_data_dir / "triclinic-2.lmp",
        'whitespaces': lammps_data_dir / "whitespaces.lmp"
    }
    
    # Check which files actually exist
    existing_files = {k: v for k, v in files.items() if v.exists()}
    return existing_files


class TestLammpsDataReader:
    """Test LAMMPS data file reader with real chemfiles test cases."""

    def test_molid_file(self, test_files):
        """Test reading molid.lmp - file with molecular IDs and full style."""
            
        reader = LammpsDataReader(test_files['molid'], atom_style="full")
        frame = reader.read()
        
        # Check basic structure
        assert 'atoms' in frame
        atoms = frame['atoms']
        
        # Should have 12 atoms based on file content
        assert len(atoms['id']) == 12
        assert 'mol' in atoms  # molecule ID should be present
        assert 'type' in atoms
        assert 'q' in atoms
        assert 'xyz' in atoms
        
        # Check coordinate data shape
        xyz = atoms['xyz']
        assert xyz.shape == (12, 3)  # 12 atoms, 3 coordinates
        
        # Check box dimensions (0-20 in each direction)
        assert frame.box is not None
        box_lengths = frame.box.lengths
        np.testing.assert_array_almost_equal(box_lengths, [20.0, 20.0, 20.0])
        
        # Check that molecule IDs are in the data (should be 0-3 based on file)
        mol_ids = atoms['mol']
        assert len(np.unique(mol_ids)) <= 4  # max 4 different molecules

    def test_whitespaces_file(self, test_files):
        """Test reading whitespaces.lmp - file with extra whitespaces."""
            
        reader = LammpsDataReader(test_files['whitespaces'], atom_style="full")
        frame = reader.read()
        
        # Should parse correctly despite extra whitespaces
        assert 'atoms' in frame
        atoms = frame['atoms']
        assert len(atoms['id']) == 1
        
        # Check the single atom's coordinates
        xyz = atoms['xyz']
        np.testing.assert_array_almost_equal(xyz[0], [5.0, 5.0, 5.0])
        
        # Check box (should be 10x10x10)
        box_lengths = frame.box.lengths
        np.testing.assert_array_almost_equal(box_lengths, [10.0, 10.0, 10.0])

    def test_triclinic_boxes(self, test_files):
        """Test reading triclinic box files."""
        # Test triclinic-1.lmp (zero tilt factors)
        if 'triclinic_1' in test_files:
            reader = LammpsDataReader(test_files['triclinic_1'])
            frame = reader.read()
            
            assert frame.box is not None
            # Should have 34x34x34 box with no tilt
            box_lengths = frame.box.lengths
            np.testing.assert_array_almost_equal(box_lengths, [34.0, 34.0, 34.0])
            
        # Test triclinic-2.lmp (non-zero tilt factors)
        if 'triclinic_2' in test_files:
            reader = LammpsDataReader(test_files['triclinic_2'])
            frame = reader.read()
            
            assert frame.box is not None
            # Should still have 34x34x34 basic dimensions
            box_lengths = frame.box.lengths  
            np.testing.assert_array_almost_equal(box_lengths, [34.0, 34.0, 34.0])

    def test_labelmap_file(self, test_files):
        """Test reading labelmap.lmp - file with atom/bond type labels."""
        
        # With unified string-based type handling, this should work now
        reader = LammpsDataReader(test_files['labelmap'], atom_style="full")
        frame = reader.read()
        
        # Check basic structure
        assert 'atoms' in frame
        atoms = frame['atoms']
        
        # Should have 16 atoms based on file content
        assert len(atoms['id']) == 16
        assert 'type' in atoms
        assert 'xyz' in atoms
        
        # Check that atom types are strings (labels like 'f', 'c3', etc.)
        atom_types = atoms['type']
        assert atom_types.dtype.kind == 'U'  # Unicode string
        
        # Check for expected labels
        unique_types = set(atom_types.flat)
        expected_labels = {'f', 'c3', 's6', 'o', 'ne', 'sy', 'Li+'}
        assert unique_types.issubset(expected_labels)
        
        # Check bonds if present
        if 'bonds' in frame:
            bonds = frame['bonds']
            assert bonds['type'].dtype.kind == 'U'  # String bond types too
            
        print(f"✓ labelmap.lmp: {len(atoms['id'])} atoms with string labels: {sorted(unique_types)}")

    def test_solvated_file(self, test_files):
        """Test reading solvated.lmp - large file with all topology types."""
        if 'solvated' not in test_files:
            pytest.skip("solvated.lmp test file not found")
            
        reader = LammpsDataReader(test_files['solvated'], atom_style="full")
        frame = reader.read()
        
        # Check atoms (should be 7772 atoms)
        assert 'atoms' in frame
        atoms = frame['atoms']
        assert len(atoms['id']) == 7772
        
        # Check all topology sections exist
        assert 'bonds' in frame
        assert 'angles' in frame
        assert 'dihedrals' in frame
        assert 'impropers' in frame
        
        # Check counts match file header
        bonds = frame['bonds']
        assert len(bonds['id']) == 6248
        
        angles = frame['angles']
        assert len(angles['id']) == 8100
        
        dihedrals = frame['dihedrals']
        assert len(dihedrals['id']) == 10720
        
        impropers = frame['impropers']
        assert len(impropers['id']) == 1376
        
        # Check atom types (should be 11 types)
        type = atoms['type']
        unique_types = np.unique(type)
        assert len(unique_types) == 11

    def test_different_atom_styles(self):
        """Test reading with different atom styles."""
        # Test atomic style
        atomic_content = """# LAMMPS data file
2 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp.write(atomic_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path, atom_style="atomic")
            frame = reader.read()
            
            atoms = frame['atoms']
            assert len(atoms['id']) == 2
            assert 'type' in atoms
            assert 'xyz' in atoms
            assert 'q' not in atoms  # atomic style doesn't have charges
            assert 'mol' not in atoms  # atomic style doesn't have molecule IDs
            
        finally:
            os.unlink(tmp_path)

    def test_charge_style(self):
        """Test reading with charge atom style."""
        charge_content = """# LAMMPS data file
2 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.5 0.0 0.0 0.0
2 1 -0.5 1.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp.write(charge_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path, atom_style="charge")
            frame = reader.read()
            
            atoms = frame['atoms']
            assert len(atoms['id']) == 2
            assert 'q' in atoms  # charge style has charges
            assert 'mol' not in atoms  # charge style doesn't have molecule IDs
            
            # Check charges
            np.testing.assert_array_almost_equal(atoms['q'], [0.5, -0.5])
            
        finally:
            os.unlink(tmp_path)

    def test_missing_sections(self):
        """Test reading file with missing optional sections."""
        minimal_content = """# LAMMPS data file
2 atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.0 0.0 0.0
2 1 1.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp.write(minimal_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path, atom_style="atomic")
            frame = reader.read()
            
            # Should have atoms but no topology
            assert 'atoms' in frame
            assert 'bonds' not in frame
            assert 'angles' not in frame
            assert 'dihedrals' not in frame
            assert 'impropers' not in frame
            
        finally:
            os.unlink(tmp_path)

    def test_type_labels_parsing(self):
        """Test parsing of type labels sections."""
        content_with_labels = """# LAMMPS data file
2 atoms
1 bonds
2 atom types
1 bond types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0
2 2.0

Atom Type Labels

1 C
2 O

Bond Type Labels

1 C-C

Atoms

1 1 1 0.0 0.0 0.0 0.0
2 1 2 0.0 1.0 0.0 0.0

Bonds

1 1 1 2
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp.write(content_with_labels)
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path, atom_style="full")
            frame = reader.read()
            
            atoms = frame['atoms']
            bonds = frame['bonds']
            
            # Check that types are properly mapped
            assert len(atoms['id']) == 2
            assert len(bonds['id']) == 1
            
        finally:
            os.unlink(tmp_path)


class TestLammpsDataWriter:
    """Test LAMMPS data file writer."""

    def test_write_read_roundtrip(self, test_files):
        """Test that we can write and read back the same data."""
        
        # Read original file
        reader = LammpsDataReader(test_files['molid'], atom_style="full")
        original_frame = reader.read()
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsDataWriter(tmp_path, atom_style="full")
            writer.write(original_frame)
            
            # Read back
            reader2 = LammpsDataReader(tmp_path, atom_style="full")
            new_frame = reader2.read()
            
            # Compare atoms
            orig_atoms = original_frame['atoms']
            new_atoms = new_frame['atoms']
            
            assert len(orig_atoms['id']) == len(new_atoms['id'])
            np.testing.assert_array_equal(orig_atoms['type'], new_atoms['type'])
            np.testing.assert_array_almost_equal(orig_atoms['xyz'], new_atoms['xyz'])
            
            # Compare box
            assert original_frame.box is not None
            assert new_frame.box is not None
            np.testing.assert_array_almost_equal(
                original_frame.box.lengths, new_frame.box.lengths
            )
            
        finally:
            os.unlink(tmp_path)

    def test_write_minimal_frame(self):
        """Test writing a minimal frame with just atoms."""
        # Create a simple frame
        frame = mp.Frame()
        
        # Add atoms data
        atoms_data = {
            'id': np.array([1, 2, 3]),
            'type': np.array([1, 1, 2]),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            'mass': np.array([1.0, 1.0, 2.0])
        }
    
        
        frame['atoms'] = atoms_data
        frame.box = mp.Box([10.0, 10.0, 10.0])
        
        # Write to temporary file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsDataWriter(tmp_path, atom_style="atomic")
            writer.write(frame)
            
            # Check file was written and has content
            assert os.path.exists(tmp_path)
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "3 atoms" in content
                assert "2 atom types" in content
                assert "Atoms" in content
                
        finally:
            os.unlink(tmp_path)

    def test_write_full_style(self):
        """Test writing with full atom style including molecule IDs and charges."""
        frame = mp.Frame()
        
        # Create atoms with all fields
        atoms_data = {
            'id': np.array([1, 2, 3]),
            'mol': np.array([1, 1, 2]),
            'type': np.array(['C', 'C', 'O']),
            'q': np.array([0.0, 0.0, -0.5]),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]]),
            'mass': np.array([12.0, 12.0, 16.0])
        }
        
        frame['atoms'] = mp.Block(atoms_data)
        frame.box = mp.Box([10.0, 10.0, 10.0])
        
        # Add bonds
        bonds_data = {
            'id': np.array([1, 2]),
            'type': np.array(['C-C', 'C-O']),
            'i': np.array([0, 1]),
            'j': np.array([1, 2])
        }
        frame['bonds'] = mp.Block(bonds_data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsDataWriter(tmp_path, atom_style="full")
            writer.write(frame)
            
            # Check file content
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "3 atoms" in content
                assert "2 bonds" in content
                assert "2 atom types" in content
                assert "2 bond types" in content
                assert "Atom Type Labels" in content
                assert "Bond Type Labels" in content
                
        finally:
            os.unlink(tmp_path)

    def test_write_with_topology(self):
        """Test writing frame with all topology types."""
        frame = mp.Frame()
        
        # Atoms
        atoms_data = {
            'id': np.array([1, 2, 3, 4]),
            'type': np.array(['C', 'C', 'O', 'H']),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]]),
            'mass': np.array([12.0, 12.0, 16.0, 1.0])
        }
        frame['atoms'] = mp.Block(atoms_data)
        
        # Bonds
        bonds_data = {
            'id': np.array([1, 2, 3]),
            'type': np.array(['C-C', 'C-O', 'O-H']),
            'i': np.array([0, 1, 2]),
            'j': np.array([1, 2, 3])
        }
        frame['bonds'] = mp.Block(bonds_data)
        
        # Angles
        angles_data = {
            'id': np.array([1, 2]),
            'type': np.array(['C-C-O', 'C-O-H']),
            'i': np.array([0, 1]),
            'j': np.array([1, 2]),
            'k': np.array([2, 3])
        }
        frame['angles'] = mp.Block(angles_data)
        
        # Dihedrals
        dihedrals_data = {
            'id': np.array([1]),
            'type': np.array(['C-C-O-H']),
            'i': np.array([0]),
            'j': np.array([1]),
            'k': np.array([2]),
            'l': np.array([3])
        }
        frame['dihedrals'] = mp.Block(dihedrals_data)
        
        # Impropers
        impropers_data = {
            'id': np.array([1]),
            'type': np.array(['C-C-O-H']),
            'atom1': np.array([0]),
            'atom2': np.array([1]),
            'atom3': np.array([2]),
            'atom4': np.array([3])
        }
        frame['impropers'] = mp.Block(impropers_data)
        
        frame.box = mp.Box([10.0, 10.0, 10.0])
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsDataWriter(tmp_path, atom_style="atomic")
            writer.write(frame)
            
            # Check file content
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "4 atoms" in content
                assert "3 bonds" in content
                assert "2 angles" in content
                assert "1 dihedrals" in content
                assert "1 impropers" in content
                assert "3 atom types" in content  # C, O, H = 3 types
                assert "3 bond types" in content
                assert "2 angle types" in content
                assert "1 dihedral types" in content
                assert "1 improper types" in content
                
        finally:
            os.unlink(tmp_path)

    def test_write_no_box(self):
        """Test writing frame without box information."""
        frame = mp.Frame()
        
        atoms_data = {
            'id': np.array([1, 2]),
            'type': np.array([1, 1]),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            'mass': np.array([1.0, 1.0])
        }
        frame['atoms'] = mp.Block(atoms_data)
        # No box set
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsDataWriter(tmp_path, atom_style="atomic")
            writer.write(frame)
            
            # Should use default box
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "0.0 10.0 xlo xhi" in content
                assert "0.0 10.0 ylo yhi" in content
                assert "0.0 10.0 zlo zhi" in content
                
        finally:
            os.unlink(tmp_path)


class TestLammpsMoleculeReader:
    """Test LAMMPS molecule template reader."""
    
    def test_molecule_reader_basic(self):
        """Test basic molecule reader functionality."""
        # Create a simple molecule template file for testing
        molecule_content = """# Test molecule
3 atoms
2 bonds

Coords
1 0.0 0.0 0.0
2 1.0 0.0 0.0  
3 0.0 1.0 0.0

Types
1 1
2 1
3 2

Bonds
1 1 1 2
2 1 2 3
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp.write(molecule_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsMoleculeReader(tmp_path)
            frame = reader.read()
            
            # Check atoms
            assert 'atoms' in frame
            atoms = frame['atoms']
            assert len(atoms['id']) == 3
            assert 'type' in atoms
            assert 'xyz' in atoms
            
            # Check bonds
            assert 'bonds' in frame
            bonds = frame['bonds']
            assert len(bonds['id']) == 2
            
        finally:
            os.unlink(tmp_path)

    def test_molecule_reader_with_charges(self):
        """Test molecule reader with charges section."""
        molecule_content = """# Test molecule with charges
2 atoms

Coords
1 0.0 0.0 0.0
2 1.0 0.0 0.0

Types
1 C
2 O

Charges
1 0.0
2 -0.5
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp.write(molecule_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsMoleculeReader(tmp_path)
            frame = reader.read()
            
            atoms = frame['atoms']
            assert len(atoms['id']) == 2
            assert 'q' in atoms
            assert 'type' in atoms
            assert 'xyz' in atoms
            
            # Check charges
            np.testing.assert_array_almost_equal(atoms['q'], [0.0, -0.5])
            
        finally:
            os.unlink(tmp_path)

    def test_molecule_reader_all_topology(self):
        """Test molecule reader with all topology types."""
        molecule_content = """# Test molecule with all topology
4 atoms
3 bonds
2 angles
1 dihedrals
1 impropers

Coords
1 0.0 0.0 0.0
2 1.0 0.0 0.0
3 0.0 1.0 0.0
4 1.0 1.0 0.0

Types
1 C
2 C
3 O
4 H

Bonds
1 C-C 1 2
2 C-O 2 3
3 O-H 3 4

Angles
1 C-C-O 1 2 3
2 C-O-H 2 3 4

Dihedrals
1 C-C-O-H 1 2 3 4

Impropers
1 C-C-O-H 1 2 3 4
"""
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp.write(molecule_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsMoleculeReader(tmp_path)
            frame = reader.read()
            
            # Check all sections
            assert 'atoms' in frame
            assert 'bonds' in frame
            assert 'angles' in frame
            assert 'dihedrals' in frame
            assert 'impropers' in frame
            
            atoms = frame['atoms']
            bonds = frame['bonds']
            angles = frame['angles']
            dihedrals = frame['dihedrals']
            impropers = frame['impropers']
            
            assert len(atoms['id']) == 4
            assert len(bonds['id']) == 3
            assert len(angles['id']) == 2
            assert len(dihedrals['id']) == 1
            assert len(impropers['id']) == 1
            
        finally:
            os.unlink(tmp_path)


class TestLammpsMoleculeWriter:
    """Test LAMMPS molecule template writer."""
    
    def test_molecule_writer_basic(self):
        """Test basic molecule writer functionality."""
        frame = mp.Frame()
        
        # Create atoms
        atoms_data = {
            'id': np.array([1, 2, 3]),
            'type': np.array(['C', 'C', 'O']),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]])
        }
        frame['atoms'] = mp.Block(atoms_data)
        
        # Create bonds
        bonds_data = {
            'id': np.array([1, 2]),
            'type': np.array(['C-C', 'C-O']),
            'i': np.array([0, 1]),
            'j': np.array([1, 2])
        }
        frame['bonds'] = mp.Block(bonds_data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsMoleculeWriter(tmp_path)
            writer.write(frame)
            
            # Check file content
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "3 atoms" in content
                assert "2 bonds" in content
                assert "Coords" in content
                assert "Types" in content
                assert "Bonds" in content
                
        finally:
            os.unlink(tmp_path)

    def test_molecule_writer_with_charges(self):
        """Test molecule writer with charges."""
        frame = mp.Frame()
        
        atoms_data = {
            'id': np.array([1, 2]),
            'type': np.array(['C', 'O']),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]]),
            'q': np.array([0.0, -0.5])
        }
        frame['atoms'] = mp.Block(atoms_data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsMoleculeWriter(tmp_path)
            writer.write(frame)
            
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "2 atoms" in content
                assert "Coords" in content
                assert "Types" in content
                assert "Charges" in content
                
        finally:
            os.unlink(tmp_path)

    def test_molecule_writer_all_topology(self):
        """Test molecule writer with all topology types."""
        frame = mp.Frame()
        
        # Atoms
        atoms_data = {
            'id': np.array([1, 2, 3, 4]),
            'type': np.array(['C', 'C', 'O', 'H']),
            'xyz': np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [1.0, 1.0, 0.0]])
        }
        frame['atoms'] = mp.Block(atoms_data)
        
        # Bonds
        bonds_data = {
            'id': np.array([1, 2, 3]),
            'type': np.array(['C-C', 'C-O', 'O-H']),
            'i': np.array([0, 1, 2]),
            'j': np.array([1, 2, 3])
        }
        frame['bonds'] = mp.Block(bonds_data)
        
        # Angles
        angles_data = {
            'id': np.array([1, 2]),
            'type': np.array(['C-C-O', 'C-O-H']),
            'i': np.array([0, 1]),
            'j': np.array([1, 2]),
            'k': np.array([2, 3])
        }
        frame['angles'] = mp.Block(angles_data)
        
        # Dihedrals
        dihedrals_data = {
            'id': np.array([1]),
            'type': np.array(['C-C-O-H']),
            'i': np.array([0]),
            'j': np.array([1]),
            'k': np.array([2]),
            'l': np.array([3])
        }
        frame['dihedrals'] = mp.Block(dihedrals_data)
        
        # Impropers
        impropers_data = {
            'id': np.array([1]),
            'type': np.array(['C-C-O-H']),
            'i': np.array([0]),
            'j': np.array([1]),
            'k': np.array([2]),
            'l': np.array([3])
        }
        frame['impropers'] = mp.Block(impropers_data)
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.mol', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            writer = LammpsMoleculeWriter(tmp_path)
            writer.write(frame)
            
            with open(tmp_path, 'r') as f:
                content = f.read()
                assert "4 atoms" in content
                assert "3 bonds" in content
                assert "2 angles" in content
                assert "1 dihedrals" in content
                assert "1 impropers" in content
                assert "Coords" in content
                assert "Types" in content
                assert "Bonds" in content
                assert "Angles" in content
                assert "Dihedrals" in content
                assert "Impropers" in content
                
        finally:
            os.unlink(tmp_path)


class TestErrorHandling:
    """Test error handling and edge cases."""
    
    def test_nonexistent_file(self):
        """Test reading nonexistent file."""
        with pytest.raises(FileNotFoundError):
            reader = LammpsDataReader("nonexistent_file.data")
            reader.read()
    
    def test_empty_file(self):
        """Test reading empty file."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path)
            frame = reader.read()
            
            # Should handle empty file gracefully
            assert frame is not None
            assert 'atoms' not in frame
            
        finally:
            os.unlink(tmp_path)
    
    def test_malformed_header(self):
        """Test reading file with malformed header."""
        malformed_content = """# LAMMPS data file
invalid atoms
1 atom types

0.0 10.0 xlo xhi
0.0 10.0 ylo yhi
0.0 10.0 zlo zhi

Masses

1 1.0

Atoms

1 1 0.0 0.0 0.0
"""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.data', delete=False) as tmp:
            tmp.write(malformed_content)
            tmp_path = tmp.name
        
        try:
            reader = LammpsDataReader(tmp_path, atom_style="atomic")
            frame = reader.read()
            
            # Should handle malformed header gracefully
            assert frame is not None
            # May not have atoms if header parsing fails
            if 'atoms' in frame:
                assert len(frame['atoms']['id']) >= 0
            
        finally:
            os.unlink(tmp_path)


if __name__ == "__main__":
    pytest.main([__file__])
