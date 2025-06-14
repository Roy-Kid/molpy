import pytest
import molpy as mp

class TestEntity:

    @pytest.fixture(scope='class')
    def entity(self):
        return mp.Entity(name='C', type='atom')
    
    def test_copy(self, entity):
        entity_copy = entity.copy()
        assert entity_copy["name"] == entity["name"]
        assert entity_copy is not entity

    def test_to_dict(self, entity):
        entity_dict = entity.to_dict()
        assert entity_dict["name"] == entity["name"]
        assert entity_dict["type"] == entity["type"]

    def test_hash(self, entity):
        
        assert hash(entity)
        entity_copy = entity.copy()
        assert hash(entity_copy) != hash(entity)
    

class TestAtom:

    @pytest.fixture
    def atom(self):
        return mp.Atom(name='C', xyz=[0, 0, 0])

    def test_compare(self, atom):
        atom1 = mp.Atom(name='C', xyz=[0, 0, 0])
        assert atom != atom1, "Different instances should have different hashes"

    def test_hash(self, atom):
        atom1 = mp.Atom(name='C', xyz=[0, 0, 0])
        assert getattr(atom, '__hash__')
        assert hash(atom) != hash(atom1), "Different instances should have different hashes"

    def test_copy(self, atom):
        atom_copy = atom.copy()
        assert atom_copy is not atom
        assert atom_copy["name"] == atom["name"]
        assert atom_copy["xyz"] == atom["xyz"]

class TestBond:

    @pytest.fixture
    def bond(self):
        return mp.Bond(
            itom=mp.Atom(name='C'),
            jtom=mp.Atom(name='H'),
            order=1
        )
    
    def test_copy(self, bond):

        bond_copy = bond.copy()
        assert bond_copy.itom is bond.itom
        assert bond_copy.jtom is bond.jtom
        assert bond_copy is not bond


class TestStruct:

    @pytest.fixture(scope='class')
    def Methane(self):

        class Methane(mp.Struct):

            def __init__(self):
                super().__init__()
                C = self.def_atom(name="C", xyz=(0.1, 0, -0.07))
                H1 = self.def_atom(name="H1", xyz=(-0.1, 0, -0.07))
                H2 = self.def_atom(name="H2", xyz=(0., 0.1, 0.07))
                H3 = self.def_atom(name="H3", xyz=(0., -0.1, 0.07))
                H4 = self.def_atom(name="H4", xyz=(0., 0, 0.1))
                self.def_bond(C, H1)
                self.def_bond(C, H2)
                self.def_bond(C, H3)
                self.def_bond(C, H4)

        return Methane
    
    @pytest.fixture(scope='class', name="methane")
    def simple_struct(self, Methane):
        return Methane()
    
    @pytest.fixture(scope='class', name="Ethane")
    def Ethane(self):
        class Ethane(mp.Struct):
            def __init__(self):
                super().__init__()
                ch3_1 = mp.builder.CH3()
                ch3_1.move(-ch3_1["atoms"][0].xyz)
                ch3_2 = mp.builder.CH3()
                ch3_2.move(-ch3_2["atoms"][0].xyz)
                ch3_2.rotate(180, axis=(0, 1, 0))
                self.add_struct(ch3_1)
                self.add_struct(ch3_2)
                self.def_bond(ch3_1["atoms"][0], ch3_2["atoms"][0])

        return Ethane
    
    @pytest.fixture(scope='class', name="ethane")
    def hierarchical_struct(self, Ethane):
        return Ethane()
    
    def test_get_atoms(self, methane):
        atoms = methane.atoms
        assert len(atoms) == 5
        assert all(isinstance(atom, mp.Atom) for atom in atoms)
        assert all(atom.name in ["C", "H1", "H2", "H3", "H4"] for atom in atoms)

    def test_get_bonds(self, methane):
        bonds = methane.bonds
        assert len(bonds) == 4
        assert all(isinstance(bond, mp.Bond) for bond in bonds)

    def test_calc_angles_dihedrals(self, ethane):
        
        topo = ethane.get_topology()
        angles = topo.angles
        dihedrals = topo.dihedrals
        assert len(ethane.angles) == 0  # angles have to be add manually
        assert len(ethane.dihedrals) == 0
        assert len(angles) == 12
        assert len(dihedrals) == 9

    def test_copy_simple_struct(self, methane):
        meth = methane.clone()
        atoms = meth.atoms
        assert len(atoms) == 5
        assert all(isinstance(atom, mp.Atom) for atom in atoms)
        assert all(not a1 is a2 for a1, a2 in zip(methane.atoms, atoms))

        bonds = meth.bonds
        assert len(bonds) == 4
        # test itom jton in bond is the reference to new atom
        assert all(not b1.itom is b2.itom for b1, b2 in zip(methane.bonds, bonds))
