import glob
import itertools as it
import os
import molpy as mp
import pytest

pytest.skip("OPLS typifier not yet ported to xarray", allow_module_level=True)


class TestOPLS:
    @pytest.fixture(autouse=True)
    def initdir(self, tmpdir):
        tmpdir.chdir()



    @pytest.fixture(autouse=True)
    def OPLS_TESTFILES_DIR(self, test_data_path):
        return os.path.join(test_data_path, "opls_validation")
    
    @pytest.fixture(autouse=True)
    def top_files(self, test_data_path):
        return glob.glob(os.path.join(test_data_path, "opls_validation", "*/*.top"))
    
    @pytest.fixture(autouse=True)
    def mol2_files(self, test_data_path):
        return glob.glob(os.path.join(test_data_path, "opls_validation", "*/*.mol2"))
    
    @pytest.fixture(autouse=True)
    def implemented_tests_path(self, test_data_path):
        return os.path.join(test_data_path, "opls_validation", "implemented_opls_tests.txt")
    
    @pytest.fixture(autouse=True)
    def correctly_implemented(self, implemented_tests_path):
        with open(implemented_tests_path) as f:
            return [line.strip() for line in f]
        
    @pytest.fixture(scope="session")
    def oplsaa(self):
        return forcefields.load_OPLSAA()


    # Please update this file if you implement atom typing for a test case.
    # You can automatically update the files by running the below function
    # `find_correctly_implemented`.
    

    def find_correctly_implemented(self):

        with open(self.implemented_tests_path, "a") as fh:
            for mol_path in it.chain(self.top_files, self.mol2_files):
                _, mol_file = os.path.split(mol_path)
                mol_name, ext = os.path.splitext(mol_file)
                try:
                    self.test_atomtyping(mol_name)
                except Exception as e:
                    print(e)
                    continue
                else:
                    if mol_name not in self.correctly_implemented:
                        fh.write("{}\n".format(mol_name))

    def test_opls_metadata(self, oplsaa):
        assert oplsaa.name == "OPLS-AA"
        assert oplsaa.version == "0.0.3"
        assert oplsaa.combining_rule == "geometric"

    def test_atomtyping(self, mol_name, oplsaa, OPLS_TESTFILES_DIR, correctly_implemented):
        for mol_name in correctly_implemented:
            files = glob.glob(os.path.join(OPLS_TESTFILES_DIR, mol_name, "*"))
            for mol_file in files:
                _, ext = os.path.splitext(mol_file)
                if ext == ".top":
                    top_filename = "{}.top".format(mol_name)
                    gro_filename = "{}.gro".format(mol_name)
                    top_path = os.path.join(OPLS_TESTFILES_DIR, mol_name, top_filename)
                    gro_path = os.path.join(OPLS_TESTFILES_DIR, mol_name, gro_filename)
                    structure = pmd.load_file(top_path, xyz=gro_path, parametrize=False)
                elif ext == ".mol2":
                    frame = mp.Frame()
                    mol2_path = os.path.join(OPLS_TESTFILES_DIR, mol_name, mol_file)
                    structure = (
                        mp.io.read_mol2(mol2_path, frame)
                        .to_struct()
                        .get_topology(attrs=["name", "number"])
                    )
            atomtype(structure, oplsaa)

    def test_full_parametrization(self, oplsaa, OPLS_TESTFILES_DIR):
        top = os.path.join(OPLS_TESTFILES_DIR, "benzene/benzene.top")
        gro = os.path.join(OPLS_TESTFILES_DIR, "benzene/benzene.gro")
        structure = pmd.load_file(top, xyz=gro)
        parametrized = oplsaa.apply(structure)

        assert sum((1 for at in parametrized.atoms if at.type == "opls_145")) == 6
        assert sum((1 for at in parametrized.atoms if at.type == "opls_146")) == 6
        assert len(parametrized.bonds) == 12
        assert all(x.type for x in parametrized.bonds)
        assert len(parametrized.angles) == 18
        assert all(x.type for x in parametrized.angles)
        assert len(parametrized.rb_torsions) == 24
        assert all(x.type for x in parametrized.dihedrals)
        assert parametrized.combining_rule == "geometric"

    def test_improper_in_structure(self, OPLS_TESTFILES_DIR):
        files_with_impropers = [
            ("o-xylene", 6),
            ("nitroethane", 1),
            ("fluorobenzene", 6),
            ("N-methylformamide", 2),
            ("formamide", 2),
            ("toluene", 6),
            ("3-methylphenol", 6),
            ("4-methylphenol", 6),
            ("2-methylphenol", 6),
            ("nitrobenzene", 7),
            ("dimethylformamide", 2),
            ("nitromethane", 1),
            ("NN-dimethylformamide", 2),
            ("124-trimethylbenzene", 6),
            ("phenol", 6),
            ("ethylbenzene", 6),
            ("NN-dimethylacetamide", 2),
            ("13-difluorobenzene", 6),
        ]  # found in the "impropers" sections of molecule_name.top
        for molecule, n_impropers in files_with_impropers:
            top = os.path.join(OPLS_TESTFILES_DIR, molecule + "/" + molecule + ".top")
            gro = os.path.join(OPLS_TESTFILES_DIR, molecule + "/" + molecule + ".gro")
            structure = pmd.load_file(top, xyz=gro)
            impropers = []
            [
                impropers.append(dihedral)
                for dihedral in structure.dihedrals
                if dihedral.improper
            ]
            assert len(impropers) == n_impropers


if __name__ == "__main__":
    TestOPLS().find_correctly_implemented()
