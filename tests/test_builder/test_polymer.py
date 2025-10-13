import pytest

import molpy as mp


class TestCoarseGrainPolymerBuilder:

    class Bead(mp.Struct):

        def __init__(self, name, xyz):
            super().__init__(name=name, xyz=xyz)
            self["xyz"] = xyz

    def test_linear(self):

        seq = ["A", "B", "B", "A"]
        bead_mapping = {
            "A": mp.Spatial(self.Bead(name="A", xyz=[0, 0, 0]), anchor=0),
            "B": mp.Spatial(self.Bead(name="B", xyz=[0, 0, 0]), anchor=0),
        }
        # random walk path
        path = [
            [0, 0, 0],
            [1, 0, 0],
            [1, 1, 0],
            [1, 1, 1],
        ]
        builder = mp.builder.PolymerBuilder(bead_mapping, style="CoarseGrain")
        polymer = builder.linear(seq, path)
        assert len(polymer["atoms"]) == 4

    def test_build(self):

        # A - B - B - B - A
        #     |       |
        #     C       C

        seq = ["A", "B", "B", "B", "A", "C", "C"]
        bead_mapping = {
            "A": mp.Spatial(self.Bead(name="A", xyz=[0, 0, 0])),
            "B": mp.Spatial(self.Bead(name="B", xyz=[0, 0, 0])),
            "C": mp.Spatial(self.Bead(name="C", xyz=[0, 0, 0])),
        }
        connect = [
            (0, 1),
            (1, 2),
            (2, 3),
            (3, 4),  # main chain
            (1, 5),  # branch 1
            (3, 6),  # branch 2
        ]
        path = [
            [0, 0, 0],  # A
            [1, 0, 0],  # B
            [2, 0, 0],  # B
            [3, 0, 0],  # B
            [4, 0, 0],  # A
            [1, 1, 0],  # C
            [3, -1, 0],  # C
        ]
        builder = mp.builder.PolymerBuilder(bead_mapping, style="CoarseGrain")
        polymer = builder.build(seq, connect, path)
        assert len(polymer) == 7


class TestAtomisticPolymerBuilder:

    class Methyl(mp.Atomistic, mp.Spatial):

        def __post_init__(self, **props):
            C1 = self.atoms.add(mp.Atom(name="C1", xyz=(0.0, 0.0, 0.0)))
            H11 = self.atoms.add(mp.Atom(name="H11", xyz=(0.0, 0.0, 1.1)))
            H12 = self.atoms.add(mp.Atom(name="H12", xyz=(1.0, 0.0, 0.0)))
            H13 = self.atoms.add(mp.Atom(name="H13", xyz=(0.0, 1.0, 0.0)))
            self.bonds.add(mp.Bond(atom1=C1, atom2=H11))
            self.bonds.add(mp.Bond(atom1=C1, atom2=H12))
            self.bonds.add(mp.Bond(atom1=C1, atom2=H13))
            return super().__post_init__()

    class Methylene(mp.Atomistic, mp.Spatial):

        def __post_init__(self, **props):
            C2 = self.atoms.add(mp.Atom(name="C2", xyz=(0.0, 0.0, 0.0)))
            H21 = self.atoms.add(mp.Atom(name="H21", xyz=(0.0, 0.0, 1.1)))
            H22 = self.atoms.add(mp.Atom(name="H22", xyz=(1.0, 0.0, 0.0)))
            self.bonds.add(mp.Bond(atom1=C2, atom2=H21))
            self.bonds.add(mp.Bond(atom1=C2, atom2=H22))
            return super().__post_init__()

    def test_xyz(self):
        methyl = self.Methyl(name="Methyl")
        assert methyl.xyz.shape == (4, 3)  # 1C + 3H

    def test_linear(self):

        methyl = self.Methyl(name="Methyl")
        MethylMonomer = mp.builder.Monomer(
            methyl, head=methyl.atoms[0], tail=methyl.atoms[0]
        )
        methylene = self.Methylene(name="Methylene")
        MethyleneMonomer = mp.builder.Monomer(
            methylene, head=methylene.atoms[0], tail=methylene.atoms[0]
        )

        seq = ["Methyl", "Methylene", "Methylene", "Methyl"]
        monomer_mapping = {
            "Methyl": MethylMonomer,
            "Methylene": MethyleneMonomer,
        }
        path = [
            [0, 0, 0],
            [1.5, 0, 0],
            [3.0, 0, 0],
            [4.5, 0, 0],
        ]
        builder = mp.builder.PolymerBuilder(monomer_mapping, style="Atomistic")
        polymer = builder.linear(seq, path)
        assert len(polymer["atoms"]) == 14
        assert len(polymer["bonds"]) == 13
