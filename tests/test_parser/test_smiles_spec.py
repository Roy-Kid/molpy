import pytest

from molpy.parser.smiles import (
    SmilesParser,
    BondStereo,
    BondOrder,
    MolGraphIR,
    parse_smiles,
)


@pytest.fixture
def parser() -> SmilesParser:
    return SmilesParser()


# --- Basic connectivity ---

def test_single_carbon_atom(parser: SmilesParser):
    ir = parser.parse("C")
    assert len(ir.atoms) == 1
    atom = ir.atoms[0]
    assert atom.element == "C"
    assert atom.aromatic is False
    assert ir.total_charge == 0
    assert ir.has_disconnected is False


def test_single_oxygen_atom(parser: SmilesParser):
    ir = parser.parse("O")
    assert len(ir.atoms) == 1
    assert ir.atoms[0].element == "O"


def test_ethane_single_bond(parser: SmilesParser):
    ir = parser.parse("CC")
    assert len(ir.atoms) == 2
    assert [a.element for a in ir.atoms] == ["C", "C"]
    assert len(ir.bonds) == 1
    bond = ir.bonds[0]
    assert bond.order == BondOrder.SINGLE
    assert {bond.u, bond.v} == {0, 1}


def test_ethanol_three_atoms(parser: SmilesParser):
    ir = parser.parse("CCO")
    assert [a.element for a in ir.atoms] == ["C", "C", "O"]
    assert len(ir.bonds) == 2
    assert all(b.order == BondOrder.SINGLE for b in ir.bonds)


def test_carbon_dioxide_double_bonds(parser: SmilesParser):
    ir = parser.parse("O=C=O")
    assert [a.element for a in ir.atoms] == ["O", "C", "O"]
    assert len(ir.bonds) == 2
    assert all(b.order == BondOrder.DOUBLE for b in ir.bonds)


def test_complex_molecule_aspirin_like(parser: SmilesParser):
    ir = parser.parse("CC(=O)OC1=CC=CC=C1C(=O)O")
    assert len(ir.atoms) > 10
    assert ir.total_charge == 0
    assert ir.has_disconnected is False


# --- Branching behaviour ---

def test_branch_double_then_single(parser: SmilesParser):
    ir = parser.parse("C(=O)O")
    doubles = sum(1 for b in ir.bonds if b.order == BondOrder.DOUBLE)
    singles = sum(1 for b in ir.bonds if b.order == BondOrder.SINGLE)
    assert doubles == 1
    assert singles == 1


def test_branching_multiple_substituents(parser: SmilesParser):
    ir = parser.parse("C(Cl)(Br)I")
    assert len(ir.atoms) == 4
    central = ir.atoms[0].id
    neighbours = [b for b in ir.bonds if central in (b.u, b.v)]
    assert len(neighbours) == 3
    assert all(b.order == BondOrder.SINGLE for b in neighbours)
    assert sorted(a.element for a in ir.atoms[1:]) == ["Br", "Cl", "I"]


def test_nested_branching(parser: SmilesParser):
    ir = parser.parse("CC(=O)(O)C")
    central = ir.atoms[1].id
    assert sum(1 for b in ir.bonds if central in (b.u, b.v)) == 4


# --- Ring closures ---

def test_three_member_ring_closure(parser: SmilesParser):
    ir = parser.parse("C1CC1")
    assert len(ir.atoms) == 3
    assert len(ir.bonds) == 3


def test_explicit_ring_bond_symbol_on_closure(parser: SmilesParser):
    ir = parser.parse("C1CC=1")
    singles = sum(1 for b in ir.bonds if b.order == BondOrder.SINGLE)
    doubles = sum(1 for b in ir.bonds if b.order == BondOrder.DOUBLE)
    assert singles == 2
    assert doubles == 1


def test_explicit_ring_bond_symbol_on_open(parser: SmilesParser):
    ir = parser.parse("C=1CC1")
    doubles = sum(1 for b in ir.bonds if b.order == BondOrder.DOUBLE)
    singles = sum(1 for b in ir.bonds if b.order == BondOrder.SINGLE)
    assert doubles == 1
    assert singles == 2


def test_two_digit_ring_indices(parser: SmilesParser):
    ir = parser.parse("C%12CC%12")
    assert [a.element for a in ir.atoms] == ["C", "C", "C"]
    assert len(ir.bonds) == 3


def test_percent_ninety_nine_ring_index(parser: SmilesParser):
    ir = parser.parse("C%99CC%99")
    assert len(ir.atoms) == 3
    assert len(ir.bonds) == 3


def test_fused_rings(parser: SmilesParser):
    ir = parser.parse("C1CCC2CCCCC12")
    assert len(ir.atoms) == 9
    # Two rings share two atoms; ensure bond count reflects fused structure
    assert len(ir.bonds) == 10


# --- Aromatic systems ---

def test_aromatic_ring(parser: SmilesParser):
    ir = parser.parse("c1ccccc1")
    assert len(ir.atoms) == 6
    assert all(atom.aromatic for atom in ir.atoms)
    assert len(ir.bonds) == 6
    assert all(bond.order == BondOrder.AROMATIC for bond in ir.bonds)


def test_aromatic_bond_two_atoms(parser: SmilesParser):
    ir = parser.parse("c:c")
    assert len(ir.atoms) == 2
    assert all(atom.aromatic for atom in ir.atoms)
    assert ir.bonds[0].order == BondOrder.AROMATIC


def test_aromatic_hetero_ring(parser: SmilesParser):
    ir = parser.parse("c1ncccc1")
    assert len(ir.atoms) == 6
    assert any(atom.element == "n" for atom in ir.atoms)
    assert sum(1 for atom in ir.atoms if atom.aromatic) == 6


def test_aromatic_bracket_protonated(parser: SmilesParser):
    ir = parser.parse("[nH+]")
    atom = ir.atoms[0]
    assert atom.element == "n"
    assert atom.aromatic is True
    assert atom.h_count_explicit == 1
    assert atom.charge == 1


# --- Bond types and stereochemistry ---

def test_double_bond(parser: SmilesParser):
    ir = parser.parse("C=C")
    assert ir.bonds[0].order == BondOrder.DOUBLE


def test_triple_bond(parser: SmilesParser):
    ir = parser.parse("C#C")
    assert ir.bonds[0].order == BondOrder.TRIPLE


def test_stereo_bonds_slashes(parser: SmilesParser):
    ir = parser.parse("C/C=C\\C")
    assert len(ir.atoms) == 4
    assert len(ir.bonds) == 3
    assert ir.bonds[0].stereo == BondStereo.UP
    assert ir.bonds[-1].stereo == BondStereo.DOWN


def test_stereo_bonds_mixed(parser: SmilesParser):
    ir = parser.parse("C\\C=C/C")
    assert len(ir.bonds) == 3
    assert ir.bonds[0].stereo == BondStereo.DOWN
    assert ir.bonds[-1].stereo == BondStereo.UP


# --- Bracket atom features ---

def test_bracket_atom_charge(parser: SmilesParser):
    ir = parser.parse("[NH4+]")
    atom = ir.atoms[0]
    assert atom.charge == 1
    assert atom.h_count_explicit == 4
    assert ir.total_charge == 1


def test_bracket_atom_isotope(parser: SmilesParser):
    ir = parser.parse("[13C]")
    assert ir.atoms[0].isotope == 13


def test_bracket_atom_chiral_simple(parser: SmilesParser):
    ir = parser.parse("[C@H](O)N")
    atom = ir.atoms[0]
    assert atom.chiral == "@"
    # Ensure explicit hydrogen recorded even with additional branches
    assert atom.h_count_explicit == 1


def test_bracket_atom_chiral_extended(parser: SmilesParser):
    ir = parser.parse("[C@TH1](F)(Cl)(Br)I")
    atom = ir.atoms[0]
    assert atom.chiral == "@TH1"
    # 1 central atom + four substituents
    assert sorted(a.element for a in ir.atoms) == ["Br", "C", "Cl", "F", "I"]


def test_multiple_charges(parser: SmilesParser):
    ir1 = parser.parse("[N++]")
    assert ir1.atoms[0].charge == 2
    ir2 = parser.parse("[O-2]")
    assert ir2.atoms[0].charge == -2


def test_atom_class_mapping(parser: SmilesParser):
    ir = parser.parse("[C:12]")
    assert ir.atoms[0].atom_map == 12


def test_bracket_wildcard(parser: SmilesParser):
    ir = parser.parse("[*]")
    assert ir.atoms[0].element == "*"


def test_bracket_atom_properties(parser: SmilesParser):
    ir = parser.parse("[13CH3+]")
    atom = ir.atoms[0]
    assert atom.element == "C"
    assert atom.isotope == 13
    assert atom.h_count_explicit == 3
    assert atom.charge == 1


def test_halogen_organic_subset(parser: SmilesParser):
    ir = parser.parse("CClBr")
    assert [a.element for a in ir.atoms] == ["C", "Cl", "Br"]


# --- Disconnected components ---

def test_disconnected_components(parser: SmilesParser):
    ir = parser.parse("C.C")
    assert ir.has_disconnected is True
    assert len(ir.atoms) == 2
    assert len(ir.bonds) == 0


def test_disconnected_with_charged_ions(parser: SmilesParser):
    ir = parser.parse("[Na+].[Cl-]")
    assert ir.has_disconnected is True
    assert {atom.element for atom in ir.atoms} == {"Na", "Cl"}
    assert ir.total_charge == 0
    charges = sorted(atom.charge for atom in ir.atoms)
    assert charges == [-1, 1]


# --- Convenience API ---

def test_parse_smiles_function():
    ir = parse_smiles("C")
    assert isinstance(ir, MolGraphIR)
    assert len(ir.atoms) == 1
    assert ir.atoms[0].element == "C"
