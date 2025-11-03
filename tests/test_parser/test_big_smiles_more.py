import pytest

from molpy.parser.big_smiles import GBigSmilesParser, NonCovalentDescriptor


@pytest.fixture
def parser():
    return GBigSmilesParser()


def test_multiple_stochastic_objects(parser: GBigSmilesParser):
    ir = parser.parse("[*]{[<]C[>]}{[<]O[>]}[*]")
    assert len(ir.stochastic_objects) == 2
    assert ir.stochastic_objects[0].repeat_units == ["C"]
    assert ir.stochastic_objects[1].repeat_units == ["O"]


def test_terminal_symbols_with_indices(parser: GBigSmilesParser):
    ir = parser.parse("[*]{[>12]CC[<3]}[*]")
    s = ir.stochastic_objects[0]
    assert s.left_symbol == ">"
    assert s.right_symbol == "<"
    assert s.repeat_units == ["CC"]


def test_ladder_descriptor_with_index(parser: GBigSmilesParser):
    ir = parser.parse("[*]{[>[<]2]C[<]}[*]")
    s = ir.stochastic_objects[0]
    assert s.left_symbol == ">"
    assert s.right_symbol == "<"
    assert s.repeat_units == ["C"]


def test_end_group_and_monomer_list(parser: GBigSmilesParser):
    ir = parser.parse("[*]{[<]C,O;N,O[>]}[*]")
    s = ir.stochastic_objects[0]
    # Two monomers before ';', two after
    assert s.repeat_units == ["C", "O", "N", "O"]


def test_stochastic_distribution_variants(parser: GBigSmilesParser):
    cases = [
        "[*]{[<]C[>]}|flory_schulz(0.1,)|[*]",
        "[*]{[<]C[>]}|schulz_zimm(1, 2)|[*]",
        "[*]{[<]C[>]}|gauss(1.0, 0.2)|[*]",
        "[*]{[<]C[>]}|uniform(0, 10)|[*]",
        "[*]{[<]C[>]}|log_normal(1, 0.25)|[*]",
        "[*]{[<]C[>]}|poisson(3,)|[*]",
    ]
    for s in cases:
        ir = parser.parse(s)
        assert len(ir.stochastic_objects) == 1
        assert ir.stochastic_objects[0].repeat_units == ["C"]
        assert ir.stochastic_objects[0].distribution is not None


def test_non_covalent_multiple_attributes_and_index_expr(parser: GBigSmilesParser):
    # Context with expression and two attributes
    ir = parser.parse("[*]{[>:|!3&5,Type=HB,Strength=weak|]C[<]}[*]")
    assert len(ir.non_covalent_descriptors) == 1
    d = ir.non_covalent_descriptors[0]
    assert d == NonCovalentDescriptor(
        symbol=">",
        index=None,
        expression="!3&5",
        attributes={"Type": "HB", "Strength": "weak"},
    )
    s = ir.stochastic_objects[0]
    assert s.left_symbol == ">"
    assert s.right_symbol == "<"


def test_terminal_bond_generation_with_spaces(parser: GBigSmilesParser):
    # Generation in terminal descriptor shouldn't break symbol extraction
    ir = parser.parse("[*]{[> | 1 | ]C[<]}[*]")
    s = ir.stochastic_objects[0]
    assert s.left_symbol == ">"
    assert s.right_symbol == "<"


def test_dot_generation_at_molecule_end(parser: GBigSmilesParser):
    # ".|10|" after a molecule should parse without affecting stochastic summary
    ir = parser.parse("[*]{[<]C[>]}[*].|10|")
    assert len(ir.stochastic_objects) == 1
    assert ir.stochastic_objects[0].repeat_units == ["C"]


def test_repeat_unit_can_be_big_smiles_molecule(parser: GBigSmilesParser):
    # A repeat unit that itself contains a simple smiles + another stochastic block
    text = "[*]{[<]C{[<]O[>]}[>]}[*]"
    ir = parser.parse(text)
    s = ir.stochastic_objects[0]
    # The outer repeat unit is parsed as a big_smiles_molecule; we record the text slice
    assert any("{[<]O[>]}" in unit for unit in s.repeat_units)


def test_fragment_definition_after_molecule(parser: GBigSmilesParser):
    # Fragment definition at the end should be accepted by the grammar
    text = "[*]{[<]C[>]}[*] . {#Frag=[*]O[*]}"
    ir = parser.parse(text)
    assert len(ir.stochastic_objects) == 1
    assert ir.stochastic_objects[0].repeat_units == ["C"]
