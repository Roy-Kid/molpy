from molpy.parser.big_smiles import GBigSmilesParser, NonCovalentDescriptor


def test_stochastic_object_summary_basic():
    parser = GBigSmilesParser()
    ir = parser.parse("[*]CCC{[<]CC[>]}[*]")

    assert len(ir.stochastic_objects) == 1
    summary = ir.stochastic_objects[0]

    assert summary.left_symbol == "<"
    assert summary.right_symbol == ">"
    assert summary.repeat_units == ["CC"]
    assert summary.has_end_group is False
    assert summary.distribution is None
    assert ir.non_covalent_descriptors == []


def test_non_covalent_descriptor_with_attributes():
    parser = GBigSmilesParser()
    ir = parser.parse("[*]{[>:|1,Type=HB|]C[<]}[*]")

    assert len(ir.non_covalent_descriptors) == 1
    descriptor = ir.non_covalent_descriptors[0]
    assert descriptor == NonCovalentDescriptor(
        symbol=">",
        index=None,
        expression="1",
        attributes={"Type": "HB"},
    )

    summary = ir.stochastic_objects[0]
    assert summary.left_symbol == ">"
    assert summary.right_symbol == "<"
    assert summary.repeat_units == ["C"]
