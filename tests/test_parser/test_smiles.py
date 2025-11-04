import pytest

from molpy.core.aa import Atom
from molpy.core.aa.base import Atomistic
from molpy.core.wrappers.monomer import Monomer
from molpy.parser.convert import (
    GBigSmilesResult,
    ir_to_atomistic,
    ir_to_bigsmiles,
    ir_to_gbigsmiles,
)
from molpy.parser.smiles_ir import IRBigConnector, IRBigIndexExpr, IRBigStochasticBlock
from molpy.parser.smiles import parse_smiles, parse_bigsmiles, parse_gbigsmiles


class TestSmilesCore:
    def test_simple_chain(self) -> None:
        doc = parse_smiles("CCO")
        assert doc.features == {"smiles"}
        structure = ir_to_atomistic(doc)
        assert isinstance(structure, Atomistic)
        assert len(structure.atoms) == 3
        assert len(structure.bonds) == 2
        assert set(structure.symbols) == {"C", "C", "O"}
        for atom in structure.atoms:
            assert atom.get("span") is not None

    def test_branch_and_ring(self) -> None:
        doc = parse_smiles("C1CCCCC1")
        assert "smiles" in doc.features
        structure = ir_to_atomistic(doc)
        assert len(structure.atoms) == 6
        assert len(structure.bonds) == 6
        first = structure.atoms[0]
        last = structure.atoms[-1]
        # ensure ring closure present
        assert any(
            (bond.itom is first and bond.jtom is last)
            or (bond.itom is last and bond.jtom is first)
            for bond in structure.bonds
        )

    def test_bracket_props(self) -> None:
        doc = parse_smiles("[13C@H-:12]")
        structure = ir_to_atomistic(doc)
        atom = structure.atoms[0]
        assert atom.get("kind") == "bracket"
        assert atom.get("isotope") == 13
        assert atom.get("chiral") == "@"
        assert atom.get("hcount") == 1
        assert atom.get("charge") == -1
        assert atom.get("element") == "C"


class TestBigSmiles:
    def test_stochastic_basic(self) -> None:
        doc = parse_bigsmiles("{[$]C,C[$]}")
        assert {"smiles", "bigsmiles"} <= doc.features
        monomers, sequence = ir_to_bigsmiles(doc, sequence_length=6, seed=5)
        assert set(monomers.keys()) == {"A", "B"}
        assert len(sequence) == 6
        assert set(sequence).issubset({"A", "B"})
        for key, monomer in monomers.items():
            assert isinstance(monomer, Monomer)
            assembly = monomer.unwrap()
            atoms = list(assembly.entities.bucket(Atom))
            assert len(atoms) == 1
            symbol = atoms[0].get("symbol")
            if key == "A":
                assert symbol == "C"
            else:
                assert symbol == "C"

    def test_terminal_and_units(self) -> None:
        doc = parse_bigsmiles("{[<]C;O[>]}")
        monomers, sequence = ir_to_bigsmiles(doc, sequence_length=4, seed=7)
        assert {"A", "B"} <= set(monomers.keys())
        assert len(sequence) == 4
        assert set(sequence).issubset(set(monomers.keys()))
        unit_atom = list(monomers["A"].unwrap().entities.bucket(Atom))[0]
        assert unit_atom.get("symbol") == "C"
        end_atom = list(monomers["B"].unwrap().entities.bucket(Atom))[0]
        assert end_atom.get("symbol") == "O"

    def test_fragment_decl_def(self) -> None:
        doc = parse_bigsmiles("[#R].{#R=C.C}")
        assert "bigsmiles" in doc.features
        monomers, sequence = ir_to_bigsmiles(doc, sequence_length=3, seed=11)
        assert "R" in monomers
        atoms = list(monomers["R"].unwrap().entities.bucket(Atom))
        assert len(atoms) == 2
        assert len(sequence) == 3
        assert set(sequence) == {"R"}


class TestGBigSmiles:
    def test_non_covalent_descriptor(self) -> None:
        doc = parse_gbigsmiles("C[::1|!1,mode=hb|]C")
        assert {"gbigsmiles", "bigsmiles"} <= doc.features
        descriptor = next(part for part in doc.molecules[0].parts if isinstance(part, IRBigConnector))
        assert descriptor.connector == ":"
        assert descriptor.index == 1
        assert descriptor.noncovalent is not None
        assert descriptor.noncovalent.label == 1
        ctx = descriptor.noncovalent.context
        assert ctx is not None
        assert ctx.expr.unary == "!"
        assert ctx.expr.left == 1
        assert ctx.kv == {"mode": "hb"}
        assert descriptor.span is not None
        result = ir_to_gbigsmiles(doc, seed=3)
        assert isinstance(result, GBigSmilesResult)
        assert not result.monomers
        assert len(result.molecules) == 1
        molecule = result.molecules[0]
        assert molecule.system_size is None
        assert not molecule.sequences
        assert len(molecule.bond_descriptors) == 1
        summary = molecule.bond_descriptors[0]
        assert summary is not None
        assert summary["connector"] == ":"
        assert summary["index"] == 1
        noncovalent = summary["noncovalent"]
        assert noncovalent["label"] == 1
        ctx_summary = noncovalent["context"]
        assert ctx_summary["expr"]["unary"] == "!"
        assert ctx_summary["expr"]["left"] == 1
        assert ctx_summary["kv"] == {"mode": "hb"}

    def test_index_expr_nested(self) -> None:
        doc = parse_gbigsmiles("C[::|(!1~2)&3|]C")
        descriptor = next(part for part in doc.molecules[0].parts if isinstance(part, IRBigConnector))
        expr = descriptor.noncovalent.context.expr
        assert expr.op == "&"
        assert expr.right == 3
        assert isinstance(expr.left, IRBigIndexExpr)
        left = expr.left
        assert left.op == "~"
        assert isinstance(left.left, IRBigIndexExpr)
        assert left.left.unary == "!"
        assert left.left.left == 1
        assert left.right == 2
        result = ir_to_gbigsmiles(doc, seed=5)
        molecule = result.molecules[0]
        summary = molecule.bond_descriptors[0]
        noncovalent = summary["noncovalent"]
        expr_summary = noncovalent["context"]["expr"]
        assert expr_summary["op"] == "&"
        assert expr_summary["right"] == 3
        left_summary = expr_summary["left"]
        assert left_summary["op"] == "~"
        assert left_summary["left"]["unary"] == "!"
        assert left_summary["left"]["left"] == 1
        assert left_summary["right"] == 2

    def test_distribution(self) -> None:
        doc = parse_gbigsmiles("{[$]C[$]}|flory_schulz(0.5)|")
        obj = next(part for part in doc.molecules[0].parts if isinstance(part, IRBigStochasticBlock))
        assert obj.distribution is not None
        assert obj.distribution.kind == "flory_schulz"
        assert obj.distribution.p1 == 0.5
        assert obj.distribution.span is not None
        result = ir_to_gbigsmiles(doc, seed=11)
        assert isinstance(result, GBigSmilesResult)
        assert set(result.monomers.keys()) == {"A"}
        molecule = result.molecules[0]
        assert molecule.system_size is None
        assert len(molecule.sequences) == 1
        sequence = molecule.sequences[0]
        assert sequence.distribution is not None
        assert sequence.distribution.kind == "flory_schulz"
        assert sequence.unit_labels == ["A"]
        assert sequence.unit_weights == [1.0]
        assert sequence.target_total is None
        assert sequence.chains
        chain = sequence.chains[0]
        assert chain.degree >= 1
        assert chain.labels == ["A"] * chain.degree

    def test_distribution1(self) -> None:
        doc = parse_gbigsmiles("CCOC(=O)C(C)(C){[>][<]CC([>])c1ccccc1, [<]CC([>])C(=O)OC [<]}|schulz_zimm(1500, 1400)|[Br].|5e5|")
        result = ir_to_gbigsmiles(doc, seed=17)
        molecule = result.molecules[0]
        assert molecule.system_size == 500000
        assert len(molecule.sequences) == 1
        sequence = molecule.sequences[0]
        assert sequence.distribution is not None
        assert sequence.distribution.kind == "schulz_zimm"
        assert set(sequence.unit_labels) == {"A", "B"}
        total_degree = sum(chain.degree for chain in sequence.chains)
        assert total_degree == 500000
        labels_used = {label for chain in sequence.chains for label in chain.labels}
        assert labels_used <= set(sequence.unit_labels)
