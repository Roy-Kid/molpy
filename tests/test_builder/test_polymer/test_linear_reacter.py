"""
Tests for linear polymer builder with Reacter-based assembly.

Tests the simplified linear() function that uses ReacterConnector
for chemical reaction-based monomer assembly.
"""

import pytest
from molpy.builder.polymer import linear, ReacterConnector
from molpy.reacter import Reacter
from molpy.reacter.selectors import port_anchor_selector, remove_one_H
from molpy.reacter.transformers import make_single_bond
from molpy.core.wrappers.monomer import Monomer
from molpy.core.atomistic import Atomistic
from molpy.parser.smiles import SmilesParser, bigsmilesir_to_monomer
from molpy.adapter.rdkit_adapter import atomistic_to_mol, mol_to_atomistic
from molpy.io import read_xml_forcefield
from molpy.typifier.atomistic import OplsAtomisticTypifier
from rdkit import Chem
from rdkit.Chem import AllChem


@pytest.fixture
def smiles_parser():
    """Create SMILES parser."""
    return SmilesParser()


@pytest.fixture
def default_reacter():
    """Create default C-C single bond reacter."""
    return Reacter(
        name="C-C",
        anchor_left=port_anchor_selector,
        anchor_right=port_anchor_selector,
        leaving_left=remove_one_H,
        leaving_right=remove_one_H,
        bond_maker=make_single_bond,
    )


@pytest.fixture
def opls_typifier():
    """Create OPLS typifier."""
    ff = read_xml_forcefield("oplsaa.xml")
    return OplsAtomisticTypifier(ff, skip_atom_typing=False)


def create_3d_monomer_from_smiles(smiles: str, parser, label: str = "A"):
    """Helper to create 3D monomer from SMILES."""
    # Parse SMILES
    ir = parser.parse_bigsmiles(smiles)
    monomer = bigsmilesir_to_monomer(ir)
    
    # Convert to RDKit and add hydrogens
    mol = atomistic_to_mol(monomer.unwrap())
    Chem.SanitizeMol(mol)
    molH = Chem.AddHs(mol)
    
    # Generate 3D coordinates
    AllChem.EmbedMolecule(molH, randomSeed=42)
    AllChem.MMFFOptimizeMolecule(molH, maxIters=200)
    
    # Convert back to atomistic
    atomistic_3d = mol_to_atomistic(molH)
    
    # Rebuild monomer with port mapping
    monomer_3d = Monomer(atomistic_3d)
    orig_atoms = list(monomer.unwrap().atoms)
    atoms_3d = list(atomistic_3d.atoms)
    
    for port_name, port in monomer.ports.items():
        old_idx = orig_atoms.index(port.target)
        new_target = atoms_3d[old_idx]
        monomer_3d.define_port(port_name, new_target)
    
    return monomer_3d


class TestLinearBasic:
    """Basic linear assembly tests."""
    
    def test_simple_dimer(self, smiles_parser, default_reacter, opls_typifier):
        """Test simple AB dimer assembly."""
        # Create monomers
        library = {
            "A": create_3d_monomer_from_smiles("CCCCO[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC(C[*:2])O[*:3]", smiles_parser),
        }
        
        # Create connector
        connector = ReacterConnector(
            default=default_reacter,
            port_map={
                ('A', 'B'): ('port_1', 'port_2'),
            },
        )
        
        # Build polymer
        polymer = linear(
            sequence="AB",
            library=library,
            connector=connector,
            typifier=opls_typifier,
            auto_retypify=True,
        )
        
        # Verify structure
        product = polymer.unwrap()
        atoms = list(product.atoms)
        bonds = list(product.bonds)
        
        assert len(atoms) == 25  # A(15) + B(12) - 2(leaving groups)
        assert len(bonds) == 24  # bonds preserved + 1 new bond
        
        # Verify all atoms are typed
        atom_types = [atom.get('type', 'NONE') for atom in atoms]
        typed = [t for t in atom_types if t != 'NONE']
        assert len(typed) == len(atoms), "All atoms should be typed"
        
        # Verify all bonds are typed
        bond_types = [bond.data.get('type', 'NONE') for bond in bonds]
        typed_bonds = [t for t in bond_types if t != 'NONE']
        assert len(typed_bonds) == len(bonds), "All bonds should be typed"
    
    def test_trimer(self, smiles_parser, default_reacter, opls_typifier):
        """Test ABC trimer assembly."""
        # Create monomers
        library = {
            "A": create_3d_monomer_from_smiles("CCCCO[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC(C[*:2])O[*:3]", smiles_parser),
            "C": create_3d_monomer_from_smiles("CCC(C[*:4])O[*:5]", smiles_parser),
        }
        
        # Create connector with port mappings
        connector = ReacterConnector(
            default=default_reacter,
            port_map={
                ('A', 'B'): ('port_1', 'port_2'),
                ('B', 'C'): ('port_3', 'port_4'),
            },
        )
        
        # Build polymer
        polymer = linear(
            sequence="ABC",
            library=library,
            connector=connector,
            typifier=opls_typifier,
            auto_retypify=True,
        )
        
        # Verify structure
        product = polymer.unwrap()
        atoms = list(product.atoms)
        bonds = list(product.bonds)
        
        # A(15) + B(12) + C(15) - 4(leaving groups for 2 connections)
        assert len(atoms) == 38
        # Original bonds + 2 new bonds
        assert len(bonds) == 37
        
        # Verify typing
        atom_types = [atom.get('type', 'NONE') for atom in atoms]
        typed = [t for t in atom_types if t != 'NONE']
        assert len(typed) == len(atoms)


class TestLinearRetypify:
    """Test retypification after assembly."""
    
    def test_retypify_angles_dihedrals(self, smiles_parser, default_reacter, opls_typifier):
        """Test that angles and dihedrals are created and typed."""
        # Create monomers (without pre-typing)
        library = {
            "A": create_3d_monomer_from_smiles("CCCCO[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC(C[*:2])O[*:3]", smiles_parser),
        }
        
        connector = ReacterConnector(
            default=default_reacter,
            port_map={('A', 'B'): ('port_1', 'port_2')},
        )
        
        # Build with retypify
        polymer = linear(
            sequence="AB",
            library=library,
            connector=connector,
            typifier=opls_typifier,
            auto_retypify=True,
        )
        
        product = polymer.unwrap()
        angles = list(product.angles)
        dihedrals = list(product.dihedrals)
        
        # Verify angles and dihedrals exist
        assert len(angles) > 0, "Angles should be created"
        assert len(dihedrals) > 0, "Dihedrals should be created"
        
        # Verify they are typed
        angle_types = [angle.data.get('type', 'NONE') for angle in angles]
        typed_angles = [t for t in angle_types if t != 'NONE']
        assert len(typed_angles) == len(angles), "All angles should be typed"
        
        dihedral_types = [dih.data.get('type', 'NONE') for dih in dihedrals]
        typed_dihedrals = [t for t in dihedral_types if t != 'NONE']
        assert len(typed_dihedrals) == len(dihedrals), "All dihedrals should be typed"
    
    def test_no_retypify(self, smiles_parser, default_reacter):
        """Test assembly without retypification."""
        library = {
            "A": create_3d_monomer_from_smiles("CCCCO[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC(C[*:2])O[*:3]", smiles_parser),
        }
        
        connector = ReacterConnector(
            default=default_reacter,
            port_map={('A', 'B'): ('port_1', 'port_2')},
        )
        
        # Build without retypify
        polymer = linear(
            sequence="AB",
            library=library,
            connector=connector,
            typifier=None,
            auto_retypify=False,
        )
        
        product = polymer.unwrap()
        atoms = list(product.atoms)
        
        # Atoms should not have types
        atom_types = [atom.get('type', 'NONE') for atom in atoms]
        untyped = [t for t in atom_types if t == 'NONE']
        assert len(untyped) == len(atoms), "Atoms should be untyped"


class TestLinearComplex:
    """Test complex sequences."""
    
    def test_repeated_monomers(self, smiles_parser, default_reacter, opls_typifier):
        """Test sequence with repeated monomers (e.g., ABCBD)."""
        library = {
            "A": create_3d_monomer_from_smiles("CCCCO[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC(C[*:2])O[*:3]", smiles_parser),
            "C": create_3d_monomer_from_smiles("CCC(C[*:4])O[*:5]", smiles_parser),
            "D": create_3d_monomer_from_smiles("CCC(C[*:6])O[*:7]", smiles_parser),
        }
        
        connector = ReacterConnector(
            default=default_reacter,
            port_map={
                ('A', 'B'): ('port_1', 'port_2'),
                ('B', 'C'): ('port_3', 'port_4'),
                ('C', 'B'): ('port_5', 'port_2'),
                ('B', 'D'): ('port_3', 'port_6'),
            },
        )
        
        # Build ABCBD
        polymer = linear(
            sequence="ABCBD",
            library=library,
            connector=connector,
            typifier=opls_typifier,
            auto_retypify=True,
        )
        
        product = polymer.unwrap()
        atoms = list(product.atoms)
        bonds = list(product.bonds)
        
        # Verify structure
        assert len(atoms) == 61  # A(15)+B(12)+C(15)+B(12)+D(15) - 8(leaving)
        assert len(bonds) == 60
        
        # Verify all typed
        atom_types = [atom.get('type', 'NONE') for atom in atoms]
        typed = [t for t in atom_types if t != 'NONE']
        assert len(typed) == len(atoms)
    
    def test_homopolymer(self, smiles_parser, default_reacter, opls_typifier):
        """Test homopolymer (AAA)."""
        library = {
            "A": create_3d_monomer_from_smiles("CC(C[*:1])O[*:2]", smiles_parser),
        }
        
        connector = ReacterConnector(
            default=default_reacter,
            port_map={
                ('A', 'A'): ('port_2', 'port_1'),
            },
        )
        
        # Build AAA
        polymer = linear(
            sequence="AAA",
            library=library,
            connector=connector,
            typifier=opls_typifier,
            auto_retypify=True,
        )
        
        product = polymer.unwrap()
        atoms = list(product.atoms)
        bonds = list(product.bonds)
        
        # 3 monomers, 2 connections
        assert len(atoms) == 12 * 3 - 4  # 3 copies - 4 leaving groups
        assert len(bonds) == 11 * 3 - 4 + 2  # bonds + 2 new connections


class TestLinearErrors:
    """Test error handling."""
    
    def test_empty_sequence(self, default_reacter):
        """Test error on empty sequence."""
        from molpy.builder.errors import SequenceError
        
        with pytest.raises(SequenceError, match="at least 2"):
            linear(
                sequence="",
                library={},
                connector=ReacterConnector(
                    default=default_reacter,
                    port_map={},
                ),
            )
    
    def test_missing_monomer(self, default_reacter, smiles_parser):
        """Test error when monomer not in library."""
        from molpy.builder.errors import SequenceError
        
        library = {
            "A": create_3d_monomer_from_smiles("CC[*:1]", smiles_parser),
        }
        
        with pytest.raises(SequenceError, match="not found"):
            linear(
                sequence="AB",  # B not in library
                library=library,
                connector=ReacterConnector(
                    default=default_reacter,
                    port_map={},
                ),
            )
    
    def test_no_port_mapping(self, default_reacter, smiles_parser):
        """Test error when port mapping not found."""
        
        library = {
            "A": create_3d_monomer_from_smiles("CC[*:1]", smiles_parser),
            "B": create_3d_monomer_from_smiles("CC[*:2]", smiles_parser),
        }
        
        # Empty port_map - should fail to find mapping
        connector = ReacterConnector(
            default=default_reacter,
            port_map={},
        )
        
        # ValueError is raised by ReacterConnector when no mapping found
        with pytest.raises(ValueError, match="No port mapping"):
            linear(
                sequence="AB",
                library=library,
                connector=connector,
            )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
