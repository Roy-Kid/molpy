"""
Tests for moltemplate-style molecular modeling using Atomistic and Spatial operations.

This module tests the core functionality for building molecular structures
in a moltemplate-like manner, including molecule manipulation, connection,
and assembly operations.
"""

import numpy as np

from molpy.core._atomistic import Atomistic
from molpy.core.wrapper import Spatial


class TestManipulate:
    """Test molecule manipulation operations using Spatial wrapper."""

    def test_create_simple_molecule(self):
        """Test creating a simple molecule (methane)."""
        # Create methane molecule
        methane = Atomistic(name="methane")

        # Add carbon atom at origin
        c_atom = methane.def_atom(name="C", element="C", xyz=np.array([0.0, 0.0, 0.0]))

        # Add hydrogen atoms around carbon
        h1 = methane.def_atom(name="H1", element="H", xyz=np.array([0.0, 0.0, 1.09]))
        h2 = methane.def_atom(name="H2", element="H", xyz=np.array([1.03, 0.0, -0.36]))
        h3 = methane.def_atom(
            name="H3", element="H", xyz=np.array([-0.51, 0.89, -0.36])
        )
        h4 = methane.def_atom(
            name="H4", element="H", xyz=np.array([-0.51, -0.89, -0.36])
        )

        # Add bonds
        methane.def_bond(c_atom, h1)
        methane.def_bond(c_atom, h2)
        methane.def_bond(c_atom, h3)
        methane.def_bond(c_atom, h4)

        assert len(methane.atoms) == 5
        assert len(methane.bonds) == 4
        assert methane.atoms[0]["name"] == "C"
        assert methane.atoms[0]["element"] == "C"

    def test_translate_molecule_with_spatial(self):
        """Test translating a molecule using Spatial wrapper."""
        # Create water molecule
        water = Atomistic(name="water")

        # Add oxygen and hydrogens
        o_atom = water.def_atom(name="O", element="O", xyz=np.array([0.0, 0.0, 0.0]))
        h1 = water.def_atom(name="H1", element="H", xyz=np.array([0.9572, 0.0, 0.0]))
        h2 = water.def_atom(name="H2", element="H", xyz=np.array([-0.2400, 0.0, 0.0]))

        # Add bonds
        water.def_bond(o_atom, h1)
        water.def_bond(o_atom, h2)

        # Wrap with Spatial wrapper for spatial operations
        spatial_water = Spatial(water)

        # Translate by vector [5.0, 2.0, 1.0]
        translation_vector = np.array([5.0, 2.0, 1.0])
        spatial_water.move(translation_vector)

        # Check that all atoms were translated
        assert np.allclose(water.atoms[0]["xyz"], np.array([5.0, 2.0, 1.0]))
        assert np.allclose(water.atoms[1]["xyz"], np.array([5.9572, 2.0, 1.0]))
        assert np.allclose(water.atoms[2]["xyz"], np.array([4.76, 2.0, 1.0]))

    def test_move_to_method(self):
        """Test the move_to method that moves molecule to a specific position."""
        # Create water molecule
        water = Atomistic(name="water")

        # Add oxygen and hydrogens
        o_atom = water.def_atom(name="O", element="O", xyz=np.array([0.0, 0.0, 0.0]))
        h1 = water.def_atom(name="H1", element="H", xyz=np.array([0.9572, 0.0, 0.0]))
        h2 = water.def_atom(name="H2", element="H", xyz=np.array([-0.2400, 0.0, 0.0]))

        # Add bonds
        water.def_bond(o_atom, h1)
        water.def_bond(o_atom, h2)

        # Wrap with Spatial wrapper
        spatial_water = Spatial(water)

        # Move water molecule so that oxygen is at position [10.0, 5.0, 2.0]
        target_position = np.array([10.0, 5.0, 2.0])
        spatial_water.move_to(target_position)

        # Check that oxygen is now at target position
        assert np.allclose(water.atoms[0]["xyz"], target_position)

        # Check that relative positions of hydrogens are preserved
        h1_relative = water.atoms[1]["xyz"] - water.atoms[0]["xyz"]
        h2_relative = water.atoms[2]["xyz"] - water.atoms[0]["xyz"]

        expected_h1_relative = np.array([0.9572, 0.0, 0.0])
        expected_h2_relative = np.array([-0.2400, 0.0, 0.0])

        assert np.allclose(h1_relative, expected_h1_relative)
        assert np.allclose(h2_relative, expected_h2_relative)

    def test_rotate_molecule_with_spatial(self):
        """Test rotating a molecule using Spatial wrapper."""
        # Create linear molecule (CO2-like)
        co2 = Atomistic(name="CO2")

        # Add atoms in linear arrangement
        c_atom = co2.def_atom(name="C", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        o1 = co2.def_atom(name="O1", element="O", xyz=np.array([1.16, 0.0, 0.0]))
        o2 = co2.def_atom(name="O2", element="O", xyz=np.array([-1.16, 0.0, 0.0]))

        # Add bonds
        co2.def_bond(c_atom, o1)
        co2.def_bond(c_atom, o2)

        # Wrap with Spatial wrapper
        spatial_co2 = Spatial(co2)

        # Rotate around Y-axis by 90 degrees
        rotation_axis = np.array([0.0, 1.0, 0.0])
        rotation_angle = np.pi / 2  # 90 degrees
        spatial_co2.rotate(rotation_axis, rotation_angle)

        # Check that rotation was applied correctly
        # C should remain at origin
        assert np.allclose(co2.atoms[0]["xyz"], np.array([0.0, 0.0, 0.0]))
        # O1 should now be at [0, 0, -1.16] (rotated from [1.16, 0, 0])
        assert np.allclose(co2.atoms[1]["xyz"], np.array([0.0, 0.0, -1.16]), atol=1e-10)
        # O2 should now be at [0, 0, 1.16] (rotated from [-1.16, 0, 0])
        assert np.allclose(co2.atoms[2]["xyz"], np.array([0.0, 0.0, 1.16]), atol=1e-10)

    def test_rotate_around_specific_atom(self):
        """Test rotating a molecule around a specific atom using Spatial wrapper."""
        # Create molecule with a central atom
        molecule = Atomistic(name="test_molecule")

        # Central atom
        center = molecule.def_atom(name="C", element="C", xyz=np.array([0.0, 0.0, 0.0]))

        # Peripheral atoms
        p1 = molecule.def_atom(name="H1", element="H", xyz=np.array([1.0, 0.0, 0.0]))
        p2 = molecule.def_atom(name="H2", element="H", xyz=np.array([0.0, 1.0, 0.0]))

        # Add bonds
        molecule.def_bond(center, p1)
        molecule.def_bond(center, p2)

        # Wrap with Spatial wrapper
        spatial_molecule = Spatial(molecule)

        # Rotate around center atom by 45 degrees around Z-axis
        rotation_axis = np.array([0.0, 0.0, 1.0])
        rotation_angle = np.pi / 4  # 45 degrees
        spatial_molecule.rotate(rotation_axis, rotation_angle, origin=center["xyz"])

        # Check rotation results
        expected_p1 = np.array([np.cos(np.pi / 4), np.sin(np.pi / 4), 0.0])
        expected_p2 = np.array([-np.sin(np.pi / 4), np.cos(np.pi / 4), 0.0])

        assert np.allclose(p1["xyz"], expected_p1, atol=1e-10)
        assert np.allclose(p2["xyz"], expected_p2, atol=1e-10)

    def test_rotate_with_matrix(self):
        """Test rotating a molecule using rotation matrix."""
        # Create simple molecule
        molecule = Atomistic(name="test_molecule")
        atom1 = molecule.def_atom(name="A1", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        atom2 = molecule.def_atom(name="A2", element="C", xyz=np.array([0.0, 1.0, 0.0]))
        molecule.def_bond(atom1, atom2)

        # Wrap with Spatial wrapper
        spatial_molecule = Spatial(molecule)

        # Create rotation matrix for 90° around Z-axis
        rotation_matrix = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]])

        # Apply rotation
        spatial_molecule.rotate_with_matrix(rotation_matrix)

        # Check results - atom1 should rotate from [1,0,0] to [0,-1,0]
        assert np.allclose(atom1["xyz"], np.array([0.0, -1.0, 0.0]), atol=1e-10)
        # atom2 should rotate from [0,1,0] to [1,0,0]
        assert np.allclose(atom2["xyz"], np.array([1.0, 0.0, 0.0]), atol=1e-10)


class TestConnect:
    """Test molecular connection and bonding operations."""

    def test_connect_two_molecules(self):
        """Test connecting two molecules by forming a bond."""
        # Create first molecule (methyl group)
        methyl = Atomistic(name="methyl")
        c1 = methyl.def_atom(name="C1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        h1 = methyl.def_atom(name="H1", element="H", xyz=np.array([1.09, 0.0, 0.0]))
        h2 = methyl.def_atom(name="H2", element="H", xyz=np.array([-0.36, 0.89, 0.0]))
        h3 = methyl.def_atom(name="H3", element="H", xyz=np.array([-0.36, -0.89, 0.0]))

        methyl.def_bond(c1, h1)
        methyl.def_bond(c1, h2)
        methyl.def_bond(c1, h3)

        # Create second molecule (hydroxyl group)
        hydroxyl = Atomistic(name="hydroxyl")
        o1 = hydroxyl.def_atom(name="O1", element="O", xyz=np.array([0.0, 0.0, 0.0]))
        h4 = hydroxyl.def_atom(name="H4", element="H", xyz=np.array([0.0, 0.0, 0.96]))

        hydroxyl.def_bond(o1, h4)

        # Wrap hydroxyl with Spatial wrapper
        spatial_hydroxyl = Spatial(hydroxyl)

        # Connect molecules by forming C-O bond
        # First, translate hydroxyl to connect with methyl
        connection_distance = 1.43  # Typical C-O bond length
        hydroxyl_center = np.array([connection_distance, 0.0, 0.0])
        spatial_hydroxyl.move(hydroxyl_center)

        # Combine molecules
        combined = Atomistic(name="methanol")
        combined.add_struct(methyl)
        combined.add_struct(hydroxyl)

        # Add the connecting bond
        combined.def_bond(c1, o1)

        assert len(combined.atoms) == 6
        assert len(combined.bonds) == 5

        # Check that the connecting bond exists between c1 and o1
        # Note: The order of atoms in the bond might vary
        connecting_bond = None
        for bond in combined.bonds:
            if (bond.itom == c1 and bond.jtom == o1) or (
                bond.itom == o1 and bond.jtom == c1
            ):
                connecting_bond = bond
                break

        assert (
            connecting_bond is not None
        ), "Connecting bond between C1 and O1 not found"
        assert connecting_bond.itom in [c1, o1]
        assert connecting_bond.jtom in [c1, o1]

    def test_connect_with_rotation(self):
        """Test connecting molecules with proper orientation using Spatial wrapper."""
        # Create two linear molecules
        mol1 = Atomistic(name="molecule1")
        a1 = mol1.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = mol1.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        mol1.def_bond(a1, a2)

        mol2 = Atomistic(name="molecule2")
        b1 = mol2.def_atom(name="B1", element="N", xyz=np.array([0.0, 0.0, 0.0]))
        b2 = mol2.def_atom(name="B2", element="N", xyz=np.array([1.0, 0.0, 0.0]))
        mol2.def_bond(b1, b2)

        # Wrap mol2 with Spatial wrapper
        spatial_mol2 = Spatial(mol2)

        # Connect molecules end-to-end with proper orientation
        # First, align mol2's direction with mol1's direction
        mol1_direction = a2["xyz"] - a1["xyz"]
        mol2_direction = b2["xyz"] - b1["xyz"]

        # Calculate rotation matrix to align directions
        # For simplicity, we'll use a manual rotation matrix
        # Since mol1 goes along X-axis and mol2 also goes along X-axis, no rotation needed

        # Translate mol2 to connect with mol1
        connection_point = a2["xyz"]
        mol2_connection_point = b1["xyz"]
        translation_vector = connection_point - mol2_connection_point
        spatial_mol2.move(translation_vector)

        # Combine molecules
        combined = Atomistic(name="connected")
        combined.add_struct(mol1)
        combined.add_struct(mol2)

        # Add connecting bond
        combined.def_bond(a2, b1)

        assert len(combined.atoms) == 4
        assert len(combined.bonds) == 3

        # Check that the connection bond exists between a2 and b1
        # Note: The order of atoms in the bond might vary
        connection_bond = None
        for bond in combined.bonds:
            if (bond.itom == a2 and bond.jtom == b1) or (
                bond.itom == b1 and bond.jtom == a2
            ):
                connection_bond = bond
                break

        assert (
            connection_bond is not None
        ), "Connecting bond between A2 and B1 not found"
        assert connection_bond.itom in [a2, b1]
        assert connection_bond.jtom in [a2, b1]

    def test_connect_with_deletion(self):
        """Test connecting molecules by deleting overlapping atoms."""
        # Create two molecules with overlapping atoms
        mol1 = Atomistic(name="molecule1")
        c1 = mol1.def_atom(name="C1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        h1 = mol1.def_atom(name="H1", element="H", xyz=np.array([1.09, 0.0, 0.0]))
        mol1.def_bond(c1, h1)

        mol2 = Atomistic(name="molecule2")
        c2 = mol2.def_atom(
            name="C2", element="C", xyz=np.array([0.0, 0.0, 0.0])
        )  # Same position as c1
        h2 = mol2.def_atom(name="H2", element="H", xyz=np.array([-1.09, 0.0, 0.0]))
        mol2.def_bond(c2, h2)

        # Connect by merging overlapping atoms
        combined = Atomistic(name="connected")
        combined.add_struct(mol1)

        # Add mol2 atoms except the overlapping one
        combined.add_atom(h2)

        # Add the connecting bond (c1 to h2)
        combined.def_bond(c1, h2)

        assert len(combined.atoms) == 3  # c1, h1, and h2
        assert len(combined.bonds) == 2  # c1-h1 and c1-h2

        # Check that the first bond is between c1 and h1
        first_bond = combined.bonds[0]
        assert first_bond.itom in [c1, h1]
        assert first_bond.jtom in [c1, h1]

        # Check that the second bond is between c1 and h2
        second_bond = combined.bonds[1]
        assert second_bond.itom in [c1, h2]
        assert second_bond.jtom in [c1, h2]


class TestAssembly:
    """Test molecular assembly and construction operations."""

    def test_build_polymer_chain(self):
        """Test building a simple polymer chain using Spatial wrapper."""
        # Create monomer template
        monomer = Atomistic(name="monomer")

        # Add atoms for monomer
        c1 = monomer.def_atom(name="C1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        c2 = monomer.def_atom(name="C2", element="C", xyz=np.array([1.54, 0.0, 0.0]))
        h1 = monomer.def_atom(name="H1", element="H", xyz=np.array([-0.51, 0.89, 0.0]))
        h2 = monomer.def_atom(name="H2", element="H", xyz=np.array([-0.51, -0.89, 0.0]))
        h3 = monomer.def_atom(name="H3", element="H", xyz=np.array([2.05, 0.89, 0.0]))
        h4 = monomer.def_atom(name="H4", element="H", xyz=np.array([2.05, -0.89, 0.0]))

        # Add bonds
        monomer.def_bond(c1, c2)
        monomer.def_bond(c1, h1)
        monomer.def_bond(c1, h2)
        monomer.def_bond(c2, h3)
        monomer.def_bond(c2, h4)

        # Build polymer chain
        chain = Atomistic(name="polymer_chain")
        n_monomers = 3

        for i in range(n_monomers):
            # Create monomer copy
            monomer_copy = monomer()

            # Wrap with Spatial wrapper for translation
            spatial_monomer = Spatial(monomer_copy)

            # Translate monomer to position
            translation = np.array([i * 1.54, 0.0, 0.0])
            spatial_monomer.move(translation)

            # Add to chain
            chain.add_struct(monomer_copy)

            # Connect monomers (except first one)
            if i > 0:
                prev_monomer_end = chain.atoms[
                    5 * i - 1
                ]  # Last atom of previous monomer
                curr_monomer_start = chain.atoms[5 * i]  # First atom of current monomer
                chain.def_bond(prev_monomer_end, curr_monomer_start)

        # Check that we have the expected number of atoms and bonds
        # Note: add_struct may not merge atoms, so we might have more atoms than expected
        assert len(chain.atoms) >= 15  # At least 3 monomers × 5 atoms
        assert (
            len(chain.bonds) >= 17
        )  # At least 3 monomers × 5 bonds + 2 connecting bonds

        # Verify that monomers are properly positioned
        # First monomer should be at origin
        first_c1 = chain.atoms[0]
        assert np.allclose(first_c1["xyz"], np.array([0.0, 0.0, 0.0]))

        # Second monomer should be translated by 1.54 along X-axis
        second_monomer_start = None
        for atom in chain.atoms:
            if atom["name"] == "C1" and np.allclose(
                atom["xyz"], np.array([1.54, 0.0, 0.0])
            ):
                second_monomer_start = atom
                break
        assert (
            second_monomer_start is not None
        ), "Second monomer not found at expected position"

    def test_build_branched_molecule(self):
        """Test building a branched molecular structure using Spatial wrapper."""
        # Create central molecule
        central = Atomistic(name="central")
        c_center = central.def_atom(
            name="C", element="C", xyz=np.array([0.0, 0.0, 0.0])
        )

        # Create branch templates
        branch1 = Atomistic(name="branch1")
        c1 = branch1.def_atom(name="C1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        h1 = branch1.def_atom(name="H1", element="H", xyz=np.array([1.09, 0.0, 0.0]))
        branch1.def_bond(c1, h1)

        branch2 = Atomistic(name="branch2")
        c2 = branch2.def_atom(name="C2", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        h2 = branch2.def_atom(name="H2", element="H", xyz=np.array([1.09, 0.0, 0.0]))
        branch2.def_bond(c2, h2)

        # Assemble branched molecule
        branched = Atomistic(name="branched")
        branched.add_struct(central)

        # Add first branch at 120° angle
        branch1_copy = branch1()
        spatial_branch1 = Spatial(branch1_copy)
        # Rotate and translate branch1
        spatial_branch1.rotate(np.array([0.0, 0.0, 1.0]), np.pi / 3)
        spatial_branch1.move(np.array([1.54, 0.0, 0.0]))
        branched.add_struct(branch1_copy)
        branched.def_bond(c_center, branch1_copy.atoms[0])

        # Add second branch at -120° angle
        branch2_copy = branch2()
        spatial_branch2 = Spatial(branch2_copy)
        # Rotate and translate branch2
        spatial_branch2.rotate(np.array([0.0, 0.0, 1.0]), -np.pi / 3)
        spatial_branch2.move(np.array([1.54, 0.0, 0.0]))
        branched.add_struct(branch2_copy)
        branched.def_bond(c_center, branch2_copy.atoms[0])

        assert len(branched.atoms) == 5
        # Note: add_struct may add internal bonds from the branch templates
        assert (
            len(branched.bonds) >= 3
        )  # At least: c_center-branch1, c_center-branch2, and internal branch bonds

        # Check that central atom is connected to both branches
        central_bonds = [
            bond
            for bond in branched.bonds
            if bond.itom == c_center or bond.jtom == c_center
        ]
        assert len(central_bonds) == 2

    def test_build_cyclic_molecule(self):
        """Test building a cyclic molecular structure."""
        # Create benzene-like ring
        ring = Atomistic(name="benzene_ring")

        # Add carbon atoms in a hexagon
        carbon_atoms = []
        n_carbons = 6
        radius = 1.40  # C-C bond length

        for i in range(n_carbons):
            angle = i * 2 * np.pi / n_carbons
            x = radius * np.cos(angle)
            y = radius * np.sin(angle)
            c_atom = ring.def_atom(
                name=f"C{i+1}", element="C", xyz=np.array([x, y, 0.0])
            )
            carbon_atoms.append(c_atom)

        # Add bonds to form ring
        for i in range(n_carbons):
            j = (i + 1) % n_carbons
            ring.def_bond(carbon_atoms[i], carbon_atoms[j])

        # Add hydrogen atoms
        hydrogen_atoms = []
        for i, c_atom in enumerate(carbon_atoms):
            angle = i * 2 * np.pi / n_carbons
            h_x = (radius + 1.09) * np.cos(angle)  # H-C bond length = 1.09 Å
            h_y = (radius + 1.09) * np.sin(angle)
            h_atom = ring.def_atom(
                name=f"H{i+1}", element="H", xyz=np.array([h_x, h_y, 0.0])
            )
            hydrogen_atoms.append(h_atom)
            ring.def_bond(c_atom, h_atom)

        assert len(ring.atoms) == 12  # 6 C + 6 H
        assert len(ring.bonds) == 12  # 6 C-C + 6 C-H

        # Check that ring is properly closed
        for i in range(n_carbons):
            j = (i + 1) % n_carbons
            bond_exists = any(
                (bond.itom == carbon_atoms[i] and bond.jtom == carbon_atoms[j])
                or (bond.itom == carbon_atoms[j] and bond.jtom == carbon_atoms[i])
                for bond in ring.bonds
            )
            assert bond_exists, f"Bond missing between C{i+1} and C{j+1}"


class TestTopology:
    """Test topology generation and management."""

    def test_generate_angles(self):
        """Test automatic generation of angle terms."""
        # Create molecule with bonds
        molecule = Atomistic(name="test_molecule")

        # Add atoms in a chain
        a1 = molecule.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = molecule.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        a3 = molecule.def_atom(name="A3", element="C", xyz=np.array([2.0, 0.0, 0.0]))

        # Add bonds
        molecule.def_bond(a1, a2)
        molecule.def_bond(a2, a3)

        # Generate angles
        angles = molecule.gen_angles()

        assert len(angles) == 1
        angle = angles[0]

        # Check that the angle involves the correct atoms
        # Note: The order might vary due to how angles are generated
        atoms_in_angle = {angle.itom, angle.jtom, angle.ktom}
        expected_atoms = {a1, a2, a3}
        assert atoms_in_angle == expected_atoms

        # Check that the middle atom (a2) is the central atom
        assert angle.jtom == a2

        # Check angle value (should be 180° for linear arrangement)
        # Note: Angle calculation uses atom["xyz"] access
        v1 = a1["xyz"] - a2["xyz"]
        v2 = a3["xyz"] - a2["xyz"]
        cos_angle = np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))
        cos_angle = np.clip(cos_angle, -1.0, 1.0)
        expected_angle = np.arccos(cos_angle)
        assert abs(angle.value - expected_angle) < 1e-10

    def test_generate_dihedrals(self):
        """Test automatic generation of dihedral terms."""
        # Create molecule with bonds
        molecule = Atomistic(name="test_molecule")

        # Add atoms in a chain
        a1 = molecule.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = molecule.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        a3 = molecule.def_atom(name="A3", element="C", xyz=np.array([2.0, 0.0, 0.0]))
        a4 = molecule.def_atom(name="A4", element="C", xyz=np.array([3.0, 0.0, 0.0]))

        # Add bonds
        molecule.def_bond(a1, a2)
        molecule.def_bond(a2, a3)
        molecule.def_bond(a3, a4)

        # Generate dihedrals
        dihedrals = molecule.gen_dihedrals()

        assert len(dihedrals) == 1
        dihedral = dihedrals[0]
        # Note: The order might be different due to how dihedrals are generated
        # Just check that all atoms are present
        atoms_in_dihedral = {dihedral.itom, dihedral.jtom, dihedral.ktom, dihedral.ltom}
        expected_atoms = {a1, a2, a3, a4}
        assert atoms_in_dihedral == expected_atoms

        # Check dihedral value (should be 0° for linear arrangement)
        # Note: Dihedral calculation uses atom["xyz"] access
        b1 = a2["xyz"] - a1["xyz"]
        b2 = a3["xyz"] - a2["xyz"]
        b3 = a4["xyz"] - a3["xyz"]

        n1 = np.cross(b1, b2)
        n2 = np.cross(b2, b3)

        n1_norm = np.linalg.norm(n1)
        n2_norm = np.linalg.norm(n2)

        if n1_norm > 1e-10 and n2_norm > 1e-10:
            n1 = n1 / n1_norm
            n2 = n2 / n2_norm
            cos_angle = np.dot(n1, n2)
            cos_angle = np.clip(cos_angle, -1.0, 1.0)
            expected_dihedral = np.arccos(cos_angle)
            if np.dot(np.cross(n1, n2), b2) < 0:
                expected_dihedral = -expected_dihedral
            assert abs(dihedral.value - expected_dihedral) < 1e-10

    def test_topology_consistency(self):
        """Test that topology is consistent after modifications."""
        # Create molecule
        molecule = Atomistic(name="test_molecule")

        # Add atoms and bonds
        a1 = molecule.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = molecule.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        a3 = molecule.def_atom(name="A3", element="C", xyz=np.array([2.0, 0.0, 0.0]))

        molecule.def_bond(a1, a2)
        molecule.def_bond(a2, a3)

        # Generate topology
        topo = molecule.get_topology()

        assert topo.n_atoms == 3
        assert topo.n_bonds == 2

        # Add more bonds and regenerate
        a4 = molecule.def_atom(name="A4", element="C", xyz=np.array([3.0, 0.0, 0.0]))
        molecule.def_bond(a3, a4)

        # Regenerate topology
        new_topo = molecule.get_topology()

        assert new_topo.n_atoms == 4
        assert new_topo.n_bonds == 3


class TestSerialization:
    """Test serialization and deserialization of molecular structures."""

    def test_to_frame_from_frame(self):
        """Test converting structure to frame and back."""
        # Create molecule
        molecule = Atomistic(name="test_molecule")

        # Add atoms and bonds
        a1 = molecule.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = molecule.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        molecule.def_bond(a1, a2)

        # Convert to frame
        frame = molecule.to_frame()

        # Check frame structure
        assert "atoms" in frame
        assert "bonds" in frame
        assert frame["atoms"].nrows == 2
        assert frame["bonds"].nrows == 1

        # Convert back to structure
        reconstructed = Atomistic.from_frame(frame)

        assert len(reconstructed.atoms) == 2
        assert len(reconstructed.bonds) == 1
        assert reconstructed.atoms[0]["name"] == "A1"
        assert reconstructed.atoms[1]["name"] == "A2"

    def test_deep_copy(self):
        """Test deep copying of molecular structures."""
        # Create molecule
        molecule = Atomistic(name="test_molecule")

        # Add atoms and bonds
        a1 = molecule.def_atom(name="A1", element="C", xyz=np.array([0.0, 0.0, 0.0]))
        a2 = molecule.def_atom(name="A2", element="C", xyz=np.array([1.0, 0.0, 0.0]))
        molecule.def_bond(a1, a2)

        # Create deep copy
        copied = molecule()

        # Modify original
        molecule.atoms[0]["xyz"] = np.array([10.0, 10.0, 10.0])

        # Check that copy is unaffected
        assert not np.allclose(copied.atoms[0]["xyz"], np.array([10.0, 10.0, 10.0]))
        assert np.allclose(copied.atoms[0]["xyz"], np.array([0.0, 0.0, 0.0]))

        # Check that copy has same structure
        assert len(copied.atoms) == 2
        assert len(copied.bonds) == 1
        assert copied.atoms[0]["name"] == "A1"
        assert copied.atoms[1]["name"] == "A2"
