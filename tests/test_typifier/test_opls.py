"""
OPLS-AA typifier tests.

This module contains tests for OPLS atom type assignment using SmartsTypifier
and validation against reference topology data.
"""

from pathlib import Path
from typing import List, Dict

import molpy as mp
import pytest


def pytest_generate_tests(metafunc):
    """Generate tests for each molecule."""
    if "molecule_name" in metafunc.fixturenames:
        # Get molecules from test data directory
        test_data_dir = metafunc.config.getoption("--testdata-dir", default=None)
        if test_data_dir is None:
            # Use default path relative to test file
            test_data_dir = Path(__file__).parent.parent / "chemfile-testcases"
        
        opls_dir = Path(test_data_dir) / "opls"
        molecules = []
        
        if opls_dir.exists():
            for molecule_dir in opls_dir.iterdir():
                if molecule_dir.is_dir():
                    gro_file = molecule_dir / f"{molecule_dir.name}.gro"
                    top_file = molecule_dir / f"{molecule_dir.name}.top"
                    if gro_file.exists() and top_file.exists():
                        molecules.append(molecule_dir.name)
        
        metafunc.parametrize("molecule_name", sorted(molecules))


class TestOPLSTypifier:
    """Test suite for OPLS-AA force field typifier."""

    @pytest.fixture(scope="class")
    def opls_validation_dir(self, TEST_DATA_DIR) -> Path:
        """Path to OPLS validation test data directory."""
        return TEST_DATA_DIR / "opls"

    @pytest.fixture(scope="class")
    def available_molecules(self, opls_validation_dir: Path) -> List[str]:
        """Get list of molecules available for testing."""
        molecules = []
        for molecule_dir in opls_validation_dir.iterdir():
            if molecule_dir.is_dir():
                gro_file = molecule_dir / f"{molecule_dir.name}.gro"
                top_file = molecule_dir / f"{molecule_dir.name}.top"
                if gro_file.exists() and top_file.exists():
                    molecules.append(molecule_dir.name)
        return sorted(molecules)

    @pytest.fixture(scope="class")
    def opls_typifier(self) -> mp.typifier.SmartsTypifier:
        """Create OPLS-AA typifier instance."""
        # Load OPLS-AA forcefield
        frame_system = mp.io.read_xml_forcefield('oplsaa')
        ff = frame_system.forcefield

        typifier = mp.typifier.SmartsTypifier(ff)
        return typifier


    def test_opls_typification_single_molecule(self, molecule_name: str, 
                                               opls_validation_dir: Path,
                                               opls_typifier: mp.typifier.SmartsTypifier):
        """Test OPLS typification for a single molecule against reference."""
        molecule_dir = opls_validation_dir / molecule_name
        gro_file = molecule_dir / f"{molecule_name}.gro"
        top_file = molecule_dir / f"{molecule_name}.top"
        
        # Read structure using molpy
        frame = mp.io.read_gro(gro_file)
        
        # Convert Frame to Atomistic structure for typification
        atomistic_structure = mp.Atomistic.from_frame(frame)
        
        # Perform typification using SmartsTypifier
        typified_frame = opls_typifier.typify(atomistic_structure)
        
        # Extract reference atom types from topology file
        reference_types = self._extract_reference_types(top_file)
        
        # Extract typified atom types
        typified_types = self._extract_typified_types(typified_frame)
        
        # Compare results
        self._compare_typification_results(typified_types, reference_types, molecule_name)

    def _extract_reference_types(self, top_file: Path) -> Dict[int, str]:
        """Extract reference atom types from topology file."""
        reference_types = {}
        
        with open(top_file, 'r') as f:
            lines = f.readlines()
        
        in_atoms_section = False
        atom_index = 1
        
        for line in lines:
            line = line.strip()
            
            # Check for atoms section
            if line.startswith('[ atoms ]'):
                in_atoms_section = True
                continue
            
            # Check for next section
            if in_atoms_section and line.startswith('['):
                break
            
            # Parse atom line
            if in_atoms_section and line and not line.startswith(';'):
                parts = line.split()
                if len(parts) >= 2:
                    atom_type = parts[1]  # Second column is atom type
                    reference_types[atom_index] = atom_type
                    atom_index += 1
        
        return reference_types

    def _extract_typified_types(self, typified_structure: mp.Atomistic) -> Dict[int, str]:
        """Extract typified atom types from atomistic structure."""
        typified_types = {}
        
        if "atoms" in typified_structure:
            atoms = typified_structure["atoms"]
            for i, atom in enumerate(atoms):
                if "type" in atom:
                    typified_types[i + 1] = atom["type"]
                else:
                    # If no type assigned, use None or empty string
                    typified_types[i + 1] = None
        
        return typified_types

    def _compare_typification_results(self, typified_types: Dict[int, str], 
                                    reference_types: Dict[int, str], 
                                    molecule_name: str):
        """Compare typification results with reference."""
        # Check atom count consistency
        assert len(typified_types) == len(reference_types), \
            f"Atom count mismatch for {molecule_name}: typified={len(typified_types)}, reference={len(reference_types)}"
        
        # Compare each atom
        mismatches = []
        for atom_id in reference_types:
            reference_type = reference_types[atom_id]
            typified_type = typified_types.get(atom_id, "UNKNOWN")
            
            if reference_type != typified_type:
                mismatches.append({
                    'atom_id': atom_id,
                    'reference': reference_type,
                    'typified': typified_type
                })
        
        # Report results
        if mismatches:
            mismatch_rate = len(mismatches) / len(reference_types)
            error_msg = f"Typification mismatches for {molecule_name} ({mismatch_rate:.2%} errors):\n"
            for mismatch in mismatches[:10]:  # Show first 10 mismatches
                error_msg += f"  Atom {mismatch['atom_id']}: "
                error_msg += f"expected {mismatch['reference']}, got {mismatch['typified']}\n"
            if len(mismatches) > 10:
                error_msg += f"  ... and {len(mismatches) - 10} more mismatches\n"
            pytest.fail(error_msg)
        
        # Success message
        print(f"✓ {molecule_name}: All {len(reference_types)} atoms correctly typified")

    def test_opls_typification_statistics(self, available_molecules: List[str], 
                                        opls_validation_dir: Path,
                                        opls_typifier: mp.typifier.SmartsTypifier):
        """Collect statistics on OPLS typification performance."""
        if not available_molecules:
            pytest.skip("No OPLS test molecules available")
        
        total_molecules = len(available_molecules)
        successful_molecules = 0
        total_atoms = 0
        successful_atoms = 0
        failed_molecules = []
        
        for molecule_name in available_molecules:
            try:
                molecule_dir = opls_validation_dir / molecule_name
                gro_file = molecule_dir / f"{molecule_name}.gro"
                top_file = molecule_dir / f"{molecule_name}.top"
                
                # Read structure and perform typification
                frame = mp.io.read_gro(gro_file)
                typified_frame = opls_typifier.typify(frame)
                
                # Extract reference and typified atom types
                reference_types = self._extract_reference_types(top_file)
                typified_types = self._extract_typified_types(typified_frame)
                
                # Count matches
                molecule_atoms = len(reference_types)
                molecule_matches = 0
                
                for atom_id in reference_types:
                    reference_type = reference_types[atom_id]
                    typified_type = typified_types.get(atom_id, "UNKNOWN")
                    
                    if reference_type == typified_type:
                        molecule_matches += 1
                
                total_atoms += molecule_atoms
                successful_atoms += molecule_matches
                
                if molecule_matches == molecule_atoms:
                    successful_molecules += 1
                else:
                    failed_molecules.append({
                        'name': molecule_name,
                        'success_rate': molecule_matches / molecule_atoms,
                        'total_atoms': molecule_atoms,
                        'matched_atoms': molecule_matches
                    })
                    
            except Exception as e:
                failed_molecules.append({
                    'name': molecule_name,
                    'error': str(e),
                    'success_rate': 0.0
                })
        
        # Print statistics
        print(f"\n=== OPLS Typification Statistics ===")
        print(f"Total molecules tested: {total_molecules}")
        print(f"Successful molecules: {successful_molecules}")
        print(f"Molecule success rate: {successful_molecules/total_molecules:.2%}")
        print(f"Total atoms tested: {total_atoms}")
        print(f"Successful atoms: {successful_atoms}")
        print(f"Atom success rate: {successful_atoms/total_atoms:.2%}")
        
        if failed_molecules:
            print(f"\nFailed molecules:")
            for failed in failed_molecules[:10]:  # Show first 10 failures
                if 'error' in failed:
                    print(f"  {failed['name']}: ERROR - {failed['error']}")
                else:
                    print(f"  {failed['name']}: {failed['success_rate']:.2%} "
                          f"({failed['matched_atoms']}/{failed['total_atoms']})")
        
        # Assert reasonable success rate
        assert successful_atoms / total_atoms > 0.5, \
            f"Atom success rate too low: {successful_atoms/total_atoms:.2%}"