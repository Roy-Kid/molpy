"""
Example usage of Compute.map() for calculating radius of gyration.
"""

import numpy as np
from molpy.analysis.base import Compute, ComputeContext, get_selected_entities, get_entities_from_frame
from molpy.analysis.selection.expression import SelectResult, ExpressionSelection
from molpy.core.frame import Frame


class GyrationCompute(Compute):
    """
    Compute radius of gyration for a single molecule.
    
    Expects the context to contain 'entity_id' identifying the molecule.
    """

    def __init__(self, name: str = "gyration"):
        super().__init__(name)

    def compute(self, context: ComputeContext) -> ComputeContext:
        """
        Calculate radius of gyration for the specified molecule.
        
        Args:
            context: Context containing entity_id for the molecule
            
        Returns:
            Context with gyration radius in result['gyration_radius']
        """
        # Get the molecule ID from context
        if 'entity_id' not in context.result:
            raise ValueError("No entity_id found in context")
        
        mol_id = context.result['entity_id']
        
        # Get atoms belonging to this molecule
        frame = context.frame
        
        # Assuming we have molecule information in the frame
        if 'atoms' not in frame:
            raise ValueError("No atoms data found in frame")
        
        atoms = frame['atoms']
        
        # Filter atoms belonging to this molecule
        # Assuming we have a 'mol_id' field in atoms that identifies which molecule each atom belongs to
        if 'mol_id' in atoms:
            atom_mask = atoms['mol_id'] == mol_id
            atom_positions = atoms['xyz'][atom_mask]
        else:
            # Fallback: assume all atoms belong to the molecule (for single molecule systems)
            atom_positions = atoms['xyz']
        
        if len(atom_positions) == 0:
            raise ValueError(f"No atoms found for molecule {mol_id}")
        
        # Calculate center of mass (assuming equal masses for simplicity)
        center_of_mass = np.mean(atom_positions, axis=0)
        
        # Calculate distances from center of mass
        distances = atom_positions - center_of_mass
        squared_distances = np.sum(distances**2, axis=1)
        
        # Radius of gyration
        gyration_radius = np.sqrt(np.mean(squared_distances))
        
        # Store result
        context.result['gyration_radius'] = gyration_radius
        context.result['center_of_mass'] = center_of_mass
        context.result['n_atoms'] = len(atom_positions)
        
        return context


class MockMolecularFrame(Frame):
    """Mock frame with molecular data for testing."""
    
    def __init__(self, n_molecules: int = 3, atoms_per_mol: int = 10):
        super().__init__()
        
        total_atoms = n_molecules * atoms_per_mol
        
        # Create atoms with random positions
        self["atoms"] = {
            "xyz": np.random.random((total_atoms, 3)) * 10,  # Random positions in 10x10x10 box
            "id": np.arange(total_atoms),
            "mol_id": np.repeat(np.arange(n_molecules), atoms_per_mol)  # Molecule assignments
        }
        
        # Create molecules
        self["molecules"] = {
            "id": np.arange(n_molecules),
            "name": [f"MOL_{i}" for i in range(n_molecules)]
        }


def demonstrate_mapped_gyration():
    """Demonstrate mapped radius of gyration calculation."""
    
    # Create a mock frame with 3 molecules
    frame = MockMolecularFrame(n_molecules=3, atoms_per_mol=10)
    context = ComputeContext.attach_frame(frame)
    
    print("=== Mapped Radius of Gyration Demo ===")
    print(f"Frame contains {len(frame['molecules']['id'])} molecules")
    print(f"Total atoms: {len(frame['atoms']['id'])}")
    
    # Create selection for molecules (select first 2)
    def select_first_two(ctx: ComputeContext) -> np.ndarray:
        n_mols = len(ctx.frame["molecules"]["id"])
        mask = np.zeros(n_mols, dtype=bool)
        mask[:2] = True  # Select first 2 molecules
        return mask
    
    # Create selection operation
    selection = ExpressionSelection("first_two_molecules", select_first_two, select_field="molecules")
    selected_context = selection.compute(context)
    
    # Type assertion for selected context
    assert isinstance(selected_context, SelectResult)
    print(f"Selected molecules: {selected_context.selected}")
    
    # Create gyration calculation
    gyration_calc = GyrationCompute("gyration")
    
    # Method 1: Map using selected entities
    selected_ids = get_selected_entities(selected_context)
    if selected_ids is not None:
        mapped_result = gyration_calc.map(context, entity_ids=selected_ids, selection_field="molecules")
        
        print("\n=== Results (Using Selected Entities) ===")
        molecules_results = mapped_result.result["molecules_results"]
        mapped_molecules = mapped_result.result["mapped_molecules"]
        
        print(f"Number of computed molecules: {len(molecules_results)}")
        
        for mol_id in mapped_molecules:
            mol_result = molecules_results[mol_id]
            gyration = mol_result.get('gyration_radius', 'N/A')
            n_atoms = mol_result.get('n_atoms', 'N/A')
            com = mol_result.get('center_of_mass', 'N/A')
            
            print(f"Molecule {mol_id}:")
            print(f"  - Radius of gyration: {gyration:.3f}")
            print(f"  - Number of atoms: {n_atoms}")
            print(f"  - Center of mass: {com}")
            print()
    
    # Method 2: Map all molecules (without any selection)
    print("\n=== Map All Molecules (No Selection) ===")
    all_mapped = gyration_calc.map(context, selection_field="molecules")
    all_results = all_mapped.result["molecules_results"]
    
    for mol_id, mol_result in all_results.items():
        gyration = mol_result['gyration_radius']
        print(f"Molecule {mol_id} gyration radius: {gyration:.3f} Å")
    
    # Method 3: Map specific molecules by manual ID specification
    print("\n=== Map Specific Molecules (Manual IDs) ===")
    specific_ids = np.array([0, 2])  # Only molecules 0 and 2
    specific_mapped = gyration_calc.map(context, entity_ids=specific_ids, selection_field="molecules")
    specific_results = specific_mapped.result["molecules_results"]
    
    for mol_id, mol_result in specific_results.items():
        gyration = mol_result['gyration_radius']
        print(f"Molecule {mol_id} gyration radius: {gyration:.3f} Å")
    
    # Method 4: Using utility function to get all entities from frame
    print("\n=== Using Utility Functions ===")
    all_mol_ids = get_entities_from_frame(context, "molecules")
    if all_mol_ids is not None:
        print(f"All molecule IDs from frame: {all_mol_ids}")
        # Could map over subset
        subset_ids = all_mol_ids[::2]  # Every other molecule
        subset_mapped = gyration_calc.map(context, entity_ids=subset_ids, selection_field="molecules")
        subset_results = subset_mapped.result["molecules_results"]
        
        print(f"Mapped subset {subset_ids}:")
        for mol_id, mol_result in subset_results.items():
            gyration = mol_result['gyration_radius']
            print(f"  Molecule {mol_id}: {gyration:.3f} Å")


if __name__ == "__main__":
    demonstrate_mapped_gyration()
