import numpy as np

from molpy.core import Atom, Atomistic, Bond, Entities, Struct, Wrapper


class Monomer(Wrapper):

    def __post_init__(self, **props):

        self.head = None
        self.tail = None

        return super().__post_init__(**props)


class Polymer(Wrapper):

    def __post_init__(self, **props):
        self.atoms = []
        self.bonds = []
        return super().__post_init__(**props)

    def __len__(self):
        """Return the number of monomers in the polymer."""
        try:
            return len(self["atoms"])
        except (KeyError, TypeError):
            return 0


class PolymerBuilder:

    def __init__(self, monomer_mapping, style="CoarseGrain"):
        self.monomer_mapping = monomer_mapping
        self.style = style

    def linear(self, sequence, path):
        """
        Build a linear polymer from a sequence of monomers, a mapping of monomer types to bead structures,
        and a path defining the 3D coordinates of the polymer.

        Args:
            sequence (list of str): List of monomer types (e.g., ['A', 'B', 'C', 'A']).
            path (list of list of float): List of 3D coordinates for each bead in the polymer.

        Returns:
            Polymer: The constructed polymer object.
        """
        # For linear polymer, connections are None (will be auto-generated as 0-1, 1-2, etc.)
        polymer = self.build(sequence, None, path)
        return polymer

    def build(self, sequence, connections, path):
        """
        Build polymer from sequence, connections, and coordinates.

        Args:
            sequence: List of monomer type names
            connections: List of (i, j) tuples defining bonds between monomers
                        If None, assumes linear chain (0-1, 1-2, 2-3, ...)
            path: List of coordinates for each monomer

        Returns:
            tuple: (monomers, connections)
        """
        if len(sequence) != len(path):
            raise ValueError("Sequence length must match path length.")

        path = np.array(path)

        # If no connections provided, assume linear chain
        if connections is None:
            connections = [(i, i + 1) for i in range(len(sequence) - 1)]

        # Calculate orientation for each monomer based on its connections
        orientations = self._calc_orientations_from_connections(
            path, connections, len(sequence)
        )

        monomers = []
        for i, (monomer_type, coord, orient) in enumerate(
            zip(sequence, path, orientations)
        ):
            if monomer_type not in self.monomer_mapping:
                raise ValueError(
                    f"Monomer type '{monomer_type}' not found in monomer mapping."
                )

            bead = self.monomer_mapping[monomer_type].copy()
            bead.move_to(coord)

            # Only apply rotation if we have a valid orientation
            if orient is not None and np.linalg.norm(orient) > 1e-12:
                # Determine current direction based on monomer structure
                if (
                    hasattr(bead, "head")
                    and hasattr(bead, "tail")
                    and bead.head is not None
                    and bead.tail is not None
                ):
                    # For atomistic monomers, use head-tail direction
                    try:
                        current_dir = bead.tail.xyz - bead.head.xyz
                        if np.linalg.norm(current_dir) > 1e-12:
                            bead.rotate_to(orient, current_dir)
                    except:
                        # If can't determine direction, use default x-axis
                        bead.rotate_to(orient, np.array([1.0, 0.0, 0.0]))
                else:
                    # For coarse-grained beads, assume initial direction is x-axis
                    bead.rotate_to(orient, np.array([1.0, 0.0, 0.0]))

            monomers.append(bead)

        connects = []
        if connections:
            for conn in connections:
                i, j = conn
                if i < 0 or i >= len(monomers) or j < 0 or j >= len(monomers):
                    raise ValueError(f"Connection indices {conn} out of range.")
                connects.append((i, j))

        # Pack into Polymer object
        polymer = self._pack(monomers, connects)
        return polymer

    def _pack(self, monomers, connects):
        """
        Pack monomers into a polymer structure.

        For CoarseGrain: Returns Polymer(Struct) with beads as 'atoms'
        For Atomistic: Returns Atomistic(Struct) with merged atoms and bonds
        """

        if self.style == "CoarseGrain":
            # CoarseGrain polymer: beads are the "atoms"
            struct = Struct()
            struct["atoms"] = monomers  # beads
            struct["monomers"] = monomers
            struct["connects"] = connects
            return Polymer(struct)

        elif self.style == "Atomistic":
            # Atomistic polymer: collect all atoms and bonds from monomers
            struct = Struct()
            struct["atoms"] = Entities[Atom]()
            struct["bonds"] = Entities[Bond]()
            struct["angles"] = Entities()
            struct["dihedrals"] = Entities()

            # Collect all atoms and bonds from each monomer
            for monomer in monomers:
                for atom in monomer.atoms:
                    struct["atoms"].add(atom)
                for bond in monomer.bonds:
                    struct["bonds"].add(bond)
                # TODO: angles, dihedrals if needed

            # Add inter-monomer bonds based on connections
            if connects:
                for i, j in connects:
                    # Get the tail of monomer i and head of monomer j
                    tail_atom = monomers[i].tail
                    head_atom = monomers[j].head

                    # Create bond between them
                    inter_bond = Bond(atom1=tail_atom, atom2=head_atom)
                    struct["bonds"].add(inter_bond)

            return Atomistic(struct)

        else:
            raise ValueError(f"Unknown polymer style: {self.style}")

    def _calc_orientations_from_connections(self, path, connections, n_monomers):
        """
        Calculate orientation for each monomer based on connection graph.

        Strategy:
        - For each monomer, find all connected neighbors from connections
        - Calculate direction vectors to all neighbors
        - Use the first connection as primary orientation
        - If multiple connections, use average direction

        This allows:
        - Linear chain: each bead points to next bead
        - Branch point: bead points toward average of connected beads
        - Branch start: bead points back to main chain

        Args:
            path: numpy array of coordinates (n_monomers, 3)
            connections: list of (i, j) tuples
            n_monomers: total number of monomers

        Returns:
            list of orientation vectors (one per monomer)
        """
        path = np.array(path)
        orientations = []

        for i in range(n_monomers):
            # Find all neighbors connected to monomer i
            neighbors = []
            for conn in connections:
                if conn[0] == i:
                    neighbors.append(conn[1])
                elif conn[1] == i:
                    neighbors.append(conn[0])

            if neighbors:
                # Calculate direction to the first neighbor (primary direction)
                # For linear chains, this is the "forward" direction
                # For branches, this points to the main chain or first branch
                primary_neighbor = neighbors[0]
                primary_direction = path[primary_neighbor] - path[i]

                # If multiple neighbors, consider their average
                if len(neighbors) > 1:
                    directions = []
                    for neighbor in neighbors:
                        direction = path[neighbor] - path[i]
                        norm = np.linalg.norm(direction)
                        if norm > 1e-12:
                            directions.append(direction / norm)

                    if directions:
                        # Use average direction
                        avg_direction = np.mean(directions, axis=0)
                        norm = np.linalg.norm(avg_direction)
                        if norm > 1e-12:
                            orientations.append(avg_direction / norm)
                        else:
                            orientations.append(None)
                    else:
                        orientations.append(None)
                else:
                    # Single neighbor - use that direction
                    norm = np.linalg.norm(primary_direction)
                    if norm > 1e-12:
                        orientations.append(primary_direction / norm)
                    else:
                        orientations.append(None)
            else:
                # No connections - no rotation needed
                orientations.append(None)

        return orientations

    def _tagent_path(self, path):
        """
        Calculate tangent vectors along a path.

        For linear polymers, returns the direction from each point to the next.
        The last point uses the same direction as the second-to-last segment.
        """
        path = np.array(path)
        diff = np.diff(path, axis=0)

        # Normalize differences
        norms = np.linalg.norm(diff, axis=1, keepdims=True)
        norms = np.where(norms < 1e-12, 1.0, norms)  # Avoid division by zero
        udiff = diff / norms

        # For the last point, use the same direction as the last segment
        orientations = list(udiff)
        if len(orientations) > 0:
            orientations.append(orientations[-1])

        return orientations
