import logging
from dataclasses import dataclass, field
from pathlib import Path

from molpy.core.wrapper import Wrapper

logger = logging.getLogger(__name__)


@dataclass
class ReactionTemplate:
    """
    Defines a reaction site with anchor atoms and edge connections.
    """

    init_atom: int  # Main anchor atom index
    edge_atom: int  # Atoms that will be deleted/modified


class ReactantWrapper(Wrapper):
    """
    Wrapper for Atomistic with reaction-specific functionality.
    Manages reaction sites, anchor points, and template extraction.
    """

    def __init__(self, struct, **kwargs):
        super().__init__(struct, **kwargs)

    def __call__(self, **kwargs):
        return self._wrapped(**kwargs)

    def extract_fragment(self, template: ReactionTemplate) -> "ReactantWrapper":
        """
        Extract molecular fragment between two anchor points.

        Args:
            start: Starting atom index
            end: Ending atom index

        Returns:
            New ReactantWrapper with extracted fragment
        """
        topology = self.get_topology()
        start = template.edge_atom
        end = template.init_atom
        main_chain, branches = get_main_chain_and_branches(topology, start, end)
        template_atoms = set(main_chain + branches)
        extracted_struct = self.get_substruct(template_atoms)

        return ReactantWrapper(extracted_struct)


def get_main_chain_and_branches(g, start: int, end: int):
    """
    Extract main chain and branches from topology graph.

    Args:
        g: Topology graph
        start: Starting atom index
        end: Ending atom index

    Returns:
        Tuple of (main_chain_atoms, branch_atoms)
    """
    from collections import deque

    def get_branch_subtree(g, start: int, excluded: set):
        visited = set()
        queue = deque([start])
        while queue:
            v = queue.popleft()
            if v in visited or v in excluded:
                continue
            visited.add(v)
            for n in g.neighbors(v):
                if n not in visited and n not in excluded:
                    queue.append(n)
        return visited

    # 1. Compute shortest path between start and end
    path = g.get_shortest_paths(start, to=end, output="vpath")[0]
    main_chain = set(path)

    # 2. Initialize visited with all main chain atoms
    visited = set(path)
    branch_atoms = set()

    # 3. Start from second atom in the path (exclude start atom's branches)
    for v in path[1:]:
        for n in g.neighbors(v):
            if n not in visited:
                # Collect entire connected component starting from n
                subtree = get_branch_subtree(g, n, main_chain)
                branch_atoms.update([x for x in subtree if x not in main_chain])
                visited.update(subtree)

    return list(path), list(branch_atoms)


class ReacterBuilder:
    """
    LAMMPS-compatible reaction template builder for automated polymer synthesis.

    This class automates the generation of pre- and post-reaction molecular templates
    for use with LAMMPS fix bond/react functionality.
    """

    def __init__(
        self,
        workdir: Path,
        typifier,
    ):
        """
        Initialize the reaction template builder.

        Args:
            template1: First reaction template
            template2: Second reaction template
            bond_changes: List of bond formation specifications
            atoms_to_delete: List of atoms to delete during reaction
            build_dir: Directory for output files
            conda_env: Conda environment for AmberTools
        """
        self.workdir = workdir
        self.typifier = typifier
        self.workdir.mkdir(parents=True, exist_ok=True)

    def export_lammps_template(self, struct, name: str) -> Path:
        """
        Export structure as LAMMPS molecule template using built-in writer.

        Args:
            struct: Structure to export
            name: Template name

        Returns:
            Path to generated template file
        """
        logger.info(f"Exporting LAMMPS template for {name}")

        template_path = self.workdir / f"{name}.template"

        # Convert to frame format
        frame = struct.to_frame()
        # Use molpy's built-in LAMMPS molecule writer
        import molpy.io as mp_io

        mp_io.write_lammps_molecule(template_path, frame)
        return template_path

    def generate_mapping(
        self,
        name,
        pre_struct: ReactantWrapper,
        post_struct: ReactantWrapper,
        inits: list[int],
        edges: list[int],
    ) -> dict:
        """
        Generate LAMMPS mapping file for bond/react.

        Args:
            pre_struct: Pre-reaction structure
            post_struct: Post-reaction structure
            inits: List of initiator atom IDs
            edges: List of edge atom IDs

        Returns:
            Dictionary containing mapping information and file path
        """
        logger.info("Generating LAMMPS mapping file")

        # Get atoms from both structures
        pre_atoms = list(pre_struct.atoms)
        post_atoms = list(post_struct.atoms)

        # Generate 1:1 atom correspondence
        equivalences = []
        min_atoms = min(len(pre_atoms), len(post_atoms))

        for i in range(min_atoms):
            pre_id = pre_atoms[i].get("id", i + 1)
            post_id = post_atoms[i].get("id", i + 1)
            equivalences.append((pre_id, post_id))

        # Use provided initiator and edge atoms
        initiator_ids = inits if inits else [1, 2]  # Default to first two atoms
        edge_ids = edges if edges else []

        # Generate mapping file content
        mapping_content = self._generate_mapping_content(
            equivalences, initiator_ids, edge_ids
        )

        # Write mapping file
        mapping_file = self.workdir / f"{name}.map"
        with open(mapping_file, "w") as f:
            f.write(mapping_content)

        return {
            "mapping_file": mapping_file,
            "equivalences": equivalences,
            "initiator_ids": initiator_ids,
            "edge_ids": edge_ids,
        }

    def _generate_mapping_content(self, equivalences, initiator_ids, edge_ids) -> str:
        """Generate the content for the LAMMPS mapping file."""
        lines = []
        lines.append("# LAMMPS bond/react mapping file")
        lines.append("")

        # Write equivalences
        lines.append("# Atom equivalences (pre -> post)")
        for pre_id, post_id in equivalences:
            lines.append(f"{pre_id} {post_id}")

        lines.append("")
        lines.append("# Initiator atoms")
        for atom_id in initiator_ids:
            lines.append(f"init {atom_id}")

        lines.append("")
        lines.append("# Edge atoms")
        for atom_id in edge_ids:
            lines.append(f"edge {atom_id}")

        return "\n".join(lines)

    def build(self, name, pre_struct, post_struct, inits=[], edges=[]):
        """Build the reaction system."""
        # Generate mapping
        mapping_info = self.generate_mapping(
            name, pre_struct, post_struct, inits, edges
        )

        # Export templates
        pre_template = self.export_lammps_template(pre_struct, f"{name}_pre")
        post_template = self.export_lammps_template(post_struct, f"{name}_post")

        return {
            "mapping": mapping_info,
            "pre_template": pre_template,
            "post_template": post_template,
        }
