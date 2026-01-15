"""
AmberPolymerBuilder: Build polymers using AmberTools backend.

This module provides a polymer builder that uses the AmberTools suite
(antechamber, parmchk2, prepgen, tleap) for construction. The API is
consistent with PolymerBuilder.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

from molpy.core.atomistic import Atomistic
from molpy.io.readers import read_amber_prmtop
from molpy.parser.smiles import parse_cgsmiles
from molpy.parser.smiles.cgsmiles_ir import CGSmilesGraphIR, CGSmilesIR

from .types import AmberBuildResult


@dataclass
class AmberPolymerBuilderConfig:
    """Configuration for AmberPolymerBuilder.
    
    Attributes:
        force_field: Force field to use (gaff or gaff2).
        charge_method: Charge method for antechamber.
        work_dir: Directory for intermediate files.
        keep_intermediates: Whether to keep intermediate files after build.
        env: Conda environment name or path for AmberTools (e.g., "AmberTools25").
        env_manager: Environment manager type ("conda" for conda environments).
    """
    
    force_field: Literal["gaff", "gaff2"] = "gaff2"
    charge_method: str = "bcc"
    work_dir: Path = field(default_factory=lambda: Path("amber_work"))
    keep_intermediates: bool = True
    env: str | Path | None = None
    env_manager: str | None = None


class AmberPolymerBuilder:
    """Build polymers from CGSmiles notation using AmberTools backend.
    
    This builder parses CGSmiles strings and constructs polymers using
    AmberTools (antechamber, parmchk2, prepgen, tleap). The API is consistent
    with PolymerBuilder.
    
    Internally, the builder:
    1. Prepares each monomer type (antechamber → parmchk2 → prepgen)
    2. Generates HEAD/CHAIN/TAIL residue variants based on port annotations
    3. Translates CGSmiles to tleap sequence command
    4. Runs tleap to build the polymer
    5. Returns Frame + ForceField from prmtop/inpcrd
    
    Example:
        >>> # Load monomer and mark ports
        >>> eo_monomer = read_pdb("PEO_initial.pdb")
        >>> eo_monomer.atoms[0]["port"] = "<"  # Head port
        >>> eo_monomer.atoms[6]["port"] = ">"  # Tail port
        >>> 
        >>> # Build polymer
        >>> builder = AmberPolymerBuilder(
        ...     library={"EO": eo_monomer},
        ...     config=AmberPolymerBuilderConfig(force_field="gaff2"),
        ... )
        >>> result = builder.build("{[#EO]|10}")
        >>> frame = result.frame
        >>> forcefield = result.forcefield
    """
    
    def __init__(
        self,
        library: Mapping[str, Atomistic],
        config: AmberPolymerBuilderConfig | None = None,
    ):
        """Initialize the polymer builder.
        
        Args:
            library: Mapping from CGSmiles labels to Atomistic monomer structures.
                Each Atomistic must have port annotations on atoms (port="<" for
                head, port=">" for tail).
            config: Builder configuration.
        """
        self.library = library
        self.config = config or AmberPolymerBuilderConfig()
        
        # Internal state
        self._prepared_monomers: dict[str, _PreparedMonomer] = {}
    
    def build(self, cgsmiles: str) -> AmberBuildResult:
        """Build a polymer from a CGSmiles string.
        
        Args:
            cgsmiles: CGSmiles notation string (e.g., "{[#EO]|10}")
            
        Returns:
            AmberBuildResult containing Frame, ForceField, and file paths.
            
        Raises:
            ValueError: If CGSmiles is invalid or labels not in library.
        """
        # Parse CGSmiles
        ir = parse_cgsmiles(cgsmiles)
        
        # Validate
        self._validate_ir(ir)
        
        # Prepare all monomers (antechamber → parmchk2 → prepgen)
        self._prepare_monomers(ir.base_graph)
        
        # Generate and run tleap
        result = self._build_with_tleap(ir.base_graph, output_prefix="polymer")
        
        return result
    
    def _validate_ir(self, ir: CGSmilesIR) -> None:
        """Validate CGSmiles IR."""
        graph = ir.base_graph
        
        if not graph.nodes:
            raise ValueError("CGSmiles graph is empty")
        
        # Check all labels exist in library
        missing_labels = set()
        for node in graph.nodes:
            if node.label not in self.library:
                missing_labels.add(node.label)
        
        if missing_labels:
            available = list(self.library.keys())
            raise ValueError(
                f"Labels {sorted(missing_labels)} not found in library. "
                f"Available labels: {available}"
            )
        
        # Validate port annotations on each monomer
        for label, monomer in self.library.items():
            head_port = None
            tail_port = None
            for atom in monomer.atoms:
                port = atom.get("port")
                if port == "<":
                    head_port = atom
                elif port == ">":
                    tail_port = atom
            
            if head_port is None or tail_port is None:
                raise ValueError(
                    f"Monomer '{label}' must have port annotations: "
                    f"port='<' for head, port='>' for tail. "
                    f"Found: head={head_port is not None}, tail={tail_port is not None}"
                )
    
    def _prepare_monomers(self, graph: CGSmilesGraphIR) -> None:
        """Prepare all monomer types used in the graph.
        
        For each monomer type:
        1. Write structure to mol2/pdb
        2. Run antechamber for atom typing and charges
        3. Run parmchk2 for missing parameters
        4. Run prepgen to create HEAD/CHAIN/TAIL variants
        """
        from molpy.wrapper.antechamber import AntechamberWrapper
        from molpy.wrapper.prepgen import Parmchk2Wrapper, PrepgenWrapper, write_prepgen_control_file
        
        # Common wrapper kwargs for environment
        wrapper_kwargs: dict = {"workdir": None}  # Will be set per-monomer
        if self.config.env is not None:
            wrapper_kwargs["env"] = self.config.env
            wrapper_kwargs["env_manager"] = self.config.env_manager
        
        # Get unique labels
        labels = {node.label for node in graph.nodes}
        
        work_dir = self.config.work_dir
        work_dir.mkdir(parents=True, exist_ok=True)
        
        for label in labels:
            if label in self._prepared_monomers:
                continue  # Already prepared
            
            monomer = self.library[label]
            monomer_dir = work_dir / "monomers" / label
            monomer_dir.mkdir(parents=True, exist_ok=True)
            
            # Find port atoms
            head_atom_name = None
            tail_atom_name = None
            for atom in monomer.atoms:
                port = atom.get("port")
                if port == "<":
                    head_atom_name = atom["name"]
                elif port == ">":
                    tail_atom_name = atom["name"]
            
            # Step 1: Write monomer to PDB
            input_pdb = monomer_dir / f"{label}.pdb"
            self._write_atomistic_pdb(monomer, input_pdb)
            
            # Step 2: Run antechamber
            ac_file = monomer_dir / f"{label}.ac"
            mol2_file = monomer_dir / f"{label}.mol2"
            
            antechamber = AntechamberWrapper(
                name="antechamber",
                workdir=monomer_dir,
                env=self.config.env,
                env_manager=self.config.env_manager,
            )
            antechamber.atomtype_assign(
                input_file=input_pdb,
                output_file=mol2_file,
                input_format="pdb",
                output_format="mol2",
                charge_method=self.config.charge_method,
                atom_type=self.config.force_field,
            )
            
            # Also generate .ac file for prepgen
            antechamber.atomtype_assign(
                input_file=input_pdb,
                output_file=ac_file,
                input_format="pdb",
                output_format="ac",
                charge_method=self.config.charge_method,
                atom_type=self.config.force_field,
            )
            
            # Step 3: Run parmchk2
            frcmod_file = monomer_dir / f"{label}.frcmod"
            parmchk2 = Parmchk2Wrapper(
                name="parmchk2",
                workdir=monomer_dir,
                env=self.config.env,
                env_manager=self.config.env_manager,
            )
            parmchk2.generate_parameters(
                input_file=mol2_file,
                output_file=frcmod_file,
            )
            
            # Step 4: Generate prepgen control files and run prepgen
            # HEAD variant (chain start): only TAIL connection active
            head_ctrl = monomer_dir / f"{label}.head"
            head_prepi = monomer_dir / f"H{label}.prepi"
            write_prepgen_control_file(
                head_ctrl,
                variant="head",
                head_name=None,  # No HEAD connection for chain start
                tail_name=tail_atom_name,
                omit_names=self._get_omit_names(monomer, "<"),
            )
            
            # CHAIN variant (middle): both connections
            chain_ctrl = monomer_dir / f"{label}.chain"
            chain_prepi = monomer_dir / f"{label}.prepi"
            write_prepgen_control_file(
                chain_ctrl,
                variant="chain",
                head_name=head_atom_name,
                tail_name=tail_atom_name,
            )
            
            # TAIL variant (chain end): only HEAD connection active
            tail_ctrl = monomer_dir / f"{label}.tail"
            tail_prepi = monomer_dir / f"T{label}.prepi"
            write_prepgen_control_file(
                tail_ctrl,
                variant="tail",
                head_name=head_atom_name,
                tail_name=None,  # No TAIL connection for chain end
                omit_names=self._get_omit_names(monomer, ">"),
            )
            
            # Run prepgen for each variant
            prepgen = PrepgenWrapper(
                name="prepgen",
                workdir=monomer_dir,
                env=self.config.env,
                env_manager=self.config.env_manager,
            )
            
            prepgen.generate_residue(
                input_file=ac_file,
                output_file=head_prepi,
                control_file=head_ctrl,
                residue_name=f"H{label[:2].upper()}",
            )
            
            prepgen.generate_residue(
                input_file=ac_file,
                output_file=chain_prepi,
                control_file=chain_ctrl,
                residue_name=label[:3].upper(),
            )
            
            prepgen.generate_residue(
                input_file=ac_file,
                output_file=tail_prepi,
                control_file=tail_ctrl,
                residue_name=f"T{label[:2].upper()}",
            )
            
            # Store prepared monomer info
            self._prepared_monomers[label] = _PreparedMonomer(
                label=label,
                frcmod_file=frcmod_file,
                head_prepi=head_prepi,
                chain_prepi=chain_prepi,
                tail_prepi=tail_prepi,
                head_resname=f"H{label[:2].upper()}",
                chain_resname=label[:3].upper(),
                tail_resname=f"T{label[:2].upper()}",
            )
    
    def _get_omit_names(self, monomer: Atomistic, port: str) -> list[str]:
        """Get atom names to omit for a port (leaving group atoms).
        
        For simplicity, we omit hydrogen atoms bonded to the port atom.
        """
        omit = []
        port_atom = None
        for atom in monomer.atoms:
            if atom["port"] == port:
                port_atom = atom
                break
        
        # Find atoms bonded to port atom that are hydrogens (leaving group)
        for bond in monomer.bonds:
            atom_i, atom_j = bond.atoms
            
            if atom_i == port_atom and atom_j["element"] == "H":
                omit.append(atom_j["name"])
            elif atom_j == port_atom and atom_i["element"] == "H":
                omit.append(atom_i["name"])
        
        return omit
    
    def _write_atomistic_pdb(self, atomistic: Atomistic, path: Path) -> None:
        """Write Atomistic to PDB format."""
        from molpy.io.writers import write_pdb
        frame = atomistic.to_frame()
        write_pdb(path, frame)
    
    def _build_with_tleap(
        self,
        graph: CGSmilesGraphIR,
        output_prefix: str,
    ) -> AmberBuildResult:
        """Generate tleap script and run tleap to build polymer."""
        from molpy.wrapper.tleap import TLeapWrapper
        
        work_dir = self.config.work_dir
        output_dir = work_dir / "output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate tleap script
        script_lines = []
        
        # Source force field
        leaprc = f"leaprc.{self.config.force_field}"
        script_lines.append(f"source {leaprc}")
        script_lines.append("")
        
        # Load frcmod and prepi files for each monomer type
        used_labels = {node.label for node in graph.nodes}
        for label in sorted(used_labels):
            prep = self._prepared_monomers[label]
            script_lines.append(f"loadamberparams {prep.frcmod_file}")
            script_lines.append(f"loadamberprep {prep.head_prepi}")
            script_lines.append(f"loadamberprep {prep.chain_prepi}")
            script_lines.append(f"loadamberprep {prep.tail_prepi}")
        script_lines.append("")
        
        # Build sequence
        sequence = self._build_sequence(graph)
        script_lines.append(f"mol = sequence {{{sequence}}}")
        script_lines.append("")
        
        # Output files
        prmtop = output_dir / f"{output_prefix}.prmtop"
        inpcrd = output_dir / f"{output_prefix}.inpcrd"
        pdb = output_dir / f"{output_prefix}.pdb"
        
        script_lines.append(f"savepdb mol {pdb}")
        script_lines.append(f"saveamberparm mol {prmtop} {inpcrd}")
        script_lines.append("quit")
        
        script = "\n".join(script_lines)
        
        # Write script
        script_path = work_dir / f"{output_prefix}.in"
        script_path.write_text(script)
        
        # Run tleap
        tleap = TLeapWrapper(
            name="tleap",
            workdir=work_dir,
            env=self.config.env,
            env_manager=self.config.env_manager,
        )
        tleap.run_from_script(script)
        
        # Load results
        frame, forcefield = read_amber_prmtop(prmtop, inpcrd)
        
        return AmberBuildResult(
            frame=frame,
            forcefield=forcefield,
            prmtop_path=prmtop,
            inpcrd_path=inpcrd,
            pdb_path=pdb if pdb.exists() else None,
            monomer_count=len(graph.nodes),
            cgsmiles=None,
        )
    
    def _build_sequence(self, graph: CGSmilesGraphIR) -> str:
        """Build tleap sequence from CGSmiles graph.
        
        Assigns HEAD/CHAIN/TAIL variants based on position:
        - First node: HEAD variant
        - Middle nodes: CHAIN variant  
        - Last node: TAIL variant
        """
        nodes = graph.nodes
        n = len(nodes)
        
        if n == 0:
            raise ValueError("Empty graph")
        
        residue_names = []
        for i, node in enumerate(nodes):
            prep = self._prepared_monomers[node.label]
            
            if n == 1:
                # Single monomer - use chain (or could use head)
                residue_names.append(prep.chain_resname)
            elif i == 0:
                # First monomer: HEAD variant
                residue_names.append(prep.head_resname)
            elif i == n - 1:
                # Last monomer: TAIL variant
                residue_names.append(prep.tail_resname)
            else:
                # Middle: CHAIN variant
                residue_names.append(prep.chain_resname)
        
        return " ".join(residue_names)


@dataclass
class _PreparedMonomer:
    """Internal state for a prepared monomer."""
    
    label: str
    frcmod_file: Path
    head_prepi: Path
    chain_prepi: Path
    tail_prepi: Path
    head_resname: str
    chain_resname: str
    tail_resname: str
