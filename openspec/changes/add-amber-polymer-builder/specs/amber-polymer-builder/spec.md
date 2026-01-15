# Capability: amber-polymer-builder

Build polymers from CGSmiles notation using AmberTools (tleap) as the backend.

## ADDED Requirements

### Requirement: AmberPolymerBuilder construction

The system SHALL provide an `AmberPolymerBuilder` class that accepts:
- A library mapping CGSmiles labels to port-marked `Atomistic` structures
- A working directory path for intermediate files
- An optional force field specification (gaff or gaff2, default gaff2)

#### Scenario: Create builder with monomer library
- **GIVEN** a dictionary `library = {"EO": eo_monomer}` where `eo_monomer` is an `Atomistic` with ports marked
- **AND** a `workdir = Path("./amber_work")`
- **WHEN** `builder = AmberPolymerBuilder(library=library, workdir=workdir)`
- **THEN** `builder` is created successfully
- **AND** `builder.library` contains the provided monomers
- **AND** `builder.workdir` equals the provided path

#### Scenario: Reject library with unmarked ports
- **GIVEN** a monomer `bad_monomer` without any port markers
- **WHEN** `AmberPolymerBuilder(library={"X": bad_monomer}, workdir=...)`
- **THEN** a `ValueError` is raised indicating missing port markers

---

### Requirement: Build linear homopolymer

The system SHALL support building linear homopolymers from CGSmiles repeat notation.

#### Scenario: Build PEO chain with repeat operator
- **GIVEN** a builder with `library = {"EO": eo_monomer}` (EO has head port `"<"` and tail port `">"`)
- **WHEN** `result = builder.build("{[#EO]|10}")`
- **THEN** `result.frame` is a `Frame` containing the polymer coordinates
- **AND** `result.force_field` is a `ForceField` with GAFF parameters
- **AND** `result.monomer_count == 10`
- **AND** intermediate files exist in `builder.workdir`

#### Scenario: Build short chain without repeat
- **GIVEN** a builder with `library = {"A": monomer_a, "B": monomer_b}`
- **WHEN** `result = builder.build("{[#A][#B][#A]}")`
- **THEN** `result.monomer_count == 3`
- **AND** the polymer contains sequence A-B-A

---

### Requirement: Automatic monomer preparation

The system SHALL automatically prepare monomer residue templates using AmberTools wrappers.

#### Scenario: Monomer preparation workflow
- **GIVEN** a builder with a monomer that has not been prepared
- **WHEN** `builder.build(...)` is called
- **THEN** the builder internally:
  1. Writes the monomer as PDB
  2. Runs antechamber for atom typing
  3. Runs parmchk2 for missing parameters
  4. Generates prepgen control files from port markers
  5. Runs prepgen to create HEAD, CHAIN, TAIL residue templates

#### Scenario: Cached monomer preparation
- **GIVEN** a builder that has already prepared monomer "EO"
- **WHEN** `builder.build("{[#EO]|5}")` is called again
- **THEN** the builder reuses existing `.prepi` files without regeneration

---

### Requirement: Result access

The system SHALL return an `AmberBuildResult` containing:
- `frame`: The polymer as a `Frame` object
- `force_field`: The `ForceField` with all parameters
- `prmtop_path`: Path to the generated `.prmtop` file
- `inpcrd_path`: Path to the generated `.inpcrd` file
- `tleap_script`: The generated tleap script (for debugging)
- `monomer_count`: Total number of monomers in the polymer

#### Scenario: Access all result fields
- **GIVEN** a successful `result = builder.build("{[#EO]|5}")`
- **THEN** `result.frame` is a `Frame` instance
- **AND** `result.force_field` is a `ForceField` instance
- **AND** `result.prmtop_path.exists()` is True
- **AND** `result.inpcrd_path.exists()` is True
- **AND** `result.tleap_script` is a non-empty string containing "sequence"
- **AND** `result.monomer_count == 5`

---

### Requirement: Port-to-residue mapping

The system SHALL translate molpy port markers to prepgen control file entries based on sequence position.

#### Scenario: Port mapping based on sequence position
- **GIVEN** a monomer with ports marked using molpy convention (e.g., `port="1"` and `port="2"`)
- **WHEN** the monomer appears at the START of a sequence
- **THEN** the builder generates a HEAD variant: only the right-side port becomes `TAIL_NAME`

#### Scenario: Middle monomer mapping
- **GIVEN** a monomer with two ports
- **WHEN** the monomer appears in the MIDDLE of a sequence
- **THEN** the builder generates a CHAIN variant: left port → `HEAD_NAME`, right port → `TAIL_NAME`

#### Scenario: Terminal monomer mapping
- **GIVEN** a monomer with two ports
- **WHEN** the monomer appears at the END of a sequence
- **THEN** the builder generates a TAIL variant: only the left-side port becomes `HEAD_NAME`

---

### Requirement: tleap script generation

The system SHALL generate valid tleap scripts from CGSmiles notation.

#### Scenario: Linear sequence script
- **GIVEN** CGSmiles `"{[#A][#B][#C]}"`
- **AND** prepared templates `HA.prepi`, `A.prepi`, `B.prepi`, `C.prepi`, `TC.prepi`
- **WHEN** the script is generated
- **THEN** it contains `loadamberprep` commands for all templates
- **AND** it contains `mol = sequence {HA B TC}`

#### Scenario: Repeat expansion in script
- **GIVEN** CGSmiles `"{[#X]|4}"`
- **WHEN** the script is generated
- **THEN** the sequence command is `mol = sequence {HX X X TX}`

---

### Requirement: Error handling

The system SHALL provide clear error messages for common failure modes.

#### Scenario: Antechamber failure
- **GIVEN** a malformed monomer structure
- **WHEN** antechamber fails during preparation
- **THEN** a `MonomerPreparationError` is raised
- **AND** the error message includes relevant output from antechamber/sqm

#### Scenario: Unknown label in CGSmiles
- **GIVEN** a builder with `library = {"A": ...}`
- **WHEN** `builder.build("{[#X]}")` is called with unknown label "X"
- **THEN** a `ValueError` is raised indicating "X" not found in library
