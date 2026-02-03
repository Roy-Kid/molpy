## ADDED Requirements

### Requirement: GAFF Atom Type Assignment
The system SHALL provide a `GaffAtomTypifier` that assigns GAFF atom types to atoms in an `Atomistic` structure using SMARTS pattern matching via the `LayeredTypingEngine`.

#### Scenario: Typify ethanol
- **WHEN** typifying an ethanol molecule with `GaffAtomTypifier`
- **THEN** the methyl carbons SHALL be assigned type `c3`
- **AND** the oxygen SHALL be assigned type `oh`
- **AND** the hydroxyl hydrogen SHALL be assigned type `ho`
- **AND** the aliphatic hydrogens SHALL be assigned type `hc`

#### Scenario: Typify benzene
- **WHEN** typifying a benzene molecule with `GaffAtomTypifier`
- **THEN** all carbons SHALL be assigned type `ca`
- **AND** all hydrogens SHALL be assigned type `ha`

### Requirement: GAFF Bond Parameter Assignment
The system SHALL provide a `GaffBondTypifier` that assigns harmonic bond parameters based on GAFF atom types, supporting wildcard (`X`) class matching for generic parameters.

#### Scenario: Assign bond parameters
- **WHEN** bond typing an ethanol molecule with assigned GAFF atom types
- **THEN** the C-C bond SHALL receive parameters matching `c3-c3` from GAFF
- **AND** the C-O bond SHALL receive parameters matching `c3-oh` from GAFF

#### Scenario: Wildcard bond matching
- **WHEN** no specific bond parameter exists for an atom type pair
- **THEN** the typifier SHALL fall back to wildcard parameters if available

### Requirement: GAFF Angle Parameter Assignment
The system SHALL provide a `GaffAngleTypifier` that assigns harmonic angle parameters based on GAFF atom types, supporting wildcard class matching.

#### Scenario: Assign angle parameters
- **WHEN** angle typing an ethanol molecule with assigned GAFF atom types
- **THEN** the C-C-O angle SHALL receive parameters matching `c3-c3-oh` from GAFF

### Requirement: GAFF Dihedral Parameter Assignment
The system SHALL provide a `GaffDihedralTypifier` that assigns periodic torsion parameters based on GAFF atom types, supporting wildcard class matching and multi-term dihedrals.

#### Scenario: Assign dihedral parameters
- **WHEN** dihedral typing an ethanol molecule with assigned GAFF atom types
- **THEN** the H-C-C-O dihedral SHALL receive periodic torsion parameters from GAFF

#### Scenario: Wildcard dihedral matching
- **WHEN** a specific 4-type dihedral parameter is not found
- **THEN** the typifier SHALL search for parameters with `X` wildcards (e.g., `X-c3-c3-X`)

### Requirement: GAFF Improper Dihedral Assignment
The system SHALL provide a `GaffImproperTypifier` that assigns improper dihedral parameters for planar centers (sp2 atoms with 3 substituents).

#### Scenario: Assign improper parameters
- **WHEN** improper typing a molecule containing a carbonyl group
- **THEN** the carbonyl carbon SHALL receive improper dihedral parameters to maintain planarity

### Requirement: GAFF Atomistic Typifier Orchestrator
The system SHALL provide a `GaffAtomisticTypifier` that orchestrates the full GAFF typing pipeline: atom typing → pair typing → bond typing → angle typing → dihedral typing → improper typing.

#### Scenario: Full GAFF typing pipeline
- **WHEN** calling `GaffAtomisticTypifier.typify(structure)`
- **THEN** all atoms SHALL have GAFF atom types assigned
- **AND** all bonds SHALL have harmonic parameters
- **AND** all angles SHALL have harmonic parameters
- **AND** all proper dihedrals SHALL have periodic torsion parameters
- **AND** all improper dihedrals SHALL have parameters where applicable

### Requirement: Shared Typifier Base Classes
The system SHALL provide `ForceFieldTypifier`, `ForceFieldAtomTypifier`, `ForceFieldBondTypifier`, `ForceFieldAngleTypifier`, and `ForceFieldDihedralTypifier` base classes that encapsulate the common SMARTS-based typing workflow. Both `OplsAtomisticTypifier` and `GaffAtomisticTypifier` SHALL extend these bases.

#### Scenario: OPLS still works after refactoring
- **WHEN** typifying a molecule with `OplsAtomisticTypifier` after refactoring
- **THEN** all results SHALL be identical to the pre-refactoring behavior
- **AND** all existing OPLS tests SHALL pass without modification

#### Scenario: No backward-compatibility code
- **WHEN** inspecting the typifier codebase after refactoring
- **THEN** there SHALL be no duplicated logic between OPLS and GAFF implementations
- **AND** there SHALL be no deprecated shims, re-exports, or backward-compatibility wrappers

## MODIFIED Requirements

## REMOVED Requirements
