## ADDED Requirements

### Requirement: Ring Connectivity Primitive
The SMARTS parser SHALL support the `x<n>` primitive for ring connectivity (number of ring bonds on an atom). When `<n>` is omitted, `x` SHALL match any atom that has at least one ring bond.

#### Scenario: Parse ring connectivity
- **WHEN** parsing the SMARTS pattern `[x2]`
- **THEN** the IR SHALL contain an `AtomPrimitiveIR` with type `ring_connectivity` and value `2`

#### Scenario: Match ring connectivity
- **WHEN** matching `[x2]` against a molecule graph
- **THEN** atoms with exactly 2 ring bonds SHALL match

### Requirement: Aliphatic Primitive
The SMARTS parser SHALL support the `A` primitive for aliphatic atoms (non-aromatic).

#### Scenario: Parse aliphatic atom
- **WHEN** parsing the SMARTS pattern `[A]`
- **THEN** the IR SHALL contain an `AtomPrimitiveIR` with type `aliphatic`

#### Scenario: Match aliphatic vs aromatic
- **WHEN** matching `[A]` against benzene
- **THEN** no ring carbon SHALL match (they are aromatic, not aliphatic)

### Requirement: Isotope Mass Prefix
The SMARTS parser SHALL support isotope mass as an integer prefix inside brackets, e.g., `[12C]`, `[2H]`.

#### Scenario: Parse isotope mass
- **WHEN** parsing the SMARTS pattern `[2H]`
- **THEN** the IR SHALL contain an `AtomPrimitiveIR` with type `isotope` and value `2`, combined with a symbol primitive for `H`

### Requirement: Chirality Primitives
The SMARTS parser SHALL support `@` (anticlockwise) and `@@` (clockwise) chirality primitives inside brackets.

#### Scenario: Parse chirality
- **WHEN** parsing the SMARTS pattern `[C@@]`
- **THEN** the IR SHALL contain an `AtomPrimitiveIR` with type `chirality` and value `@@`

### Requirement: Extended Charge Syntax
The SMARTS parser SHALL support charge notations: `+`, `-`, `+<n>`, `-<n>`, `++` (= +2), `--` (= -2).

#### Scenario: Parse double-plus charge
- **WHEN** parsing the SMARTS pattern `[Fe++]`
- **THEN** the IR SHALL produce a charge primitive with value `+2`

#### Scenario: Parse numeric charge
- **WHEN** parsing the SMARTS pattern `[O-2]`
- **THEN** the IR SHALL produce a charge primitive with value `-2`

### Requirement: GAFF SMARTS Compatibility
The SMARTS parser SHALL successfully parse all atom type SMARTS patterns used in the GAFF force field definition file (`gaff.prm`).

#### Scenario: Parse all GAFF patterns
- **WHEN** loading the 571 SMARTS patterns from `gaff.prm`
- **THEN** all patterns SHALL parse without error
- **AND** each parsed IR SHALL round-trip back to a semantically equivalent SMARTS string

## MODIFIED Requirements

## REMOVED Requirements
