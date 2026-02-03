## ADDED Requirements

### Requirement: GAFF XML Forcefield Data
The system SHALL provide a `gaff.xml` forcefield file in `src/molpy/data/forcefield/` containing all GAFF parameters converted from the AMBER `gaff.dat` and `gaff.prm` source files.

#### Scenario: Load GAFF forcefield
- **WHEN** calling `get_forcefield_path("gaff.xml")`
- **THEN** a valid path to the GAFF XML file SHALL be returned
- **AND** the XML SHALL parse as a valid `ForceField` object

#### Scenario: GAFF atom types with SMARTS
- **WHEN** reading `<AtomTypes>` from `gaff.xml`
- **THEN** each atom type SHALL have a `def` attribute containing a SMARTS pattern
- **AND** each atom type SHALL have `name`, `class`, `element`, `mass`, and `desc` attributes
- **AND** atom type priorities SHALL reflect the ordering from `gaff.prm` (later rules = higher priority)

#### Scenario: GAFF combining rule
- **WHEN** reading the `<ForceField>` root element
- **THEN** `combining_rule` SHALL be `"arithmetic"` (Lorentz-Berthelot)

### Requirement: Periodic Torsion Force Parameters
The forcefield XML schema SHALL support `<PeriodicTorsionForce>` elements containing periodic (Fourier) dihedral parameters with attributes: `class1`, `class2`, `class3`, `class4`, `periodicity`, `phase`, `k`.

#### Scenario: Read GAFF dihedral parameters
- **WHEN** parsing `<PeriodicTorsionForce>` from `gaff.xml`
- **THEN** dihedral entries SHALL contain periodicity (integer), phase (degrees), and force constant k (kJ/mol)
- **AND** wildcard class `X` SHALL be supported for generic dihedrals

### Requirement: Periodic Improper Force Parameters
The forcefield XML schema SHALL support `<PeriodicImproperForce>` elements for improper dihedral parameters with the same attribute structure as periodic torsions.

#### Scenario: Read GAFF improper parameters
- **WHEN** parsing `<PeriodicImproperForce>` from `gaff.xml`
- **THEN** improper entries SHALL contain periodicity, phase, and force constant
- **AND** the central atom SHALL be identifiable by its position in the class attributes

### Requirement: Unit Conversion Correctness
All parameters in `gaff.xml` SHALL be stored in molpy's internal unit system (kJ/mol for energy, nm for length, radians for angles where applicable).

#### Scenario: Bond parameter units
- **WHEN** reading a bond parameter from `gaff.xml`
- **THEN** `length` SHALL be in nm (original Å ÷ 10)
- **AND** `k` SHALL be in kJ/mol/nm² (original kcal/mol/Å² × 4.184 × 100)

#### Scenario: VDW parameter units
- **WHEN** reading a nonbonded parameter from `gaff.xml`
- **THEN** `sigma` SHALL be in nm (converted from AMBER R* via σ = R* × 2^(5/6) / 10)
- **AND** `epsilon` SHALL be in kJ/mol (original kcal/mol × 4.184)

## MODIFIED Requirements

## REMOVED Requirements
