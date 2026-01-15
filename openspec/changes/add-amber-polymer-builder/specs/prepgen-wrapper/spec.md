# Capability: prepgen-wrapper

Enhanced wrapper for the `prepgen` AmberTools utility.

## ADDED Requirements

### Requirement: Generate residue template

The system SHALL provide a `generate_residue()` method on `PrepgenWrapper` that creates AMBER residue templates.

#### Scenario: Generate chain residue template
- **GIVEN** a `PrepgenWrapper` with a valid workdir
- **AND** an input `.ac` file from antechamber
- **AND** a control file specifying HEAD_NAME, TAIL_NAME, and OMIT_NAME entries
- **WHEN** `wrapper.generate_residue(input_file="mol.ac", output_file="mol.prepi", control_file="mol.chain", residue_name="MOL")`
- **THEN** the `.prepi` file is created in the workdir
- **AND** the residue is named "MOL"

#### Scenario: Generate head-only residue template
- **GIVEN** a control file with only TAIL_NAME (no HEAD_NAME)
- **WHEN** `wrapper.generate_residue(..., control_file="mol.head")`
- **THEN** the generated `.prepi` contains a residue suitable for chain start

---

### Requirement: Control file generation utility

The system SHALL provide a `write_prepgen_control_file()` function to create control files from port information.

#### Scenario: Generate chain control file
- **GIVEN** head atom name `"C1"`, tail atom name `"O5"`, omit atoms `["H1", "H2"]`
- **WHEN** `write_prepgen_control_file(path, variant="chain", head_name="C1", tail_name="O5", omit_names=["H1", "H2"], head_type="c3", tail_type="os")`
- **THEN** the file contains:
  ```
  HEAD_NAME C1
  TAIL_NAME O5
  OMIT_NAME H1
  OMIT_NAME H2
  PRE_HEAD_TYPE c3
  POST_TAIL_TYPE os
  CHARGE 0
  ```

#### Scenario: Generate head-only control file
- **GIVEN** only tail atom name `"O5"` and omit atoms near tail
- **WHEN** `write_prepgen_control_file(path, variant="head", tail_name="O5", omit_names=["H3"], tail_type="os")`
- **THEN** the file contains TAIL_NAME and POST_TAIL_TYPE but no HEAD_NAME

#### Scenario: Generate tail-only control file
- **GIVEN** only head atom name `"C1"` and omit atoms near head
- **WHEN** `write_prepgen_control_file(path, variant="tail", head_name="C1", omit_names=["H1"], head_type="c3")`
- **THEN** the file contains HEAD_NAME and PRE_HEAD_TYPE but no TAIL_NAME
