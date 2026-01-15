# Tasks: Add AmberPolymerBuilder

## Prerequisites
- [x] Verify AmberTools wrappers (antechamber, parmchk2, prepgen, tleap) work correctly
- [x] Ensure `read_amber_prmtop` returns correct `Frame` + `ForceField`

---

## 1. Enhance PrepgenWrapper API
**Goal**: Add high-level method for generating residue templates with HEAD/TAIL/CHAIN control files.

- [x] Add `PrepgenWrapper.generate_residue()` method that accepts:
  - Input `.ac` or `.mol2` file
  - Output `.prepi` file
  - Control file (`.chain`, `.head`, or `.tail`)
  - Residue name
- [x] Add utility function `write_prepgen_control_file()` to generate control files from port info
- [x] Add tests for `PrepgenWrapper.generate_residue()` *(covered in wrapper tests)*

**Validation**: Unit test that generates a `.prepi` file from a simple molecule ✅

**Note**: Implementation in `src/molpy/wrapper/prepgen.py`

---

## 2. Create MonomerPreparer class
**Goal**: Encapsulate the antechamber → parmchk2 → prepgen workflow for a single monomer.

- [x] ~~Create `MonomerPreparer` class~~ *Integrated directly into `AmberPolymerBuilder._prepare_monomers()` method*
- [x] Implement `prepare()` method:
  - Write `Atomistic` to PDB
  - Run antechamber for atom typing
  - Run parmchk2 for missing parameters
  - Generate control files from port markers
  - Run prepgen for HEAD, CHAIN, TAIL variants
- [x] Return `_PreparedMonomer` dataclass with paths to all generated files
- [x] Add tests for monomer preparation (requires AmberTools, marked as external)

**Validation**: Integration test verifies .prepi files exist ✅

**Note**: Merged into `AmberPolymerBuilder` in `src/molpy/builder/ambertools/amber_builder.py`

---

## 3. Create TLeapScriptGenerator
**Goal**: Translate CGSmiles IR to tleap script.

- [x] ~~Create `TLeapScriptGenerator` class~~ *Integrated directly into `AmberPolymerBuilder._build_with_tleap()` method*
- [x] Implement script generation:
  - Accept `CGSmilesGraphIR` and monomer template paths
  - Output tleap script string
- [x] Handle linear sequences: first=HEAD, middle=CHAIN, last=TAIL
- [x] Handle repeat operator `|N`: expand to N monomers
- [x] Add tests for script generation (pure Python, no AmberTools needed)

**Validation**: Unit test in `tests/test_builder/test_amber_polymer_builder.py` verifies correct sequence for `{[#EO]|5}` ✅

**Note**: Merged into `AmberPolymerBuilder` in `src/molpy/builder/ambertools/amber_builder.py`

---

## 4. Create AmberBuildResult dataclass
**Goal**: Define result type consistent with PolymerBuildResult.

- [x] Create `AmberBuildResult` in `src/molpy/builder/ambertools/types.py`
- [x] Include: `frame`, `forcefield`, `prmtop_path`, `inpcrd_path`, `pdb_path`, `monomer_count`, `cgsmiles`
- [x] Add docstrings explaining each field

**Validation**: Import and instantiate in test ✅ (`TestAmberBuildResult::test_fields`)

---

## 5. Create AmberPolymerBuilder main class
**Goal**: Orchestrate the entire workflow.

- [x] Create `AmberPolymerBuilder` in `src/molpy/builder/ambertools/amber_builder.py`
- [x] Implement `__init__()`:
  - Accept `library`, `config` (includes workdir, force_field gaff/gaff2)
  - Create wrapper instances internally (on-demand)
- [x] Implement `build(cgsmiles: str) -> AmberBuildResult`:
  - Parse CGSmiles using `parse_cgsmiles()`
  - Prepare monomers (with caching)
  - Generate tleap script
  - Execute tleap via wrapper
  - Read results using `read_amber_prmtop()`
- [x] Add `_prepare_monomers()` internal method with caching
- [x] Add `_validate_ir()` to check port markers and label existence

**Validation**: Integration tests pass ✅

---

## 6. Add exports and documentation
**Goal**: Make AmberPolymerBuilder accessible and documented.

- [x] Export from `src/molpy/builder/ambertools/__init__.py`
- [x] Add docstrings with usage examples
- [x] Update `src/molpy/builder/__init__.py` to export `AmberPolymerBuilder`, `AmberPolymerBuilderConfig`

**Validation**: `from molpy.builder.ambertools import AmberPolymerBuilder` works ✅

---

## 7. Integration tests
**Goal**: End-to-end tests with real AmberTools.

- [x] Test config validation
- [x] Test sequence generation for linear homopolymer
- [x] Test error handling: invalid CGSmiles, missing port markers, missing labels
- [x] Mark all as `external` (require AmberTools) in `tests/conftest.py`

**Validation**: 9 tests pass in `tests/test_builder/test_amber_polymer_builder.py` ✅

---

## 8. (Future) Branch support
**Goal**: Handle branched polymers. **Deferred until linear case works.**

- [ ] Extend sequence generation for branch syntax
- [ ] Generate `bond` commands for branch connections
- [ ] Add tests for branched structures

**Validation**: Test `{[#A]([#B])[#C]}` builds correctly

---

## Scope for Initial Implementation

Focus on the PEO example:
- ✅ Linear homopolymer: `{[#EO]|N}`
- ✅ Single monomer type with two ports
- ✅ Generate HEAD, CHAIN, TAIL variants
- ⏳ Branch support is explicitly deferred

---

## Implementation Notes

### Module Location Change
The original proposal specified `src/molpy/builder/polymer/`. Implementation uses `src/molpy/builder/ambertools/` to better reflect the AmberTools-specific nature of this builder.

### Merged Classes
- `MonomerPreparer` merged into `AmberPolymerBuilder._prepare_monomers()`
- `TLeapScriptGenerator` merged into `AmberPolymerBuilder._build_with_tleap()`

This simplification reduces file count while maintaining the same functionality.

### Port Convention
- `port="<"` for head connection point
- `port=">"` for tail connection point
- Direct dictionary access (`atom["port"]`) - missing port raises KeyError
