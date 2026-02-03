## 1. SMARTS Parser Enhancement
- [ ] 1.1 Add `x<n>` (ring connectivity) primitive to grammar and transformer
- [ ] 1.2 Add `A` (aliphatic) primitive to grammar and transformer
- [ ] 1.3 Add isotope mass prefix support to grammar (`[12C]`, `[2H]`)
- [ ] 1.4 Add chirality (`@`, `@@`) support to grammar
- [ ] 1.5 Fix charge syntax to support `++`, `--`, `+2`, `-2` forms
- [ ] 1.6 Update `SmartsIR` types enum for new primitives (`ring_connectivity`, `aliphatic`, `isotope`, `chirality`)
- [ ] 1.7 Update `SMARTSGraph._atom_primitive_matches()` to evaluate new primitives against molecule graph
- [ ] 1.8 Update `adapter.py` `build_mol_graph()` to populate new vertex attributes needed by new primitives (ring connectivity, aliphatic flag)
- [ ] 1.9 Validate all 571 GAFF SMARTS patterns from `gaff.prm` parse without error
- [ ] 1.10 Write unit tests for new SMARTS primitives (grammar parsing + IR generation + round-trip)

## 2. GAFF XML Data Conversion
- [ ] 2.1 Write conversion script `scripts/convert_gaff.py` that reads `gaff.prm` + `gaff.dat` and outputs `gaff.xml`
- [ ] 2.2 Parse `gaff.prm`: extract atom type name, SMARTS pattern, description for each of 571 rules
- [ ] 2.3 Parse `gaff.dat`: extract atom masses, bond parameters, angle parameters, dihedral parameters, improper parameters, VDW parameters
- [ ] 2.4 Implement unit conversions (AMBER → molpy units: kcal→kJ, Å→nm, R*→σ)
- [ ] 2.5 Generate `<AtomTypes>` section with SMARTS `def`, computed priority from rule order
- [ ] 2.6 Generate `<HarmonicBondForce>` section
- [ ] 2.7 Generate `<HarmonicAngleForce>` section
- [ ] 2.8 Generate `<PeriodicTorsionForce>` section (new XML element for GAFF dihedrals)
- [ ] 2.9 Generate `<PeriodicImproperForce>` section (new XML element)
- [ ] 2.10 Generate `<NonbondedForce>` section (LJ parameters)
- [ ] 2.11 Write `gaff.xml` to `src/molpy/data/forcefield/gaff.xml`
- [ ] 2.12 Write tests verifying a sample of converted parameters match expected values

## 3. Forcefield XML Reader Extension
- [ ] 3.1 Extend XML reader to parse `<PeriodicTorsionForce>` element
- [ ] 3.2 Extend XML reader to parse `<PeriodicImproperForce>` element
- [ ] 3.3 Add `PeriodicTorsionType` and `ImproperType` to core forcefield types if not present
- [ ] 3.4 Write tests for reading `gaff.xml` metadata, atom types, and all parameter sections

## 4. Typifier Refactoring
- [ ] 4.1 Extract `ForceFieldAtomTypifier` base class from `OplsAtomTypifier` with shared SMARTS-based typing logic
- [ ] 4.2 Extract `ForceFieldBondTypifier` base from `OplsBondTypifier`
- [ ] 4.3 Extract `ForceFieldAngleTypifier` base from `OplsAngleTypifier`
- [ ] 4.4 Extract `ForceFieldDihedralTypifier` base from `OplsDihedralTypifier`
- [ ] 4.5 Extract `ForceFieldTypifier` base orchestrator from `OplsAtomisticTypifier`
- [ ] 4.6 Rewrite `OplsAtomisticTypifier` and sub-typifiers as subclasses of new bases
- [ ] 4.7 Verify all existing OPLS tests pass after refactoring
- [ ] 4.8 Remove any dead code or backward-compatibility shims

## 5. GaffTypifier Implementation
- [ ] 5.1 Implement `GaffAtomTypifier(ForceFieldAtomTypifier)` — priority from rule order, no type references
- [ ] 5.2 Implement `GaffBondTypifier(ForceFieldBondTypifier)` — wildcard `X` class matching
- [ ] 5.3 Implement `GaffAngleTypifier(ForceFieldAngleTypifier)` — wildcard `X` class matching
- [ ] 5.4 Implement `GaffDihedralTypifier(ForceFieldDihedralTypifier)` — periodic torsions, wildcard matching
- [ ] 5.5 Implement `GaffImproperTypifier` — periodic improper matching
- [ ] 5.6 Implement `GaffAtomisticTypifier(ForceFieldTypifier)` — orchestrates all GAFF sub-typifiers
- [ ] 5.7 Write unit tests for each GAFF sub-typifier with small molecules
- [ ] 5.8 Write integration test: typify ethanol, acetic acid, benzene with GAFF and verify all atom types match expected GAFF assignments

## Dependencies
- Task group 1 (SMARTS parser) has no dependencies; can start immediately
- Task group 2 (XML conversion) depends on 1.9 (all GAFF patterns must parse)
- Task group 3 (XML reader) can run in parallel with group 2
- Task group 4 (refactoring) has no dependencies; can start immediately
- Task group 5 (GaffTypifier) depends on groups 2, 3, and 4
