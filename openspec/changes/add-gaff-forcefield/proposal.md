# Change: Add GAFF Force Field Support

## Why
molpy currently only supports OPLS-AA for automatic atom typing. GAFF (Generalized Amber Force Field) is widely used in the AMBER ecosystem for organic molecules. Adding GAFF support requires: (1) converting GAFF parameters to molpy's XML format, (2) extending the SMARTS parser to handle GAFF's SMARTS patterns (which use features not yet supported), and (3) implementing a `GaffTypifier` that follows the same layered architecture as `OplsAtomisticTypifier` but adapted for GAFF's typing rules.

## What Changes

### 1. GAFF Forcefield XML Data
- Convert `gaff.prm` (571 SMARTS atom type definitions) and `gaff.dat` (bond, angle, dihedral, improper, VDW parameters) into `src/molpy/data/forcefield/gaff.xml` following the same XML schema as `oplsaa.xml`
- GAFF uses class-based parameter lookup (e.g., atom types `c3`, `ca`, `ct` map to classes for bond/angle/dihedral matching)
- GAFF uses `combining_rule="arithmetic"` (Lorentz-Berthelot) vs OPLS-AA's `"geometric"`

### 2. SMARTS Parser Enhancement
- **BREAKING**: Refactor grammar to properly support the full Daylight SMARTS specification
- Add missing primitives: `x<n>` (ring connectivity), `A` (aliphatic), isotope mass, chirality (`@`, `@@`)
- Fix charge syntax: support `++`, `--`, `+2`, `-2` forms
- Add ring bond primitive `@` in bond context (distinct from chirality `@` in atom context)
- Add directional bond unspecified forms `/?`, `\?`
- Ensure atom list notation `[C,N,O]` inside brackets works correctly with existing OR logic
- Ensure implicit AND between adjacent primitives in bracket atoms (e.g., `[#6X4]` = atomic number 6 AND connectivity 4)

### 3. GaffTypifier
- Implement `GaffAtomTypifier` following the same `LayeredTypingEngine` pattern as OPLS
- Implement `GaffBondTypifier`, `GaffAngleTypifier`, `GaffDihedralTypifier`, `GaffImproperTypifier`
- Implement `GaffAtomisticTypifier` as the top-level orchestrator
- **BREAKING**: Refactor shared typifier infrastructure to eliminate duplication between OPLS and GAFF implementations; extract common base classes

### 4. Refactoring (no backward compatibility)
- Unify typifier base classes so OPLS and GAFF share parameter-lookup logic
- Refactor `atomistic.py` to extract force-field-agnostic logic into base classes
- Clean up any dead code paths

## Impact
- Affected specs: new capabilities `gaff-forcefield`, `smarts-parser`, `gaff-typifier`
- Affected code:
  - `src/molpy/parser/grammar/smarts.lark` — grammar extension
  - `src/molpy/parser/smarts.py` — transformer updates for new primitives
  - `src/molpy/typifier/atomistic.py` — refactor to shared base + GAFF/OPLS specifics
  - `src/molpy/typifier/graph.py` — support new SMARTS primitives in matching
  - `src/molpy/data/forcefield/` — new `gaff.xml`
  - `tests/` — new test files for all three capabilities
