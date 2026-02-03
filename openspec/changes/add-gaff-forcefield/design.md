## Context

molpy's typifier system currently supports OPLS-AA only. The architecture is well-structured with a layered typing engine, SMARTS graph matching, and dependency analysis. Adding GAFF requires extending this system without duplicating logic. The SMARTS parser needs to support more of the Daylight specification since GAFF's `gaff.prm` uses patterns like `[#6X3]`, `[#1X1][C]([F,Cl,Br,I])`, and ring/aromatic primitives that partially overlap with existing support.

### Stakeholders
- Users building AMBER-based simulation pipelines
- Existing OPLS-AA users (refactoring must not break their workflows)

## Goals / Non-Goals

### Goals
- GAFF atom type assignment via SMARTS matching identical to antechamber's behavior
- Full bonded parameter assignment (bonds, angles, proper dihedrals, improper dihedrals)
- Nonbonded parameter assignment (LJ σ/ε, charges from atom types)
- SMARTS parser covers all primitives used in `gaff.prm` and common chemical informatics use cases
- Shared typifier infrastructure between OPLS and GAFF

### Non-Goals
- GAFF2 support (can be added later with same infrastructure)
- AM1-BCC charge calculation (requires QM; handled by external tools like antechamber)
- Automatic GAFF atom type assignment for metal-containing systems
- Full Daylight SMARTS spec compliance for reaction SMARTS (`>>`) or component-level grouping (`.`)

## Decisions

### 1. XML Format for GAFF

**Decision**: Use the same XML schema as `oplsaa.xml` with minor extensions for improper dihedrals.

**Rationale**: The existing schema already supports `AtomTypes` with SMARTS `def`, `HarmonicBondForce`, `HarmonicAngleForce`, and torsion forces. GAFF uses the same harmonic functional forms. The main addition is `PeriodicTorsionForce` (GAFF uses periodic dihedrals, not Ryckaert-Bellemans like OPLS) and `PeriodicImproperForce`.

**Format mapping**:
| GAFF section | XML element | Notes |
|---|---|---|
| Atom types (gaff.prm) | `<AtomTypes>` | SMARTS in `def`, class = type name |
| Bond params (gaff.dat) | `<HarmonicBondForce>` | k in kcal/mol/Å² → kJ/mol/nm² |
| Angle params (gaff.dat) | `<HarmonicAngleForce>` | k in kcal/mol/rad² → kJ/mol/rad² |
| Dihedral params (gaff.dat) | `<PeriodicTorsionForce>` | New element; V, phase, periodicity |
| Improper params (gaff.dat) | `<PeriodicImproperForce>` | New element |
| VDW params (gaff.dat) | `<NonbondedForce>` | R*/2 → σ, ε in kcal/mol → kJ/mol |

**Unit conversions** (AMBER → molpy internal SI-like units):
- Energy: 1 kcal/mol = 4.184 kJ/mol
- Length: Å → nm (÷10)
- Bond force constant: kcal/mol/Å² → kJ/mol/nm² (×4.184×100)
- Angle force constant: kcal/mol/rad² → kJ/mol/rad² (×4.184)
- VDW radius: R* (Å, half Rmin) → σ (nm) = R* × 2^(5/6) / 10

### 2. SMARTS Grammar Extension Strategy

**Decision**: Extend the existing Lark grammar incrementally, keeping it Earley-parsed.

**Alternatives considered**:
- Rewrite grammar in LALR: Rejected because SMARTS has inherent ambiguities (e.g., `C` as element vs atom class) that Earley handles naturally
- Use a hand-written recursive descent parser: Rejected; Lark is already working well and maintainable

**New grammar additions**:
```
atom_id: ...
       | "x" ring_connectivity?   // ring bond count
       | "A"                       // aliphatic
       | isotope atom_symbol       // isotope mass prefix
       | "@" "@"?                  // chirality
```

**GAFF-specific SMARTS patterns that drive requirements**:
- `[#6X3]` — atomic number + connectivity (already supported)
- `[#1X1][C]([F,Cl,Br,I])` — atom lists in brackets (already works via OR)
- `[SX2]=*` — wildcard with bond (supported)
- `[#7X3][CX3](=[OX1])` — nested neighbor spec (supported)
- `[#6r5]` — ring size (supported)
- No `%label` type references in GAFF (unlike OPLS) → simpler dependency graph
- No chirality patterns in GAFF → chirality can be low priority

### 3. Typifier Refactoring

**Decision**: Extract a `ForceFieldTypifier` base class from `OplsAtomisticTypifier` that handles the common workflow (atom typing → pair typing → bond typing → angle typing → dihedral typing). OPLS and GAFF subclasses override only what differs.

**Class hierarchy**:
```
ForceFieldTypifier (base)
├── OplsAtomisticTypifier
│   ├── OplsAtomTypifier (uses LayeredTypingEngine)
│   ├── OplsBondTypifier
│   ├── OplsAngleTypifier
│   └── OplsDihedralTypifier (RB torsions)
└── GaffAtomisticTypifier
    ├── GaffAtomTypifier (uses LayeredTypingEngine, same engine)
    ├── GaffBondTypifier
    ├── GaffAngleTypifier
    ├── GaffDihedralTypifier (periodic torsions)
    └── GaffImproperTypifier (periodic impropers)
```

**Key differences between OPLS and GAFF typing**:
| Aspect | OPLS-AA | GAFF |
|---|---|---|
| Atom type naming | `opls_135`, `opls_145` | `c3`, `ca`, `hc` |
| Class concept | Separate class field (e.g., `CT`, `CA`) | Class = type name (no separate class) |
| Type references | Yes (`%opls_145` in SMARTS) | No |
| Override mechanism | `overrides` attribute | Rule ordering (last match wins) |
| Torsion form | Ryckaert-Bellemans | Periodic (Fourier) |
| Improper dihedrals | Not in current XML | Periodic impropers |
| Combining rule | Geometric | Arithmetic (Lorentz-Berthelot) |

### 4. GAFF Atom Typing Priority

**Decision**: GAFF uses a "last match wins" strategy ordered by specificity. Convert this to priority values: rules later in `gaff.prm` get higher priority. The `LayeredTypingEngine` already supports priority-based conflict resolution.

Since GAFF has no `%label` type references, all patterns are level 0 (no dependencies). This means the `LayeredTypingEngine` processes everything in a single pass — simpler than OPLS.

## Risks / Trade-offs

- **Risk**: GAFF SMARTS in `gaff.prm` may use non-standard SMARTS extensions specific to OpenBabel.
  - **Mitigation**: Validate all 571 patterns parse correctly; flag and handle any non-standard patterns.

- **Risk**: Unit conversion errors in parameter conversion.
  - **Mitigation**: Write conversion tests comparing a few known molecules against AMBER reference values.

- **Risk**: Refactoring typifier base classes may break OPLS-AA typing.
  - **Mitigation**: Run all existing OPLS tests after refactoring; add regression tests.

- **Trade-off**: Not supporting GAFF2. The infrastructure will support it, but the data files and any GAFF2-specific typing rules are out of scope.

## Open Questions

1. Should `gaff.xml` include all ~571 atom types from `gaff.prm`, or only the commonly used organic subset? → **Proposed**: Include all; let users filter if needed.
2. Should the XML conversion be a one-time script or a reusable parser for AMBER `.dat`/`.prm` files? → **Proposed**: One-time script (store in `scripts/`), since GAFF updates are infrequent.
3. How to handle GAFF wildcard dihedral parameters (e.g., `X-c3-c3-X`)? → **Proposed**: Use `class="X"` as a wildcard marker, consistent with how OPLS uses `""` for wildcards.
