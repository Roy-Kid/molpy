# Design: AmberPolymerBuilder Architecture

## Overview

This document captures architectural decisions for integrating AmberTools-based polymer building into molpy while maintaining API consistency with the existing `PolymerBuilder`.

## Component Relationships

```
┌─────────────────────────────────────────────────────────────────────┐
│                        AmberPolymerBuilder                          │
│  (Orchestrates entire workflow, manages workdir, caches templates)  │
└─────────────────────────────────────────────────────────────────────┘
                                   │
           ┌───────────────────────┼───────────────────────┐
           ▼                       ▼                       ▼
┌─────────────────────┐  ┌─────────────────────┐  ┌─────────────────────┐
│  MonomerPreparer    │  │ TLeapScriptGenerator│  │   Result Reader     │
│  (per-monomer prep) │  │ (CGSmiles → tleap)  │  │ (prmtop → Frame)    │
└─────────────────────┘  └─────────────────────┘  └─────────────────────┘
           │                       │                       │
           ▼                       ▼                       ▼
┌─────────────────────┐  ┌─────────────────────┐  ┌─────────────────────┐
│ AntechamberWrapper  │  │   TLeapWrapper      │  │  read_amber_prmtop  │
│ Parmchk2Wrapper     │  │                     │  │  (existing I/O)     │
│ PrepgenWrapper      │  │                     │  │                     │
└─────────────────────┘  └─────────────────────┘  └─────────────────────┘
```

## Key Design Decisions

### 1. Monomer Library Contract

The `library` parameter follows the same contract as `PolymerBuilder`:
- Keys are CGSmiles labels (e.g., "EO", "PS")
- Values are `Atomistic` objects with ports marked on atoms

**Port convention for AmberPolymerBuilder:**
- Use molpy's standard port markers (e.g., `port="1"`, `port="2"` or `ports=["<", ">"]`)
- The builder internally translates to AmberTools terminology:
  - Port on **left side** of sequence → `HEAD_NAME` in prepgen
  - Port on **right side** of sequence → `TAIL_NAME` in prepgen
- Atoms to be removed during linking → `OMIT_NAME` entries (leaving groups)

**Translation rule:**
For a monomer with two ports, the position in the CGSmiles sequence determines which port is HEAD vs TAIL:
- First monomer in sequence: only TAIL port is active (HEAD variant)
- Middle monomers: both ports active (CHAIN variant)
- Last monomer in sequence: only HEAD port is active (TAIL variant)

### 2. Residue Template Generation Strategy

For each monomer type, we generate three variants:
- **HEAD variant** (e.g., `HPT.prepi`): Only TAIL connection point, used for chain start
- **CHAIN variant** (e.g., `PEO.prepi`): Both HEAD and TAIL connection points
- **TAIL variant** (e.g., `TPT.prepi`): Only HEAD connection point, used for chain end

Naming convention: `{label}.prepi`, `H{label}.prepi`, `T{label}.prepi`

### 3. CGSmiles to tleap Translation

| CGSmiles Pattern | tleap Command |
|------------------|---------------|
| `{[#A][#B][#C]}` | `mol = sequence {HA A B TC}` |
| `{[#A]\|5}` | `mol = sequence {HA A A A TA}` |
| `{[#A]([#B])[#C]}` | Branch handling (see below) |

**Branch handling:**
tleap's `sequence` creates linear chains. For branches, we:
1. Build main chain with placeholder
2. Use `bond` command to attach branches

Example for `{[#A]([#B])[#C]}`:
```
main = sequence {HA A TC}
branch = sequence {HB TB}
mol = combine {main branch}
bond mol.1.C3 mol.2.C1  # Connect branch point
```

### 4. Workdir Structure

```
workdir/
├── monomers/
│   ├── EO/
│   │   ├── EO.pdb           # Input structure
│   │   ├── EO.ac            # Antechamber output
│   │   ├── EO.mol2          # Typed structure
│   │   ├── EO.frcmod        # Missing parameters
│   │   ├── EO.chain         # Prepgen control (chain)
│   │   ├── EO.head          # Prepgen control (head)
│   │   ├── EO.tail          # Prepgen control (tail)
│   │   ├── EO.prepi         # Chain residue template
│   │   ├── HEO.prepi        # Head residue template
│   │   └── TEO.prepi        # Tail residue template
│   └── PS/
│       └── ...
├── tleap.in                  # Generated script
├── tleap.out                 # tleap output log
├── polymer.prmtop            # Final topology
└── polymer.inpcrd            # Final coordinates
```

### 5. Caching Strategy

To avoid redundant computation:
- Check if `.prepi` files exist before regenerating
- Store monomer preparation state in builder
- Option: `force_rebuild=True` to regenerate all

### 6. Error Handling

| Stage | Failure Mode | Handling |
|-------|--------------|----------|
| antechamber | Charge calculation fails | Raise `MonomerPreparationError` with sqm.out contents |
| parmchk2 | Missing parameters | Log warning, proceed (tleap may still work) |
| prepgen | Invalid port mapping | Raise `PortMappingError` with details |
| tleap | Build fails | Raise `TLeapBuildError` with tleap.out contents |

### 7. Result Type

```python
@dataclass
class AmberBuildResult:
    """Result of AmberPolymerBuilder.build()."""
    
    frame: Frame                    # Coordinates and topology
    force_field: ForceField         # Parameters
    prmtop_path: Path               # Path to prmtop file
    inpcrd_path: Path               # Path to inpcrd file
    tleap_script: str               # Generated script (for debugging)
    monomer_count: int              # Total monomers in polymer
    connection_count: int           # Number of bonds formed
```

## Trade-offs Considered

### Alternative A: Reuse Reacter with tleap backend
- **Rejected**: Reacter assumes step-by-step reactions; tleap builds all at once
- tleap's `sequence` command is more efficient than iterative bonding

### Alternative B: Subclass PolymerBuilder
- **Rejected**: Different internal workflow (no Reacter, no incremental building)
- Would require overriding most methods, defeating purpose of inheritance

### Alternative C: Separate monomer prep from building
- **Considered**: User calls `prepare_monomer()` separately, then `build()`
- **Decision**: Internal preparation is cleaner; user only provides Atomistic with ports

## Interface Consistency

To maintain consistency with `PolymerBuilder`:

| PolymerBuilder | AmberPolymerBuilder |
|----------------|---------------------|
| `library: Mapping[str, Atomistic]` | Same |
| `build(cgsmiles: str)` | Same |
| `PolymerBuildResult.polymer` | `AmberBuildResult.frame` |
| Uses `ReacterConnector` | Uses `TLeapWrapper` internally |

## Future Extensions

1. **Branched polymers**: Full support for tleap's `bond` command (deferred)
2. **Multi-chain systems**: `combine` multiple sequences
3. **Mixed force fields**: GAFF + ff19SB for protein-polymer systems
4. **Copolymer support**: Alternating/random sequences
