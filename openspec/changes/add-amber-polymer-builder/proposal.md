# Change: Add AmberPolymerBuilder for tleap-based polymer construction

## Why

Currently, `PolymerBuilder` uses the `Reacter` system for step-by-step chemical reactions to build polymers. However, when using AmberTools as a backend, we cannot perform incremental reactions—tleap builds the entire polymer in one shot via its `sequence` command. We need a dedicated builder that:

1. Translates CGSmiles notation to tleap's sequence/branch commands
2. Wraps the entire AmberTools workflow (antechamber → parmchk2 → prepgen → tleap)
3. Keeps API consistency with `PolymerBuilder` while leveraging tleap's native polymer building

## What Changes

### New Components

- **`AmberPolymerBuilder`** (`src/molpy/builder/polymer/amber_builder.py`): A polymer builder that uses tleap as backend instead of Reacter
- **`TLeapScriptGenerator`** (`src/molpy/builder/polymer/tleap_script.py`): Translates CGSmiles IR to tleap script commands
- **`PrepgenWrapper`** improvements: Enhanced API for generating residue templates (.prepi) with HEAD/TAIL/CHAIN definitions
- **Monomer preparation utilities**: Functions to generate prepgen control files from port-marked Atomistic

### Workflow

```
User provides:
├── Monomer library: Dict[str, Atomistic] (with ports marked)
├── CGSmiles string: e.g., "{[#EO]|10}"
└── Workdir for intermediate files

AmberPolymerBuilder orchestrates:
├── Step 1: For each monomer type, prepare residue templates
│   ├── Write monomer as PDB (internal)
│   ├── antechamber → assign atom types and charges
│   ├── parmchk2 → generate missing parameters
│   └── prepgen → generate .prepi files (HEAD, CHAIN, TAIL variants)
│
├── Step 2: Parse CGSmiles and generate tleap script
│   ├── Translate CGSmiles graph to tleap sequence command
│   └── Handle terminal residues (first=HEAD, middle=CHAIN, last=TAIL)
│
├── Step 3: Execute tleap
│   └── Run tleap script → prmtop + inpcrd
│
└── Step 4: Read back results
    └── read_amber_prmtop → Frame + ForceField
```

### API Design

```python
from molpy.builder.polymer import AmberPolymerBuilder
from molpy.wrapper import TLeapWrapper

# 1. Prepare monomers (ports already marked externally)
eo_monomer = ...  # Atomistic with port="<" on head, port=">" on tail

# 2. Create builder
builder = AmberPolymerBuilder(
    library={"EO": eo_monomer},
    workdir=Path("./amber_work"),
    force_field="gaff2",  # or "gaff"
)

# 3. Build polymer (consistent with PolymerBuilder.build())
result = builder.build("{[#EO]|10}")

# 4. Access results
frame = result.frame        # Frame with coordinates
ff = result.force_field     # ForceField with parameters
```

### Key Design Decisions

1. **No Reacter involvement**: tleap handles all bonding internally via `sequence` command
2. **Port semantics preserved**: Ports on Atomistic map to HEAD_NAME/TAIL_NAME in prepgen control files
3. **All intermediate files in workdir**: User never touches .ac, .prepi, .frcmod files directly
4. **Consistent result type**: Returns a result object similar to `PolymerBuildResult` but with `Frame` + `ForceField`

## Impact

- **Affected code**: 
  - `src/molpy/builder/polymer/` (new files)
  - `src/molpy/wrapper/prepgen.py` (enhanced API)
  - `src/molpy/wrapper/tleap.py` (possibly minor additions)
- **New dependencies**: None (AmberTools is external, already wrapped)
- **Breaking changes**: None (additive change)
