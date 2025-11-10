# SMARTS Matcher Implementation

A **minimal, deterministic** SMARTS matcher for atomtyping using predicate-based pattern matching and VF2 subgraph isomorphism.

## Overview

This implementation provides:
1. **Predicate-based pattern matching** - No SMARTS parser dependency in matcher
2. **Deterministic conflict resolution** - Total ordering: priority → specificity → size → definition order
3. **Programmatic pattern construction** - Builder API for tests and programmatic use
4. **VF2 subgraph isomorphism** - Using igraph's efficient implementation
5. **Explain mode** - JSON debugging output showing all candidates and resolution

## Architecture

### Core Components

```
src/molpy/typifier/
├── predicates.py    # Predicate factories (is_element, bond_order, etc.)
├── builder.py       # PatternBuilder for programmatic construction
├── graph.py         # SMARTSGraph (refactored to support predicates)
├── matcher.py       # SmartsMatcher with VF2 and conflict resolution
└── adapter.py       # Molecule → igraph.Graph conversion
```

### Data Flow

```
Molecule (Atomistic)
    ↓
build_mol_graph() → igraph.Graph with attributes
    ↓                  (element, is_aromatic, charge, degree, etc.)
    ↓
SmartsMatcher.find_candidates()
    ↓
    ├─ For each pattern:
    │   ├─ VF2 subgraph isomorphism
    │   ├─ Create Candidate objects
    │   └─ Respect target_vertices
    ↓
List[Candidate]
    ↓
SmartsMatcher.resolve()
    ↓
    ├─ Group by atom_id
    ├─ Sort by (priority, score, size, def_order)
    └─ Pick winner
    ↓
Dict[atom_id, atomtype]
```

## Usage

### 1. Building Patterns Programmatically

```python
from molpy.typifier.builder import PatternBuilder, quick_pattern
from molpy.typifier.predicates import is_element, charge, bond_order

# Quick single-atom pattern
pattern1 = quick_pattern("t_O_neg", "O", priority=0, charge=-1)

# Multi-atom pattern with bonds
builder = PatternBuilder("t_CO_single", priority=0)
c = builder.add_vertex(preds=[is_element("C")])
o = builder.add_vertex(preds=[is_element("O")])
builder.add_edge(c, o, preds=[bond_order(1)])
pattern2 = builder.build()

# Target-vertex pattern (only type specific atoms)
builder = PatternBuilder("t_hydroxyl_O", priority=5)
c = builder.add_vertex(preds=[is_element("C")])
o = builder.add_vertex(preds=[is_element("O")])
builder.add_edge(c, o, preds=[bond_order(1)])
builder.set_target_vertices([o])  # Only type the oxygen
pattern3 = builder.build()
```

### 2. Matching Molecules

```python
from molpy.core.atomistic import Atomistic
from molpy.typifier.matcher import SmartsMatcher
from molpy.typifier.adapter import build_mol_graph

# Create molecule
mol = Atomistic()
c1 = mol.add_atom(element="C", charge=0)
o = mol.add_atom(element="O", charge=-1)
mol.add_bond(c1, o, order=1)

# Build graph
graph, vs_to_atomid, atomid_to_vs = build_mol_graph(mol)

# Match patterns
patterns = [pattern1, pattern2, pattern3]
matcher = SmartsMatcher(patterns)

# Find candidates
candidates = matcher.find_candidates(graph, vs_to_atomid)

# Resolve conflicts
result = matcher.resolve(candidates)
print(result)  # {atom_id: "atomtype"}
```

### 3. Explain Mode (Debugging)

```python
# Get detailed explanation
explain = matcher.explain(candidates)

# Export as JSON
import json
print(json.dumps(explain, indent=2))
```

Output:
```json
{
  "atom_id": {
    "winner": "t_O_neg",
    "candidates": [
      {
        "atomtype": "t_O_neg",
        "priority": 0,
        "score": 1,
        "pattern_size": [1, 0],
        "rank": 1
      }
    ]
  }
}
```

## Predicate System

### Available Predicates

**Vertex (Atom) Predicates:**
- `is_element(symbol: str)` - Element matching (weight=0, baseline)
- `is_aromatic(value: bool)` - Aromaticity (weight=2)
- `charge(value: int)` - Formal charge (weight=1)
- `degree(value: int)` - Connectivity degree (weight=1)
- `hyb(value: int)` - Hybridization (1=sp, 2=sp2, 3=sp3) (weight=1)
- `in_ring(value: bool)` - Ring membership (weight=2)
- `ring_size(size: int)` - Specific ring size (weight=3)
- `ring_count(count: int)` - Number of rings (weight=2)
- `custom_vertex(func, name, weight)` - Custom predicate (weight=4)

**Edge (Bond) Predicates:**
- `bond_order(order: int | str)` - Bond order (1, 2, 3, ":") (weight=3)
- `bond_aromatic(value: bool)` - Aromatic bond (weight=2)
- `bond_in_ring(value: bool)` - Bond in ring (weight=2)

### Specificity Scoring

Patterns are scored based on predicate weights:
```
score = Σ(vertex_predicate_weights) + Σ(edge_predicate_weights)
```

Higher scores indicate more specific patterns.

## Conflict Resolution

Deterministic total ordering (highest wins):
1. **Priority** (user-defined, higher first)
2. **Specificity score** (computed from predicates)
3. **Pattern size** (# vertices + # edges, larger first)
4. **Definition order** (later-defined wins)

## Molecule Adapter

The adapter converts `Atomistic` structures to `igraph.Graph`:

### Vertex Attributes Set:
- `element`: str (e.g., "C", "N", "O")
- `number`: int (atomic number)
- `is_aromatic`: bool
- `charge`: int
- `degree`: int (computed)
- `hyb`: int | None (1=sp, 2=sp2, 3=sp3)
- `in_ring`: bool (computed)
- `cycles`: set[tuple] (computed)

### Edge Attributes Set:
- `order`: int | str (1, 2, 3, ":")
- `is_aromatic`: bool
- `is_in_ring`: bool (computed)

## Testing

All tests are pure graph tests with no external chemistry libraries:

```bash
# Run all tests
pytest tests/test_typifier/test_matcher_*.py -v

# Run specific test suite
pytest tests/test_typifier/test_matcher_basic.py -v        # Tests 1-4
pytest tests/test_typifier/test_matcher_intermediate.py -v # Tests 5-8
pytest tests/test_typifier/test_matcher_advanced.py -v     # Tests 9-12
```

### Test Coverage

1. **Single-atom typing** - Element baseline
2. **Charge constraint** - Specificity beats plain
3. **Bond order** - C-O vs C=O
4. **Aromatic vs aliphatic** - Benzene ring
5. **Ring membership** - Cyclohexane vs hexane
6. **Hybridization/degree** - sp2 vs sp3
7. **Branch constraint** - Isopropyl central carbon
8. **Target-vertex indexing** - C-O pattern typing only O
9. **Overlapping matches** - Toluene with multiple patterns
10. **Priority overrides** - High priority beats specificity
11. **Edge-only specificity** - Formaldehyde vs methanol
12. **Performance** - 200-atom alkane chain

**All 13 tests pass** (including 4 existing + 9 new)

## Example

See `examples/matcher_example.py` for a complete working example with ethanol atomtyping:

```bash
python examples/matcher_example.py
```

Output demonstrates:
- Pattern creation
- Molecule construction
- Candidate finding
- Conflict resolution
- JSON explain output

## Performance

- 200-atom linear alkane: < 0.2s (including ring detection)
- Deterministic results (same input → same output)
- Efficient VF2 implementation from igraph

## Future Enhancements

- [ ] Support for recursive SMARTS ($(...)patterns)
- [ ] More sophisticated ring detection (SSSR)
- [ ] Parallel pattern matching for large pattern sets
- [ ] Pattern caching and reuse
- [ ] Integration with SMARTS string parser (already exists)
- [ ] Support for atom/bond labels (:1, :2, etc.)

## API Reference

### PatternBuilder

```python
class PatternBuilder:
    def __init__(self, atomtype_name: str, priority: int = 0, source: str = "")
    def add_vertex(self, preds: list[VertexPredicate]) -> int
    def add_edge(self, u: int, v: int, preds: list[EdgePredicate]) -> None
    def set_target_vertices(self, vertices: list[int]) -> None
    def build(self) -> SMARTSGraph
```

### SmartsMatcher

```python
class SmartsMatcher:
    def __init__(self, patterns: list[SMARTSGraph], scoring: ScoringPolicy | None = None)
    def find_candidates(self, mol_graph: Graph, vs_to_atomid: dict[int, int]) -> list[Candidate]
    def resolve(self, candidates: list[Candidate]) -> dict[int, str]
    def explain(self, candidates: list[Candidate]) -> dict[int, Any]
```

### Candidate

```python
@dataclass
class Candidate:
    atom_id: int
    atomtype: str
    source: str
    priority: int
    score: int
    pattern_size: tuple[int, int]
    definition_order: int
```

## License

Part of the molpy project.
