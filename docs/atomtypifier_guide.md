# OplsAtomTypifier - SMARTS-based Atom Typing

## Overview

`OplsAtomTypifier` is a wrapper around the SMARTS matcher that enables atom typing for molecular structures using pattern-based matching. It assigns atom types from a `ForceField` to atoms in an `Atomistic` structure based on SMARTS patterns.

## Architecture

```
ForceField (with AtomTypes)
    ↓
OplsAtomTypifier._extract_patterns()
    ↓
List[SMARTSGraph] patterns
    ↓
SmartsMatcher
    ↓
find_candidates() → resolve()
    ↓
atom.data["type"] = atomtype
```

## Key Components

### 1. OplsAtomTypifier

The main class that:
- Extracts or creates SMARTS patterns from a `ForceField`
- Uses `SmartsMatcher` for pattern matching
- Applies results to `Atomistic` structures

### 2. Pattern Extraction

The `_extract_patterns()` method converts `AtomType` objects into `SMARTSGraph` patterns. There are two approaches:

**Approach A: Store patterns in ForceField**
```python
# Create patterns manually
patterns = [pattern1, pattern2, pattern3]

# Attach to forcefield
ff._atom_patterns = patterns

# Custom typifier uses them
class CustomOplsAtomTypifier(OplsAtomTypifier):
    def _extract_patterns(self):
        if hasattr(self.ff, '_atom_patterns'):
            return self.ff._atom_patterns
        return super()._extract_patterns()
```

**Approach B: Extract from AtomType metadata**
```python
# Store SMARTS or pattern info in AtomType
at = AtomType("CT", element="C", smarts="[C;X4]", priority=0)

# Typifier extracts patterns from metadata
def _extract_patterns(self):
    patterns = []
    for at in self.ff.get_types(AtomType):
        smarts = at.get("smarts")
        if smarts:
            pattern = parse_smarts(smarts, at.name, at.get("priority", 0))
            patterns.append(pattern)
    return patterns
```

## Usage Examples

### Basic Usage

```python
from molpy.core.atomistic import Atomistic
from molpy.core.forcefield import ForceField, AtomType, Style
from molpy.typifier.atomistic import OplsAtomTypifier
from molpy.typifier.builder import quick_pattern

# Create forcefield
ff = ForceField(name="my_ff")

# Define atom types
at_c = AtomType("CT", element="C", priority=0, mass=12.01, charge=0.0)
at_o = AtomType("OH", element="O", priority=5, mass=16.00, charge=-0.5)

# Add to forcefield
atom_style = Style("atom")
atom_style.types.add(at_c)
atom_style.types.add(at_o)
ff.styles.add(atom_style)

# Create patterns (manual approach)
patterns = [
    quick_pattern("CT", "C", priority=0),
    quick_pattern("OH", "O", priority=5)
]
ff._atom_patterns = patterns

# Create custom typifier
class MyTypifier(OplsAtomTypifier):
    def _extract_patterns(self):
        if hasattr(self.ff, '_atom_patterns'):
            return self.ff._atom_patterns
        return super()._extract_patterns()

# Create molecule
mol = Atomistic()
c = mol.add_atom(element="C", is_aromatic=False, charge=0, hyb=3)
o = mol.add_atom(element="O", is_aromatic=False, charge=0)
mol.add_bond(c, o, order=1)

# Apply typifier
typifier = MyTypifier(ff)
typifier.typify(mol)

# Check results
print(c.get("type"))  # "CT"
print(o.get("type"))  # "OH"
```

### Advanced Pattern Matching

```python
from molpy.typifier.builder import PatternBuilder
from molpy.typifier.predicates import is_element, bond_order, degree

# Create complex pattern for alcohol carbon
pb = PatternBuilder("C_OH", priority=2, source="alcohol_pattern")
v_c = pb.add_vertex(preds=[is_element("C"), degree(4)])
v_o = pb.add_vertex(preds=[is_element("O")])
pb.add_edge(v_c, v_o, preds=[bond_order(1)])
pb.set_target_vertices([v_c])  # Only type the carbon
pattern = pb.build()

# Add to forcefield patterns
ff._atom_patterns.append(pattern)
```

### Full Atomistic Typing

To type an entire structure (atoms, bonds, angles, dihedrals):

```python
from molpy.typifier.atomistic import OplsAtomisticTypifier

# Create full typifier
full_typifier = OplsAtomisticTypifier(ff)

# This will type atoms, bonds, angles, and dihedrals
full_typifier.typify(mol)
```

## Priority Resolution

When multiple patterns match the same atom, conflicts are resolved using:

1. **Priority**: Higher priority patterns win
2. **Specificity Score**: More specific patterns win
3. **Pattern Size**: Larger patterns win
4. **Definition Order**: Later defined patterns win

Example:
```python
# Generic carbon (priority=0)
p1 = quick_pattern("CT", "C", priority=0)

# Alcohol carbon (priority=2)
p2 = quick_pattern("C_OH", "C", priority=2)

# If both match, C_OH wins due to higher priority
```

## Implementation Notes

### Current Limitations

1. **Pattern Extraction**: The default `_extract_patterns()` implementation is minimal. You need to either:
   - Subclass and override `_extract_patterns()`
   - Store patterns directly in `ff._atom_patterns`
   - Extend `AtomType` to include SMARTS strings or pattern metadata

2. **SMARTS Parsing**: There's no built-in SMARTS string parser yet. Patterns must be built programmatically using `PatternBuilder`.

### Extension Points

To integrate with existing forcefield formats:

```python
class OplsAtomTypifier(OplsAtomTypifier):
    def _extract_patterns(self):
        """Extract patterns from OPLS forcefield format."""
        patterns = []
        
        for at in self.ff.get_types(AtomType):
            # Parse OPLS-specific metadata
            smarts = at.get("smarts") or at.get("definition")
            priority = at.get("priority", 0)
            
            if smarts:
                # Parse SMARTS string (requires implementation)
                pattern = self._parse_smarts(smarts, at.name, priority)
                patterns.append(pattern)
            else:
                # Fallback to element-only matching
                element = at.get("element")
                if element:
                    pattern = quick_pattern(at.name, element, priority)
                    patterns.append(pattern)
        
        return patterns
    
    def _parse_smarts(self, smarts, name, priority):
        """Parse SMARTS string into SMARTSGraph pattern."""
        # TODO: Implement SMARTS parser
        raise NotImplementedError("SMARTS parsing not yet implemented")
```

## Testing

See `tests/test_typifier/test_atom_typifier.py` for comprehensive test examples.

Run tests:
```bash
pytest tests/test_typifier/test_atom_typifier.py -v
```

## See Also

- `matcher_example.py` - SMARTS matcher basics
- `atom_typifier_example.py` - Complete usage example
- `SmartsMatcher` - Core matching engine
- `PatternBuilder` - Pattern construction utilities
