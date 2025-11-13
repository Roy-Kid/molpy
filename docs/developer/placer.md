# Placer System for Polymer Assembly

## Overview

The **Placer** system is responsible for positioning monomers during polymer assembly. It determines:

1. **Separation distance** between monomers (via `Separator`)
2. **Orientation** of monomers (via `Orienter`)

This ensures that monomers are spatially arranged before chemical reactions occur, resulting in more realistic polymer structures.

## Architecture

### Components

```
Placer
├── Separator     # Calculates distance between monomers
│   └── VdWSeparator    # Based on van der Waals radii
└── Orienter      # Determines monomer orientation
    └── LinearOrienter  # Linear alignment of monomers
```

### VdWSeparator

Calculates separation distance based on the **van der Waals radii** of port anchor atoms.

**Formula**:
```
separation = vdw_left + vdw_right + buffer
```

**Parameters**:
- `buffer` (float): Additional spacing in Angstroms (default: 0.0)

**VdW radii** are defined in `molpy.core.element.Element`:
- H: 1.20 Å
- C: 1.70 Å
- N: 1.55 Å
- O: 1.52 Å
- (and more...)

### LinearOrienter

Arranges monomers in a **linear fashion** by:

1. Calculating the direction vector from each port
2. Aligning the next monomer's port axis with the previous monomer's port axis
3. Positioning the monomer at the calculated separation distance

**Port Direction Calculation**:
- **Primary method**: Use the bond direction (from bonded neighbor to port anchor)
- **Fallback**: Use direction from monomer centroid to port anchor

## Usage

### Basic Usage

```python
from molpy.builder.placer import create_vdw_linear_placer
from molpy.builder.polymer import linear

# Create a VdW-based placer with 0.5 Å buffer
placer = create_vdw_linear_placer(buffer=0.5)

# Use in linear polymer assembly
polymer = linear(
    sequence="ABC",
    library=monomers,
    connector=connector,
    placer=placer,  # Apply placer!
)
```

### Advanced Usage

```python
from molpy.builder.placer import Placer, VdWSeparator, LinearOrienter

# Create custom placer
separator = VdWSeparator(buffer=1.0)  # 1.0 Å buffer
orienter = LinearOrienter()

placer = Placer(separator=separator, orienter=orienter)

# Use in assembly
polymer = linear(
    sequence="ABCD",
    library=monomers,
    connector=connector,
    placer=placer,
)
```

### Manual Usage (Outside linear() Function)

```python
# Get monomers
monomer_A = monomers["A"]
monomer_B = monomers["B"].copy()

# Get ports
port_A = monomer_A.ports["out"]
port_B = monomer_B.ports["in"]

# Apply placer
placer.place_monomer(monomer_A, monomer_B, port_A, port_B)

# Now monomer_B is positioned relative to monomer_A
# Its coordinates have been modified in-place
```

## How It Works

### Step-by-Step Process

1. **Calculate Separation**:
   - Get VdW radii of port anchor atoms (from `Element`)
   - Add buffer distance
   - Total separation = `vdw_left + vdw_right + buffer`

2. **Determine Port Directions**:
   - For each port, find the bonded neighbor
   - Direction = vector from neighbor to port anchor (normalized)
   - If no neighbor found, use direction from centroid

3. **Calculate Transformation**:
   - Target position = `left_anchor_pos + left_direction * separation`
   - Rotation = align `right_direction` with `-left_direction` (face each other)
   - Uses Rodrigues' rotation formula

4. **Apply Transformation**:
   - Rotate all atoms in right monomer around its port anchor
   - Translate all atoms to final position

### Example

```
Before placer:
  A: ... -C-O*         (port_1 at O)
  B:      *C-C- ...    (port_2 at C)
  (overlapping or random positions)

After placer (separation = vdw_O + vdw_C + buffer = 1.52 + 1.70 + 0.5 = 3.72 Å):
  A: ... -C-O*                        *C-C- ...
                 <-- 3.72 Å -->
  (properly separated, aligned in line)
```

## Integration with linear()

The `linear()` function applies the placer **before each chemical reaction**:

```python
def linear(..., placer=None):
    ...
    for step in range(len(sequence) - 1):
        # Copy next monomer
        right_monomer = library[right_label].copy()
        
        # Select ports
        left_port, right_port = connector.select_ports(...)
        
        # Position right monomer (if placer provided)
        if placer is not None:
            placer.place_monomer(
                current_polymer,
                right_monomer,
                left_port,
                right_port,
            )
        
        # Execute reaction
        product = connector.connect(...)
```

**Benefits**:
- More extended polymer structures
- Avoids atom overlap before reaction
- Better initial geometry for MMFF optimization
- More realistic linear polymer conformations

## Custom Placers

You can create custom placers by implementing the protocols:

### Custom Separator

```python
class CustomSeparator:
    def get_separation(self, left_monomer, right_monomer, left_port, right_port) -> float:
        # Custom logic
        return 5.0  # Fixed 5.0 Å separation
```

### Custom Orienter

```python
class CustomOrienter:
    def get_orientation(self, left_monomer, right_monomer, left_port, right_port, separation):
        # Custom logic: calculate translation and rotation
        translation = np.array([separation, 0, 0])  # Simple x-axis translation
        rotation = np.eye(3)  # No rotation
        return translation, rotation
```

### Use Custom Placer

```python
custom_placer = Placer(
    separator=CustomSeparator(),
    orienter=CustomOrienter(),
)

polymer = linear(..., placer=custom_placer)
```

## Testing

See `test_placer_integration.py` for a complete example:

```bash
python test_placer_integration.py
```

Expected output:
```
Bounding box change:
  Without placer: 4.36 × 3.45 × 3.06 Å
  With placer:    7.13 × 4.83 × 5.04 Å

Max extent:
  Without placer: 4.36 Å
  With placer:    7.13 Å
  Difference:     2.78 Å

✅ SUCCESS: Placer created more extended structure!
```

## Future Enhancements

Potential improvements:
1. **Helical/Spiral Orienter**: For twisted polymer structures
2. **Branching Placer**: For branched polymers
3. **Collision Detection**: Check for atom overlap
4. **Energy-based Placement**: Use force field to find optimal positions
5. **Solvent-aware Placement**: Account for solvent molecules
6. **Periodic Boundary Placer**: For bulk polymer systems

## API Reference

### `create_vdw_linear_placer(buffer=0.0) -> Placer`

Factory function to create a VdW-based linear placer.

**Parameters**:
- `buffer` (float): Additional spacing in Angstroms

**Returns**:
- `Placer` instance configured with `VdWSeparator` and `LinearOrienter`

### `Placer(separator, orienter)`

Main placer class.

**Methods**:
- `place_monomer(left_monomer, right_monomer, left_port, right_port) -> None`

### `VdWSeparator(buffer=0.0)`

Van der Waals radius-based separator.

**Methods**:
- `get_separation(left_monomer, right_monomer, left_port, right_port) -> float`

### `LinearOrienter()`

Linear alignment orienter.

**Methods**:
- `get_orientation(left_monomer, right_monomer, left_port, right_port, separation) -> tuple[np.ndarray, np.ndarray]`

## Notes

- Placer modifies monomer coordinates **in-place**
- VdW radii are taken from `Element` database
- Rotation uses **Rodrigues' formula** for numerical stability
- Port direction defaults to centroid-to-anchor if no bonds found
- Works with any monomer that has 3D coordinates (`pos` or `xyz` keys)
