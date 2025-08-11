# Region & Geometric Operations

Geometric regions define spatial boundaries for molecular analysis. This tutorial teaches you how to work with MolPy's Region classes to create spatial filters, analyze molecular distributions, and perform geometric operations on molecular data.

## Understanding Geometric Regions

Geometric regions define spatial volumes in three-dimensional space where molecular analysis can be performed. They can represent simple shapes like cubes and spheres, or complex combinations created through boolean operations. Regions are essential for spatial analysis, confinement studies, and defining analysis volumes.

MolPy's Region system provides a unified interface for geometric operations. Each region implements the `MaskPredicate` interface, making it compatible with the selection system. Regions can be combined using logical operators (`&`, `|`, `~`) to create complex spatial filters. This design enables sophisticated spatial analysis workflows.

Regions are fundamental to understanding spatial organization in molecular systems. They help analyze molecular distributions, study confinement effects, and define analysis boundaries for complex molecular assemblies.

## Basic Region Types

### Box Regions

Box regions define rectangular volumes with specified dimensions and origin. They're useful for defining analysis boundaries and studying confinement effects.

```python
import molpy as mp
import numpy as np

# Create a box region
box_region = mp.BoxRegion(
    lengths=[10.0, 8.0, 6.0],  # Length in x, y, z directions
    origin=[0.0, 0.0, 0.0],    # Origin point
    coord_field="xyz"           # Coordinate field name
)

# Create a cube region (special case of box)
cube_region = mp.Cube(
    edge=5.0,                   # Edge length
    origin=[2.0, 2.0, 2.0],    # Center of cube
    coord_field="xyz"
)

# Check region properties
print(f"Box region: {box_region}")
print(f"Box bounds: {box_region.bounds}")
print(f"Cube region: {cube_region}")
print(f"Cube bounds: {cube_region.bounds}")
```

Box regions are defined by their dimensions and origin point. The `bounds` property provides the axis-aligned bounding box, which is useful for understanding the spatial extent of the region.

### Sphere Regions

Sphere regions define spherical volumes with specified radius and center. They're useful for studying molecular environments around specific points.

```python
# Create a sphere region
sphere_region = mp.SphereRegion(
    radius=3.0,                 # Sphere radius
    center=[5.0, 5.0, 5.0],    # Sphere center
    coord_field="xyz"
)

# Check sphere properties
print(f"Sphere region: {sphere_region}")
print(f"Sphere bounds: {sphere_region.bounds}")
print(f"Sphere center: {sphere_region.center}")
print(f"Sphere radius: {sphere_region.radius}")
```

Sphere regions are defined by their radius and center point. They're particularly useful for studying molecular solvation shells, binding sites, and local environments.

## Working with Regions

### Coordinate Field Specification

Regions work with different coordinate field names, making them flexible for various data formats.

```python
# Create molecular data with different coordinate field names
molecular_data = {
    'x': [0.0, 1.0, 2.0, 3.0, 4.0],
    'y': [0.0, 1.0, 2.0, 3.0, 4.0],
    'z': [0.0, 1.0, 2.0, 3.0, 4.0],
    'element': ['C', 'H', 'O', 'N', 'S']
}

# Create block with xyz coordinates
molecule = mp.Block(molecular_data)

# Create regions with different coordinate field specifications
region_xyz = mp.Cube(3.0, origin=[1.0, 1.0, 1.0], coord_field="xyz")
region_individual = mp.Cube(3.0, origin=[1.0, 1.0, 1.0], coord_field=["x", "y", "z"])

# Test both regions
print(f"Region with xyz field: {region_xyz}")
print(f"Region with individual fields: {region_individual}")
```

The coordinate field specification determines how regions extract spatial information from molecular data. This flexibility enables regions to work with various data formats.

### Region Testing

Regions can test whether points lie inside their boundaries, providing boolean masks for spatial filtering.

```python
# Test points against regions
test_points = np.array([
    [0.0, 0.0, 0.0],  # Inside box, outside sphere
    [5.0, 5.0, 5.0],  # Inside sphere, outside box
    [1.5, 1.5, 1.5],  # Inside both
    [10.0, 10.0, 10.0] # Outside both
])

# Test against different regions
box_mask = box_region.isin(test_points)
sphere_mask = sphere_region.isin(test_points)
cube_mask = cube_region.isin(test_points)

print("Point testing results:")
for i, point in enumerate(test_points):
    print(f"Point {i} {point}:")
    print(f"  In box: {box_mask[i]}")
    print(f"  In sphere: {sphere_mask[i]}")
    print(f"  In cube: {cube_mask[i]}")
```

This testing capability enables spatial filtering of molecular data, allowing you to analyze specific regions of your system.

## Region Combinations

### Boolean Operations

Regions can be combined using logical operators to create complex spatial filters.

```python
# Create complex regions through combination
# Region inside box but outside sphere
box_not_sphere = box_region & ~sphere_region

# Region inside either cube or sphere
cube_or_sphere = cube_region | sphere_region

# Region inside all three
all_regions = box_region & cube_region & sphere_region

print("Combined regions:")
print(f"Box not sphere: {box_not_sphere}")
print(f"Cube or sphere: {cube_or_sphere}")
print(f"All regions: {all_regions}")

# Test combined regions
box_not_sphere_mask = box_not_sphere.isin(test_points)
cube_or_sphere_mask = cube_or_sphere.isin(test_points)
all_regions_mask = all_regions.isin(test_points)

print("\nCombined region testing:")
for i, point in enumerate(test_points):
    print(f"Point {i} {point}:")
    print(f"  Box not sphere: {box_not_sphere_mask[i]}")
    print(f"  Cube or sphere: {cube_or_sphere_mask[i]}")
    print(f"  All regions: {all_regions_mask[i]}")
```

Boolean combinations enable sophisticated spatial analysis. You can create regions that represent complex spatial relationships and use them for targeted molecular analysis.

### Complex Spatial Filters

Complex regions can be built by combining multiple simple regions in various ways.

```python
def create_complex_region():
    """Create a complex spatial region for analysis."""
    # Define basic regions
    central_cube = mp.Cube(4.0, origin=[3.0, 3.0, 3.0])
    outer_shell = mp.SphereRegion(6.0, center=[5.0, 5.0, 5.0])
    exclusion_box = mp.BoxRegion([2.0, 2.0, 2.0], origin=[4.0, 4.0, 4.0])

    # Create complex region: outer shell minus central cube and exclusion box
    complex_region = outer_shell & ~central_cube & ~exclusion_box

    return complex_region

# Create and test complex region
complex_region = create_complex_region()
print(f"Complex region: {complex_region}")
print(f"Complex region bounds: {complex_region.bounds}")

# Test against test points
complex_mask = complex_region.isin(test_points)
for i, point in enumerate(test_points):
    print(f"Point {i} {point}: In complex region: {complex_mask[i]}")
```

This demonstrates how complex spatial filters can be built from simple components. Such regions are useful for studying specific molecular environments or defining complex analysis volumes.

## Region Applications

### Molecular Confinement Analysis

Regions can be used to study how molecules behave under spatial confinement.

```python
def analyze_confinement(molecule, region):
    """Analyze molecular properties within a confined region."""
    # Get atoms inside the region
    region_mask = region(molecule)
    confined_atoms = molecule[region_mask]

    if confined_atoms.nrows == 0:
        return None

    # Analyze confined atoms
    analysis = {
        'n_atoms': confined_atoms.nrows,
        'elements': list(set(confined_atoms['element'])),
        'center': np.mean(np.column_stack([
            confined_atoms['x'],
            confined_atoms['y'],
            confined_atoms['z']
        ]), axis=0)
    }

    return analysis

# Analyze confinement in different regions
regions_to_test = [box_region, sphere_region, cube_region]

print("Confinement analysis:")
for i, region in enumerate(regions_to_test):
    result = analyze_confinement(molecule, region)
    if result:
        print(f"Region {i}: {result}")
    else:
        print(f"Region {i}: No atoms inside")
```

This analysis shows how regions can be used to study molecular confinement effects, which are important for understanding molecular behavior in restricted environments.

### Spatial Distribution Analysis

Regions enable analysis of molecular spatial distributions and clustering.

```python
def analyze_spatial_distribution(molecule, regions):
    """Analyze molecular distribution across multiple regions."""
    distribution = {}

    for i, region in enumerate(regions):
        region_mask = region(molecule)
        atoms_in_region = molecule[region_mask]

        distribution[f'region_{i}'] = {
            'n_atoms': atoms_in_region.nrows,
            'density': atoms_in_region.nrows / region.volume if hasattr(region, 'volume') else 'N/A'
        }

    return distribution

# Analyze distribution across regions
regions = [box_region, sphere_region, cube_region]
distribution = analyze_spatial_distribution(molecule, regions)

print("Spatial distribution analysis:")
for region_name, data in distribution.items():
    print(f"{region_name}: {data}")
```

This analysis provides insights into how molecules are distributed across different spatial regions, which is crucial for understanding molecular organization and interactions.

## Region Best Practices

### Coordinate System Consistency

Ensure that regions use the same coordinate system as your molecular data.

```python
def validate_coordinate_system(region, molecule):
    """Validate that region coordinates match molecule coordinates."""
    # Check if coordinate fields exist
    if hasattr(region, 'coord_field'):
        if isinstance(region.coord_field, str):
            if region.coord_field == "xyz":
                # Check if xyz field exists or if individual fields exist
                if "xyz" not in molecule and not all(f in molecule for f in ["x", "y", "z"]):
                    print(f"Warning: Region uses 'xyz' field but molecule doesn't have it")
            elif region.coord_field in molecule:
                print(f"Region uses field: {region.coord_field}")
            else:
                print(f"Warning: Region field {region.coord_field} not found in molecule")
        elif isinstance(region.coord_field, list):
            if all(f in molecule for f in region.coord_field):
                print(f"Region uses fields: {region.coord_field}")
            else:
                print(f"Warning: Some region fields not found in molecule")

# Validate coordinate systems
validate_coordinate_system(box_region, molecule)
```

### Region Bounds

Always check region bounds to ensure they're appropriate for your analysis.

```python
def check_region_bounds(region, molecule):
    """Check if region bounds are appropriate for the molecule."""
    # Get molecule bounds
    molecule_coords = np.column_stack([
        molecule['x'], molecule['y'], molecule['z']
    ])
    molecule_min = np.min(molecule_coords, axis=0)
    molecule_max = np.max(molecule_coords, axis=0)

    # Get region bounds
    region_bounds = region.bounds
    region_min = region_bounds[0]
    region_max = region_bounds[1]

    # Check overlap
    overlap = np.all((region_max >= molecule_min) & (region_min <= molecule_max))

    if overlap:
        print("✓ Region overlaps with molecule")
    else:
        print("✗ Region does not overlap with molecule")

    return overlap

# Check region bounds
check_region_bounds(box_region, molecule)
check_region_bounds(sphere_region, molecule)
```

## Summary

This tutorial covered the fundamental concepts of MolPy's Region system. You learned how to create basic geometric regions, combine them using boolean operations, and apply them to molecular analysis workflows.

Geometric regions provide powerful tools for spatial analysis in molecular systems. They enable confinement studies, spatial distribution analysis, and complex spatial filtering. The Region system integrates seamlessly with MolPy's selection framework, making it easy to create sophisticated spatial analysis workflows.

### Next Steps

Continue your MolPy journey by exploring selection and filtering, mastering system organization for complex assemblies, understanding trajectory analysis, and learning about advanced analysis techniques.

Understanding geometric regions is essential for spatial molecular analysis in MolPy. The Region system provides the foundation for sophisticated spatial filtering and analysis workflows.
