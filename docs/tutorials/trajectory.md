# Trajectory Analysis

Trajectories represent sequences of molecular frames that evolve over time. This tutorial teaches you how to work with MolPy's Trajectory class to analyze molecular dynamics data, process time series, and extract meaningful information from simulation results.

## Understanding Trajectories

Trajectories are collections of molecular frames that represent the evolution of a system over time. Each frame contains the atomic positions, velocities, and other properties at a specific time point. Trajectories are essential for studying molecular dynamics, analyzing structural changes, and understanding time-dependent phenomena.

MolPy's Trajectory class provides a unified interface for working with trajectory data. It supports different frame providers, including lists of frames and generators, and offers efficient access to individual frames and frame subsets. The class also provides tools for trajectory manipulation, analysis, and splitting.

Trajectories enable you to study how molecular systems evolve, analyze structural dynamics, and extract time-dependent properties. They're fundamental to molecular dynamics analysis and structural biology studies.

## Creating and Working with Trajectories

### Basic Trajectory Creation

Trajectories can be created from various frame sources, including lists of frames and frame generators.

```python
import molpy as mp
import numpy as np

# Create a simple water molecule frame
def create_water_frame(time, offset=0.0):
    """Create a water molecule frame at a specific time."""
    water_atoms = mp.Block({
        'x': [0.0 + offset, 0.9572 + offset, -0.2400 + offset],
        'y': [0.0, 0.0, 0.0],
        'z': [0.0, 0.0, 0.0],
        'element': ['O', 'H', 'H'],
        'type': [1, 2, 2]
    })

    frame = mp.Frame({'atoms': water_atoms})
    frame.metadata = {'time': time, 'step': int(time * 1000)}
    return frame

# Create a list of frames for a trajectory
frames = []
for i in range(10):
    time = i * 0.001  # 1 ps intervals
    offset = i * 0.01  # Small movement over time
    frame = create_water_frame(time, offset)
    frames.append(frame)

# Create trajectory from frame list
trajectory = mp.Trajectory(frames)
print(f"Trajectory length: {len(trajectory)} frames")
print(f"Time range: {frames[0].metadata['time']:.3f} to {frames[-1].metadata['time']:.3f} ps")
```

This creates a simple trajectory where the water molecule moves slightly over time. The trajectory maintains the sequence of frames and provides access to their properties.

### Trajectory Access and Iteration

Trajectories support standard Python iteration and indexing operations, making them easy to work with.

```python
# Iterate through all frames
print("Trajectory analysis:")
for i, frame in enumerate(trajectory):
    time = frame.metadata['time']
    oxygen_pos = frame['atoms']['x'][0]
    print(f"Frame {i}: time={time:.3f} ps, O_x={oxygen_pos:.3f} Å")

# Access specific frames
first_frame = trajectory[0]
last_frame = trajectory[-1]
middle_frame = trajectory[4]

print(f"\nFirst frame time: {first_frame.metadata['time']:.3f} ps")
print(f"Last frame time: {last_frame.metadata['time']:.3f} ps")
print(f"Middle frame time: {middle_frame.metadata['time']:.3f} ps")

# Slice trajectory
early_trajectory = trajectory[0:5]
late_trajectory = trajectory[5:10]

print(f"Early trajectory: {len(early_trajectory)} frames")
print(f"Late trajectory: {len(late_trajectory)} frames")
```

This demonstrates the basic trajectory operations. You can iterate through frames, access specific frames by index, and create trajectory subsets using slicing.

## Trajectory Analysis

### Time Series Analysis

Trajectories enable analysis of how molecular properties change over time.

```python
def analyze_trajectory_properties(trajectory):
    """Analyze how molecular properties change over time."""
    times = []
    oxygen_positions = []
    molecular_sizes = []

    for frame in trajectory:
        # Extract time
        time = frame.metadata['time']
        times.append(time)

        # Extract oxygen position
        oxygen_x = frame['atoms']['x'][0]
        oxygen_positions.append(oxygen_x)

        # Calculate molecular size (distance between atoms)
        positions = np.column_stack([
            frame['atoms']['x'],
            frame['atoms']['y'],
            frame['atoms']['z']
        ])

        # Calculate O-H distances
        oh_distances = []
        for i in range(1, 3):  # Hydrogen atoms
            distance = np.linalg.norm(positions[0] - positions[i])
            oh_distances.append(distance)

        avg_oh_distance = np.mean(oh_distances)
        molecular_sizes.append(avg_oh_distance)

    return times, oxygen_positions, molecular_sizes

# Analyze the trajectory
times, ox_pos, mol_sizes = analyze_trajectory_properties(trajectory)

print("Trajectory analysis results:")
for i, (time, pos, size) in enumerate(zip(times, ox_pos, mol_sizes)):
    print(f"Frame {i}: t={time:.3f} ps, O_x={pos:.3f} Å, O-H={size:.3f} Å")
```

This analysis shows how the water molecule's position and geometry change over time. Such time series analysis is fundamental to understanding molecular dynamics.

### Statistical Analysis

Trajectories enable statistical analysis of molecular properties across multiple time points.

```python
def calculate_trajectory_statistics(trajectory):
    """Calculate statistical properties across the trajectory."""
    all_oxygen_positions = []
    all_oh_distances = []

    for frame in trajectory:
        # Oxygen position
        oxygen_x = frame['atoms']['x'][0]
        all_oxygen_positions.append(oxygen_x)

        # O-H distances
        positions = np.column_stack([
            frame['atoms']['x'],
            frame['atoms']['y'],
            frame['atoms']['z']
        ])

        for i in range(1, 3):
            distance = np.linalg.norm(positions[0] - positions[i])
            all_oh_distances.append(distance)

    # Calculate statistics
    ox_pos_array = np.array(all_oxygen_positions)
    oh_dist_array = np.array(all_oh_distances)

    stats = {
        'oxygen_x_mean': np.mean(ox_pos_array),
        'oxygen_x_std': np.std(ox_pos_array),
        'oxygen_x_range': np.ptp(ox_pos_array),
        'oh_distance_mean': np.mean(oh_dist_array),
        'oh_distance_std': np.std(oh_dist_array),
        'oh_distance_range': np.ptp(oh_dist_array)
    }

    return stats

# Calculate and display statistics
stats = calculate_trajectory_statistics(trajectory)
print("Trajectory statistics:")
for key, value in stats.items():
    print(f"  {key}: {value:.3f}")
```

Statistical analysis provides insights into the average behavior and variability of your molecular system over time.

## Trajectory Manipulation

### Frame Mapping

Trajectories support mapping operations that apply functions to each frame, enabling batch processing and analysis.

```python
def add_center_of_mass(frame):
    """Add center of mass to frame metadata."""
    positions = np.column_stack([
        frame['atoms']['x'],
        frame['atoms']['y'],
        frame['atoms']['z']
    ])

    # Calculate center of mass (simplified - equal masses)
    com = np.mean(positions, axis=0)
    frame.metadata['center_of_mass'] = com
    return frame

# Apply mapping to trajectory
trajectory_with_com = trajectory.map(add_center_of_mass)

# Check results
for i, frame in enumerate(trajectory_with_com):
    com = frame.metadata['center_of_mass']
    print(f"Frame {i}: COM = {com}")
```

Mapping operations are useful for adding derived properties to frames or transforming frame data in consistent ways.

### Trajectory Splitting

Trajectories can be split into smaller segments for analysis or processing. This is useful for parallel processing or focused analysis of specific time periods.

```python
# Split trajectory by frame interval
frame_splitter = mp.TrajectorySplitter(trajectory)

# Split every 3 frames
frame_splits = frame_splitter.split_frames(3)
print(f"Frame-based splits: {len(frame_splits)} segments")
for i, split in enumerate(frame_splits):
    print(f"  Split {i}: {len(split)} frames")

# Split by time interval (if time metadata is available)
time_splits = frame_splitter.split_time(0.003)  # 3 ps intervals
print(f"Time-based splits: {len(time_splits)} segments")
for i, split in enumerate(time_splits):
    if len(split) > 0:
        start_time = split[0].metadata['time']
        end_time = split[-1].metadata['time']
        print(f"  Split {i}: {start_time:.3f} to {end_time:.3f} ps")
```

Splitting enables focused analysis of specific time periods and can improve computational efficiency for large trajectories.

## Advanced Trajectory Operations

### Custom Splitting Strategies

MolPy provides flexible splitting strategies that can be customized for specific analysis needs.

```python
class CustomTimeStrategy(mp.SplitStrategy):
    """Custom strategy for splitting based on time windows."""

    def __init__(self, window_size, overlap=0.0):
        self.window_size = window_size
        self.overlap = overlap

    def get_split_indices(self, trajectory):
        """Get indices for time-window based splitting."""
        split_indices = []
        current_time = 0.0

        for i, frame in enumerate(trajectory):
            frame_time = frame.metadata['time']

            if frame_time >= current_time:
                split_indices.append(i)
                current_time = frame_time + self.window_size - self.overlap

        return split_indices

# Use custom splitting strategy
custom_strategy = CustomTimeStrategy(window_size=0.005, overlap=0.001)
custom_splits = frame_splitter.split(custom_strategy)

print(f"Custom time-window splits: {len(custom_splits)} segments")
for i, split in enumerate(custom_splits):
    if len(split) > 0:
        time_range = split[-1].metadata['time'] - split[0].metadata['time']
        print(f"  Split {i}: {len(split)} frames, {time_range:.3f} ps duration")
```

Custom strategies enable sophisticated trajectory analysis workflows tailored to specific research questions.

### Generator-Based Trajectories

For very large trajectories, generator-based approaches provide memory-efficient processing.

```python
def frame_generator():
    """Generate frames on demand."""
    for i in range(20):
        time = i * 0.001
        offset = i * 0.01
        frame = create_water_frame(time, offset)
        yield frame

# Create generator-based trajectory
gen_trajectory = mp.Trajectory(mp.FrameGenerator(frame_generator()))

# Process frames without loading all into memory
print("Processing generator-based trajectory:")
for i, frame in enumerate(gen_trajectory):
    if i < 5:  # Only show first few
        time = frame.metadata['time']
        oxygen_x = frame['atoms']['x'][0]
        print(f"  Frame {i}: t={time:.3f} ps, O_x={oxygen_x:.3f} Å")
    else:
        break

print(f"Processed trajectory with generator approach")
```

Generator-based trajectories are essential for handling very large datasets that don't fit in memory.

## Trajectory Best Practices

### Memory Management

For large trajectories, use appropriate strategies to manage memory usage.

```python
def process_large_trajectory(trajectory, chunk_size=100):
    """Process large trajectory in chunks to manage memory."""
    results = []

    for i in range(0, len(trajectory), chunk_size):
        chunk = trajectory[i:i+chunk_size]
        chunk_results = []

        for frame in chunk:
            # Process individual frame
            oxygen_x = frame['atoms']['x'][0]
            chunk_results.append(oxygen_x)

        results.extend(chunk_results)
        print(f"Processed chunk {i//chunk_size + 1}")

    return results

# Example usage for large trajectory
# large_results = process_large_trajectory(large_trajectory, chunk_size=50)
```

### Metadata Consistency

Ensure that all frames in your trajectory have consistent metadata structure.

```python
def validate_trajectory_metadata(trajectory):
    """Validate that all frames have consistent metadata."""
    required_keys = ['time', 'step']
    metadata_keys = set()

    for frame in trajectory:
        if hasattr(frame, 'metadata') and frame.metadata:
            metadata_keys.update(frame.metadata.keys())

    missing_keys = set(required_keys) - metadata_keys
    if missing_keys:
        print(f"Warning: Missing metadata keys: {missing_keys}")

    return len(missing_keys) == 0

# Validate trajectory
is_valid = validate_trajectory_metadata(trajectory)
print(f"Trajectory metadata valid: {is_valid}")
```

## Summary

This tutorial covered the fundamental concepts of MolPy's Trajectory class. You learned how to create and work with trajectories, analyze time-dependent properties, and manipulate trajectory data for various analysis workflows.

Trajectories are essential for studying molecular dynamics and time-dependent phenomena. The Trajectory class provides efficient access to frame sequences, supports various analysis operations, and enables sophisticated trajectory manipulation. Understanding trajectory handling is crucial for molecular dynamics analysis and structural biology studies.

### Next Steps

Continue your MolPy journey by understanding region-based operations, exploring selection and filtering, mastering system organization for complex assemblies, and learning about advanced analysis techniques.

Understanding trajectory analysis is essential for molecular dynamics and time-dependent studies in MolPy. The Trajectory class provides the foundation for all time-series molecular analysis workflows.
