# C3DExport

**A Python library for converting C3D biomechanics motion capture files to OpenSim-compatible formats.**

C3DExport converts C3D files (commonly used in biomechanics motion capture) into OpenSim-compatible `.trc` (marker trajectory) and `.mot` (ground reaction force) files. It supports Kistler force plates with both 8-channel and 6-channel configurations.

---

## Table of Contents

- [Features](#features)
- [Project Structure](#project-structure)
- [Installation](#installation)
- [Quick Start](#quick-start)
- [Module: `C3DExport.core`](#module-c3dexportcore)
- [Module: `C3DExport.io`](#module-c3dexportio)
- [Module: `C3DExport.utils`](#module-c3dexportutils)
- [Output File Formats](#output-file-formats)
- [Dependencies](#dependencies)
- [License](#license)

---

## Features

- **Three conversion pipelines** for different force plate configurations
- **Automatic stance phase detection** (threshold-based and peak-based)
- **COP (Center of Pressure) anomaly detection and correction**
- **Coordinate system transforms**: Kistler plate-local -> Lab global -> OpenSim Y-up
- **Butterworth low-pass filtering** for marker and force data
- **Polyphase resampling** to synchronize analog and marker data rates
- **Batch processing** support via helper scripts in `Transform/`

---

## Project Structure

```
c3dexport/
├── C3DExport/                 # Core library package
│   ├── __init__.py            # Package entry point, exports main converters
│   ├── core.py                # Main conversion pipelines
│   ├── io.py                  # OpenSim file I/O (.trc and .mot writers)
│   └── utils.py               # Utility functions (rotation, filtering, detection)
│
├── Transform/                 # Helper scripts and batch processing
│   ├── c3d_to_opensim.py      # Standalone 8-channel converter
│   ├── c3d_to_opensim_v2.py   # Standalone 8-channel converter (v2)
│   ├── c3d6_to_opensim.py     # Standalone 6-channel converter
│   ├── batch_process.py       # Batch processing script
│   ├── batch_process_v2.py    # Batch processing script (v2)
│   ├── batch_c3d6_to_opensim.py # Batch 6-channel processing
│   ├── draw_picture_check_grf.py # GRF visualization
│   ├── draw_check_v2.py       # Visualization (v2)
│   ├── transform_utils.py     # Additional transform utilities
│   └── README.md              # Transform module documentation
│
└── README.md                  # This file
```

---

## Installation

### Prerequisites

- Python 3.8+
- Required packages:

```bash
pip install numpy scipy ezc3d pandas
```

### Install

Clone or copy the `C3DExport` folder into your project directory, then import:

```python
from C3DExport import convert_c3d, convert_c3d_v2, convert_c3d6
```

---

## Quick Start

```python
from C3DExport import convert_c3d

# Convert an 8-channel Kistler C3D file
result = convert_c3d("data/trial.c3d", "output/")

print(result['trc_path'])      # -> "output/trial.trc"
print(result['mot_path'])      # -> "output/trial.mot"
print(result['csv_path'])      # -> "output/cut_records.csv"
print(result['stance_info'])   # -> stance phase metadata dict
```

### Using V2 (Peak Stance + COP Correction)

```python
from C3DExport import convert_c3d_v2

result = convert_c3d_v2("data/trial.c3d", "output/",
                         cop_middle_ratio=0.3,
                         cop_rate_multiplier=2.0,
                         cop_jump_threshold=0.03)
```

### Using 6-Channel Converter

```python
from C3DExport import convert_c3d6

result = convert_c3d6("data/trial.c3d", "output/")
```

---

## Module: `C3DExport.core`

The core module provides the main conversion pipelines. All three converters follow the same general flow:

1. **Read C3D** - Parse marker and force plate data from the C3D file
2. **Resample** - Synchronize analog (force) data to marker frame rate
3. **Force Computation** - Convert raw force plate channels to Fx/Fy/Fz/COP/Tz
4. **Coordinate Transforms** - Plate-local -> Lab -> OpenSim Y-up
5. **Filtering** - Apply Butterworth low-pass filter to markers
6. **Stance Detection** - Identify heel-strike and toe-off events
7. **Output** - Write `.trc`, `.mot`, and CSV records

### Main Converter Functions

#### `convert_c3d(c3d_path, output_dir, ...)`

Convert **8-channel Kistler** C3D to OpenSim files using **threshold-based** stance detection.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `c3d_path` | `str` | required | Path to the input `.c3d` file |
| `output_dir` | `str` | required | Directory for output files |
| `marker_cutoff` | `float` | `6.0` | Marker low-pass filter cutoff frequency (Hz) |
| `force_cutoff` | `float` | `50.0` | Force low-pass filter cutoff frequency (Hz) |
| `stance_threshold` | `float` | `30.0` | Vertical force threshold for stance detection (N) |
| `stance_pad_frames` | `int` | `25` | Padding frames before/after stance phase |
| `verbose` | `bool` | `True` | Print progress messages to stdout |

**Returns:** `dict` with keys:
- `trc_path` (`str`): Path to the output `.trc` file
- `mot_path` (`str`): Path to the output `.mot` file
- `csv_path` (`str`): Path to the cut records CSV file
- `stance_info` (`dict`): Stance phase metadata (heel-strike frame, toe-off frame, cut range, etc.)

---

#### `convert_c3d_v2(c3d_path, output_dir, ...)`

Convert **8-channel Kistler** C3D to OpenSim files using **peak-based** stance detection with **COP slope correction**.

Additional parameters beyond `convert_c3d`:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `cop_middle_ratio` | `float` | `0.3` | Middle section ratio for COP slope fitting (30%-70%) |
| `cop_rate_multiplier` | `float` | `2.0` | Outlier detection rate multiplier |
| `cop_jump_threshold` | `float` | `0.03` | Frame-to-frame COP jump threshold (m) |

**Returns:** `dict` with same keys as `convert_c3d`, plus additional COP correction info in `stance_info`:
- `peak_frame` (`int`): Frame index of peak vertical force
- `peak_value` (`float`): Peak vertical force value (N)
- `copx_slope` (`float`): Fitted COPx slope
- `copz_slope` (`float`): Fitted COPz slope
- `copx_outlier_count` (`int`): Number of COPx outlier frames corrected
- `copz_outlier_count` (`int`): Number of COPz outlier frames corrected

---

#### `convert_c3d6(c3d_path, output_dir, ...)`

Convert **6-channel** force plate C3D to OpenSim files.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `c3d_path` | `str` | required | Path to the input `.c3d` file with 6-channel data |
| `output_dir` | `str` | required | Directory for output files |
| `marker_cutoff` | `float` | `6.0` | Marker low-pass filter cutoff frequency (Hz) |
| `force_cutoff` | `float` | `50.0` | Force low-pass filter cutoff frequency (Hz) |
| `stance_threshold` | `float` | `30.0` | Vertical force threshold for stance detection (N) |
| `stance_pad_frames` | `int` | `25` | Padding frames before/after stance phase |
| `verbose` | `bool` | `True` | Print progress messages to stdout |

**Returns:** `dict` with same keys as `convert_c3d`.

### Internal Helper Functions

These are private functions used internally by the converters:

| Function | Description |
|----------|-------------|
| `_read_c3d(c3d_path)` | Parse C3D file, extract marker/force data, plate geometry |
| `_resample_analog(...)` | Resample analog signals to match marker frame rate |
| `_compute_type2_8ch(...)` | Compute Type 2 force data from 8-channel Kistler raw data |
| `_compute_type2_6ch(...)` | Compute Type 2 force data from 6-channel raw data |
| `_transform_forces(...)` | Transform forces: plate-local -> lab -> OpenSim |
| `_rotate_markers(...)` | Rotate marker coordinates to OpenSim system |
| `_filter_markers(...)` | Apply Butterworth low-pass filter to markers |
| `_build_structured_plates_v1(...)` | Build structured plate data with COP threshold filtering |
| `_detect_active_plate(...)` | Detect which force plate has the strongest vertical force |
| `_cut_and_zero_padding(...)` | Cut data to stance range and zero-pad outside stance |
| `_markers_to_flat(...)` | Reshape markers from (3, N, T) to (T, N*3) for .trc output |
| `_build_mot_plates(...)` | Convert structured plates to format expected by `write_mot` |
| `_write_csv_record(...)` | Append stance info to CSV record file |

---

## Module: `C3DExport.io`

File I/O functions for writing OpenSim-compatible output files.

### Functions

#### `write_trc(filepath, marker_data, marker_labels, frame_rate, units='mm', orig_start_frame=1)`

Write marker trajectory data to an OpenSim `.trc` file.

| Parameter | Type | Description |
|-----------|------|-------------|
| `filepath` | `str` | Output file path |
| `marker_data` | `ndarray (n_frames, n_markers*3)` | Marker XYZ coordinates, interleaved |
| `marker_labels` | `list[str]` | Marker names (e.g., `['R.ASIS', 'L.ASIS', ...]`) |
| `frame_rate` | `float` | Sampling rate in Hz |
| `units` | `str` | `'mm'` or `'m'` (default: `'mm'`) |
| `orig_start_frame` | `int` | Original first frame number (default: `1`) |

The `.trc` file format includes:
- Header with file metadata (data rate, frame count, marker count, units)
- Column headers with marker names and axis labels (`X1, Y1, Z1, X2, Y2, Z2, ...`)
- Data rows with frame number, timestamp, and XYZ coordinates for each marker

---

#### `write_mot(filepath, force_data_per_plate, n_plates, frame_rate, filename=None)`

Write ground reaction force data to an OpenSim `.mot` file.

| Parameter | Type | Description |
|-----------|------|-------------|
| `filepath` | `str` | Output file path |
| `force_data_per_plate` | `list[dict]` | Per-plate data, each dict has keys: `force` (N,3), `cop` (N,3), `torque` (N,3) |
| `n_plates` | `int` | Number of force plates |
| `frame_rate` | `float` | Sampling rate in Hz |
| `filename` | `str` or `None` | Header filename (default: basename of filepath) |

The `.mot` file column layout:
```
time
ground_force_vx  ground_force_vy  ground_force_vz    (FP1)
ground_force_px  ground_force_py  ground_force_pz    (FP1 COP)
1_ground_force_vx ...                                (FP2, if present)
ground_torque_x  ground_torque_y  ground_torque_z    (FP1)
1_ground_torque_x ...                                (FP2, if present)
```

---

## Module: `C3DExport.utils`

Core utility functions organized into five categories.

### Rotation Matrices

#### `rotation_matrix(axis, angle_deg)`

Generate a 3x3 rotation matrix following the right-hand rule.

| Parameter | Type | Description |
|-----------|------|-------------|
| `axis` | `str` | Rotation axis: `'X'`, `'Y'`, or `'Z'` |
| `angle_deg` | `float` | Rotation angle in degrees |

**Returns:** `ndarray (3, 3)` - Rotation matrix

---

#### `chain_rotations(*rotations)`

Chain multiple rotations together (applied left-to-right).

| Parameter | Type | Description |
|-----------|------|-------------|
| `*rotations` | `tuples` | Sequence of `(axis, angle_deg)` tuples |

**Returns:** `ndarray (3, 3)` - Combined rotation matrix

**Example:**
```python
R = chain_rotations(('X', -90), ('Y', 90))
```

---

#### `apply_rotation(data_3xN, R)`

Apply a 3x3 rotation matrix to (3, N) data.

| Parameter | Type | Description |
|-----------|------|-------------|
| `data_3xN` | `ndarray (3, N)` | Input data (3 rows x N columns) |
| `R` | `ndarray (3, 3)` | Rotation matrix |

**Returns:** `ndarray (3, N)` - Rotated data

---

### Kistler Force-Plate Computation

#### `compute_kistler_channel8(channels_8, a, b, az0)`

Compute Type 2 force-plate data from **8-channel Kistler Type 3** raw data.

| Parameter | Type | Description |
|-----------|------|-------------|
| `channels_8` | `ndarray (8, N)` | Raw channels: `[fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]` |
| `a` | `float` | Sensor offset in local-Y (AP direction), mm |
| `b` | `float` | Sensor offset in local-X (ML direction), mm |
| `az0` | `float` | Top-plate offset, mm (typically negative) |

**Returns:** `dict` with keys:
- `Fx`, `Fy`, `Fz` (`ndarray`): Ground reaction forces (N)
- `ax`, `ay` (`ndarray`): COP coordinates from plate center (mm)
- `Tz` (`ndarray`): Free vertical moment (N*mm)

**Algorithm:**
1. Sum raw channels to get Fx_raw, Fy_raw, Fz_raw
2. Compute moments Mx, My, Mz from sensor geometry
3. Transfer moments to plate surface using `az0` offset
4. Compute COP: `ax = -Myp/Fz`, `ay = Mxp/Fz` (valid only when |Fz| >= 20N)
5. Compute free vertical moment: `Tz = Mz - Fy*ax + Fx*ay`
6. Negate to get ground reaction force convention

---

#### `compute_kistler_channel6(channels_6, az0)`

Compute COP and free moment from **6-channel** force-plate data.

| Parameter | Type | Description |
|-----------|------|-------------|
| `channels_6` | `ndarray (6, N)` | Raw channels: `[Fx, Fy, Fz, Mx, My, Mz]` |
| `az0` | `float` | Top-plate offset, mm |

**Returns:** `dict` with same keys as `compute_kistler_channel8`

---

### Coordinate Transforms

#### `plate_local_to_lab(type2, corners)`

Convert force data from Kistler plate-local to lab global coordinates.

| Coordinate System | X | Y | Z |
|-------------------|---|---|---|
| Kistler local | right+ | posterior+ | down+ |
| Lab global | forward+ | left+ | up+ |

| Parameter | Type | Description |
|-----------|------|-------------|
| `type2` | `dict` | Output of `compute_kistler_channel8` or `compute_kistler_channel6` |
| `corners` | `ndarray (3, 4)` | Force plate corners from C3D `FORCE_PLATFORM:CORNERS` |

**Returns:** `dict` with keys: `Fx`, `Fy`, `Fz` (N), `COPx`, `COPy`, `COPz` (mm), `Tz` (N*mm)

---

#### `lab_to_opensim_force(data_lab)`

Convert force data from lab global to OpenSim coordinate system.

| Coordinate System | X | Y | Z |
|-------------------|---|---|---|
| Lab global | forward | left | up |
| OpenSim | forward | up | right |

| Parameter | Type | Description |
|-----------|------|-------------|
| `data_lab` | `dict` | Force data in lab coordinates |

**Returns:** `dict` with same keys, values in OpenSim coordinates

---

### Filtering & Resampling

#### `butter_lowpass_filter(data, cutoff, fs, order=4)`

4th-order Butterworth low-pass filter with zero-phase implementation (via `filtfilt`).

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `data` | `ndarray (N,)` | required | Input signal |
| `cutoff` | `float` | required | Cutoff frequency (Hz) |
| `fs` | `float` | required | Sampling frequency (Hz) |
| `order` | `int` | `4` | Filter order |

**Returns:** `ndarray (N,)` - Filtered signal

---

#### `resample_to_target_rate(data_1d, src_rate, tgt_rate)`

Resample a 1D signal from source rate to target rate using polyphase resampling.

| Parameter | Type | Description |
|-----------|------|-------------|
| `data_1d` | `ndarray (N,)` | Input signal |
| `src_rate` | `float` | Source sampling rate (Hz) |
| `tgt_rate` | `float` | Target sampling rate (Hz) |

**Returns:** `ndarray (M,)` - Resampled signal

---

### Stance Phase Detection

#### `detect_stance_phase(vertical_force, threshold=30.0, pad_frames=25)`

Detect stance phase using the **threshold method**.

**Algorithm:** Find the first and last frames where |vertical_force| > threshold, then add padding frames.

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vertical_force` | `ndarray (N,)` | required | Vertical GRF (positive = up in OpenSim Y) |
| `threshold` | `float` | `30.0` | Force threshold (N) |
| `pad_frames` | `int` | `25` | Extra frames before/after stance |

**Returns:** `(start_idx, end_idx, hs_idx, toe_idx)` - All 0-indexed frame indices

**Raises:** `ValueError` if no stance phase is detected

---

#### `detect_stance_phase_from_peak(vertical_force, threshold=30.0, pad_frames=25)`

Detect stance phase using the **peak-based method**.

**Algorithm:**
1. Find the peak (maximum absolute vertical force)
2. Walk left from peak until force < threshold -> heel-strike
3. Walk right from peak until force < threshold -> toe-off
4. Add padding frames

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `vertical_force` | `ndarray (N,)` | required | Vertical GRF (OpenSim Y, positive = up) |
| `threshold` | `float` | `30.0` | Force threshold (N) |
| `pad_frames` | `int` | `25` | Padding frames |

**Returns:** `(start_idx, end_idx, hs_idx, toe_idx, peak_idx, peak_val)` - All 0-indexed

**Raises:** `ValueError` if peak value <= threshold

---

### COP Anomaly Detection & Correction

#### `detect_cop_anomalies(copx, copy, time, force_vertical, threshold=20.0, jump_threshold=0.03)`

Detect COPx/COPy anomaly frames during stance phase.

**Dual detection strategy:**
1. Force < threshold but COP != 0 -> anomaly
2. Frame-to-frame COP jump > jump_threshold (m) -> anomaly

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `copx`, `copy` | `ndarray (N,)` | required | COP coordinates in meters |
| `time` | `ndarray (N,)` | required | Time vector in seconds |
| `force_vertical` | `ndarray (N,)` | required | Vertical force in N |
| `threshold` | `float` | `20.0` | Force threshold (N) |
| `jump_threshold` | `float` | `0.03` | COP jump threshold (m) |

**Returns:**
- `copx_anomaly` (`ndarray (N,) bool`): Anomaly flags for COPx
- `copy_anomaly` (`ndarray (N,) bool`): Anomaly flags for COPy
- `info` (`dict`): Anomaly statistics

---

#### `correct_cop_slope(copx, copy, time, force_vertical, threshold=20.0, middle_ratio=0.3, rate_multiplier=2.0)`

Correct COP outliers using middle-section slope fitting.

**Algorithm:**
1. Identify stance phase (force >= threshold)
2. Fit linear slope to the middle section (e.g., 30%-70% of stance)
3. Walk outward from middle; if frame-to-frame delta > rate_multiplier * |k*dt|, replace with slope-predicted value
4. Non-stance COP values set to zero

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `copx`, `copy` | `ndarray (N,)` | required | COP coordinates in meters |
| `time` | `ndarray (N,)` | required | Time vector in seconds |
| `force_vertical` | `ndarray (N,)` | required | Vertical force in N |
| `threshold` | `float` | `20.0` | Force threshold (N) |
| `middle_ratio` | `float` | `0.3` | Middle section ratio for fitting |
| `rate_multiplier` | `float` | `2.0` | Outlier rate multiplier |

**Returns:**
- `copx_corrected` (`ndarray`): Corrected COPx
- `copy_corrected` (`ndarray`): Corrected COPy
- `info` (`dict`): Slope info and outlier counts

---

## Output File Formats

### `.trc` File (Marker Trajectory)

Tab-separated file with:
- **Line 1:** File type identifier (`PathFileType 4 (X/Y/Z)`)
- **Line 2:** Header field names
- **Line 3:** Header values (data rate, frame count, marker count, units, etc.)
- **Line 4:** Marker names (each spanning 3 columns)
- **Line 5:** Axis sub-headers (`X1, Y1, Z1, X2, Y2, Z2, ...`)
- **Line 6:** Blank line
- **Data rows:** `Frame#  Time  X1  Y1  Z1  X2  Y2  Z2  ...`

### `.mot` File (Ground Reaction Force)

Tab-separated file with:
- **Header:** filename, version, nRows, nColumns, inDegrees, endheader
- **Column labels:** `time`, `ground_force_vx/vy/vz`, `ground_force_px/py/pz`, `ground_torque_x/y/z`
- **Data rows:** Time-series of forces (N), COP (m), and torques (N*m)

### `cut_records.csv` (Stance Metadata)

CSV file recording stance phase information for each processed trial:
- `trial`: Trial name
- `heel_strike_frame`: Heel-strike frame number
- `toe_off_frame`: Toe-off frame number
- `cut_start_frame`, `cut_end_frame`: Cut range
- `total_frames`: Number of frames in cut region

---

## Dependencies

| Package | Version | Purpose |
|---------|---------|---------|
| `numpy` | >= 1.20 | Array operations, linear algebra |
| `scipy` | >= 1.7 | Butterworth filter, polyphase resampling |
| `ezc3d` | >= 1.5 | C3D file parsing |
| `pandas` | >= 1.3 | CSV record writing |

Install all dependencies:

```bash
pip install numpy scipy ezc3d pandas
```

---

## License

This project is for research and educational purposes.
