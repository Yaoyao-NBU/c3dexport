"""Quick diagnostic: trace data shapes through the pipeline."""
import os, sys
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from C3DExport.core import _read_c3d, _resample_analog, _compute_forces, _transform_forces
from C3DExport.utils import plate_local_to_lab, lab_to_opensim_force

c3d_path = os.path.join(os.path.dirname(__file__), 'Sample10', 'TYPE-2.C3d')

if not os.path.exists(c3d_path):
    # try uppercase
    c3d_path = os.path.join(os.path.dirname(__file__), 'Sample10', 'TYPE-2.C3D')

print(f"Reading: {c3d_path}")
data = _read_c3d(c3d_path)

print(f"\n--- Raw data ---")
print(f"  analog_data shape : {data['analog_data'].shape}")
print(f"  channels shape    : {data['channels'].shape}")
print(f"  origin shape      : {data['origin'].shape}")
print(f"  corners shape     : {data['corners'].shape}")
print(f"  point_rate        : {data['point_rate']}")
print(f"  analog_rate       : {data['analog_rate']}")
print(f"  n_point_frames    : {data['n_point_frames']}")
print(f"  fp_types          : {data['fp_types']}")

print(f"\n--- After resampling ---")
analog_data, analog_rate = _resample_analog(
    data['analog_data'], data['analog_rate'],
    data['point_rate'], data['n_point_frames']
)
print(f"  resampled shape   : {analog_data.shape}")
print(f"  returned rate     : {analog_rate}")

print(f"\n--- After force computation ---")
plate_type2 = _compute_forces(
    analog_data, data['channels'], data['origin'],
    data['n_plates'], data['fp_types']
)
for pi, p in enumerate(plate_type2):
    for k, v in p.items():
        print(f"  plate[{pi}]['{k}'] shape: {v.shape}")

print(f"\n--- After coordinate transforms ---")
plate_os, R_lab2os = _transform_forces(
    plate_type2, data['corners'], data['n_plates']
)
for pi, p in enumerate(plate_os):
    for k, v in p.items():
        print(f"  plate_os[{pi}]['{k}'] shape: {v.shape}")

print(f"\n--- Restructured (force/cop/torque) ---")
for pi in range(data['n_plates']):
    d = plate_os[pi]
    n = len(d['Fx'])
    force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])
    cop    = np.column_stack([d['COPx'] / 1000.0, d['COPy'] / 1000.0, d['COPz'] / 1000.0])
    torque = np.column_stack([np.zeros(n), d['Tz'] / 1000.0, np.zeros(n)])
    print(f"  plate[{pi}] force shape: {force.shape}")
    print(f"  plate[{pi}] cop   shape: {cop.shape}")
    print(f"  plate[{pi}] torque shape: {torque.shape}")

print("\nAll shapes OK!")
