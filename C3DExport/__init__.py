"""
C3DExport -- C3D to OpenSim Conversion Library
===============================================

Convert C3D biomechanics motion capture files to OpenSim-compatible
.trc (marker) and .mot (ground reaction force) files.

Quick Start
-----------
    from C3DExport import convert_c3d

    result = convert_c3d("data/trial.c3d", "output/")
    print(result['trc_path'])
    print(result['mot_path'])

Available Converters
--------------------
    convert_c3d     -- 8-channel Kistler (threshold stance detection)
    convert_c3d_v2  -- 8-channel Kistler (peak stance + COP correction)
    convert_c3d6    -- 6-channel force data

Utility Functions (from C3DExport.utils)
----------------------------------------
    rotation_matrix, apply_rotation, chain_rotations
    compute_kistler_channel8, compute_kistler_channel6
    plate_local_to_lab, lab_to_opensim_force
    butter_lowpass_filter, resample_to_target_rate
    detect_stance_phase, detect_stance_phase_from_peak
    detect_cop_anomalies, correct_cop_slope

I/O Functions (from C3DExport.io)
---------------------------------
    write_trc, write_mot
"""

from .core import convert_c3d, convert_c3d_v2, convert_c3d6

__all__ = [
    'convert_c3d',
    'convert_c3d_v2',
    'convert_c3d6',
]

__version__ = '1.0.0'
