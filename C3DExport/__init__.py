"""
C3DExport -- C3D to OpenSim Conversion Library
===============================================

Convert C3D biomechanics motion capture files to OpenSim-compatible
.trc (marker) and .mot (ground reaction force) files.

Quick Start
-----------
    from C3DExport import convert_c3d_type3

    result = convert_c3d_type3("data/trial.c3d", "output/")
    print(result['trc_path'])
    print(result['mot_path'])

Available Converters
--------------------
    convert_c3d_type1 -- 6-channel AMTI Type 1 (direct COP measurement)
    convert_c3d_type2 -- 6-channel AMTI Type 2 (moment-derived COP)
    convert_c3d_type3 -- 8-channel Kistler Type 3 (threshold stance detection)

Utility Functions (from C3DExport.utils)
----------------------------------------
    rotation_matrix, apply_rotation, chain_rotations
    compute_forceplate_type1, compute_forceplate_type2, compute_forceplate_type3
    plate_local_to_lab, lab_to_opensim_force
    butter_lowpass_filter, resample_to_target_rate
    detect_stance_phase, detect_stance_phase_from_peak
    detect_cop_anomalies, correct_cop_slope

I/O Functions (from C3DExport.io)
---------------------------------
    write_trc, write_mot
"""

from .core import convert_c3d_type1, convert_c3d_type2, convert_c3d_type3

__all__ = [
    'convert_c3d_type1',
    'convert_c3d_type2',
    'convert_c3d_type3',
]

__version__ = '1.0.0'
