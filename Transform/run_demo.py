"""
Demo: Convert S15T1V11.c3d → OpenSim .trc + .mot
=================================================
Reads the sample C3D file from the Data folder and outputs
OpenSim-compatible files to the Transform/output folder.
"""

import os
import sys

# Ensure this script's directory is on the path for local imports
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from c3d_to_opensim import process_c3d

# ── Paths ────────────────────────────────────────────────────────────────────
C3D_FILE   = os.path.join(SCRIPT_DIR, '..', 'Data', 'S15T1V11.c3d')
OUTPUT_DIR = os.path.join(SCRIPT_DIR, 'output')

# ── Parameters ───────────────────────────────────────────────────────────────
MARKER_CUTOFF     = 6.0     # Hz — Butterworth low-pass for markers
FORCE_CUTOFF      = 50.0    # Hz — Butterworth low-pass for force data
STANCE_THRESHOLD  = 30.0    # N  — vertical-force threshold for stance detection
STANCE_PAD_FRAMES = 25      # frames to pad before / after stance

# ── Run ──────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    result = process_c3d(
        c3d_path          = C3D_FILE,
        output_dir        = OUTPUT_DIR,
        marker_cutoff     = MARKER_CUTOFF,
        force_cutoff      = FORCE_CUTOFF,
        stance_threshold  = STANCE_THRESHOLD,
        stance_pad_frames = STANCE_PAD_FRAMES,
    )

    print("\n" + "-" * 40)
    print("Output files:")
    for key in ('trc_path', 'mot_path', 'csv_path'):
        print(f"   {key}: {result[key]}")

    print("\nStance Info:")
    for key, val in result['stance_info'].items():
        print(f"   {key}: {val}")
