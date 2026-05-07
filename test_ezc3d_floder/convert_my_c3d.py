"""
Batch C3D -> OpenSim converter.

Edit the paths below, then run:
    python convert_my_c3d.py
"""

import os
import sys

# === EDIT THESE PATHS ===
C3D_DIR   = r"E:\Python_Learn\c3dexport\test_ezc3d_floder\Sample10"       # folder containing .c3d files
OUTPUT_DIR = r"E:\Python_Learn\c3dexport\test_ezc3d_floder\output"          # where .trc / .mot / record.csv go
# =========================

# Add parent directory to import C3DExport
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from C3DExport import convert_c3d_type1, convert_c3d_type2, convert_c3d_type3
from C3DExport.core import _read_c3d


def auto_convert(c3d_path, output_dir, **kwargs):
    """Read fp_types from C3D and call the matching converter."""
    data = _read_c3d(c3d_path)
    fp_types = data['fp_types']

    if not fp_types:
        print(f"  [SKIP] No force plate type found: {c3d_path}")
        return None

    fp_type = fp_types[0]
    if fp_type == 1:
        print(f"  -> Type 1 (AMTI 6-ch, direct COP)")
        return convert_c3d_type1(c3d_path, output_dir, **kwargs)
    elif fp_type == 2:
        print(f"  -> Type 2 (AMTI 6-ch, moment-derived COP)")
        return convert_c3d_type2(c3d_path, output_dir, **kwargs)
    elif fp_type == 3:
        print(f"  -> Type 3 (Kistler 8-ch)")
        return convert_c3d_type3(c3d_path, output_dir, **kwargs)
    else:
        print(f"  [SKIP] Unsupported force plate type {fp_type}: {c3d_path}")
        return None


def main():
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    c3d_files = sorted([
        f for f in os.listdir(C3D_DIR) if f.lower().endswith('.c3d')
    ])

    if not c3d_files:
        print(f"No .c3d files found in {C3D_DIR}")
        return

    print(f"Found {len(c3d_files)} C3D file(s) in {C3D_DIR}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("=" * 60)

    for i, fname in enumerate(c3d_files, 1):
        c3d_path = os.path.join(C3D_DIR, fname)
        print(f"\n[{i}/{len(c3d_files)}] {fname}")

        try:
            result = auto_convert(c3d_path, OUTPUT_DIR)
            if result:
                print(f"  .trc: {result['trc_path']}")
                print(f"  .mot: {result['mot_path']}")
        except Exception as e:
            print(f"  [ERROR] {e}")

    print(f"\n{'=' * 60}")
    print("Done!")


if __name__ == '__main__':
    main()
