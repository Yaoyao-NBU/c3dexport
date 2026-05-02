"""
C3DExport.io
============
File I/O functions for writing OpenSim .trc (markers) and .mot (forces) files.
"""

import os
import numpy as np


def write_trc(filepath, marker_data, marker_labels, frame_rate,
              units='mm', orig_start_frame=1):
    """Write marker data to OpenSim .trc file.

    Parameters
    ----------
    filepath : str
    marker_data : ndarray (n_frames, n_markers * 3) -- X Y Z interleaved
    marker_labels : list[str] -- marker names
    frame_rate : float -- Hz
    units : str -- 'mm' or 'm'
    orig_start_frame : int -- original first frame number
    """
    n_frames  = marker_data.shape[0]
    n_markers = len(marker_labels)
    filename  = os.path.basename(filepath)

    os.makedirs(os.path.dirname(filepath) or '.', exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        # Line 1: file type identifier
        f.write(f"PathFileType\t4\t(X/Y/Z)\t{filename}\n")

        # Line 2: header field names
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\t"
                "OrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")

        # Line 3: header field values
        f.write(f"{frame_rate:.0f}\t{frame_rate:.0f}\t{n_frames}\t{n_markers}\t"
                f"{units}\t{frame_rate:.0f}\t{orig_start_frame}\t{n_frames}\n")

        # Line 4: marker name row (each name spans 3 columns)
        name_parts = ["Frame#", "Time"]
        for lbl in marker_labels:
            name_parts.extend([lbl, "", ""])
        f.write("\t".join(name_parts) + "\n")

        # Line 5: axis sub-headers
        sub_parts = ["", ""]
        for i in range(n_markers):
            idx = i + 1
            sub_parts.extend([f"X{idx}", f"Y{idx}", f"Z{idx}"])
        f.write("\t".join(sub_parts) + "\n")

        # Line 6: blank line
        f.write("\n")

        # Data rows
        for i in range(n_frames):
            frame_num = i + 1
            time_val  = i / frame_rate
            row = [f"{frame_num}", f"{time_val:.6f}"]
            for j in range(marker_data.shape[1]):
                row.append(f"{marker_data[i, j]:.6f}")
            f.write("\t".join(row) + "\n")


def write_mot(filepath, force_data_per_plate, n_plates, frame_rate,
              filename=None):
    """Write force-plate data to OpenSim .mot file.

    Column layout:
      time
      [FP1 force 3] [FP1 COP 3]  [FP2 force 3] [FP2 COP 3] ...
      [FP1 torque 3]            [FP2 torque 3] ...

    Parameters
    ----------
    filepath : str
    force_data_per_plate : list[dict]
        Each dict has keys 'force' (n, 3), 'cop' (n, 3), 'torque' (n, 3)
    n_plates : int
    frame_rate : float -- Hz
    filename : str or None -- header filename (default: basename)
    """
    n_frames = force_data_per_plate[0]['force'].shape[0]
    fname    = filename or os.path.basename(filepath)

    os.makedirs(os.path.dirname(filepath) or '.', exist_ok=True)

    # Build column labels
    labels = ["time"]

    # Force + COP columns for all plates first
    for pi in range(n_plates):
        prefix = "" if pi == 0 else f"{pi}_"
        labels.extend([
            f"{prefix}ground_force_vx",
            f"{prefix}ground_force_vy",
            f"{prefix}ground_force_vz",
            f"{prefix}ground_force_px",
            f"{prefix}ground_force_py",
            f"{prefix}ground_force_pz",
        ])
    # Then torque columns for all plates
    for pi in range(n_plates):
        prefix = "" if pi == 0 else f"{pi}_"
        labels.extend([
            f"{prefix}ground_torque_x",
            f"{prefix}ground_torque_y",
            f"{prefix}ground_torque_z",
        ])

    n_cols = len(labels)

    with open(filepath, 'w', newline='') as f:
        # Header
        f.write(f"{fname}\n")
        f.write("version=1\n")
        f.write(f"nRows={n_frames}\n")
        f.write(f"nColumns={n_cols}\n")
        f.write("inDegrees=yes\n")
        f.write("endheader\n")

        # Column labels
        f.write("\t".join(labels) + "\n")

        # Data rows
        for i in range(n_frames):
            time_val = i / frame_rate
            row = [f"{time_val:.6f}"]

            # Force + COP for each plate
            for pi in range(n_plates):
                d = force_data_per_plate[pi]
                for v in d['force'][i, :]:
                    row.append(f"{v:.6f}")
                for v in d['cop'][i, :]:
                    row.append(f"{v:.6f}")

            # Torque for each plate
            for pi in range(n_plates):
                d = force_data_per_plate[pi]
                for v in d['torque'][i, :]:
                    row.append(f"{v:.6f}")

            f.write("\t".join(row) + "\n")
