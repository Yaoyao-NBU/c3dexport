"""
C3DExport.core
==============
Main conversion pipeline: C3D -> OpenSim .trc + .mot

Provides three converters:
  - convert_c3d_type1: 6-channel AMTI Type 1 (direct COP measurement)
  - convert_c3d_type2: 6-channel AMTI Type 2 (moment-derived COP)
  - convert_c3d_type3: 8-channel Kistler Type 3 (threshold stance detection)
"""

import os
import numpy as np
import ezc3d

from .utils import (
    rotation_matrix, apply_rotation,
    compute_forceplate_type1, compute_forceplate_type2, compute_forceplate_type3,
    plate_local_to_lab, lab_to_opensim_force,
    butter_lowpass_filter, resample_to_target_rate,
    detect_stance_phase, detect_stance_phase_from_peak,
)
from .io import write_trc, write_mot


def _strip_label_prefix(label):
    """Remove prefix before colon (e.g. 'Trial:R.ASIS' -> 'R.ASIS')."""
    return label.split(':')[-1] if ':' in label else label


def _read_c3d(c3d_path):
    """Read C3D file and return parsed data.

    Returns
    -------
    dict with keys:
      point_rate, analog_rate, first_frame,
      pt_data, marker_labels, marker_units, n_markers, n_point_frames,
      analog_data, n_analog_frames,
      n_plates, channels, origin, corners, fp_types
    """
    c = ezc3d.c3d(c3d_path)

    point_rate  = float(c['header']['points']['frame_rate'])
    analog_rate = float(c['header']['analogs']['frame_rate'])
    first_frame = c['header']['points']['first_frame']

    pt_data       = c['data']['points']
    raw_labels    = c['parameters']['POINT']['LABELS']['value']
    marker_labels = [_strip_label_prefix(lbl) for lbl in raw_labels]

    n_markers      = len(marker_labels)
    n_point_frames = pt_data.shape[2]

    try:
        marker_units = c['parameters']['POINT']['UNITS']['value'][0]
    except Exception:
        marker_units = 'mm'

    analog_data     = c['data']['analogs'][0]
    n_analog_frames = analog_data.shape[1]

    fp_params = c['parameters']['FORCE_PLATFORM']
    n_plates  = int(fp_params['USED']['value'][0])
    channels  = fp_params['CHANNEL']['value']
    origin    = fp_params['ORIGIN']['value']
    corners   = np.array(fp_params['CORNERS']['value'])

    try:
        fp_types = [int(t) for t in fp_params['TYPE']['value']]
    except KeyError:
        fp_types = []
        print(f"  [WARN] FORCE_PLATFORM TYPE parameter not found, defaulting to empty list (missing in {c3d_path})")

    return dict(
        point_rate=point_rate, analog_rate=analog_rate, first_frame=first_frame,
        pt_data=pt_data, marker_labels=marker_labels,
        marker_units=marker_units, n_markers=n_markers,
        n_point_frames=n_point_frames,
        analog_data=analog_data, n_analog_frames=n_analog_frames,
        n_plates=n_plates, channels=channels, origin=origin, corners=corners,
        fp_types=fp_types,
    )


# ============================================================================
#  Shared Helpers
# ============================================================================

def _resample_analog(analog_data, analog_rate, point_rate, n_point_frames):
    """Resample analog data to marker rate."""
    if analog_rate != point_rate:
        n_target = n_point_frames
        resampled = np.zeros((analog_data.shape[0], n_target))
        for ch in range(analog_data.shape[0]):
            rs = resample_to_target_rate(analog_data[ch, :], analog_rate, point_rate)
            resampled[ch, :] = rs[:n_target]
        return resampled, point_rate
    return analog_data, analog_rate


def _compute_forces(analog_data, channels, origin, n_plates, fp_types):
    """Dispatch force computation based on force plate type."""
    if not fp_types:
        raise ValueError("Cannot determine force plate type (fp_types is empty). "
                         "Please specify the correct converter manually.")

    fp_type = fp_types[0]
    if fp_type == 1:
        return _compute_type1_6ch(analog_data, channels, origin, n_plates)
    elif fp_type == 2:
        return _compute_type2_6ch(analog_data, channels, origin, n_plates)
    elif fp_type == 3:
        return _compute_type2_8ch(analog_data, channels, origin, n_plates)
    else:
        raise ValueError(f"Unsupported force plate type: {fp_type}. "
                         f"Supported types: 1 (AMTI 6-ch), 2/3 (Kistler 8-ch).")


def _compute_type1_6ch(analog_data, channels, origin, n_plates):
    """Compute force data from 6-channel AMTI raw data (Type 1)."""
    plate_type1 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch6    = analog_data[ch_idx, :]
        az0    = float(origin[2, pi])
        plate_type1.append(compute_forceplate_type1(ch6, az0))
    return plate_type1


def _compute_type2_8ch(analog_data, channels, origin, n_plates):
    """Compute force data from 8-channel Kistler raw data (Type 3)."""
    plate_type2 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch8    = analog_data[ch_idx, :]
        a   = float(origin[0, pi])
        b   = float(origin[1, pi])
        az0 = float(origin[2, pi])
        plate_type2.append(compute_forceplate_type3(ch8, a, b, az0))
    return plate_type2


def _compute_type2_6ch(analog_data, channels, origin, n_plates):
    """Compute force data from 6-channel raw data (Type 2)."""
    plate_type2 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch6    = analog_data[ch_idx, :]
        az0    = float(origin[2, pi])
        plate_type2.append(compute_forceplate_type2(ch6, az0))
    return plate_type2


def _transform_forces(plate_type2, corners, n_plates):
    """Transform force data: plate-local -> lab -> OpenSim."""
    R_lab2os = rotation_matrix('X', -90)
    plate_os = []
    for pi in range(n_plates):
        lab_data = plate_local_to_lab(plate_type2[pi], corners[:, :, pi])
        plate_os.append(lab_to_opensim_force(lab_data))
    return plate_os, R_lab2os


def _rotate_markers(pt_data, n_markers, R_lab2os):
    """Rotate marker data from lab to OpenSim coordinates."""
    markers_xyz = pt_data[:3, :, :]
    markers_os  = np.zeros_like(markers_xyz)
    for mi in range(n_markers):
        markers_os[:, mi, :] = apply_rotation(markers_xyz[:, mi, :], R_lab2os)
    return markers_os


def _filter_markers(markers_os, n_markers, marker_cutoff, point_rate):
    """Apply low-pass filter to marker data."""
    for mi in range(n_markers):
        for axis in range(3):
            markers_os[axis, mi, :] = butter_lowpass_filter(
                markers_os[axis, mi, :], marker_cutoff, point_rate
            )
    return markers_os


def _early_pipeline(data, marker_cutoff, force_cutoff, log):
    """Common pipeline steps: resample, compute forces, transform, filter markers.

    Returns: plate_os, markers_os, R_lab2os, analog_rate
    """
    log("\n[Step 2] Resampling analog signals ...")
    analog_data, analog_rate = _resample_analog(
        data['analog_data'], data['analog_rate'],
        data['point_rate'], data['n_point_frames']
    )

    log("\n[Step 3] Computing forces from raw channels ...")
    plate_type2 = _compute_forces(
        analog_data, data['channels'], data['origin'],
        data['n_plates'], data['fp_types']
    )

    log("\n[Step 4-5] Coordinate transforms (plate-local -> lab -> OpenSim) ...")
    plate_os, R_lab2os = _transform_forces(
        plate_type2, data['corners'], data['n_plates']
    )
    markers_os = _rotate_markers(data['pt_data'], data['n_markers'], R_lab2os)

    log(f"\n[Step 6] Marker filtering: {marker_cutoff} Hz ...")
    markers_os = _filter_markers(markers_os, data['n_markers'], marker_cutoff, data['point_rate'])

    return plate_os, markers_os, R_lab2os, analog_rate


def _build_structured_plates_v1(plate_os, n_plates, cop_threshold=20.0):
    """Build structured plate data with COP threshold and unit conversion (V1/V2 shared)."""
    structured = []
    for pi in range(n_plates):
        d = plate_os[pi]
        n = len(d['Fx'])
        valid = np.abs(d['Fy']) >= cop_threshold
        structured.append(dict(
            Fx=d['Fx'], Fy=d['Fy'], Fz=d['Fz'],
            COPx=np.where(valid, d['COPx'], 0.0) / 1000.0,
            COPy=np.where(valid, d['COPy'], 0.0),
            COPz=np.where(valid, d['COPz'], 0.0) / 1000.0,
            Tz=d['Tz'] / 1000.0,
        ))
    return structured


def _detect_active_plate(structured_plates, n_plates):
    """Detect which plate has the strongest vertical force."""
    max_fy = [np.max(np.abs(structured_plates[pi]['Fy'])) for pi in range(n_plates)]
    return int(np.argmax(max_fy))


def _cut_and_zero_padding(structured_plates, n_plates, start_idx, end_idx, hs_idx, to_idx):
    """Cut force data to stance range and zero padding frames."""
    for pi in range(n_plates):
        for key in ('Fx', 'Fy', 'Fz', 'COPx', 'COPy', 'COPz', 'Tz'):
            structured_plates[pi][key] = structured_plates[pi][key][start_idx:end_idx+1]

    rel_hs = hs_idx - start_idx
    rel_to = to_idx - start_idx

    for pi in range(n_plates):
        for key in ('Fx', 'Fy', 'Fz', 'COPx', 'COPy', 'COPz', 'Tz'):
            structured_plates[pi][key][:rel_hs] = 0
            structured_plates[pi][key][rel_to+1:] = 0

    return structured_plates, rel_hs, rel_to


def _markers_to_flat(markers_os, n_markers):
    """Reshape markers (3, n_markers, n_frames) -> (n_frames, n_markers*3)."""
    n_out = markers_os.shape[2]
    flat  = np.zeros((n_out, n_markers * 3))
    for mi in range(n_markers):
        flat[:, mi * 3]     = markers_os[0, mi, :]
        flat[:, mi * 3 + 1] = markers_os[1, mi, :]
        flat[:, mi * 3 + 2] = markers_os[2, mi, :]
    return flat


def _build_mot_plates(structured_plates, n_plates):
    """Convert structured plate dicts to format expected by write_mot."""
    mot_plates = []
    for pi in range(n_plates):
        d = structured_plates[pi]
        n = len(d['Fy'])
        force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])
        cop    = np.column_stack([d['COPx'], d['COPy'], d['COPz']])
        torque = np.column_stack([np.zeros(n), d['Tz'], np.zeros(n)])
        mot_plates.append(dict(force=force, cop=cop, torque=torque))
    return mot_plates


def _write_csv_record(csv_path, stance_info):
    """Append stance info record to CSV file."""
    import pandas as pd
    file_exists = os.path.isfile(csv_path)
    pd.DataFrame([stance_info]).to_csv(
        csv_path, mode='a', index=False,
        header=not file_exists, encoding='utf-8-sig'
    )


def _write_outputs(output_dir, trial, marker_flat, marker_labels, marker_units,
                   point_rate, mot_plates, n_plates, csv_path, stance_info, log):
    """Common output-writing steps for all converters."""
    trc_path = os.path.join(output_dir, f"{trial}.trc")
    write_trc(trc_path, marker_flat, marker_labels, point_rate,
              units=marker_units, orig_start_frame=1)

    mot_path = os.path.join(output_dir, f"{trial}.mot")
    write_mot(mot_path, mot_plates, n_plates, point_rate,
              filename=f"{trial}.mot")

    _write_csv_record(csv_path, stance_info)

    log(f"\n{'=' * 60}")
    log(f"  [DONE] .trc: {trc_path}")
    log(f"  [DONE] .mot: {mot_path}")
    log(f"{'=' * 60}")

    return trc_path, mot_path


# ============================================================================
#  Main Conversion Functions
# ============================================================================

def convert_c3d_type2(c3d_path, output_dir,
                      marker_cutoff=6.0, force_cutoff=50.0,
                      stance_threshold=30.0, stance_pad_frames=25,
                      verbose=True):
    """Convert 6-channel AMTI (Type 2) C3D to OpenSim .trc + .mot.

    Parameters
    ----------
    c3d_path : str -- input C3D file with 6-channel force data
    output_dir : str -- output directory
    marker_cutoff : float -- marker low-pass cutoff (Hz, default 6)
    force_cutoff : float -- force low-pass cutoff (Hz, default 50)
    stance_threshold : float -- stance detection force threshold (N, default 30)
    stance_pad_frames : int -- padding frames (default 25)
    verbose : bool -- print progress messages (default True)

    Returns
    -------
    dict with keys: trc_path, mot_path, csv_path, stance_info
    """
    def _log(msg):
        if verbose:
            print(msg)

    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    _log("=" * 60)
    _log(f"  C3D (6-Channel AMTI Type 2) -> OpenSim Converter")
    _log(f"  Input : {c3d_path}")
    _log(f"  Output: {output_dir}")
    _log("=" * 60)

    # Step 1: Read C3D
    _log("\n[Step 1] Reading C3D file ...")
    data = _read_c3d(c3d_path)
    _log(f"  Markers: {data['n_markers']} markers, {data['n_point_frames']} frames @ {data['point_rate']} Hz")
    _log(f"  Forces: {data['n_plates']} plate(s), {data['n_analog_frames']} frames @ {data['analog_rate']} Hz")
    _log(f"  Force plate types: {data['fp_types']}")

    # Step 2-6: Common pipeline
    plate_os, markers_os, R_lab2os, analog_rate = _early_pipeline(
        data, marker_cutoff, force_cutoff, _log
    )

    # Step 5-6: Restructure + unit conversion
    _log("\n[Step 5-6] Force array restructuring & unit conversion ...")
    structured_plates = []
    for pi in range(data['n_plates']):
        d = plate_os[pi]
        n = len(d['Fx'])
        force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])
        cop    = np.column_stack([d['COPx'] / 1000.0, d['COPy'] / 1000.0, d['COPz'] / 1000.0])
        torque = np.column_stack([np.zeros(n), d['Tz'] / 1000.0, np.zeros(n)])
        structured_plates.append(dict(force=force, cop=cop, torque=torque))

    # Step 7: Filter & resample
    _log("\n[Step 7] Filtering & resampling ...")
    for pi in range(data['n_plates']):
        for key in ('force', 'cop', 'torque'):
            arr = structured_plates[pi][key]
            for col in range(arr.shape[1]):
                if np.any(arr[:, col] != 0):
                    arr[:, col] = butter_lowpass_filter(arr[:, col], force_cutoff, data['analog_rate'])

    if data['analog_rate'] != data['point_rate']:
        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                arr = structured_plates[pi][key]
                n_target = data['n_point_frames']
                resampled = np.zeros((n_target, arr.shape[1]))
                for col in range(arr.shape[1]):
                    rs = resample_to_target_rate(arr[:, col], data['analog_rate'], data['point_rate'])
                    resampled[:, col] = rs[:n_target]
                structured_plates[pi][key] = resampled

    # Step 8: Stance detection
    _log("\n[Step 8] Stance phase detection & extraction ...")
    max_fy = [np.max(np.abs(structured_plates[pi]['force'][:, 1])) for pi in range(data['n_plates'])]
    active_plate = int(np.argmax(max_fy))
    _log(f"  Active plate: FP{active_plate+1}")

    csv_path = os.path.join(output_dir, "record.csv")

    try:
        start_idx, end_idx, hs_idx, to_idx = detect_stance_phase(
            structured_plates[active_plate]['force'][:, 1],
            threshold=stance_threshold, pad_frames=stance_pad_frames,
        )
        n_cut = end_idx - start_idx + 1
        _log(f"  Heel Strike: frame {hs_idx+1}, Toe Off: frame {to_idx+1}")
        _log(f"  Cut range: {start_idx+1}-{end_idx+1} ({n_cut} frames)")

        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                structured_plates[pi][key] = structured_plates[pi][key][start_idx:end_idx+1, :]

        rel_hs = hs_idx - start_idx
        rel_to = to_idx - start_idx
        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                structured_plates[pi][key][:rel_hs, :] = 0
                structured_plates[pi][key][rel_to+1:, :] = 0

        markers_os = markers_os[:, :, start_idx:end_idx+1]

        stance_info = dict(
            trial=trial, heel_strike_frame=hs_idx+1, toe_off_frame=to_idx+1,
            cut_start_frame=start_idx+1, cut_end_frame=end_idx+1, total_frames=n_cut,
        )
    except ValueError as e:
        _log(f"  [WARN] {e} -> using all data")
        n_cut = structured_plates[0]['force'].shape[0]
        stance_info = dict(
            trial=trial, heel_strike_frame='N/A', toe_off_frame='N/A',
            cut_start_frame=1, cut_end_frame=n_cut, total_frames=n_cut,
        )

    # Step 9: Write output files
    _log("\n[Step 9] Writing OpenSim files ...")
    marker_flat = _markers_to_flat(markers_os, data['n_markers'])
    trc_path, mot_path = _write_outputs(
        output_dir, trial, marker_flat, data['marker_labels'], data['marker_units'],
        data['point_rate'], structured_plates, data['n_plates'],
        csv_path, stance_info, _log
    )

    return dict(trc_path=trc_path, mot_path=mot_path, csv_path=csv_path, stance_info=stance_info)


def convert_c3d_type3(c3d_path, output_dir,
                      marker_cutoff=6.0, force_cutoff=50.0,
                      stance_threshold=30.0, stance_pad_frames=25,
                      verbose=True):
    """Convert 8-channel Kistler (Type 3) C3D to OpenSim .trc + .mot (threshold stance detection).

    Parameters
    ----------
    c3d_path : str -- input C3D file path
    output_dir : str -- output directory
    marker_cutoff : float -- marker low-pass cutoff (Hz, default 6)
    force_cutoff : float -- force low-pass cutoff (Hz, default 50)
    stance_threshold : float -- stance detection force threshold (N, default 30)
    stance_pad_frames : int -- padding frames before/after stance (default 25)
    verbose : bool -- print progress messages (default True)

    Returns
    -------
    dict with keys: trc_path, mot_path, csv_path, stance_info
    """
    def _log(msg):
        if verbose:
            print(msg)

    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    _log("=" * 60)
    _log(f"  C3D (8-Channel Kistler Type 3) -> OpenSim Converter")
    _log(f"  Input : {c3d_path}")
    _log(f"  Output: {output_dir}")
    _log("=" * 60)

    # Step 1: Read C3D
    _log("\n[Step 1] Reading C3D file ...")
    data = _read_c3d(c3d_path)
    _log(f"  Markers: {data['n_markers']} markers, {data['n_point_frames']} frames @ {data['point_rate']} Hz")
    _log(f"  Forces: {data['n_plates']} plate(s), {data['n_analog_frames']} frames @ {data['analog_rate']} Hz")
    _log(f"  Force plate types: {data['fp_types']}")

    # Step 2-6: Common pipeline
    plate_os, markers_os, R_lab2os, analog_rate = _early_pipeline(
        data, marker_cutoff, force_cutoff, _log
    )

    # Step 7: COP threshold + unit conversion
    _log("\n[Step 7] COP threshold & unit conversion ...")
    structured_plates = _build_structured_plates_v1(plate_os, data['n_plates'])

    # Step 8: Stance detection
    _log("\n[Step 8] Stance phase detection & extraction ...")
    active_plate = _detect_active_plate(structured_plates, data['n_plates'])
    _log(f"  Active plate: FP{active_plate+1}")

    csv_path = os.path.join(output_dir, "record.csv")

    try:
        start_idx, end_idx, hs_idx, to_idx = detect_stance_phase(
            structured_plates[active_plate]['Fy'],
            threshold=stance_threshold, pad_frames=stance_pad_frames,
        )
        n_cut = end_idx - start_idx + 1
        _log(f"  Heel Strike: frame {hs_idx+1}, Toe Off: frame {to_idx+1}")
        _log(f"  Cut range: {start_idx+1}-{end_idx+1} ({n_cut} frames)")

        structured_plates, _, _ = _cut_and_zero_padding(
            structured_plates, data['n_plates'],
            start_idx, end_idx, hs_idx, to_idx
        )
        markers_os = markers_os[:, :, start_idx:end_idx+1]

        stance_info = dict(
            trial=trial, heel_strike_frame=hs_idx+1, toe_off_frame=to_idx+1,
            cut_start_frame=start_idx+1, cut_end_frame=end_idx+1, total_frames=n_cut,
        )
    except ValueError as e:
        _log(f"  [WARN] {e} -> using all data")
        n_cut = len(structured_plates[0]['Fy'])
        stance_info = dict(
            trial=trial, heel_strike_frame='N/A', toe_off_frame='N/A',
            cut_start_frame=1, cut_end_frame=n_cut, total_frames=n_cut,
        )

    # Step 9: Write output files
    _log("\n[Step 9] Writing OpenSim files ...")
    marker_flat = _markers_to_flat(markers_os, data['n_markers'])
    mot_plates = _build_mot_plates(structured_plates, data['n_plates'])
    trc_path, mot_path = _write_outputs(
        output_dir, trial, marker_flat, data['marker_labels'], data['marker_units'],
        data['point_rate'], mot_plates, data['n_plates'],
        csv_path, stance_info, _log
    )

    return dict(trc_path=trc_path, mot_path=mot_path, csv_path=csv_path, stance_info=stance_info)


def convert_c3d_type1(c3d_path, output_dir,
                      marker_cutoff=6.0, force_cutoff=50.0,
                      stance_threshold=30.0, stance_pad_frames=25,
                      verbose=True):
    """Convert 6-channel AMTI (Type 1) C3D to OpenSim .trc + .mot.

    Type 1 force plates output forces and COP directly (Px, Py are COP
    coordinates, Mz is free vertical moment). No moment-to-COP conversion.

    Parameters
    ----------
    c3d_path : str -- input C3D file with 6-channel Type 1 force data
    output_dir : str -- output directory
    marker_cutoff : float -- marker low-pass cutoff (Hz, default 6)
    force_cutoff : float -- force low-pass cutoff (Hz, default 50)
    stance_threshold : float -- stance detection force threshold (N, default 30)
    stance_pad_frames : int -- padding frames (default 25)
    verbose : bool -- print progress messages (default True)

    Returns
    -------
    dict with keys: trc_path, mot_path, csv_path, stance_info
    """
    def _log(msg):
        if verbose:
            print(msg)

    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    _log("=" * 60)
    _log(f"  C3D (6-Channel AMTI Type 1) -> OpenSim Converter")
    _log(f"  Input : {c3d_path}")
    _log(f"  Output: {output_dir}")
    _log("=" * 60)

    # Step 1: Read C3D
    _log("\n[Step 1] Reading C3D file ...")
    data = _read_c3d(c3d_path)
    _log(f"  Markers: {data['n_markers']} markers, {data['n_point_frames']} frames @ {data['point_rate']} Hz")
    _log(f"  Forces: {data['n_plates']} plate(s), {data['n_analog_frames']} frames @ {data['analog_rate']} Hz")
    _log(f"  Force plate types: {data['fp_types']}")

    # Step 2-6: Common pipeline
    plate_os, markers_os, R_lab2os, analog_rate = _early_pipeline(
        data, marker_cutoff, force_cutoff, _log
    )

    # Step 5-6: Restructure + unit conversion
    _log("\n[Step 5-6] Force array restructuring & unit conversion ...")
    structured_plates = []
    for pi in range(data['n_plates']):
        d = plate_os[pi]
        n = len(d['Fx'])
        force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])
        cop    = np.column_stack([d['COPx'] / 1000.0, d['COPy'] / 1000.0, d['COPz'] / 1000.0])
        torque = np.column_stack([np.zeros(n), d['Tz'] / 1000.0, np.zeros(n)])
        structured_plates.append(dict(force=force, cop=cop, torque=torque))

    # Step 7: Filter & resample
    _log("\n[Step 7] Filtering & resampling ...")
    for pi in range(data['n_plates']):
        for key in ('force', 'cop', 'torque'):
            arr = structured_plates[pi][key]
            for col in range(arr.shape[1]):
                if np.any(arr[:, col] != 0):
                    arr[:, col] = butter_lowpass_filter(arr[:, col], force_cutoff, data['analog_rate'])

    if data['analog_rate'] != data['point_rate']:
        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                arr = structured_plates[pi][key]
                n_target = data['n_point_frames']
                resampled = np.zeros((n_target, arr.shape[1]))
                for col in range(arr.shape[1]):
                    rs = resample_to_target_rate(arr[:, col], data['analog_rate'], data['point_rate'])
                    resampled[:, col] = rs[:n_target]
                structured_plates[pi][key] = resampled

    # Step 8: Stance detection
    _log("\n[Step 8] Stance phase detection & extraction ...")
    max_fy = [np.max(np.abs(structured_plates[pi]['force'][:, 1])) for pi in range(data['n_plates'])]
    active_plate = int(np.argmax(max_fy))
    _log(f"  Active plate: FP{active_plate+1}")

    csv_path = os.path.join(output_dir, "record.csv")

    try:
        start_idx, end_idx, hs_idx, to_idx = detect_stance_phase(
            structured_plates[active_plate]['force'][:, 1],
            threshold=stance_threshold, pad_frames=stance_pad_frames,
        )
        n_cut = end_idx - start_idx + 1
        _log(f"  Heel Strike: frame {hs_idx+1}, Toe Off: frame {to_idx+1}")
        _log(f"  Cut range: {start_idx+1}-{end_idx+1} ({n_cut} frames)")

        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                structured_plates[pi][key] = structured_plates[pi][key][start_idx:end_idx+1, :]

        rel_hs = hs_idx - start_idx
        rel_to = to_idx - start_idx
        for pi in range(data['n_plates']):
            for key in ('force', 'cop', 'torque'):
                structured_plates[pi][key][:rel_hs, :] = 0
                structured_plates[pi][key][rel_to+1:, :] = 0

        markers_os = markers_os[:, :, start_idx:end_idx+1]

        stance_info = dict(
            trial=trial, heel_strike_frame=hs_idx+1, toe_off_frame=to_idx+1,
            cut_start_frame=start_idx+1, cut_end_frame=end_idx+1, total_frames=n_cut,
        )
    except ValueError as e:
        _log(f"  [WARN] {e} -> using all data")
        n_cut = structured_plates[0]['force'].shape[0]
        stance_info = dict(
            trial=trial, heel_strike_frame='N/A', toe_off_frame='N/A',
            cut_start_frame=1, cut_end_frame=n_cut, total_frames=n_cut,
        )

    # Step 9: Write output files
    _log("\n[Step 9] Writing OpenSim files ...")
    marker_flat = _markers_to_flat(markers_os, data['n_markers'])
    trc_path, mot_path = _write_outputs(
        output_dir, trial, marker_flat, data['marker_labels'], data['marker_units'],
        data['point_rate'], structured_plates, data['n_plates'],
        csv_path, stance_info, _log
    )

    return dict(trc_path=trc_path, mot_path=mot_path, csv_path=csv_path, stance_info=stance_info)
