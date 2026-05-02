"""
C3DExport.utils
===============
Core utility functions for C3D to OpenSim conversion.

Provides:
  - Rotation matrices and coordinate transforms
  - Kistler force-plate computation (8-channel and 6-channel)
  - Butterworth low-pass filtering and polyphase resampling
  - Stance phase detection (threshold-based and peak-based)
  - COP anomaly detection and slope correction
"""

import numpy as np
from math import gcd
from scipy.signal import butter, filtfilt, resample_poly


# ============================================================================
#  Rotation Matrices
# ============================================================================

def rotation_matrix(axis, angle_deg):
    """Generate a 3x3 rotation matrix (right-hand rule).

    Parameters
    ----------
    axis : str -- 'X', 'Y', or 'Z'
    angle_deg : float -- rotation angle in degrees

    Returns
    -------
    R : ndarray (3, 3)
    """
    theta = np.radians(angle_deg)
    c, s = np.cos(theta), np.sin(theta)

    if axis.upper() == 'X':
        return np.array([[1, 0,  0],
                         [0, c, -s],
                         [0, s,  c]])
    elif axis.upper() == 'Y':
        return np.array([[ c, 0, s],
                         [ 0, 1, 0],
                         [-s, 0, c]])
    elif axis.upper() == 'Z':
        return np.array([[c, -s, 0],
                         [s,  c, 0],
                         [0,  0, 1]])
    else:
        raise ValueError(f"Invalid axis '{axis}'. Use 'X', 'Y', or 'Z'.")


def chain_rotations(*rotations):
    """Chain multiple rotations. e.g. chain_rotations(('X', -90), ('Y', 90))

    Parameters
    ----------
    *rotations : tuples of (axis, angle_deg)

    Returns
    -------
    R : ndarray (3, 3) -- combined rotation matrix (applied left-to-right)
    """
    R = np.eye(3)
    for axis, angle in rotations:
        R = rotation_matrix(axis, angle) @ R
    return R


def apply_rotation(data_3xN, R):
    """Apply a 3x3 rotation matrix to (3, N) data.

    Parameters
    ----------
    data_3xN : ndarray (3, N)
    R : ndarray (3, 3)

    Returns
    -------
    rotated : ndarray (3, N)
    """
    return R @ data_3xN


# ============================================================================
#  Kistler Force-Plate Computation
# ============================================================================

FZ_THRESHOLD = 20.0  # N -- COP / Tz set to 0 when |Fz| below this


def compute_kistler_channel8(channels_8, a, b, az0):
    """Compute Type 2 force-plate data from 8-channel Kistler Type 3 raw data.

    Parameters
    ----------
    channels_8 : ndarray (8, N)
        [fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]
    a   : float -- sensor offset in local-Y (AP / walking), mm
    b   : float -- sensor offset in local-X (ML), mm
    az0 : float -- top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz (N, ground reaction),
          ax, ay (mm, COP from plate centre),
          Tz (N*mm, free vertical moment)
    """
    fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4 = channels_8

    Fx_raw = fx12 + fx34
    Fy_raw = fy14 + fy23
    Fz_raw = fz1 + fz2 + fz3 + fz4

    Mx = b * (fz1 + fz2 - fz3 - fz4)
    My = a * (-fz1 + fz2 + fz3 - fz4)
    Mz = b * (-fx12 + fx34) + a * (fy14 - fy23)

    # Transfer moments to plate surface
    Mxp = Mx + Fy_raw * az0
    Myp = My - Fx_raw * az0

    # COP (valid only when |Fz| >= threshold)
    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    with np.errstate(invalid='ignore', divide='ignore'):
        ax = np.where(valid, -Myp / Fz_raw, 0.0)
        ay = np.where(valid,  Mxp / Fz_raw, 0.0)

    # Free vertical moment
    Tz_raw = Mz - Fy_raw * ax + Fx_raw * ay
    Tz = np.where(valid, -Tz_raw, 0.0)

    # Negate to get ground reaction force (from ground ON the person)
    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


def compute_kistler_channel6(channels_6, az0):
    """Compute COP and free moment from 6-channel force-plate data.

    Parameters
    ----------
    channels_6 : ndarray (6, N)
        [Fx, Fy, Fz, Mx, My, Mz] in Kistler local coords
    az0 : float -- top-plate offset, mm

    Returns
    -------
    dict  Fx, Fy, Fz (N), ax, ay (mm), Tz (N*mm)
    """
    Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz = channels_6

    Mxp = Mx + Fy_raw * az0
    Myp = My - Fx_raw * az0

    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    with np.errstate(invalid='ignore', divide='ignore'):
        ax = np.where(valid, -Myp / Fz_raw, 0.0)
        ay = np.where(valid,  Mxp / Fz_raw, 0.0)

    Tz_raw = Mz - Fy_raw * ax + Fx_raw * ay
    Tz = np.where(valid, -Tz_raw, 0.0)

    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


# ============================================================================
#  Coordinate Transforms
# ============================================================================

def plate_local_to_lab(type2, corners):
    """Convert force data from Kistler plate-local to lab global.

    Kistler local:  X=right+,    Y=posterior+, Z=down+
    Lab global:     X=forward+,  Y=left+,      Z=up+

    Parameters
    ----------
    type2 : dict -- output of compute_kistler_channel8/6
    corners : ndarray (3, 4) -- FORCE_PLATFORM:CORNERS for this plate

    Returns
    -------
    dict  Fx, Fy, Fz (N), COPx, COPy, COPz (mm), Tz (N*mm)
    """
    plate_cx = np.mean(corners[0, :])
    plate_cy = np.mean(corners[1, :])

    Fx_lab = type2['Fy']
    Fy_lab = type2['Fx']
    Fz_lab = type2['Fz']

    COPx_lab = plate_cx - type2['ay']
    COPy_lab = plate_cy - type2['ax']
    COPz_lab = np.zeros_like(COPx_lab)

    Tz_lab = -type2['Tz']

    return dict(Fx=Fx_lab, Fy=Fy_lab, Fz=Fz_lab,
                COPx=COPx_lab, COPy=COPy_lab, COPz=COPz_lab,
                Tz=Tz_lab)


def lab_to_opensim_force(data_lab):
    """Convert force data from lab global to OpenSim (Y-up).

    Lab:     X=forward, Y=left,  Z=up
    OpenSim: X=forward, Y=up,    Z=right

    Parameters
    ----------
    data_lab : dict -- Fx, Fy, Fz, COPx, COPy, COPz, Tz (lab global)

    Returns
    -------
    dict -- same keys, values in OpenSim coordinate system
    """
    return dict(
        Fx  =  data_lab['Fx'],
        Fy  =  data_lab['Fz'],
        Fz  = -data_lab['Fy'],
        COPx =  data_lab['COPx'],
        COPy =  data_lab['COPz'],
        COPz = -data_lab['COPy'],
        Tz   =  data_lab['Tz'],
    )


# ============================================================================
#  Filtering & Resampling
# ============================================================================

def butter_lowpass_filter(data, cutoff, fs, order=4):
    """4th-order Butterworth low-pass filter (zero-phase via filtfilt).

    Parameters
    ----------
    data : ndarray (N,)
    cutoff : float -- cutoff frequency (Hz)
    fs : float -- sampling frequency (Hz)
    order : int -- filter order (default 4)

    Returns
    -------
    filtered : ndarray (N,)
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)


def resample_to_target_rate(data_1d, src_rate, tgt_rate):
    """Resample a 1-D signal from src_rate to tgt_rate (polyphase).

    Parameters
    ----------
    data_1d : ndarray (N,)
    src_rate : float -- source sampling rate (Hz)
    tgt_rate : float -- target sampling rate (Hz)

    Returns
    -------
    resampled : ndarray (M,)
    """
    up   = int(tgt_rate)
    down = int(src_rate)
    g    = gcd(up, down)
    up  //= g
    down //= g
    return resample_poly(data_1d, up, down)


# ============================================================================
#  Stance Phase Detection
# ============================================================================

def detect_stance_phase(vertical_force, threshold=30.0, pad_frames=25):
    """Detect stance phase from vertical ground reaction force (threshold method).

    Parameters
    ----------
    vertical_force : ndarray (N,) -- vertical GRF (positive = up in OpenSim Y)
    threshold : float -- force threshold (N, default 30)
    pad_frames : int -- extra frames before/after (default 25)

    Returns
    -------
    start_idx, end_idx, hs_idx, to_idx : int (0-indexed)
    """
    contact = np.abs(vertical_force) > threshold

    if not np.any(contact):
        raise ValueError("No stance phase detected (force never exceeds threshold)")

    contact_indices = np.where(contact)[0]
    hs_idx = int(contact_indices[0])
    to_idx = int(contact_indices[-1])

    start_idx = max(0, hs_idx - pad_frames)
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)

    return start_idx, end_idx, hs_idx, to_idx


def detect_stance_phase_from_peak(vertical_force, threshold=30.0, pad_frames=25):
    """Detect stance phase using peak-based method.

    Algorithm:
      1. Find peak (max absolute vertical force)
      2. Walk left from peak until force < threshold -> heel-strike
      3. Walk right from peak until force < threshold -> toe-off
      4. Add padding frames; zero out padding region

    Parameters
    ----------
    vertical_force : ndarray (N,) -- vertical GRF (OpenSim Y, positive = up)
    threshold : float -- force threshold (N, default 30)
    pad_frames : int -- padding frames (default 25)

    Returns
    -------
    start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val
    """
    abs_force = np.abs(vertical_force)

    peak_idx = int(np.argmax(abs_force))
    peak_val = float(abs_force[peak_idx])

    if peak_val <= threshold:
        raise ValueError(
            f"Peak value ({peak_val:.2f} N) <= threshold ({threshold} N), "
            "cannot detect stance"
        )

    hs_idx = 0
    for i in range(peak_idx, -1, -1):
        if abs_force[i] < threshold:
            hs_idx = i + 1
            break

    to_idx = len(vertical_force) - 1
    for i in range(peak_idx, len(vertical_force)):
        if abs_force[i] < threshold:
            to_idx = i - 1
            break

    start_idx = max(0, hs_idx - pad_frames)
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)

    return start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val


# ============================================================================
#  COP Anomaly Detection & Slope Correction
# ============================================================================

def detect_cop_anomalies(copx, copy, time, force_vertical,
                         threshold=20.0, jump_threshold=0.03):
    """Detect COPx/COPy anomaly frames (stance phase only).

    Dual detection:
      1. Force < threshold but COP != 0 -> anomaly
      2. Frame-to-frame COP jump > jump_threshold (m) -> anomaly

    Parameters
    ----------
    copx, copy : ndarray (N,) -- COP in m
    time : ndarray (N,) -- time in s
    force_vertical : ndarray (N,) -- vertical force in N
    threshold : float -- force threshold (N, default 20)
    jump_threshold : float -- COP jump threshold (m, default 0.03)

    Returns
    -------
    copx_anomaly, copy_anomaly : ndarray (N,) bool
    info : dict
    """
    n = len(copx)
    copx_anomaly = np.zeros(n, dtype=bool)
    copy_anomaly = np.zeros(n, dtype=bool)

    invalid_force = np.abs(force_vertical) < threshold
    copx_nonzero = np.abs(copx) > 1e-6
    copy_nonzero = np.abs(copy) > 1e-6
    copx_anomaly |= (invalid_force & copx_nonzero)
    copy_anomaly |= (invalid_force & copy_nonzero)

    stance_mask = np.abs(force_vertical) >= threshold
    if np.sum(stance_mask) > 1:
        stance_indices = np.where(stance_mask)[0]
        for arr, anomaly in [(copx, copx_anomaly), (copy, copy_anomaly)]:
            diffs = np.abs(np.diff(arr[stance_mask]))
            jump_frames = diffs > jump_threshold
            for k in range(len(jump_frames)):
                if jump_frames[k]:
                    anomaly[stance_indices[k]] = True
                    anomaly[stance_indices[k + 1]] = True

    info = {
        'copx_invalid_count': int(np.sum(invalid_force & copx_nonzero)),
        'copy_invalid_count': int(np.sum(invalid_force & copy_nonzero)),
        'copx_jump_count': int(np.sum(copx_anomaly) - np.sum(invalid_force & copx_nonzero)),
        'copy_jump_count': int(np.sum(copy_anomaly) - np.sum(invalid_force & copy_nonzero)),
        'copx_total_anomaly': int(np.sum(copx_anomaly)),
        'copy_total_anomaly': int(np.sum(copy_anomaly)),
    }

    return copx_anomaly, copy_anomaly, info


def correct_cop_slope(copx, copy, time, force_vertical,
                      threshold=20.0, middle_ratio=0.3, rate_multiplier=2.0):
    """Correct COPx/COPy outliers using middle-section slope fitting.

    Algorithm:
      1. Identify stance phase (force >= threshold)
      2. Fit linear slope from middle section (e.g. 30%-70%)
      3. Walk outward from middle; if frame-to-frame delta > rate_multiplier * |k*dt|, replace
      4. Non-stance COP set to zero

    Parameters
    ----------
    copx, copy : ndarray (N,) -- COP in m
    time : ndarray (N,) -- time in s
    force_vertical : ndarray (N,) -- vertical force in N
    threshold : float -- force threshold (N, default 20)
    middle_ratio : float -- middle section ratio (default 0.3)
    rate_multiplier : float -- outlier rate multiplier (default 2.0)

    Returns
    -------
    copx_corrected, copy_corrected : ndarray (N,)
    info : dict
    """
    copx_corrected = copx.copy()
    copy_corrected = copy.copy()

    invalid_mask = np.abs(force_vertical) < threshold
    copx_corrected[invalid_mask] = 0.0
    copy_corrected[invalid_mask] = 0.0

    stance_mask = ~invalid_mask
    stance_indices = np.where(stance_mask)[0]
    n_stance = len(stance_indices)

    info = {
        'copx_slope': 0.0, 'copy_slope': 0.0,
        'copx_outlier_count': 0, 'copy_outlier_count': 0,
    }

    if n_stance < 10:
        return copx_corrected, copy_corrected, info

    margin = int(n_stance * middle_ratio / 2)
    mid_start = margin
    mid_end = n_stance - margin
    if mid_end - mid_start < 5:
        return copx_corrected, copy_corrected, info

    mid_indices = np.arange(mid_start, mid_end)
    t_stance = time[stance_mask]

    for col_name, data in [('copx', copx_corrected), ('copy', copy_corrected)]:
        stance_vals = data[stance_mask].copy()

        mid_vals = stance_vals[mid_indices]
        k, b = np.polyfit(t_stance[mid_indices], mid_vals, 1)

        dt = np.median(np.diff(t_stance)) if len(t_stance) > 1 else 1.0 / 200
        expected_delta = abs(k * dt)
        outlier_delta = expected_delta * rate_multiplier

        corrected = stance_vals.copy()

        for i in range(mid_start - 1, -1, -1):
            if abs(corrected[i] - corrected[i + 1]) > outlier_delta:
                corrected[i] = corrected[i + 1] - k * dt
                info[f'{col_name}_outlier_count'] += 1

        for i in range(mid_end, n_stance):
            if abs(corrected[i] - corrected[i - 1]) > outlier_delta:
                corrected[i] = corrected[i - 1] + k * dt
                info[f'{col_name}_outlier_count'] += 1

        data[stance_mask] = corrected
        info[f'{col_name}_slope'] = k

    return copx_corrected, copy_corrected, info
