"""
Transform Utilities for C3D → OpenSim Conversion
=================================================
Provides coordinate rotation, force computation, filtering, resampling,
stance detection, and .trc / .mot file I/O functions.

Coordinate Systems:
  Kistler plate-local:  X=right+,    Y=posterior+, Z=down+
  Lab global:           X=forward+,  Y=left+,      Z=up+
  OpenSim:              X=forward+,  Y=up+,        Z=right+
"""

import os
import numpy as np
from math import gcd
from scipy.signal import butter, filtfilt, resample_poly


# ═══════════════════════════════════════════════════════════════════════════════
#  Rotation Matrices
# ═══════════════════════════════════════════════════════════════════════════════

def rotation_matrix(axis, angle_deg):
    """
    Generate a 3×3 rotation matrix (right-hand rule).

    Parameters
    ----------
    axis : str — 'X', 'Y', or 'Z'
    angle_deg : float — rotation angle in degrees

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
    """
    Chain multiple rotations.  e.g. chain_rotations(('X', -90), ('Y', 90))

    Parameters
    ----------
    *rotations : tuples of (axis, angle_deg)

    Returns
    -------
    R : ndarray (3, 3) — combined rotation matrix (applied left-to-right)
    """
    R = np.eye(3)
    for axis, angle in rotations:
        R = rotation_matrix(axis, angle) @ R
    return R


def apply_rotation(data_3xN, R):
    """
    Apply a 3×3 rotation matrix to (3, N) data.

    Parameters
    ----------
    data_3xN : ndarray (3, N)
    R : ndarray (3, 3)

    Returns
    -------
    rotated : ndarray (3, N)
    """
    return R @ data_3xN


# ═══════════════════════════════════════════════════════════════════════════════
#  Kistler Type 3 → Type 2  Force-Plate Computation
# ═══════════════════════════════════════════════════════════════════════════════

FZ_THRESHOLD = 20.0   # N — COP / Tz set to 0 when |Fz| below this


def compute_kistler_channel8(channels_8, a, b, az0):
    """
    Compute Type 2 force-plate data from 8-channel Kistler Type 3 raw data.

    Parameters
    ----------
    channels_8 : ndarray (8, N)
        [fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]
    a   : float — sensor offset in local-Y (AP / walking), mm
    b   : float — sensor offset in local-X (ML), mm
    az0 : float — top-plate offset, mm (typically negative)

    Returns
    -------
    dict  Fx, Fy, Fz  (N, ground reaction),
          ax, ay      (mm, COP from plate centre),
          Tz          (N·mm, free vertical moment)
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

    # COP (valid only when |Fz| ≥ threshold)
    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    with np.errstate(invalid='ignore', divide='ignore'):
        ax = np.where(valid, -Myp / Fz_raw, 0.0)
        ay = np.where(valid,  Mxp / Fz_raw, 0.0)

    # Free vertical moment
    Tz_raw = Mz - Fy_raw * ax + Fx_raw * ay
    Tz = np.where(valid, -Tz_raw, 0.0)

    # Negate to get ground reaction force (from ground ON the person)
    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


def compute_kistler_channel6(channels_6 , az0):
    """
    Compute COP and free moment from 6-channel force-plate data (Fx, Fy, Fz, Mx, My, Mz).

    The input data is already in Kistler plate-local coordinates:
      X=right+, Y=posterior+, Z=down+

    Parameters
    ----------
    channels_6 : ndarray (6, N)
        [Fx, Fy, Fz, Mx, My, Mz] — forces and moments in Kistler local coords

    Returns
    -------
    dict  Fx, Fy, Fz  (N, ground reaction forces),
          ax, ay      (N, mm, COP from plate centre),
          Tz          (N, N·mm, free vertical moment)
    """
    Fx_raw, Fy_raw, Fz_raw, Mx, My, Mz = channels_6

    Mxp = Mx + Fy_raw * az0
    Myp = My - Fx_raw * az0
    # COP (valid only when |Fz| ≥ threshold)
    # Note: COP is computed from the given moments, assuming they are already
    # at the plate surface or at the measurement reference point
    valid = np.abs(Fz_raw) >= FZ_THRESHOLD
    with np.errstate(invalid='ignore', divide='ignore'):
        ax = np.where(valid, -Myp / Fz_raw, 0.0)
        ay = np.where(valid,  Mxp / Fz_raw, 0.0)

    # Free vertical moment (moment about the COP)
    Tz_raw = Mz - Fy_raw * ax + Fx_raw * ay
    Tz = np.where(valid, -Tz_raw, 0.0)

    # Negate to get ground reaction force (from ground ON the person)
    return dict(Fx=-Fx_raw, Fy=-Fy_raw, Fz=-Fz_raw, ax=ax, ay=ay, Tz=Tz)


# ═══════════════════════════════════════════════════════════════════════════════
#  Coordinate Transforms
# ═══════════════════════════════════════════════════════════════════════════════

def plate_local_to_lab(type2, corners):
    """
    Convert Type 2 force data from Kistler plate-local → lab global.

    Kistler local:  X=right+,     Y=posterior+,  Z=down+
    Lab global:     X=forward+,   Y=left+,       Z=up+

    Mapping (ground reaction, double negation on axis flip):
        GRF_Lab_Fx =  type2['Fy']         (forward)
        GRF_Lab_Fy =  type2['Fx']         (left)
        GRF_Lab_Fz =  type2['Fz']         (up)
        COP_Lab_X  = plate_cx − ay        (forward, mm)
        COP_Lab_Y  = plate_cy − ax        (left, mm)
        COP_Lab_Z  = 0                    (ground level)
        Tz_lab     = −Tz_kistler          (about Z-up)

    Parameters
    ----------
    type2 : dict — output of compute_kistler_type2()
    corners : ndarray (3, 4) — FORCE_PLATFORM:CORNERS for this plate

    Returns
    -------
    dict  Fx, Fy, Fz (N),  COPx, COPy, COPz (mm),  Tz (N·mm)
    """
    plate_cx = np.mean(corners[0, :])   # lab X (forward)
    plate_cy = np.mean(corners[1, :])   # lab Y (left)

    # Forces — ground reaction in lab global
    Fx_lab = type2['Fy']    # forward
    Fy_lab = type2['Fx']    # left
    Fz_lab = type2['Fz']    # up

    # COP in lab global (mm)
    COPx_lab = plate_cx - type2['ay']    # forward
    COPy_lab = plate_cy - type2['ax']    # left
    COPz_lab = np.zeros_like(COPx_lab)   # ground level

    # Free moment: Kistler Z is down, lab Z is up → negate
    Tz_lab = -type2['Tz']

    return dict(Fx=Fx_lab, Fy=Fy_lab, Fz=Fz_lab,
                COPx=COPx_lab, COPy=COPy_lab, COPz=COPz_lab,
                Tz=Tz_lab)


def lab_to_opensim_force(data_lab):
    """
    Convert force data from lab global → OpenSim (Y-up).

    Lab:     X=forward, Y=left,  Z=up
    OpenSim: X=forward, Y=up,    Z=right

    Mapping:  X_os = X_lab,  Y_os = Z_lab,  Z_os = −Y_lab

    Parameters
    ----------
    data_lab : dict — Fx, Fy, Fz, COPx, COPy, COPz, Tz  (lab global)

    Returns
    -------
    dict — same keys, values in OpenSim coordinate system
    """
    return dict(
        Fx  =  data_lab['Fx'],       # forward
        Fy  =  data_lab['Fz'],       # up   (lab Z → OS Y)
        Fz  = -data_lab['Fy'],       # right (lab −Y → OS Z)
        COPx =  data_lab['COPx'],    # forward
        COPy =  data_lab['COPz'],    # up = 0  (ground)
        COPz = -data_lab['COPy'],    # right
        Tz   =  data_lab['Tz'],      # free vertical moment (stays about vert axis)
    )


# ═══════════════════════════════════════════════════════════════════════════════
#  Filtering & Resampling
# ═══════════════════════════════════════════════════════════════════════════════

def butter_lowpass_filter(data, cutoff, fs, order=4):
    """
    4th-order Butterworth low-pass filter (zero-phase via filtfilt).

    Parameters
    ----------
    data   : ndarray (N,)
    cutoff : float — cutoff frequency (Hz)
    fs     : float — sampling frequency (Hz)
    order  : int   — filter order (default 4)

    Returns
    -------
    filtered : ndarray (N,)
    """
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return filtfilt(b, a, data)


def resample_to_target_rate(data_1d, src_rate, tgt_rate):
    """
    Resample a 1-D signal from src_rate to tgt_rate (polyphase).

    Parameters
    ----------
    data_1d  : ndarray (N,)
    src_rate : float — source sampling rate (Hz)
    tgt_rate : float — target sampling rate (Hz)

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


# ═══════════════════════════════════════════════════════════════════════════════
#  Stance Phase Detection
# ═══════════════════════════════════════════════════════════════════════════════

def detect_stance_phase(vertical_force, threshold=30.0, pad_frames=25):
    """
    Detect stance phase from vertical ground reaction force.

    Parameters
    ----------
    vertical_force : ndarray (N,) — vertical GRF (positive = up in OpenSim Y)
    threshold      : float        — force threshold (N, default 30)
    pad_frames     : int          — extra frames before / after (default 25)

    Returns
    -------
    start_idx : int — 0-indexed start (with padding)
    end_idx   : int — 0-indexed end   (with padding)
    hs_idx    : int — 0-indexed heel-strike frame
    to_idx    : int — 0-indexed toe-off frame
    """
    contact = np.abs(vertical_force) > threshold

    if not np.any(contact):
        raise ValueError("No stance phase detected (force never exceeds threshold)")

    contact_indices = np.where(contact)[0]
    hs_idx = int(contact_indices[0])
    to_idx = int(contact_indices[-1])

    start_idx = max(0, hs_idx - pad_frames)   #边界保护作用
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)  #边界保护作用

    return start_idx, end_idx, hs_idx, to_idx


def detect_stance_phase_from_peak(vertical_force, threshold=30.0, pad_frames=25):
    """
    从峰值向两侧遍历检测stance阶段（适用于单力台踩踏场景）

    算法:
      1. 找到垂直力最大值(峰值)位置
      2. 从峰值向左遍历，找到第一个低于阈值的帧 → 定义为 load (heel-strike)
      3. 从峰值向右遍历，找到第一个低于阈值的帧 → 定义为 off (toe-off)
      4. 在 load/off 基础上添加 padding 帧，padding 区域的力数据手动归零

    Parameters
    ----------
    vertical_force : ndarray (N,) — 垂直GRF (OpenSim Y方向, 正=向上)
    threshold      : float        — 力阈值 (N, 默认30)
    pad_frames     : int          — 前后补偿帧数 (默认25)

    Returns
    -------
    start_idx : int — 截取起始索引 (含padding, 0-indexed)
    end_idx   : int — 截取结束索引 (含padding, 0-indexed)
    hs_idx    : int — heel-strike 帧索引 (0-indexed, 即 load 点)
    to_idx    : int — toe-off 帧索引 (0-indexed, 即 off 点)
    peak_idx  : int — 峰值帧索引 (0-indexed)
    peak_val  : float — 峰值力大小 (N)
    """
    abs_force = np.abs(vertical_force)

    # Step1: 找峰值
    peak_idx = int(np.argmax(abs_force))
    peak_val = float(abs_force[peak_idx])

    if peak_val <= threshold:
        raise ValueError(
            f"Peak value ({peak_val:.2f} N) <= threshold ({threshold} N), "
            "无法检测stance"
        )

    # Step2: 从峰值向左搜索，找到第一个 < threshold 的帧 → hs_idx
    hs_idx = 0  # 默认从最左开始
    for i in range(peak_idx, -1, -1):
        if abs_force[i] < threshold:
            hs_idx = i + 1  # 第一个低于阈值的帧的下一帧即为load点
            break

    # Step3: 从峰值向右搜索，找到第一个 < threshold 的帧 → to_idx
    to_idx = len(vertical_force) - 1  # 默认到最右
    for i in range(peak_idx, len(vertical_force)):
        if abs_force[i] < threshold:
            to_idx = i - 1  # 第一个低于阈值的帧的上一帧即为off点
            break

    # Step4: 添加padding
    start_idx = max(0, hs_idx - pad_frames)
    end_idx   = min(len(vertical_force) - 1, to_idx + pad_frames)

    return start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val


def detect_cop_anomalies(copx, copy, time, force_vertical,
                         threshold=20.0, jump_threshold=0.03):
    """
    检测 COPx / COPy 的异常帧（仅在stance阶段检测）

    双重检测机制:
      1. 有效性检测 — 垂直力 < threshold 时，COP应被置零，非零则标记异常
      2. 帧间跳变检测 — 相邻帧COP变化超过 jump_threshold (m) 则标记异常

    Parameters
    ----------
    copx           : ndarray (N,) — COP前后方向 (m)
    copy           : ndarray (N,) — COP左右方向 (m)
    time           : ndarray (N,) — 时间序列 (s)
    force_vertical : ndarray (N,) — 垂直力 (N)
    threshold      : float        — 力阈值，低于此值COP应归零 (N, 默认20)
    jump_threshold : float        — 帧间跳变阈值 (m, 默认0.03 即30mm)

    Returns
    -------
    copx_anomaly : ndarray (N,) bool — COPx异常帧标记
    copy_anomaly : ndarray (N,) bool — COPy异常帧标记
    info         : dict — 统计信息
    """
    n = len(copx)
    copx_anomaly = np.zeros(n, dtype=bool)
    copy_anomaly = np.zeros(n, dtype=bool)

    # 检测1: 力低于阈值但COP非零 → 异常
    invalid_force = np.abs(force_vertical) < threshold
    copx_nonzero = np.abs(copx) > 1e-6
    copy_nonzero = np.abs(copy) > 1e-6
    copx_anomaly |= (invalid_force & copx_nonzero)
    copy_anomaly |= (invalid_force & copy_nonzero)

    # 检测2: 帧间跳变异常（只在stance阶段内检测）
    stance_mask = np.abs(force_vertical) >= threshold
    if np.sum(stance_mask) > 1:
        stance_indices = np.where(stance_mask)[0]
        for arr, anomaly in [(copx, copx_anomaly), (copy, copy_anomaly)]:
            diffs = np.abs(np.diff(arr[stance_mask]))
            jump_frames = diffs > jump_threshold
            # 前后帧都标记
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
    """
    基于中间部分斜率纠正COPx/COPy异常值

    算法（参考 Data_ProcessFunction.process_cop_outliers_slope）:
      1. 识别stance阶段（力 >= threshold）
      2. 取stance中间部分（如30%~70%），线性拟合得到通用斜率 k
      3. 从中间向两端遍历，若帧间变化率 > rate_multiplier * |k*dt|，视为异常
      4. 异常帧用 上一正常值 + k*dt 替换（保持线性趋势连续性）
      5. 非stance阶段（力 < threshold）的COP置零

    Parameters
    ----------
    copx           : ndarray (N,) — COP前后方向 (m)
    copy           : ndarray (N,) — COP左右方向 (m)
    time           : ndarray (N,) — 时间序列 (s)
    force_vertical : ndarray (N,) — 垂直力 (N)
    threshold      : float        — 力阈值 (N, 默认20)
    middle_ratio   : float        — 中间部分占比 (默认0.3)
    rate_multiplier: float        — 异常变化率倍数 (默认2.0)

    Returns
    -------
    copx_corrected : ndarray (N,) — 纠正后的COPx
    copy_corrected : ndarray (N,) — 纠正后的COPy
    info           : dict — 纠正统计信息
    """
    copx_corrected = copx.copy()
    copy_corrected = copy.copy()

    # 非stance阶段COP置零
    invalid_mask = np.abs(force_vertical) < threshold
    copx_corrected[invalid_mask] = 0.0
    copy_corrected[invalid_mask] = 0.0

    # stance阶段
    stance_mask = ~invalid_mask
    stance_indices = np.where(stance_mask)[0]
    n_stance = len(stance_indices)

    info = {
        'copx_slope': 0.0, 'copy_slope': 0.0,
        'copx_outlier_count': 0, 'copy_outlier_count': 0,
    }

    if n_stance < 10:
        return copx_corrected, copy_corrected, info

    # 中间部分范围
    margin = int(n_stance * middle_ratio / 2)
    mid_start = margin
    mid_end = n_stance - margin
    if mid_end - mid_start < 5:
        return copx_corrected, copy_corrected, info

    mid_indices = np.arange(mid_start, mid_end)
    t_stance = time[stance_mask]

    for col_name, data in [('copx', copx_corrected), ('copy', copy_corrected)]:
        stance_vals = data[stance_mask].copy()

        # 用中间部分拟合线性趋势
        mid_vals = stance_vals[mid_indices]
        k, b = np.polyfit(t_stance[mid_indices], mid_vals, 1)

        # 帧间时间间隔
        dt = np.median(np.diff(t_stance)) if len(t_stance) > 1 else 1.0 / 200
        expected_delta = abs(k * dt)
        outlier_delta = expected_delta * rate_multiplier

        corrected = stance_vals.copy()

        # 向左遍历
        for i in range(mid_start - 1, -1, -1):
            if abs(corrected[i] - corrected[i + 1]) > outlier_delta:
                corrected[i] = corrected[i + 1] - k * dt
                info[f'{col_name}_outlier_count'] += 1

        # 向右遍历
        for i in range(mid_end, n_stance):
            if abs(corrected[i] - corrected[i - 1]) > outlier_delta:
                corrected[i] = corrected[i - 1] + k * dt
                info[f'{col_name}_outlier_count'] += 1

        data[stance_mask] = corrected
        info[f'{col_name}_slope'] = k

    return copx_corrected, copy_corrected, info


# ═══════════════════════════════════════════════════════════════════════════════
#  File Writers — .trc  (Markers)
# ═══════════════════════════════════════════════════════════════════════════════

def write_trc(filepath, marker_data, marker_labels, frame_rate,
              units='mm', orig_start_frame=1):
    """
    Write marker data to OpenSim .trc file.

    Parameters
    ----------
    filepath       : str
    marker_data    : ndarray (n_frames, n_markers * 3)  — X Y Z interleaved
    marker_labels  : list[str]       — marker names
    frame_rate     : float           — Hz
    units          : str             — 'mm' or 'm'
    orig_start_frame : int           — original first frame number
    """
    n_frames  = marker_data.shape[0]
    n_markers = len(marker_labels)
    filename  = os.path.basename(filepath)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    with open(filepath, 'w', newline='') as f:
        # ── Line 1: file type identifier ──
        f.write(f"PathFileType\t4\t(X/Y/Z)\t{filename}\n")

        # ── Line 2: header field names ──
        f.write("DataRate\tCameraRate\tNumFrames\tNumMarkers\tUnits\t"
                "OrigDataRate\tOrigDataStartFrame\tOrigNumFrames\n")

        # ── Line 3: header field values ──
        f.write(f"{frame_rate:.0f}\t{frame_rate:.0f}\t{n_frames}\t{n_markers}\t"
                f"{units}\t{frame_rate:.0f}\t{orig_start_frame}\t{n_frames}\n")

        # ── Line 4: marker name row (each name spans 3 columns) ──
        name_parts = ["Frame#", "Time"]
        for lbl in marker_labels:
            name_parts.extend([lbl, "", ""])
        f.write("\t".join(name_parts) + "\n")

        # ── Line 5: axis sub-headers ──
        sub_parts = ["", ""]
        for i in range(n_markers):
            idx = i + 1
            sub_parts.extend([f"X{idx}", f"Y{idx}", f"Z{idx}"])
        f.write("\t".join(sub_parts) + "\n")

        # ── Line 6: blank line ──
        f.write("\n")

        # ── Data rows ──
        for i in range(n_frames):
            frame_num = i + 1
            time_val  = i / frame_rate
            row = [f"{frame_num}", f"{time_val:.6f}"]
            for j in range(marker_data.shape[1]):
                row.append(f"{marker_data[i, j]:.6f}")
            f.write("\t".join(row) + "\n")

    print(f"  [OK] TRC -> {filepath}  ({n_frames} frames, {n_markers} markers)")


# ═══════════════════════════════════════════════════════════════════════════════
#  File Writers — .mot  (Ground Reaction Forces)
# ═══════════════════════════════════════════════════════════════════════════════

def write_mot(filepath, force_data_per_plate, n_plates, frame_rate,
              filename=None):
    """
    Write force-plate data to OpenSim .mot file.

    Column layout (matching standard OpenSim template):
      time
      [FP1 force 3] [FP1 COP 3]  [FP2 force 3] [FP2 COP 3] …
      [FP1 torque 3]            [FP2 torque 3] …

    Parameters
    ----------
    filepath             : str
    force_data_per_plate : list[dict]
        Each dict has keys 'force' (n, 3), 'cop' (n, 3), 'torque' (n, 3)
    n_plates             : int
    frame_rate           : float — Hz
    filename             : str or None — header filename (default: basename)
    """
    n_frames = force_data_per_plate[0]['force'].shape[0]
    fname    = filename or os.path.basename(filepath)

    os.makedirs(os.path.dirname(filepath), exist_ok=True)

    # ── Build column labels ──
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
        # ── Header ──
        f.write(f"{fname}\n")
        f.write("version=1\n")
        f.write(f"nRows={n_frames}\n")
        f.write(f"nColumns={n_cols}\n")
        f.write("inDegrees=yes\n")
        f.write("endheader\n")

        # ── Column labels ──
        f.write("\t".join(labels) + "\n")

        # ── Data rows ──
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

    print(f"  [OK] MOT -> {filepath}  ({n_frames} frames, {n_plates} plates, {n_cols} columns)")
