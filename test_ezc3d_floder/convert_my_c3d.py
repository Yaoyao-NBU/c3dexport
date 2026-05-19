"""
C3D -> OpenSim 转换脚本 (Type 3 Kistler, 自定义力台方向)
=========================================================
直接调用 C3DExport 库的底层函数完成转换。
力台局部坐标系 Y 轴与实验室 X 轴方向相同（绕 Z 轴旋转 90°）。
COP 值做 6 Hz 低通滤波去除异常跳变。
不做 stance 阶段检测，输出全段数据。

Usage:
  python convert_my_c3d.py
"""

import os
import sys
import glob
import numpy as np
import ezc3d

# 将项目根目录加入搜索路径，以便导入 C3DExport 库
sys.path.insert(0, r"E:\Python_Learn\c3dexport")

input_file  = r"G:\houtitui_project\ceshi\stic.c3d"
output_file = r"G:\houtitui_project\ceshi\output"
# ── 从 C3DExport 库导入底层函数 ──────────────────────────────────────────
from C3DExport.utils import (
    rotation_matrix,
    apply_rotation,
    compute_forceplate_type3,
    plate_local_to_lab,
    lab_to_opensim_force,
    butter_lowpass_filter,
    resample_to_target_rate,
)
from C3DExport.io import write_trc, write_mot


# ============================================================================
#  Local Functions
# ============================================================================

def filter_cop(cop_data, cutoff, fs):
    """对 COP 数据进行低通滤波，去除异常跳变值。

    Parameters
    ----------
    cop_data : ndarray (3, N) -- COPx, COPy, COPz
    cutoff : float -- 截止频率 (Hz)
    fs : float -- 采样率 (Hz)

    Returns
    -------
    filtered : ndarray (3, N)
    """
    filtered = np.zeros_like(cop_data)
    for axis in range(3):
        filtered[axis, :] = butter_lowpass_filter(cop_data[axis, :], cutoff, fs)
    return filtered


def strip_label_prefix(label):
    """去除标签前缀，如 'Trial:R.ASIS' -> 'R.ASIS'。"""
    return label.split(':')[-1] if ':' in label else label


# ============================================================================
#  Main Conversion Pipeline
# ============================================================================

def convert_c3d(c3d_path, output_dir):
    """将 Type 3 Kistler C3D 文件转换为 OpenSim .trc + .mot。

    力台局部坐标系 Y 轴与实验室 X 轴方向相同（绕 Z 轴旋转 90°）。

    Parameters
    ----------
    c3d_path : str -- 输入 C3D 文件路径
    output_dir : str -- 输出目录
    """
    MARKER_CUTOFF = 6.0     # Hz, marker 低通截止频率
    FORCE_CUTOFF  = 50.0    # Hz, 力信号低通截止频率
    COP_CUTOFF    = 6.0     # Hz, COP 低通截止频率

    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print(f"  C3D -> OpenSim (Type 3 Kistler, 自定义力台方向)")
    print(f"  Input : {c3d_path}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # ─────────────────────────────────────────────────────────────────────
    #  Step 1: 读取 C3D 文件
    # ─────────────────────────────────────────────────────────────────────
    print("\n[Step 1] 读取 C3D 文件 ...")

    c = ezc3d.c3d(c3d_path)

    point_rate  = float(c['header']['points']['frame_rate'])
    analog_rate = float(c['header']['analogs']['frame_rate'])

    pt_data       = c['data']['points']
    raw_labels    = c['parameters']['POINT']['LABELS']['value']
    marker_labels = [strip_label_prefix(lbl) for lbl in raw_labels]
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

    print(f"  Markers : {n_markers} markers, {n_point_frames} frames @ {point_rate} Hz")
    print(f"  Forces  : {n_plates} plate(s), {n_analog_frames} frames @ {analog_rate} Hz")
    print(f"  FP types: {fp_types}")

    # ─────────────────────────────────────────────────────────────────────
    #  Step 2: 模拟信号滤波（在原始采样率下）
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 2] 模拟信号滤波: {FORCE_CUTOFF} Hz @ {analog_rate} Hz ...")

    for ch in range(analog_data.shape[0]):
        if np.any(analog_data[ch, :] != 0):
            analog_data[ch, :] = butter_lowpass_filter(
                analog_data[ch, :], FORCE_CUTOFF, analog_rate
            )

    # ─────────────────────────────────────────────────────────────────────
    #  Step 3: 重采样到 marker 采样率
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 3] 重采样: {analog_rate} Hz -> {point_rate} Hz ...")

    if analog_rate != point_rate:
        n_target = n_point_frames
        resampled = np.zeros((analog_data.shape[0], n_target))
        for ch in range(analog_data.shape[0]):
            rs = resample_to_target_rate(analog_data[ch, :], analog_rate, point_rate)
            resampled[ch, :] = rs[:n_target]
        analog_data = resampled

    # ─────────────────────────────────────────────────────────────────────
    #  Step 4: Kistler Type 3 力台计算
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 4] Kistler Type 3 力台计算 ...")

    plate_local = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch8    = analog_data[ch_idx, :]
        a   = float(origin[0, pi])
        b   = float(origin[1, pi])
        az0 = float(origin[2, pi])
        plate_local.append(compute_forceplate_type3(ch8, a, b, az0))
        print(f"  FP{pi+1}: a={a:.1f} mm, b={b:.1f} mm, az0={az0:.1f} mm")

    # ─────────────────────────────────────────────────────────────────────
    #  Step 5: COP 低通滤波 (6 Hz，去除异常跳变)
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 5] COP 低通滤波: {COP_CUTOFF} Hz ...")

    for pi in range(n_plates):
        cop_raw = np.vstack([
            plate_local[pi]['ax'],
            plate_local[pi]['ay'],
            np.zeros_like(plate_local[pi]['ax']),
        ])
        cop_filtered = filter_cop(cop_raw, COP_CUTOFF, point_rate)
        plate_local[pi]['ax'] = cop_filtered[0, :]
        plate_local[pi]['ay'] = cop_filtered[1, :]

    # ─────────────────────────────────────────────────────────────────────
    #  Step 6: 力台坐标系 -> 实验室坐标系
    #  力台局部 Y 轴与实验室 X 轴方向相同（绕 Z 轴旋转 90°）
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 6] 力台坐标系 -> 实验室坐标系 (Y→X, 绕Z旋转90°) ...")

    # 从力台角点推导力台在实验室中的朝向
    R_plate_standard, _ = plate_local_to_lab.__wrapped__ if hasattr(plate_local_to_lab, '__wrapped__') else (None, None)

    # 力台局部 Y 轴对齐实验室 X 轴的旋转矩阵
    # Lab_X = Plate_Y, Lab_Y = Plate_X
    R_fp = np.array([[0, 1, 0],
                     [1, 0, 0],
                     [0, 0, 1]], dtype=float)

    plate_lab = []
    for pi in range(n_plates):
        P_center = np.mean(corners[:, :, pi], axis=1)  # (3,) 力台中心

        # 力: 旋转到实验室坐标系
        F_plate = np.vstack([
            plate_local[pi]['Fx'],
            plate_local[pi]['Fy'],
            plate_local[pi]['Fz'],
        ])
        F_lab = R_fp @ F_plate

        # COP: 旋转偏移量 + 平移到实验室坐标系
        cop_offset = np.vstack([
            plate_local[pi]['ax'],
            plate_local[pi]['ay'],
            np.zeros_like(plate_local[pi]['ax']),
        ])
        cop_lab = R_fp @ cop_offset + P_center[:, np.newaxis]

        # 自由力矩
        Tz_plate = np.vstack([
            np.zeros_like(plate_local[pi]['Tz']),
            np.zeros_like(plate_local[pi]['Tz']),
            plate_local[pi]['Tz'],
        ])
        Tz_lab = (R_fp @ Tz_plate)[2, :]

        plate_lab.append(dict(
            Fx=F_lab[0], Fy=F_lab[1], Fz=F_lab[2],
            COPx=cop_lab[0], COPy=cop_lab[1], COPz=cop_lab[2],
            Tz=Tz_lab,
        ))

    # ─────────────────────────────────────────────────────────────────────
    #  Step 7: 实验室坐标系 -> OpenSim 坐标系 (Y-up)
    #  Lab(X-forward, Y-left, Z-up) -> OpenSim(X-forward, Y-up, Z-right)
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 7] 实验室坐标系 -> OpenSim 坐标系 (Y-up) ...")

    R_lab2os = rotation_matrix('X', -90)

    # 力台数据: lab -> OpenSim
    plate_os = []
    for pi in range(n_plates):
        plate_os.append(lab_to_opensim_force(plate_lab[pi]))

    # Marker 数据: lab -> OpenSim
    markers_xyz = pt_data[:3, :, :]  # (3, n_markers, n_frames)
    markers_os  = np.zeros_like(markers_xyz)
    for mi in range(n_markers):
        markers_os[:, mi, :] = apply_rotation(markers_xyz[:, mi, :], R_lab2os)

    # ─────────────────────────────────────────────────────────────────────
    #  Step 8: Marker 滤波
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 8] Marker 滤波: {MARKER_CUTOFF} Hz ...")

    for mi in range(n_markers):
        for axis in range(3):
            markers_os[axis, mi, :] = butter_lowpass_filter(
                markers_os[axis, mi, :], MARKER_CUTOFF, point_rate
            )

    # ─────────────────────────────────────────────────────────────────────
    #  Step 9: 力台数据结构整理 + 单位转换
    #  COP: mm -> m,  Tz: N·mm -> N·m
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 9] 力台数据整理 & 单位转换 ...")

    structured_plates = []
    for pi in range(n_plates):
        d = plate_os[pi]
        n = len(d['Fx'])
        valid = np.abs(d['Fy']) >= 20.0  # COP 阈值

        structured_plates.append(dict(
            Fx=d['Fx'], Fy=d['Fy'], Fz=d['Fz'],
            COPx=np.where(valid, d['COPx'], 0.0) / 1000.0,
            COPy=np.where(valid, d['COPy'], 0.0),
            COPz=np.where(valid, d['COPz'], 0.0) / 1000.0,
            Tz=d['Tz'] / 1000.0,
        ))

    # ─────────────────────────────────────────────────────────────────────
    #  Step 10: 写入输出文件
    # ─────────────────────────────────────────────────────────────────────
    print(f"\n[Step 10] 写入 OpenSim 文件 ...")

    # Marker: (3, n_markers, n_frames) -> (n_frames, n_markers*3)
    n_out = markers_os.shape[2]
    marker_flat = np.zeros((n_out, n_markers * 3))
    for mi in range(n_markers):
        marker_flat[:, mi * 3]     = markers_os[0, mi, :]
        marker_flat[:, mi * 3 + 1] = markers_os[1, mi, :]
        marker_flat[:, mi * 3 + 2] = markers_os[2, mi, :]

    # .trc
    trc_path = os.path.join(output_dir, f"{trial}.trc")
    write_trc(trc_path, marker_flat, marker_labels, point_rate,
              units=marker_units, orig_start_frame=1)

    # .mot
    mot_plates = []
    for pi in range(n_plates):
        d = structured_plates[pi]
        n = len(d['Fx'])
        force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])
        cop    = np.column_stack([d['COPx'], d['COPy'], d['COPz']])
        torque = np.column_stack([np.zeros(n), d['Tz'], np.zeros(n)])
        mot_plates.append(dict(force=force, cop=cop, torque=torque))

    mot_path = os.path.join(output_dir, f"{trial}.mot")
    write_mot(mot_path, mot_plates, n_plates, point_rate,
              filename=f"{trial}.mot")

    print(f"\n{'=' * 60}")
    print(f"  [DONE] .trc: {trc_path}")
    print(f"  [DONE] .mot: {mot_path}")
    print(f"{'=' * 60}")

    return trc_path, mot_path


# ============================================================================
#  Entry Point
# ============================================================================

if __name__ == '__main__':
    os.makedirs(output_file, exist_ok=True)

    if os.path.isfile(input_file):
        c3d_files = [input_file]
    elif os.path.isdir(input_file):
        c3d_files = sorted(set(
            glob.glob(os.path.join(input_file, '*.c3d'))
            + glob.glob(os.path.join(input_file, '*.C3D'))
        ))
    else:
        print(f"[ERROR] 路径不存在: {input_file}")
        sys.exit(1)

    if not c3d_files:
        print(f"[ERROR] 在 {input_file} 中未找到 C3D 文件")
        sys.exit(1)

    print(f"找到 {len(c3d_files)} 个 C3D 文件\n")

    for c3d_path in c3d_files:
        try:
            convert_c3d(c3d_path, output_file)
        except Exception as e:
            print(f"\n[FAIL] {os.path.basename(c3d_path)}: {e}")
            import traceback
            traceback.print_exc()
