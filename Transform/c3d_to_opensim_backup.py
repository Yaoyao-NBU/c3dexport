"""
C3D → OpenSim Converter
========================
Complete pipeline to convert C3D biomechanics files to OpenSim-compatible
.trc (marker) and .mot (ground reaction force) files.

Processing Steps:
  1. Read C3D file with ezc3d
  2. Compute Kistler Type 3 → Type 2 force-plate data
  3. Transform force data: plate-local → lab global
  4. Transform all data:   lab global  → OpenSim (Y-up)
  5. Restructure force arrays for OpenSim format
  6. Unit conversion (COP: mm→m, torque: N·mm→N·m)
  7. Filter (Butterworth) & resample force data to marker rate
  8. Detect and extract stance phase
  9. Write .trc and .mot output files

Usage:
  python c3d_to_opensim.py <c3d_file> [output_dir]
"""

import os
import sys
import numpy as np
import ezc3d
import pandas as pd

from transform_utils import (
    rotation_matrix, apply_rotation,
    compute_kistler_channel8,
    compute_kistler_channel6,
    plate_local_to_lab,
    lab_to_opensim_force,
    butter_lowpass_filter,
    resample_to_target_rate,
    detect_stance_phase,
    write_trc,
    write_mot,
)

#去除前缀的函数
def _strip_label_prefix(label):
    """Remove any prefix before a colon (e.g. 'Trial:R.ASIS' → 'R.ASIS')."""
    return label.split(':')[-1] if ':' in label else label


def process_c3d(c3d_path, output_dir,
                marker_cutoff=6.0,
                force_cutoff=50.0,
                stance_threshold=30.0,
                stance_pad_frames=25):
    """
    Full pipeline: C3D → OpenSim .trc + .mot

    Parameters
    ----------
    c3d_path         : str   — input C3D file
    output_dir       : str   — output directory
    marker_cutoff    : float — low-pass cutoff for markers (Hz, default 6)
    force_cutoff     : float — low-pass cutoff for forces  (Hz, default 50)
    stance_threshold : float — vertical-force threshold for stance (N, default 30)
    stance_pad_frames: int   — padding frames before / after stance (default 25)

    Returns
    -------
    dict with output file paths and stance information
    """
    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print(f"  C3D -> OpenSim Converter")
    print(f"  Input : {c3d_path}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 1 — Read C3D file
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 1] 读取 C3D 文件 ...")

    c = ezc3d.c3d(c3d_path)

    point_rate  = float(c['header']['points']['frame_rate'])
    analog_rate = float(c['header']['analogs']['frame_rate'])
    first_frame = c['header']['points']['first_frame']

    # Marker data: (4, n_markers, n_point_frames)
    pt_data       = c['data']['points']
    raw_labels    = c['parameters']['POINT']['LABELS']['value']
    marker_labels = [_strip_label_prefix(lbl) for lbl in raw_labels]

    n_markers      = len(marker_labels)
    n_point_frames = pt_data.shape[2]

    # Marker units
    try:
        marker_units = c['parameters']['POINT']['UNITS']['value'][0]
    except Exception:
        marker_units = 'mm'

    print(f"  Markers : {n_markers} markers,  {n_point_frames} frames @ {point_rate} Hz,  units = {marker_units}")

    # Analog (force platform) data: (1, n_channels, n_analog_frames)
    analog_data    = c['data']['analogs'][0]        # (n_channels, n_analog_frames)
    n_analog_frames = analog_data.shape[1]

    # Force-platform parameters
    fp_params = c['parameters']['FORCE_PLATFORM']
    n_plates  = int(fp_params['USED']['value'][0])
    channels  = fp_params['CHANNEL']['value']
    origin    = fp_params['ORIGIN']['value']
    corners   = np.array(fp_params['CORNERS']['value'])  # (3, 4, n_plates)

    print(f"  Forces  : {n_plates} plate(s),  {n_analog_frames} frames @ {analog_rate} Hz")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 2 — Kistler Type 3 → Type 2
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 2] Kistler Type 3 -> Type 2 力台计算 ...")

    plate_type2 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch8    = analog_data[ch_idx, :]

        b   = float(origin[0, pi])
        a   = float(origin[1, pi])
        az0 = float(origin[2, pi])

        t2 = compute_kistler_channel8(ch8, a, b, az0)
        plate_type2.append(t2)
        print(f"  FP{pi+1}: a={a:.0f} mm, b={b:.0f} mm, az0={az0:.0f} mm")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 3 — Plate-local → Lab global
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 3] 力台本地坐标系 -> 实验室坐标系 ...")

    plate_lab = []
    for pi in range(n_plates):
        plate_corners = corners[:, :, pi]          # (3, 4)
        lab_data = plate_local_to_lab(plate_type2[pi], plate_corners)
        plate_lab.append(lab_data)

        cx = np.mean(plate_corners[0, :])
        cy = np.mean(plate_corners[1, :])
        print(f"  FP{pi+1}: plate centre = ({cx:.1f}, {cy:.1f}) mm")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 4 — Lab global → OpenSim (Y-up)
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 4] 实验室坐标系 -> OpenSim 坐标系 (Y-up) ...")

    # Rotation matrix: Lab(X-fwd, Y-left, Z-up) → OpenSim(X-fwd, Y-up, Z-right)
    # Equivalent to rotation about X-axis by −90°
    R_lab2os = rotation_matrix('X', -90)
    print(f"  旋转矩阵 (绕 X 轴 -90 度):")
    for row in R_lab2os:
        print(f"    [{row[0]:6.3f}  {row[1]:6.3f}  {row[2]:6.3f}]")

    # 4a — Rotate markers
    markers_xyz = pt_data[:3, :, :]                    # (3, n_markers, n_frames)
    markers_os  = np.zeros_like(markers_xyz)
    for mi in range(n_markers):
        markers_os[:, mi, :] = apply_rotation(markers_xyz[:, mi, :], R_lab2os)
    print(f"  Marker 旋转完成: {n_markers} markers")

    # 4b — Rotate / remap force data
    plate_os = []
    for pi in range(n_plates):
        os_data = lab_to_opensim_force(plate_lab[pi])
        plate_os.append(os_data)
    print(f"  力台旋转完成: {n_plates} plate(s)")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 5 + 6 — Restructure force arrays & unit conversion
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 5-6] 力台数据重构 & 单位转换 (COP->m, Torque->N*m) ...")

    # Build per-plate structured arrays
    # Each plate → dict{'force': (N,3), 'cop': (N,3), 'torque': (N,3)}
    #
    # In OpenSim (Y-up):
    #   ground_force_vy = vertical force
    #   ground_force_py = 0  (COP on ground plane)
    #   ground_torque_y = free vertical moment
    #   ground_torque_x = 0,  ground_torque_z = 0

    structured_plates = []
    for pi in range(n_plates):
        d = plate_os[pi]
        n = len(d['Fx'])

        force  = np.column_stack([d['Fx'], d['Fy'], d['Fz']])                 # (N, 3)  [N]
        cop    = np.column_stack([d['COPx'] / 1000.0,                          # mm → m
                                  d['COPy'] / 1000.0,
                                  d['COPz'] / 1000.0])                         # (N, 3)  [m]
        torque = np.column_stack([np.zeros(n),                                 # Tx = 0
                                  d['Tz'] / 1000.0,                            # N·mm → N·m (vertical)
                                  np.zeros(n)])                                 # Tz = 0

        structured_plates.append(dict(force=force, cop=cop, torque=torque))
        print(f"  FP{pi+1}: {n} frames - COP [m], Torque [N*m]")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 7 — Filter & resample
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 7] 滤波 & 重采样 ...")

    # 7a — Filter force data at analog_rate
    print(f"  力台滤波: 截止 {force_cutoff} Hz @ {analog_rate} Hz ...")
    for pi in range(n_plates):
        for key in ('force', 'cop', 'torque'):
            arr = structured_plates[pi][key]
            for col in range(arr.shape[1]):
                # Skip columns that are all-zero (e.g. torque_x, torque_z, cop_y)
                if np.any(arr[:, col] != 0):
                    arr[:, col] = butter_lowpass_filter(arr[:, col], force_cutoff, analog_rate)

    # 7b — Resample force data: analog_rate → point_rate
    if analog_rate != point_rate:
        ratio = analog_rate / point_rate
        print(f"  力台重采样: {analog_rate} Hz -> {point_rate} Hz (ratio={ratio:.1f}) ...")
        for pi in range(n_plates):
            for key in ('force', 'cop', 'torque'):
                arr = structured_plates[pi][key]
                n_target = n_point_frames
                resampled = np.zeros((n_target, arr.shape[1]))
                for col in range(arr.shape[1]):
                    rs = resample_to_target_rate(arr[:, col], analog_rate, point_rate)
                    resampled[:, col] = rs[:n_target]
                structured_plates[pi][key] = resampled

    # 7c — Filter markers at point_rate
    print(f"  Marker 滤波: 截止 {marker_cutoff} Hz @ {point_rate} Hz ...")
    for mi in range(n_markers):
        for axis in range(3):
            markers_os[axis, mi, :] = butter_lowpass_filter(
                markers_os[axis, mi, :], marker_cutoff, point_rate
            )

    print(f"  [OK] 滤波完成")

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 8 — Stance phase detection & extraction
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 8] Stance 阶段检测 & 截取 ...")

    # Auto-detect which plate has the strongest vertical contact
    # Vertical force = Fy (OpenSim Y-up) = column index 1
    max_fy_per_plate = [np.max(np.abs(structured_plates[pi]['force'][:, 1]))
                        for pi in range(n_plates)]
    active_plate = int(np.argmax(max_fy_per_plate))
    print(f"  Plate max |Fy|: {['FP{}: {:.1f} N'.format(i+1, v) for i, v in enumerate(max_fy_per_plate)]}")
    print(f"  Active plate for stance detection: FP{active_plate+1}")

    vertical_force = structured_plates[active_plate]['force'][:, 1]

    try:
        start_idx, end_idx, hs_idx, to_idx = detect_stance_phase(
            vertical_force,
            threshold=stance_threshold,
            pad_frames=stance_pad_frames,
        )
        n_cut = end_idx - start_idx + 1
        print(f"  Heel Strike : frame {hs_idx+1}  (FP{active_plate+1})")
        print(f"  Toe Off     : frame {to_idx+1}  (FP{active_plate+1})")
        print(f"  Cut range (+/-{stance_pad_frames} frames): frame {start_idx+1} - {end_idx+1}  ({n_cut} frames)")

        # Cut force
        for pi in range(n_plates):
            for key in ('force', 'cop', 'torque'):
                structured_plates[pi][key] = structured_plates[pi][key][start_idx:end_idx+1, :]

        #补偿帧归零算法
        # 计算在截取后的新数组中的相对位置
        rel_hs = hs_idx - start_idx
        rel_to = to_idx - start_idx

        for pi in range(n_plates):
            # 1. 归零 COP (ax, ay) 
            # 在 0 到 rel_hs 帧，以及 rel_to 帧之后的数据设为 0
            structured_plates[pi]['cop'][:rel_hs, :] = 0
            structured_plates[pi]['cop'][rel_to+1:, :] = 0

            # 2. 归零 自由力矩 (Tz) 
            structured_plates[pi]['torque'][:rel_hs, :] = 0
            structured_plates[pi]['torque'][rel_to+1:, :] = 0

            # 3. 归零 水平力 (Fx, Fy) 
            # 假设 force 数组的列顺序是 [Fx, Fy, Fz]
            # 我们只对第 0 列 (Fx) 和 第 1 列 (Fy) 进行归零，保留垂直力 Fz 的原始噪声或趋势
            structured_plates[pi]['force'][:rel_hs, :] = 0
            structured_plates[pi]['force'][rel_to+1:, :] = 0
            
            # 可选：如果你希望 Fz 在非接触期间也绝对为 0
            # structured_plates[pi]['force'][:rel_hs, 2] = 0
            # structured_plates[pi]['force'][rel_to+1:, 2] = 0
        
        # Cut markers
        markers_os = markers_os[:, :, start_idx:end_idx+1]

        stance_info = dict(
            trial            = trial,
            heel_strike_frame = hs_idx + 1,
            toe_off_frame     = to_idx + 1,
            cut_start_frame   = start_idx + 1,
            cut_end_frame     = end_idx + 1,
            total_frames      = n_cut,
        )
    except ValueError as e:
        print(f"  [WARN] {e}")
        print(f"  -> 使用全部数据（不截取）")
        n_cut = structured_plates[0]['force'].shape[0]
        stance_info = dict(
            trial            = trial,
            heel_strike_frame = 'N/A',
            toe_off_frame     = 'N/A',
            cut_start_frame   = 1,
            cut_end_frame     = n_cut,
            total_frames      = n_cut,
        )

    # ══════════════════════════════════════════════════════════════════════════
    #  Step 9 — Write output files
    # ══════════════════════════════════════════════════════════════════════════
    print("\n[Step 9] 写入 OpenSim 文件 ...")

    n_out_frames = markers_os.shape[2]

    # 9a — Reshape markers: (3, n_markers, n_frames) → (n_frames, n_markers*3)
    marker_flat = np.zeros((n_out_frames, n_markers * 3))
    for mi in range(n_markers):
        marker_flat[:, mi * 3]     = markers_os[0, mi, :]   # X
        marker_flat[:, mi * 3 + 1] = markers_os[1, mi, :]   # Y
        marker_flat[:, mi * 3 + 2] = markers_os[2, mi, :]   # Z

    # 9b — Write .trc
    trc_path = os.path.join(output_dir, f"{trial}.trc")
    write_trc(
        trc_path, marker_flat, marker_labels, point_rate,
        units=marker_units, orig_start_frame=1,
    )

    # 9c — Write .mot
    mot_path = os.path.join(output_dir, f"{trial}.mot")
    write_mot(
        mot_path, structured_plates, n_plates, point_rate,
        filename=f"{trial}.mot",
    )

    # 9d — Write stance-info CSV
    csv_path = os.path.join(c3d_path, f"cut_records.csv")
    file_exists = os.path.isfile(csv_path)
    pd.DataFrame([stance_info]).to_csv(
    csv_path, 
    mode='a', 
    index=False, 
    header=not file_exists,
    encoding='utf-8-sig' # 解决中文乱码逻辑（可选）
   )
    print(f"  [OK] 截取记录 -> {csv_path}")

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 60)
    print("  [DONE] C3D -> OpenSim Done!")
    print(f"  .trc : {trc_path}")
    print(f"  .mot : {mot_path}")
    print(f"  记录 : {csv_path}")
    print("=" * 60)

    return dict(
        trc_path   = trc_path,
        mot_path   = mot_path,
        csv_path   = csv_path,
        stance_info = stance_info,
    )


# ── CLI entry point ──────────────────────────────────────────────────────────

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python c3d_to_opensim.py <c3d_file> [output_dir]")
        sys.exit(1)

    c3d_file = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(c3d_file), 'output'
    )

    process_c3d(c3d_file, out_dir)
