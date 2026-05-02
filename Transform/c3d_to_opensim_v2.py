"""
C3D → OpenSim Converter V2 (Peak-based Stance Detection)
=========================================================
基于峰值检测stance的C3D转换工具，改进点：
  - 从垂直力峰值向两侧遍历检测 load/off，而非简单阈值
  - padding 帧力数据手动归零
  - 集成 COPx/COPy 异常检测与斜率纠正
  - 使用绝对路径配置，无需交互输入

Usage:
  直接修改下方配置区域，或:
  python c3d_to_opensim_v2.py <c3d_file> [output_dir]
"""

import os
import sys
import numpy as np
import ezc3d
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from transform_utils import (
    rotation_matrix, apply_rotation,
    compute_kistler_channel8,
    compute_kistler_channel6,
    plate_local_to_lab,
    lab_to_opensim_force,
    butter_lowpass_filter,
    resample_to_target_rate,
    detect_stance_phase_from_peak,
    detect_cop_anomalies,
    correct_cop_slope,
    write_trc,
    write_mot,
)


def _strip_label_prefix(label):
    """去除标记名冒号前缀 (e.g. 'Trial:R.ASIS' → 'R.ASIS')"""
    return label.split(':')[-1] if ':' in label else label


def process_c3d_v2(c3d_path, output_dir,
                   marker_cutoff=6.0,
                   force_cutoff=50.0,
                   stance_threshold=30.0,
                   stance_pad_frames=25,
                   cop_middle_ratio=0.3,
                   cop_rate_multiplier=2.0,
                   cop_jump_threshold=0.03):
    """
    完整处理流程: C3D → OpenSim .trc + .mot
    使用峰值检测stance + COP斜率纠正

    Parameters
    ----------
    c3d_path            : str   — 输入C3D文件路径
    output_dir          : str   — 输出目录
    marker_cutoff       : float — marker低通滤波截止频率 (Hz, 默认6)
    force_cutoff        : float — 力数据低通滤波截止频率 (Hz, 默认50)
    stance_threshold    : float — stance检测力阈值 (N, 默认30)
    stance_pad_frames   : int   — 前后补偿帧数 (默认25)
    cop_middle_ratio    : float — COP斜率纠正中间部分占比 (默认0.3)
    cop_rate_multiplier : float — COP变化率异常倍数 (默认2.0)
    cop_jump_threshold  : float — COP帧间跳变阈值 (m, 默认0.03)

    Returns
    -------
    dict — 输出文件路径及stance信息
    """
    trial = os.path.splitext(os.path.basename(c3d_path))[0]
    os.makedirs(output_dir, exist_ok=True)

    print("=" * 60)
    print(f"  C3D -> OpenSim V2 (Peak Stance + COP Correction)")
    print(f"  Input : {c3d_path}")
    print(f"  Output: {output_dir}")
    print("=" * 60)

    # ════════════════════════════════════════════════════════════════
    #  Step 1 — 读取C3D文件
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 1] 读取 C3D 文件 ...")

    c = ezc3d.c3d(c3d_path)

    point_rate  = float(c['header']['points']['frame_rate'])
    analog_rate = float(c['header']['analogs']['frame_rate'])

    pt_data       = c['data']['points']
    raw_labels    = c['parameters']['POINT']['LABELS']['value']
    marker_labels = [_strip_label_prefix(lbl) for lbl in raw_labels]

    n_markers      = len(marker_labels)
    n_point_frames = pt_data.shape[2]

    try:
        marker_units = c['parameters']['POINT']['UNITS']['value'][0]
    except Exception:
        marker_units = 'mm'

    print(f"  Markers: {n_markers} markers, {n_point_frames} frames @ {point_rate} Hz")

    analog_data     = c['data']['analogs'][0]
    n_analog_frames = analog_data.shape[1]

    fp_params = c['parameters']['FORCE_PLATFORM']
    n_plates  = int(fp_params['USED']['value'][0])
    channels  = fp_params['CHANNEL']['value']
    origin    = fp_params['ORIGIN']['value']
    corners   = np.array(fp_params['CORNERS']['value'])

    print(f"  Forces : {n_plates} plate(s), {n_analog_frames} frames @ {analog_rate} Hz")

    # ════════════════════════════════════════════════════════════════
    #  Step 2 — 重采样模拟数据到marker频率
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 2] 重采样模拟信号 ...")

    if analog_rate != point_rate:
        ratio = analog_rate / point_rate
        print(f"  重采样: {analog_rate} Hz -> {point_rate} Hz")
        n_target_frames = n_point_frames
        resampled_data = np.zeros((analog_data.shape[0], n_target_frames))
        for ch_idx in range(analog_data.shape[0]):
            rs = resample_to_target_rate(analog_data[ch_idx, :], analog_rate, point_rate)
            resampled_data[ch_idx, :] = rs[:n_target_frames]
        analog_data = resampled_data
        analog_rate = point_rate

    print(f"  [OK] 重采样完成")

    # ════════════════════════════════════════════════════════════════
    #  Step 3 — Kistler Type 3 → Type 2
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 3] Kistler Type 3 -> Type 2 ...")

    plate_type2 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch8    = analog_data[ch_idx, :]

        # b   = float(origin[0, pi])
        # a   = float(origin[1, pi])
        a   = float(origin[0, pi])
        b   = float(origin[1, pi])
        az0 = float(origin[2, pi])

        t2 = compute_kistler_channel8(ch8, a, b, az0)
        plate_type2.append(t2)
        print(f"  FP{pi+1}: a={a:.0f}mm, b={b:.0f}mm, az0={az0:.0f}mm")

    # ════════════════════════════════════════════════════════════════
    #  Step 4 — 力台本地 → 实验室坐标系
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 4] Plate-local -> Lab global ...")

    plate_lab = []
    for pi in range(n_plates):
        plate_corners = corners[:, :, pi]
        lab_data = plate_local_to_lab(plate_type2[pi], plate_corners)
        plate_lab.append(lab_data)

    # ════════════════════════════════════════════════════════════════
    #  Step 5 — 实验室 → OpenSim (Y-up)
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 5] Lab -> OpenSim (Y-up) ...")

    R_lab2os = rotation_matrix('X', -90)

    # 旋转marker
    markers_xyz = pt_data[:3, :, :]
    markers_os  = np.zeros_like(markers_xyz)
    for mi in range(n_markers):
        markers_os[:, mi, :] = apply_rotation(markers_xyz[:, mi, :], R_lab2os)

    # 旋转力台
    plate_os = []
    for pi in range(n_plates):
        os_data = lab_to_opensim_force(plate_lab[pi])
        plate_os.append(os_data)

    # ════════════════════════════════════════════════════════════════
    #  Step 6 — Marker滤波
    # ════════════════════════════════════════════════════════════════
    print(f"\n[Step 6] Marker滤波: {marker_cutoff} Hz ...")
    for mi in range(n_markers):
        for axis in range(3):
            markers_os[axis, mi, :] = butter_lowpass_filter(
                markers_os[axis, mi, :], marker_cutoff, point_rate
            )

    # ════════════════════════════════════════════════════════════════
    #  Step 7 — 力阈值处理 + COP单位转换
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 7] COP阈值处理 & 单位转换 ...")

    cop_threshold = 20.0
    structured_plates = []

    for pi in range(n_plates):
        d_os = plate_os[pi]
        n = len(d_os['Fx'])

        valid = np.abs(d_os['Fy']) >= cop_threshold
        COPx = np.where(valid, d_os['COPx'], 0.0)
        COPy = np.where(valid, d_os['COPy'], 0.0)
        COPz = np.where(valid, d_os['COPz'], 0.0)

        structured_plates.append(dict(
            Fx=d_os['Fx'],
            Fy=d_os['Fy'],
            Fz=d_os['Fz'],
            COPx=COPx / 1000.0,   # mm -> m
            COPy=COPy,             # 已经是0 (ground level)
            COPz=COPz / 1000.0,   # mm -> m
            Tz=d_os['Tz'] / 1000.0,  # N·mm -> N·m
        ))
        print(f"  FP{pi+1}: 阈值={cop_threshold}N, 单位转换完成")

    # ════════════════════════════════════════════════════════════════
    #  Step 8 — 峰值检测stance + 截取 + padding归零
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 8] Peak-based Stance检测 & 截取 ...")

    # 选择最强垂直力的力台（踩上去的那个）
    max_fy_per_plate = [np.max(np.abs(structured_plates[pi]['Fy']))
                        for pi in range(n_plates)]
    active_plate = int(np.argmax(max_fy_per_plate))
    print(f"  Active plate: FP{active_plate+1} (max Fy={max_fy_per_plate[active_plate]:.1f}N)")

    vertical_force = structured_plates[active_plate]['Fy']

    try:
        start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val = \
            detect_stance_phase_from_peak(
                vertical_force,
                threshold=stance_threshold,
                pad_frames=stance_pad_frames,
            )

        n_cut = end_idx - start_idx + 1
        print(f"  Peak: {peak_val:.2f}N @ frame {peak_idx+1}")
        print(f"  Load (HS): frame {hs_idx+1}")
        print(f"  Off  (TO): frame {to_idx+1}")
        print(f"  Cut range (+/-{stance_pad_frames} frames): {start_idx+1} - {end_idx+1} ({n_cut} frames)")

        # 截取所有力台数据
        for pi in range(n_plates):
            for key in ('Fx', 'Fy', 'Fz', 'COPx', 'COPy', 'COPz', 'Tz'):
                structured_plates[pi][key] = structured_plates[pi][key][start_idx:end_idx+1]

        # padding帧归零
        rel_hs = hs_idx - start_idx
        rel_to = to_idx - start_idx
        zero_keys = ('Fx', 'Fy', 'Fz', 'COPx', 'COPy', 'COPz', 'Tz')

        for pi in range(n_plates):
            for key in zero_keys:
                # load之前归零
                structured_plates[pi][key][:rel_hs] = 0
                # off之后归零
                structured_plates[pi][key][rel_to+1:] = 0

        print(f"  Padding帧归零完成 (0~{rel_hs} 和 {rel_to+1}~{n_cut-1})")

        # 截取marker
        markers_os = markers_os[:, :, start_idx:end_idx+1]

        # 生成时间列（用于COP纠正）
        time_col = np.arange(n_cut) / point_rate

        # ══════════════════════════════════════════════════════════
        #  Step 9 — COP异常检测与斜率纠正
        # ══════════════════════════════════════════════════════════
        print("\n[Step 9] COP异常检测 & 斜率纠正 ...")

        copx_anomaly, copy_anomaly, anomaly_info = detect_cop_anomalies(
            structured_plates[active_plate]['COPx'],
            structured_plates[active_plate]['COPz'],
            time_col,
            structured_plates[active_plate]['Fy'],
            threshold=cop_threshold,
            jump_threshold=cop_jump_threshold,
        )
        print(f"  FP{active_plate+1} 异常检测: "
              f"COPx={anomaly_info['copx_total_anomaly']}帧, "
              f"COPz={anomaly_info['copy_total_anomaly']}帧")

        # 斜率纠正
        copx_corr, copz_corr, slope_info = correct_cop_slope(
            structured_plates[active_plate]['COPx'],
            structured_plates[active_plate]['COPz'],
            time_col,
            structured_plates[active_plate]['Fy'],
            threshold=cop_threshold,
            middle_ratio=cop_middle_ratio,
            rate_multiplier=cop_rate_multiplier,
        )

        structured_plates[active_plate]['COPx'] = copx_corr
        structured_plates[active_plate]['COPz'] = copz_corr
        print(f"  斜率纠正: COPx纠正{slope_info['copx_outlier_count']}帧 "
              f"(k={slope_info['copx_slope']:.6f}m/s), "
              f"COPz纠正{slope_info['copy_outlier_count']}帧 "
              f"(k={slope_info['copy_slope']:.6f}m/s)")

        stance_info = dict(
            trial=trial,
            heel_strike_frame=hs_idx + 1,
            toe_off_frame=to_idx + 1,
            cut_start_frame=start_idx + 1,
            cut_end_frame=end_idx + 1,
            total_frames=n_cut,
            peak_frame=peak_idx + 1,
            peak_value=round(peak_val, 2),
            copx_slope=round(slope_info['copx_slope'], 6),
            copz_slope=round(slope_info['copy_slope'], 6),
            copx_outlier_count=slope_info['copx_outlier_count'],
            copz_outlier_count=slope_info['copy_outlier_count'],
        )
    except ValueError as e:
        print(f"  [WARN] {e} -> 使用全部数据")
        n_cut = len(structured_plates[0]['Fy'])
        stance_info = dict(
            trial=trial,
            heel_strike_frame='N/A',
            toe_off_frame='N/A',
            cut_start_frame=1,
            cut_end_frame=n_cut,
            total_frames=n_cut,
            peak_frame='N/A',
            peak_value='N/A',
            copx_slope=0,
            copz_slope=0,
            copx_outlier_count=0,
            copz_outlier_count=0,
        )

    # ════════════════════════════════════════════════════════════════
    #  Step 10 — 写入输出文件
    # ════════════════════════════════════════════════════════════════
    print("\n[Step 10] 写入 OpenSim 文件 ...")

    n_out_frames = markers_os.shape[2]

    # Marker: (3, n_markers, n_frames) → (n_frames, n_markers*3)
    marker_flat = np.zeros((n_out_frames, n_markers * 3))
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
        n = len(structured_plates[pi]['Fy'])
        force = np.column_stack([
            structured_plates[pi]['Fx'],
            structured_plates[pi]['Fy'],
            structured_plates[pi]['Fz']
        ])
        cop = np.column_stack([
            structured_plates[pi]['COPx'],
            structured_plates[pi]['COPy'],
            structured_plates[pi]['COPz']
        ])
        torque = np.column_stack([
            np.zeros(n),
            structured_plates[pi]['Tz'],
            np.zeros(n)
        ])
        mot_plates.append(dict(force=force, cop=cop, torque=torque))

    mot_path = os.path.join(output_dir, f"{trial}.mot")
    write_mot(mot_path, mot_plates, n_plates, point_rate,
              filename=f"{trial}.mot")

    # CSV记录
    csv_path = os.path.join(output_dir, "cut_records_v2.csv")
    file_exists = os.path.isfile(csv_path)
    pd.DataFrame([stance_info]).to_csv(
        csv_path, mode='a', index=False,
        header=not file_exists, encoding='utf-8-sig'
    )

    print(f"\n{'=' * 60}")
    print(f"  [DONE] {trial} -> .trc + .mot")
    print(f"  .trc: {trc_path}")
    print(f"  .mot: {mot_path}")
    print(f"{'=' * 60}")

    return dict(
        trc_path=trc_path,
        mot_path=mot_path,
        csv_path=csv_path,
        stance_info=stance_info,
    )


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print("Usage: python c3d_to_opensim_v2.py <c3d_file> [output_dir]")
        sys.exit(1)

    c3d_file = sys.argv[1]
    out_dir  = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(c3d_file), 'output_v2'
    )

    process_c3d_v2(c3d_file, out_dir)
