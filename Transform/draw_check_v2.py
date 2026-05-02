"""
V2数据可视化检查工具 — Peak Stance + COP Correction
=====================================================
批量读取C3D文件，提取GRF/COP数据，绘制对比图用于检查V2处理效果

功能:
  1. 从C3D文件提取2号力台的GRF/COP数据
  2. 使用峰值检测stance + COP斜率纠正
  3. 绘制单文件组合图 (Fy/COPx/COPz)
  4. 绘制原始 vs 纠正后对比图
  5. 按trial叠加折线图

Usage:
  修改下方配置区域，然后运行:
  python draw_check_v2.py
"""

import os
import sys
import numpy as np
import ezc3d
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rcParams
import copy

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from transform_utils import (
    compute_kistler_channel8,
    plate_local_to_lab,
    lab_to_opensim_force,
    butter_lowpass_filter,
    detect_stance_phase_from_peak,
    detect_cop_anomalies,
    correct_cop_slope,
)


# ═══════════════════════════════════════════════════════════════════
#  配置区域
# ═══════════════════════════════════════════════════════════════════

# C3D文件目录（绝对路径）
INPUT_DIR = r"E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Transform\data"

# 输出图像目录（绝对路径）
OUTPUT_DIR = r"E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Transform\output_v2\plots"

# 文件匹配模式
PATTERN = "*.c3d"

# 采样频率 (Hz)
FS = 100.0   # 根据你的C3D文件调整

# Stance检测参数
STANCE_THRESHOLD  = 30.0    # 力阈值 (N)
STANCE_PAD_FRAMES = 25      # 补偿帧数
COP_THRESHOLD     = 20.0    # COP力阈值 (N)

# COP斜率纠正参数
COP_MIDDLE_RATIO    = 0.3   # 中间部分占比
COP_RATE_MULTIPLIER = 2.0   # 异常倍数
COP_JUMP_THRESHOLD  = 0.03  # 帧间跳变阈值 (m)

# 绘图参数
DPI = 300

# ═══════════════════════════════════════════════════════════════════

# 全局样式
rcParams['font.family'] = 'Arial'
rcParams['font.size'] = 10
rcParams['lines.linewidth'] = 1.0
rcParams['figure.dpi'] = DPI
rcParams['savefig.dpi'] = DPI
rcParams['savefig.bbox'] = 'tight'


def extract_plate2_data(c3d_path, plate_idx=1):
    """
    从C3D文件提取指定力台数据，执行峰值stance检测 + COP纠正

    Parameters
    ----------
    c3d_path  : str — C3D文件路径
    plate_idx : int — 力台索引（默认1，即2号力台）

    Returns
    -------
    dict — 包含原始/纠正后的数据及stance信息
    """
    trial_name = os.path.splitext(os.path.basename(c3d_path))[0]
    c = ezc3d.c3d(c3d_path)

    point_rate  = float(c['header']['points']['frame_rate'])
    analog_rate = float(c['header']['analogs']['frame_rate'])
    analog_data = c['data']['analogs'][0]

    fp_params = c['parameters']['FORCE_PLATFORM']
    n_plates  = int(fp_params['USED']['value'][0])
    channels  = fp_params['CHANNEL']['value']
    origin    = fp_params['ORIGIN']['value']
    corners   = np.array(fp_params['CORNERS']['value'])

    # 确保力台存在
    if plate_idx >= n_plates:
        raise ValueError(f"请求力台{plate_idx+1}，但只有{n_plates}个")

    # 重采样到marker频率
    if analog_rate != point_rate:
        n_target = c['data']['points'].shape[2]
        resampled = np.zeros((analog_data.shape[0], n_target))
        from transform_utils import resample_to_target_rate
        for ch in range(analog_data.shape[0]):
            rs = resample_to_target_rate(analog_data[ch, :], analog_rate, point_rate)
            resampled[ch, :] = rs[:n_target]
        analog_data = resampled
        analog_rate = point_rate

    # Type3 → Type2
    ch_idx = channels[:, plate_idx].astype(int) - 1
    ch8 = analog_data[ch_idx, :]
    # b   = float(origin[0, plate_idx])
    # a   = float(origin[1, plate_idx])
    a   = float(origin[0, plate_idx])
    b   = float(origin[1, plate_idx])
    az0 = float(origin[2, plate_idx])
    t2 = compute_kistler_channel8(ch8, a, b, az0)

    # Lab → OpenSim
    plate_corners = corners[:, :, plate_idx]
    lab_data = plate_local_to_lab(t2, plate_corners)
    os_data = lab_to_opensim_force(lab_data)

    # COP阈值处理 + 单位转换
    valid = np.abs(os_data['Fy']) >= COP_THRESHOLD
    COPx = np.where(valid, os_data['COPx'], 0.0) / 1000.0  # mm -> m
    COPz = np.where(valid, os_data['COPz'], 0.0) / 1000.0  # mm -> m
    Fy = os_data['Fy']
    Fx = os_data['Fx']
    Fz = os_data['Fz']

    # 峰值检测stance
    start_idx, end_idx, hs_idx, to_idx, peak_idx, peak_val = \
        detect_stance_phase_from_peak(Fy, threshold=STANCE_THRESHOLD,
                                       pad_frames=STANCE_PAD_FRAMES)

    n_cut = end_idx - start_idx + 1
    time_col = np.arange(n_cut) / point_rate

    # 截取
    Fy_cut   = Fy[start_idx:end_idx+1].copy()
    Fx_cut   = Fx[start_idx:end_idx+1].copy()
    Fz_cut   = Fz[start_idx:end_idx+1].copy()
    COPx_cut = COPx[start_idx:end_idx+1].copy()
    COPz_cut = COPz[start_idx:end_idx+1].copy()

    # padding归零
    rel_hs = hs_idx - start_idx
    rel_to = to_idx - start_idx
    for arr in [Fx_cut, Fy_cut, Fz_cut, COPx_cut, COPz_cut]:
        arr[:rel_hs] = 0
        arr[rel_to+1:] = 0

    # 保存原始COP（用于对比图）
    COPx_orig = COPx_cut.copy()
    COPz_orig = COPz_cut.copy()

    # COP斜率纠正
    COPx_corr, COPz_corr, slope_info = correct_cop_slope(
        COPx_cut, COPz_cut, time_col, Fy_cut,
        threshold=COP_THRESHOLD,
        middle_ratio=COP_MIDDLE_RATIO,
        rate_multiplier=COP_RATE_MULTIPLIER,
    )

    return {
        'name': trial_name,
        'time': time_col,
        'Fy': Fy_cut,
        'COPx': COPx_corr,
        'COPz': COPz_corr,
        'COPx_orig': COPx_orig,
        'COPz_orig': COPz_orig,
        'peak_idx': peak_idx - start_idx,
        'peak_val': peak_val,
        'hs_idx': rel_hs,
        'to_idx': rel_to,
        'slope_info': slope_info,
        'fs': point_rate,
    }


def plot_single_combined(file_data, output_dir):
    """单文件组合图: Fy + COPx + COPz"""
    configs = [
        ('Fy',   'Vertical GRF (Fy)',        'Force (N)'),
        ('COPx', 'COPx (Slope Corrected)',   'Position (m)'),
        ('COPz', 'COPz (Slope Corrected)',   'Position (m)'),
    ]

    fig, axes = plt.subplots(3, 1, figsize=(12, 10), sharex=True)
    fig.suptitle(f"{file_data['name']} — V2 Processing Result",
                 fontsize=14, fontweight='bold')

    for idx, (var, title, ylabel) in enumerate(configs):
        ax = axes[idx]
        ax.plot(file_data['time'], file_data[var],
                linewidth=1.2, color='blue', alpha=0.9)
        # 标记peak位置
        if var == 'Fy':
            ax.axvline(file_data['time'][file_data['peak_idx']],
                       color='red', linestyle='--', alpha=0.5, label='Peak')
            ax.axvline(file_data['time'][file_data['hs_idx']],
                       color='green', linestyle=':', alpha=0.5, label='Load')
            ax.axvline(file_data['time'][file_data['to_idx']],
                       color='orange', linestyle=':', alpha=0.5, label='Off')
            ax.legend(fontsize=9)
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, loc='left')
        ax.grid(True, linestyle='--', alpha=0.6)

    axes[-1].set_xlabel('Time (s)', fontsize=12)
    plt.tight_layout()

    path = os.path.join(output_dir, f"{file_data['name']}_combined.png")
    plt.savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close()
    return path


def plot_single_compare(file_data, output_dir):
    """单文件对比图: 原始 vs 纠正后 COPx/COPz"""
    cop_configs = [
        ('COPx', 'COPx: Original vs Slope Corrected', 'Position (m)'),
        ('COPz', 'COPz: Original vs Slope Corrected', 'Position (m)'),
    ]

    fig, axes = plt.subplots(2, 1, figsize=(12, 8), sharex=True)
    fig.suptitle(f"{file_data['name']} — COP Slope Correction Comparison",
                 fontsize=14, fontweight='bold')

    for idx, (var, title, ylabel) in enumerate(cop_configs):
        ax = axes[idx]
        ax.plot(file_data['time'], file_data[f'{var}_orig'],
                linewidth=1.0, color='red', alpha=0.7, label='Original')
        ax.plot(file_data['time'], file_data[var],
                linewidth=1.2, color='blue', alpha=0.9, label='Corrected')
        ax.set_ylabel(ylabel, fontsize=11)
        ax.set_title(title, fontsize=12, loc='left')
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.legend(loc='best', fontsize=9)

    axes[-1].set_xlabel('Time (s)', fontsize=12)
    plt.tight_layout()

    path = os.path.join(output_dir, f"{file_data['name']}_compare.png")
    plt.savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close()
    return path


def plot_stance_detection(file_data, output_dir):
    """Stance检测可视化: 标注peak/load/off位置"""
    fig, ax = plt.subplots(figsize=(12, 5))

    ax.plot(file_data['time'], file_data['Fy'],
            linewidth=1.2, color='blue', label='Fy (Vertical GRF)')
    ax.axhline(STANCE_THRESHOLD, color='gray', linestyle='--',
               alpha=0.5, label=f'Threshold={STANCE_THRESHOLD}N')

    # 标记关键点
    peak_t = file_data['time'][file_data['peak_idx']]
    hs_t   = file_data['time'][file_data['hs_idx']]
    to_t   = file_data['time'][file_data['to_idx']]

    ax.plot(peak_t, file_data['Fy'][file_data['peak_idx']],
            'rv', markersize=10, label=f'Peak={file_data["peak_val"]:.1f}N')
    ax.plot(hs_t, file_data['Fy'][file_data['hs_idx']],
            'g^', markersize=8, label='Load')
    ax.plot(to_t, file_data['Fy'][file_data['to_idx']],
            'm^', markersize=8, label='Off')

    # 标记padding区域
    ax.axvspan(file_data['time'][0], hs_t, alpha=0.1, color='gray',
               label='Padding (zeroed)')
    ax.axvspan(to_t, file_data['time'][-1], alpha=0.1, color='gray')

    ax.set_xlabel('Time (s)', fontsize=12)
    ax.set_ylabel('Force (N)', fontsize=12)
    ax.set_title(f"{file_data['name']} — Peak-based Stance Detection",
                 fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.6)
    plt.tight_layout()

    path = os.path.join(output_dir, f"{file_data['name']}_stance.png")
    plt.savefig(path, dpi=DPI, bbox_inches='tight')
    plt.close()
    return path


def plot_trial_overlay(trial_name, files_data, output_dir):
    """按trial叠加折线图"""
    configs = [
        ('Fy',   'Vertical GRF - Fy',        'Force (N)'),
        ('COPx', 'COPx (Slope Corrected)',   'Position (m)'),
        ('COPz', 'COPz (Slope Corrected)',   'Position (m)'),
    ]

    saved = []
    for var, title, ylabel in configs:
        fig, ax = plt.subplots(figsize=(12, 6.75))
        for i, fd in enumerate(files_data):
            ax.plot(fd['time'], fd[var], label=fd['name'],
                    color=plt.cm.tab10(i % 10), linewidth=0.8, alpha=0.8)
        ax.set_title(f"{trial_name} — {title}", fontsize=14, fontweight='bold')
        ax.set_xlabel('Time (s)', fontsize=12)
        ax.set_ylabel(ylabel, fontsize=12)
        ax.grid(True, linestyle='--', alpha=0.6)
        ax.legend(loc='center left', bbox_to_anchor=(1.02, 0.5),
                  fontsize=8, framealpha=0.8)
        plt.tight_layout()

        path = os.path.join(output_dir, f'{trial_name}_{var}.png')
        plt.savefig(path, dpi=DPI, bbox_inches='tight')
        plt.close()
        saved.append(path)
    return saved


def main():
    """主函数"""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

    print("=" * 60)
    print("V2 数据可视化检查工具 (Peak Stance + COP Correction)")
    print("=" * 60)
    print(f"输入目录: {INPUT_DIR}")
    print(f"输出目录: {OUTPUT_DIR}")
    print(f"文件模式: {PATTERN}")
    print(f"Stance阈值: {STANCE_THRESHOLD}N, Padding: {STANCE_PAD_FRAMES}帧")
    print(f"COP纠正: 中间比={COP_MIDDLE_RATIO}, 倍数={COP_RATE_MULTIPLIER}x")
    print("=" * 60)

    # 查找C3D文件
    c3d_files = []
    for root, dirs, files in os.walk(INPUT_DIR):
        for filename in files:
            if fnmatch.fnmatch(filename.lower(), PATTERN.lower()):
                c3d_files.append(os.path.join(root, filename))
    c3d_files.sort()

    if not c3d_files:
        print(f"\n未找到C3D文件！")
        sys.exit(1)

    print(f"\n找到 {len(c3d_files)} 个C3D文件")

    # 提取数据
    all_data = []
    for i, c3d_path in enumerate(c3d_files, 1):
        print(f"\n[{i}/{len(c3d_files)}] {os.path.basename(c3d_path)}")
        try:
            data = extract_plate2_data(c3d_path)
            all_data.append(data)
            si = data['slope_info']
            print(f"  Peak={data['peak_val']:.1f}N, "
                  f"COPx纠正{si['copx_outlier_count']}帧, "
                  f"COPz纠正{si['copy_outlier_count']}帧")
        except Exception as e:
            print(f"  [ERROR] {e}")
            continue

    if not all_data:
        print("\n无可用数据，退出")
        sys.exit(1)

    # 绘图
    print(f"\n开始绘图 ({len(all_data)} 个文件)...")

    # 1) Stance检测图
    for fd in all_data:
        path = plot_stance_detection(fd, OUTPUT_DIR)
        print(f"  [OK] {path}")

    # 2) 单文件组合图
    for fd in all_data:
        path = plot_single_combined(fd, OUTPUT_DIR)
        print(f"  [OK] {path}")

    # 3) 单文件对比图
    for fd in all_data:
        path = plot_single_compare(fd, OUTPUT_DIR)
        print(f"  [OK] {path}")

    # 4) 叠加图（如果有多个文件）
    if len(all_data) > 1:
        paths = plot_trial_overlay("all_trials", all_data, OUTPUT_DIR)
        for p in paths:
            print(f"  [OK] {p}")

    print(f"\n{'=' * 60}")
    print(f"绘图完成! 图像保存在: {OUTPUT_DIR}")
    print(f"{'=' * 60}")


if __name__ == '__main__':
    main()
