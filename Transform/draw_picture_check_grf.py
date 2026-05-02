"""
批量读取 C3D 文件，提取 COP 和自由力矩数据并可视化
====================================================

功能：
1. 使用 transform_utils.py 中的函数读取 C3D 文件
2. 计算力台数组的数据
3. 将 GRF 数组中的 COP 和自由力矩存储到 numpy 数组
4. 批量读取 C3D 并补充数组
5. 绘制所有数据的折线图：
   - COP 前后方向 (COPx) 到一起
   - COP 左右方向 (COPz) 到一起
   - 自由力矩 (Tz) 到一起

Usage:
方法1：在代码中配置路径（推荐） - 修改 main() 函数中的配置区域
方法2：命令行参数 - python draw_picture_check_grf.py <input_dir> [output_dir] [c3d_pattern]
"""

import os
import sys
import numpy as np
import ezc3d
import fnmatch
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties

# 添加 transform_utils 模块路径
sys.path.insert(0, os.path.dirname(__file__))
from transform_utils import (
    compute_kistler_type2,
    plate_local_to_lab,
    lab_to_opensim_force,
    butter_lowpass_filter,
    detect_stance_phase,
)


def extract_grf_data_from_c3d(c3d_path):
    """
    从单个 C3D 文件提取 GRF 数据（COP 和自由力矩）

    Parameters
    ----------
    c3d_path : str
        C3D 文件路径

    Returns
    -------
    dict
        包含提取的数据：
        - cop_x: 前后方向 COP (mm)
        - cop_y: 左右方向 COP (mm)
        - free_moment: 自由力矩 Tz (N·mm)
        - trial_name: 试次名称
        - n_frames: 帧数
    """
    print(f"正在读取: {c3d_path}")

    # 读取 C3D 文件
    c = ezc3d.c3d(c3d_path)
    trial_name = os.path.splitext(os.path.basename(c3d_path))[0]

    # 获取模拟数据采样率
    analog_rate = float(c['header']['analogs']['frame_rate'])

    # 获取模拟数据
    analog_data = c['data']['analogs'][0]  # (n_channels, n_analog_frames)

    # 获取力台参数
    fp_params = c['parameters']['FORCE_PLATFORM']
    n_plates = int(fp_params['USED']['value'][0])
    channels = fp_params['CHANNEL']['value']
    origin = fp_params['ORIGIN']['value']
    corners = np.array(fp_params['CORNERS']['value'])  # (3, 4, n_plates)

    print(f"  力台数量: {n_plates}, 采样率: {analog_rate} Hz")

    # 计算每个力台的 Type 2 数据
    plate_type2 = []
    for pi in range(n_plates):
        ch_idx = channels[:, pi].astype(int) - 1
        ch8 = analog_data[ch_idx, :]

        b = float(origin[0, pi])
        a = float(origin[1, pi])
        az0 = float(origin[2, pi])

        t2 = compute_kistler_type2(ch8, a, b, az0)
        plate_type2.append(t2)

    # 转换到实验室坐标系
    plate_lab = []
    for pi in range(n_plates):
        plate_corners = corners[:, :, pi]  # (3, 4)
        lab_data = plate_local_to_lab(plate_type2[pi], plate_corners)
        plate_lab.append(lab_data)

    # 转换到 OpenSim 坐标系
    plate_os = []
    for pi in range(n_plates):
        os_data = lab_to_opensim_force(plate_lab[pi])
        plate_os.append(os_data)

    # 只使用第二个测力台（索引1），因为第一个测力台没有数据
    target_plate_idx = 1
    if target_plate_idx >= n_plates:
        raise ValueError(f"需要第二个测力台，但只有 {n_plates} 个测力台")

    d = plate_os[target_plate_idx]

    # 在 OpenSim 坐标系中，垂直力是 Fy
    vertical_force = d['Fy']

    # 使用 detect_stance_phase 检测 stance 阶段，阈值设为 35N
    try:
        start_idx, end_idx, hs_idx, to_idx = detect_stance_phase(
            vertical_force,
            threshold=35.0,
            pad_frames=25,
        )
        print(f"  Stance 检测: HS 帧={hs_idx+1}, TO 帧={to_idx+1}, 截取范围={start_idx+1}-{end_idx+1}")
    except ValueError as e:
        print(f"  [WARN] Stance 检测失败: {e}")
        print(f"  -> 使用全部数据")
        start_idx, end_idx = 0, len(vertical_force) - 1

    # 截取 stance 阶段的数据
    n_frames = end_idx - start_idx + 1

    # 提取 COP 和自由力矩数据
    # 在 OpenSim 坐标系中：
    # - COPx: 前后方向
    # - COPz: 左右方向
    # - Tz: 自由垂直力矩

    cop_x_data = [d['COPx'][start_idx:end_idx+1] / 1000.0]  # 转换为 M
    cop_y_data = [d['COPz'][start_idx:end_idx+1] / 1000.0]  # 转换为 M
    free_moment_data = [d['Tz'][start_idx:end_idx+1] / 1000.0]  # 转换为 N·M

    print(f"  提取完成: {n_frames} 帧 (stance 阶段), COP 单位=M, 自由力矩单位=N·M")

    return dict(
        cop_x=cop_x_data,
        cop_y=cop_y_data,
        free_moment=free_moment_data,
        trial_name=trial_name,
        n_frames=n_frames,
        n_plates=1,  # 只使用一个测力台
    )


def batch_extract_grf_data(input_dir, pattern="*.c3d"):
    """
    批量提取 C3D 文件的 GRF 数据（递归遍历所有子目录）

    Parameters
    ----------
    input_dir : str
        输入目录
    pattern : str
        C3D 文件模式（支持通配符，如 *.c3d 或 S*.c3d）

    Returns
    -------
    list[dict]
        每个文件的提取数据
    """
    
    # 递归遍历目录，查找所有匹配的 C3D 文件
    c3d_files = []

    print(f"正在递归搜索目录: {input_dir}")
    print(f"文件匹配模式: {pattern}")
    print("-" * 60)

    for root, dirs, files in os.walk(input_dir):
        for filename in files:
            # 检查文件是否匹配模式
            if fnmatch.fnmatch(filename.lower(), pattern.lower()):
                c3d_files.append(os.path.join(root, filename))

    c3d_files.sort()

    if len(c3d_files) == 0:
        raise ValueError(f"在 {input_dir} 及其子目录中未找到匹配 '{pattern}' 的 C3D 文件")

    print(f"\n找到 {len(c3d_files)} 个 C3D 文件")
    print("=" * 60)

    all_data = []

    for idx, c3d_file in enumerate(c3d_files, 1):
        try:
            print(f"[{idx}/{len(c3d_files)}] ", end="")
            data = extract_grf_data_from_c3d(c3d_file)
            all_data.append(data)
            print("-" * 60)
        except Exception as e:
            print(f"  [ERROR] 处理 {c3d_file} 时出错: {e}")
            import traceback
            traceback.print_exc()
            continue

    return all_data


def plot_grf_data(all_data, output_dir, show_plot=False):
    """
    绘制所有数据的折线图

    Parameters
    ----------
    all_data : list[dict]
        所有提取的数据
    output_dir : str
        输出目录
    show_plot : bool
        是否显示图形
    """
    # 设置中文字体
    font = FontProperties(fname=r"C:\Windows\Fonts\msyh.ttc", size=12)
    title_font = FontProperties(fname=r"C:\Windows\Fonts\msyh.ttc", size=14, weight='bold')

    # 创建输出目录
    os.makedirs(output_dir, exist_ok=True)

    # 确定最大帧数用于时间轴
    max_frames = max(d['n_frames'] for d in all_data)

    # 创建时间轴
    time_axis = np.arange(max_frames) / 1000.0  # 假设 1000 Hz 采样率

    # ──────────────────────────────────────────────────────────────────────
    # 图 1: COP 前后方向 (COPx)
    # ──────────────────────────────────────────────────────────────────────
    fig1, ax1 = plt.subplots(figsize=(12, 6))

    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            cop_x = data['cop_x'][pi]
            ax1.plot(time_axis[:len(cop_x)], cop_x,
                    #label=f"{trial_name}_FP{pi+1}", #为了避免图例过多，这里暂时注释掉标签
                    linewidth=1)

    ax1.set_xlabel('时间 (s)', fontproperties=font)
    ax1.set_ylabel('COP 前后方向 (m)', fontproperties=font)
    ax1.set_title('COP 前后方向 (COPx)', fontproperties=title_font)
    ax1.grid(True, alpha=0.3)

    # 放置图例
    ax1.legend(prop=font, loc='best', fontsize=8)

    output_path1 = os.path.join(output_dir, "cop_x_all.png")
    plt.savefig(output_path1, dpi=300, bbox_inches='tight')
    print(f"保存: {output_path1}")

    #if show_plot:
    #    plt.show()
    #else:
    plt.close()

    # ──────────────────────────────────────────────────────────────────────
    # 图 2: COP 左右方向 (COPz)
    # ──────────────────────────────────────────────────────────────────────
    fig2, ax2 = plt.subplots(figsize=(12, 6))

    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            cop_y = data['cop_y'][pi]
            ax2.plot(time_axis[:len(cop_y)], cop_y,
                    #label=f"{trial_name}_FP{pi+1}", #为了避免图例过多，这里暂时注释掉标签
                    linewidth=1)

    ax2.set_xlabel('时间 (s)', fontproperties=font)
    ax2.set_ylabel('COP 左右方向 (m)', fontproperties=font)
    ax2.set_title('COP 左右方向 (COPz)', fontproperties=title_font)
    ax2.grid(True, alpha=0.3)

    # 放置图例
    ax2.legend(prop=font, loc='best', fontsize=8)

    output_path2 = os.path.join(output_dir, "cop_z_all.png")
    plt.savefig(output_path2, dpi=300, bbox_inches='tight')
    print(f"保存: {output_path2}")

    #if show_plot:
    #    plt.show()
    #else:
    plt.close()

    # ──────────────────────────────────────────────────────────────────────
    # 图 3: 自由力矩 (Tz)
    # ──────────────────────────────────────────────────────────────────────
    fig3, ax3 = plt.subplots(figsize=(12, 6))

    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            free_moment = data['free_moment'][pi]
            ax3.plot(time_axis[:len(free_moment)], free_moment,
                    #label=f"{trial_name}_FP{pi+1}", #为了避免图例过多，这里暂时注释掉标签
                    linewidth=1)

    ax3.set_xlabel('时间 (s)', fontproperties=font)
    ax3.set_ylabel('自由力矩 Tz (N·m)', fontproperties=font)
    ax3.set_title('自由力矩 (Tz)', fontproperties=title_font)
    ax3.grid(True, alpha=0.3)

    # 放置图例
    ax3.legend(prop=font, loc='best', fontsize=8)

    output_path3 = os.path.join(output_dir, "free_moment_all.png")
    plt.savefig(output_path3, dpi=300, bbox_inches='tight')
    print(f"保存: {output_path3}")

    #if show_plot:
    #    plt.show()
    #else:
    plt.close()

"""
    # ──────────────────────────────────────────────────────────────────────
    # 图 4: 子图组合 (所有三个图在一张图上)
    # ──────────────────────────────────────────────────────────────────────
    fig4, (ax4a, ax4b, ax4c) = plt.subplots(3, 1, figsize=(12, 10))

    # COPx
    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            cop_x = data['cop_x'][pi]
            ax4a.plot(time_axis[:len(cop_x)], cop_x,
                      label=f"{trial_name}_FP{pi+1}",
                      linewidth=0.8)

    ax4a.set_ylabel('COPx (mm)', fontproperties=font)
    ax4a.set_title('COP 前后方向', fontproperties=title_font)
    ax4a.grid(True, alpha=0.3)
    ax4a.legend(prop=font, loc='best', fontsize=8)

    # COPz
    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            cop_y = data['cop_y'][pi]
            ax4b.plot(time_axis[:len(cop_y)], cop_y,
                      label=f"{trial_name}_FP{pi+1}",
                      linewidth=0.8)

    ax4b.set_ylabel('COPz (mm)', fontproperties=font)
    ax4b.set_title('COP 左右方向', fontproperties=title_font)
    ax4b.grid(True, alpha=0.3)
    ax4b.legend(prop=font, loc='best', fontsize=8)

    # Free Moment
    for idx, data in enumerate(all_data):
        trial_name = data['trial_name']
        for pi in range(data['n_plates']):
            free_moment = data['free_moment'][pi]
            ax4c.plot(time_axis[:len(free_moment)], free_moment,
                      label=f"{trial_name}_FP{pi+1}",
                      linewidth=0.8)

    ax4c.set_xlabel('时间 (s)', fontproperties=font)
    ax4c.set_ylabel('Tz (N·mm)', fontproperties=font)
    ax4c.set_title('自由力矩', fontproperties=title_font)
    ax4c.grid(True, alpha=0.3)
    ax4c.legend(prop=font, loc='best', fontsize=8)

    plt.tight_layout()

    output_path4 = os.path.join(output_dir, "grf_all_combined.png")
    plt.savefig(output_path4, dpi=300, bbox_inches='tight')
    print(f"保存: {output_path4}")

    if show_plot:
        plt.show()
    else:
        plt.close()
"""

def main():
    """主函数"""
    # ==================== 配置区域 ====================
    # 修改下面的路径配置，然后直接运行脚本即可

    # C3D 文件所在目录
    INPUT_DIR = r"G:\Carbon_Plate_Shoes_Data\New_Opensim_Transform\2_try\raw_data\input"

    # 输出图表目录（可选，默认为输入目录下的 plots 文件夹）
    OUTPUT_DIR = r"G:\Carbon_Plate_Shoes_Data\New_Opensim_Transform\2_try\raw_data\Picture"

    # C3D 文件匹配模式（例如：*.c3d 或 S*.c3d）
    PATTERN = "*.c3d"
    # =================================================

    # 支持命令行参数覆盖（优先级高于上面配置）
    if len(sys.argv) > 1:
        input_dir = sys.argv[1]
        output_dir = sys.argv[2] if len(sys.argv) > 2 else OUTPUT_DIR
        pattern = sys.argv[3] if len(sys.argv) > 3 else PATTERN
    else:
        input_dir = INPUT_DIR
        output_dir = OUTPUT_DIR
        pattern = PATTERN

    print("=" * 60)
    print("批量提取 C3D GRF 数据并绘图")
    print(f"输入目录: {input_dir}")
    print(f"输出目录: {output_dir}")
    print(f"文件模式: {pattern}")
    print("=" * 60)

    try:
        # 批量提取数据
        all_data = batch_extract_grf_data(input_dir, pattern)

        if len(all_data) == 0:
            print("\n[ERROR] 没有成功提取任何数据")
            sys.exit(1)

        print(f"\n成功提取 {len(all_data)} 个文件的数据")

        # 绘制图表
        print("\n开始绘制图表...")
        plot_grf_data(all_data, output_dir, show_plot=True)

        print("\n" + "=" * 60)
        print("[DONE] 处理完成！")
        print("=" * 60)

    except Exception as e:
        print(f"\n[ERROR] 处理失败: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
