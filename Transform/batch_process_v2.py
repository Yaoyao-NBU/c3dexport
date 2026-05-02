"""
Batch C3D to OpenSim Converter V2 (Peak-based Stance Detection)
================================================================
批量C3D转换工具 - 使用峰值检测stance + COP斜率纠正

改进点 (对比 batch_process.py):
  - 使用 detect_stance_phase_from_peak 替代 detect_stance_phase
  - padding帧力数据手动归零
  - 集成 COPx/COPy 异常检测与斜率纠正
  - 绝对路径配置，无需交互输入

Usage:
  直接修改下方配置区域，然后运行:
  python batch_process_v2.py
"""

import os
import sys
import glob
import traceback

from typing import List, Dict, Tuple
from concurrent.futures import ProcessPoolExecutor, as_completed

# 本地模块
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from c3d_to_opensim_v2 import process_c3d_v2


# ═══════════════════════════════════════════════════════════════════
#  配置区域 — 修改此处参数即可
# ═══════════════════════════════════════════════════════════════════

# 输入输出目录（绝对路径）
INPUT_ROOT  = r"G:\Carbon_Plate_Shoes_Data\New_Opensim_Transform\5_try\Raw_data\input"
OUTPUT_ROOT = r"G:\Carbon_Plate_Shoes_Data\New_Opensim_Transform\5_try\Raw_data\output_v2"

# 搜索选项
RECURSIVE = True       # 是否递归搜索子目录
PARALLEL  = False      # 是否启用并行处理
MAX_WORKERS = None     # 并行进程数（None=自动）

# 处理参数
PROCESS_PARAMS = {
    'marker_cutoff':       6.0,    # Marker低通滤波截止频率 (Hz)
    'force_cutoff':        50.0,   # 力数据低通滤波截止频率 (Hz)
    'stance_threshold':    30.0,   # stance检测力阈值 (N)
    'stance_pad_frames':   20,     # 前后补偿帧数
    'cop_middle_ratio':    0.3,    # COP斜率纠正中间部分占比
    'cop_rate_multiplier': 3.0,    # COP变化率异常倍数
    'cop_jump_threshold':  0.03,   # COP帧间跳变阈值 (m)
}

# ═══════════════════════════════════════════════════════════════════


def find_c3d_files(root_dir: str, recursive: bool = True) -> List[str]:
    """查找目录中所有C3D文件"""
    pattern = '**/*.c3d' if recursive else '*.c3d'
    files = glob.glob(os.path.join(root_dir, pattern), recursive=recursive)
    return sorted(files)


def get_output_path(c3d_path: str, input_root: str, output_root: str) -> str:
    """获取输出路径，保持原有文件夹层级"""
    input_root = os.path.normpath(input_root)
    output_root = os.path.normpath(output_root)
    c3d_path = os.path.normpath(c3d_path)

    c3d_dir = os.path.dirname(c3d_path)
    try:
        rel_path = os.path.relpath(c3d_dir, input_root)
    except ValueError:
        raise ValueError(f"输入目录和C3D文件不在同一驱动器: {input_root} vs {c3d_dir}")

    return os.path.join(output_root, rel_path)


def process_single_file(args: Tuple) -> Dict:
    """处理单个C3D文件（用于多进程）"""
    c3d_path, output_dir, process_params = args
    trial_name = os.path.splitext(os.path.basename(c3d_path))[0]

    result = {
        'c3d_path': c3d_path,
        'trial_name': trial_name,
        'status': 'success',
        'error': None,
        'output': None
    }

    try:
        output = process_c3d_v2(
            c3d_path=c3d_path,
            output_dir=output_dir,
            marker_cutoff=process_params.get('marker_cutoff', 6.0),
            force_cutoff=process_params.get('force_cutoff', 50.0),
            stance_threshold=process_params.get('stance_threshold', 30),
            stance_pad_frames=process_params.get('stance_pad_frames', 25),
            cop_middle_ratio=process_params.get('cop_middle_ratio', 0.3),
            cop_rate_multiplier=process_params.get('cop_rate_multiplier', 2.0),
            cop_jump_threshold=process_params.get('cop_jump_threshold', 0.03),
        )
        result['output'] = output
    except Exception as e:
        result['status'] = 'failed'
        result['error'] = str(e) + '\n' + traceback.format_exc()

    return result


def print_progress(current: int, total: int, filename: str = None):
    """打印进度"""
    bar_length = 40
    filled = int(bar_length * current / total)
    bar = '█' * filled + '-' * (bar_length - filled)
    percent = (current / total) * 100

    if filename:
        print(f"\r进度: [{bar}] {percent:.1f}% ({current}/{total}) - {filename}",
              end='', flush=True)
    else:
        print(f"\r进度: [{bar}] {percent:.1f}% ({current}/{total})",
              end='', flush=True)


def batch_process(input_root: str, output_root: str,
                  recursive: bool = True,
                  parallel: bool = False,
                  max_workers: int = None,
                  process_params: Dict = None):
    """批量处理C3D文件"""
    if process_params is None:
        process_params = {}

    print("=" * 70)
    print("     C3D -> OpenSim V2 批量转换工具 (Peak Stance + COP Correction)")
    print("=" * 70)
    print(f"输入目录:       {input_root}")
    print(f"输出目录:       {output_root}")
    print(f"递归搜索:       {'是' if recursive else '否'}")
    print(f"并行处理:       {'是' if parallel else '否'}")
    if parallel:
        print(f"最大进程数:     {max_workers if max_workers else '自动'}")
    print(f"Marker滤波:     {process_params.get('marker_cutoff', 6.0)} Hz")
    print(f"Force滤波:      {process_params.get('force_cutoff', 50.0)} Hz")
    print(f"Stance阈值:     {process_params.get('stance_threshold', 30.0)} N")
    print(f"Padding帧数:    {process_params.get('stance_pad_frames', 25)}")
    print(f"COP中间区间比:  {process_params.get('cop_middle_ratio', 0.3)}")
    print(f"COP异常倍数:    {process_params.get('cop_rate_multiplier', 2.0)}x")
    print("=" * 70)

    if not os.path.isdir(input_root):
        print(f"\n错误: 输入目录不存在 - {input_root}")
        sys.exit(1)

    os.makedirs(output_root, exist_ok=True)

    # 查找C3D文件
    print(f"\n扫描目录中...")
    c3d_files = find_c3d_files(input_root, recursive)

    if not c3d_files:
        print(f"\n未找到C3D文件！")
        return

    print(f"找到 {len(c3d_files)} 个C3D文件")
    for f in c3d_files[:5]:
        print(f"  - {os.path.relpath(f, input_root)}")
    if len(c3d_files) > 5:
        print(f"  ... 还有 {len(c3d_files) - 5} 个文件")

    # 准备任务
    tasks = []
    for c3d_path in c3d_files:
        output_dir = get_output_path(c3d_path, input_root, output_root)
        tasks.append((c3d_path, output_dir, process_params))

    # 执行
    results = []
    success_count = 0
    failed_count = 0

    print(f"\n开始处理...")
    print("-" * 70)

    if parallel:
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(process_single_file, task): task
                              for task in tasks}
            for i, future in enumerate(as_completed(future_to_task)):
                result = future.result()
                results.append(result)
                rel_path = os.path.relpath(result['c3d_path'], input_root)
                if result['status'] == 'success':
                    success_count += 1
                    print(f"[{i+1}/{len(tasks)}] OK {rel_path}")
                else:
                    failed_count += 1
                    print(f"[{i+1}/{len(tasks)}] FAIL {rel_path} - "
                          f"{result['error'][:50]}...")
    else:
        for i, task in enumerate(tasks):
            c3d_path = task[0]
            rel_path = os.path.relpath(c3d_path, input_root)
            print_progress(i, len(tasks), rel_path)

            result = process_single_file(task)
            results.append(result)

            if result['status'] == 'success':
                success_count += 1
            else:
                failed_count += 1

        print_progress(len(tasks), len(tasks))

    # 总结
    print(f"\n{'=' * 70}")
    print(f"处理完成! 总计: {len(results)}, 成功: {success_count}, 失败: {failed_count}")
    print(f"{'=' * 70}")

    if failed_count > 0:
        print("\n失败文件:")
        print("-" * 70)
        for result in results:
            if result['status'] == 'failed':
                rel_path = os.path.relpath(result['c3d_path'], input_root)
                print(f"\n文件: {rel_path}")
                print(f"错误: {result['error']}")

    if success_count > 0:
        print("\n成功输出文件（前10个）:")
        print("-" * 70)
        count = 0
        for result in results:
            if result['status'] == 'success':
                out = result['output']
                rel_path = os.path.relpath(result['c3d_path'], input_root)
                print(f"{rel_path}:")
                print(f"  .trc: {out['trc_path']}")
                print(f"  .mot: {out['mot_path']}")
                count += 1
                if count >= 10:
                    print(f"... 还有 {success_count - 10} 个")
                    break


def main():
    """直接使用配置区域参数运行"""
    batch_process(
        input_root=INPUT_ROOT,
        output_root=OUTPUT_ROOT,
        recursive=RECURSIVE,
        parallel=PARALLEL,
        max_workers=MAX_WORKERS,
        process_params=PROCESS_PARAMS,
    )


if __name__ == '__main__':
    main()
