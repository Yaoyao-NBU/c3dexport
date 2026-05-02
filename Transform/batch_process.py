"""
Batch C3D to OpenSim Converter
===============================
批量遍历文件夹中的C3D文件并转换为OpenSim格式(.trc + .mot)

功能:
  - 递归搜索指定目录中的所有.c3d文件
  - 对每个C3D文件执行转换处理
  - 输出文件保持原有文件夹层级结构
  - 支持并行处理以提高效率
  - 完善的错误处理和进度显示

使用方法:
  # 基本使用 - 输入目录和输出目录作为参数
  python batch_process.py

  # 运行后会提示输入：
  # 1. 输入目录路径（包含C3D文件）
  # 2. 输出目录路径（输出结果）
"""

import os
import sys
import glob
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import List, Dict, Tuple
import traceback

# 添加当前目录到路径以支持本地导入
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, SCRIPT_DIR)

from c3d_to_opensim import process_c3d


def find_c3d_files(root_dir: str, recursive: bool = True) -> List[str]:
    """
    查找指定目录中的所有C3D文件

    Parameters
    ----------
    root_dir : str
        根目录路径
    recursive : bool
        是否递归搜索子目录（默认True）

    Returns
    -------
    List[str] — C3D文件路径列表
    """
    pattern = '**/*.c3d' if recursive else '*.c3d'
    files = glob.glob(os.path.join(root_dir, pattern), recursive=recursive)
    return sorted(files)


def get_output_path(c3d_path: str, input_root: str, output_root: str) -> str:
    """
    获取输出目录路径，保持原有文件夹层级结构

    Parameters
    ----------
    c3d_path : str
        C3D文件完整路径
    input_root : str
        输入根目录
    output_root : str
        输出根目录

    Returns
    -------
    str — 输出目录路径（保持原有文件夹层级结构）

    Examples
    --------
    >>> input_root = "C:/Data/input"
    >>> output_root = "C:/Data/output"
    >>> c3d_path = "C:/Data/input/subject1/trial1/file.c3d"
    >>> get_output_path(c3d_path, input_root, output_root)
    'C:/Data/output/subject1/trial1'
    """
    # 规范化路径（去掉尾部斜杠）
    input_root = os.path.normpath(input_root)
    output_root = os.path.normpath(output_root)
    c3d_path = os.path.normpath(c3d_path)

    # 获取C3D文件所在目录
    c3d_dir = os.path.dirname(c3d_path)

    # 计算相对路径（保持原有文件夹结构）
    try:
        rel_path = os.path.relpath(c3d_dir, input_root)
    except ValueError:
        # 如果不在同一驱动器上，使用完整路径
        raise ValueError(f"输入目录和C3D文件不在同一驱动器上: {input_root} vs {c3d_dir}")

    # 构建输出目录
    output_dir = os.path.join(output_root, rel_path)

    return output_dir


def process_single_file(args: Tuple) -> Dict:
    """
    处理单个C3D文件（用于多进程）

    Parameters
    ----------
    args : Tuple
        (c3d_path, output_dir, process_params)

    Returns
    -------
    Dict — 处理结果字典
    """
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
        output = process_c3d(
            c3d_path=c3d_path,
            output_dir=output_dir,
            marker_cutoff=process_params.get('marker_cutoff', 6.0),
            force_cutoff=process_params.get('force_cutoff', 50.0),
            stance_threshold=process_params.get('stance_threshold', 30),
            stance_pad_frames=process_params.get('stance_pad_frames', 25),
        )
        result['output'] = output
    except Exception as e:
        result['status'] = 'failed'
        result['error'] = str(e) + '\n' + traceback.format_exc()

    return result


def print_progress(current: int, total: int, filename: str = None):
    """打印进度信息"""
    bar_length = 40
    filled = int(bar_length * current / total)
    bar = '█' * filled + '-' * (bar_length - filled)
    percent = (current / total) * 100

    if filename:
        print(f"\r进度: [{bar}] {percent:.1f}% ({current}/{total}) - {filename}", end='', flush=True)
    else:
        print(f"\r进度: [{bar}] {percent:.1f}% ({current}/{total})", end='', flush=True)


def batch_process(input_root: str,
                  output_root: str,
                  recursive: bool = True,
                  parallel: bool = False,
                  max_workers: int = None,
                  process_params: Dict = None):
    """
    批量处理C3D文件

    Parameters
    ----------
    input_root : str
        输入根目录路径
    output_root : str
        输出根目录路径
    recursive : bool
        是否递归搜索子目录
    parallel : bool
        是否启用并行处理
    max_workers : int or None
        最大工作进程数（None表示自动）
    process_params : Dict
        处理参数（marker_cutoff, force_cutoff等）
    """
    if process_params is None:
        process_params = {}

    print("=" * 70)
    print("          C3D -> OpenSim 批量转换工具")
    print("=" * 70)
    print(f"输入目录:     {input_root}")
    print(f"输出目录:     {output_root}")
    print(f"递归搜索:     {'是' if recursive else '否'}")
    print(f"并行处理:     {'是' if parallel else '否'}")
    if parallel:
        print(f"最大进程数:   {max_workers if max_workers else '自动'}")
    print(f"Marker滤波:   {process_params.get('marker_cutoff', 6.0)} Hz")
    print(f"Force滤波:    {process_params.get('force_cutoff', 50.0)} Hz")
    print("=" * 70)

    # 验证目录存在
    if not os.path.isdir(input_root):
        print(f"\n错误: 输入目录不存在 - {input_root}")
        sys.exit(1)

    # 创建输出目录
    os.makedirs(output_root, exist_ok=True)
    print(f"\n已创建输出目录: {output_root}")

    # 查找C3D文件
    print(f"\n扫描目录中...")
    c3d_files = find_c3d_files(input_root, recursive)

    if not c3d_files:
        print(f"\n未找到C3D文件！")
        return

    print(f"找到 {len(c3d_files)} 个C3D文件")
    print(f"前5个文件:")
    for f in c3d_files[:5]:
        rel_path = os.path.relpath(f, input_root)
        print(f"  - {rel_path}")
    if len(c3d_files) > 5:
        print(f"  ... 还有 {len(c3d_files) - 5} 个文件")

    # 准备处理参数
    tasks = []
    for c3d_path in c3d_files:
        output_dir = get_output_path(c3d_path, input_root, output_root)
        tasks.append((c3d_path, output_dir, process_params))

    # 执行处理
    results = []
    success_count = 0
    failed_count = 0

    print(f"\n开始处理...")
    print("-" * 70)

    if parallel:
        # 并行处理
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            future_to_task = {executor.submit(process_single_file, task): task for task in tasks}

            for i, future in enumerate(as_completed(future_to_task)):
                result = future.result()
                results.append(result)

                if result['status'] == 'success':
                    success_count += 1
                    rel_path = os.path.relpath(result['c3d_path'], input_root)
                    print(f"[{i+1}/{len(tasks)}] ✓ {rel_path}")
                else:
                    failed_count += 1
                    rel_path = os.path.relpath(result['c3d_path'], input_root)
                    print(f"[{i+1}/{len(tasks)}] ✗ {rel_path} - {result['error'][:50]}...")
    else:
        # 串行处理
        for i, task in enumerate(tasks):
            c3d_path, output_dir, _ = task
            rel_path = os.path.relpath(c3d_path, input_root)

            print_progress(i, len(tasks), rel_path)

            result = process_single_file(task)
            results.append(result)

            if result['status'] == 'success':
                success_count += 1
            else:
                failed_count += 1

        print_progress(len(tasks), len(tasks))

    # 打印总结
    print("\n" + "=" * 70)
    print("处理完成!")
    print(f"总计:   {len(results)} 个文件")
    print(f"成功:   {success_count} 个")
    print(f"失败:   {failed_count} 个")
    print("=" * 70)

    # 显示失败详情
    if failed_count > 0:
        print("\n失败文件详情:")
        print("-" * 70)
        for result in results:
            if result['status'] == 'failed':
                rel_path = os.path.relpath(result['c3d_path'], input_root)
                print(f"\n文件: {rel_path}")
                print(f"错误: {result['error']}")

    # 显示成功输出路径
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
                    print(f"... 还有 {success_count - 10} 个成功文件")
                    break


def main():
    """命令行入口 - 交互式输入"""
    print("=" * 70)
    print("          C3D -> OpenSim 批量转换工具")
    print("=" * 70)

    # 输入目录
    input_root = input("\n请输入包含C3D文件的目录路径: ").strip()
    input_root = input_root.strip('"\'')  # 去除引号

    if not os.path.isdir(input_root):
        print(f"错误: 目录不存在 - {input_root}")
        sys.exit(1)

    # 输出目录
    output_root = input("\n请输入输出目录路径: ").strip()
    output_root = output_root.strip('"\'')  # 去除引号

    # 询问是否递归
    recursive_input = input("\n是否递归搜索子目录？ (y/n，默认y): ").strip().lower()
    recursive = recursive_input != 'n'

    # 询问是否并行
    parallel_input = input("\n是否启用并行处理？ (y/n，默认n): ").strip().lower()
    parallel = parallel_input == 'y'

    max_workers = None
    if parallel:
        workers_input = input("最大进程数（留空自动）: ").strip()
        if workers_input:
            try:
                max_workers = int(workers_input)
            except ValueError:
                print("警告: 进程数格式错误，使用自动模式")

    # 构建处理参数（使用默认值）
    process_params = {
        'marker_cutoff': 6.0,
        'force_cutoff': 60,
        'stance_threshold': 35.0,
        'stance_pad_frames': 25,
    }

    # 执行批量处理
    batch_process(
        input_root=input_root,
        output_root=output_root,
        recursive=recursive,
        parallel=parallel,
        max_workers=max_workers,
        process_params=process_params,
    )


if __name__ == '__main__':
    main()
