"""
Batch C3D (6-Channel) → OpenSim Converter
=========================================
Batch process all C3D files in a directory with 6-channel force data.

Usage:
  python batch_c3d6_to_opensim.py <input_dir> [output_dir]
"""

import os
import sys
import glob
import traceback
from c3d6_to_opensim import process_c3d6


def batch_process(input_dir, output_dir=None,
                 marker_cutoff=6.0,
                 force_cutoff=50.0,
                 stance_threshold=30.0,
                 stance_pad_frames=25):
    """
    Batch process all C3D files in input_dir.

    Parameters
    ----------
    input_dir        : str   — directory containing C3D files
    output_dir       : str   — output directory (default: input_dir/output)
    marker_cutoff    : float — low-pass cutoff for markers (Hz, default 6)
    force_cutoff     : float — low-pass cutoff for forces  (Hz, default 50)
    stance_threshold : float — vertical-force threshold for stance (N, default 30)
    stance_pad_frames: int   — padding frames before / after stance (default 25)
    """
    # 查找所有 C3D 文件
    c3d_files = glob.glob(os.path.join(input_dir, '*.c3d'))
    c3d_files.extend(glob.glob(os.path.join(input_dir, '*.C3D')))

    if not c3d_files:
        print(f"[ERROR] 在目录 {input_dir} 中未找到 C39D 文件")
        return

    # 设置输出目录
    if output_dir is None:
        output_dir = os.path.join(input_dir, 'output')

    # 排序文件名
    c3d_files.sort()

    n_total = len(c3d_files)
    n_success = 0
    n_failed = 0
    failed_files = []

    print("=" * 70)
    print("  Batch C3D (6-Channel) -> OpenSim Converter")
    print(f"  输入目录: {input_dir}")
    print(f"  输出目录: {output_dir}")
    print(f"  找到 {n_total} 个 C3D 文件")
    print("=" * 70)

    for i, c3d_file in enumerate(c3d_files, 1):
        filename = os.path.basename(c3d_file)
        print(f"\n[{i}/{n_total}] 处理: {filename}")
        print("-" * 70)

        try:
            result = process_c3d6(
                c3d_file,
                output_dir,
                marker_cutoff=marker_cutoff,
                force_cutoff=force_cutoff,
                stance_threshold=stance_threshold,
                stance_pad_frames=stance_pad_frames,
            )
            n_success += 1
        except Exception as e:
            n_failed += 1
            failed_files.append((filename, str(e)))
            print(f"[ERROR] 处理失败: {filename}")
            print(f"  错误信息: {e}")
            traceback.print_exc()

    # ── Summary ──────────────────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  [批量处理完成]")
    print(f"  总文件数: {n_total}")
    print(f"  成功: {n_success}")
    print(f"  失败: {n_failed}")
    print("=" * 70)

    if failed_files:
        print("\n  失败文件列表:")
        for filename, error in failed_files:
            print(f"    - {filename}: {error}")


if __name__ == '__main__':
    print("=" * 70)
    print("  Batch C3D (6-Channel) -> OpenSim Converter")
    print("=" * 70)

    # 交互式输入
    input_dir = input("请输入 C3D 文件所在目录: ").strip()
    output_dir = input("请输入输出目录 (直接回车默认为输入目录/output): ").strip()

    if output_dir == '':
        output_dir = os.path.join(input_dir, 'output')

    print("\n" + "=" * 70)
    batch_process(input_dir, output_dir)
