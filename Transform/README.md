# C3D to OpenSim Data Transformation Tool

这是一个完整的 Python 工具包，用于将 C3D 生物力学运动捕捉数据转换为 OpenSim 兼容的 `.trc`（标记点）和 `.mot`（地面反力）文件，同时完成坐标系的转换和数据预处理。

## 📋 目录

- [功能概述](#功能概述)
- [坐标系系统](#坐标系系统)
- [文件说明](#文件说明)
- [安装依赖](#安装依赖)
- [使用方法](#使用方法)
- [处理流程详解](#处理流程详解)
- [参数说明](#参数说明)
- [输出文件](#输出文件)
- [示例](#示例)
- [注意事项](#注意事项)
- [参考文献](#参考文献)

---

## 功能概述

本工具包的主要功能包括：

1. **读取 C3D 文件**：使用 `ezc3d` 库读取运动捕捉数据
2. **坐标系统转换**：
   - Kistler 力台本地坐标系 → 实验室全局坐标系
   - 实验室坐标系 → OpenSim Y轴垂直坐标系
3. **力台数据处理**：
   - 支持 Kistler 8通道 Type3 → Type2 力学计算
   - 支持 6通道力数据（Fx, Fy, Fz, Mx, My, Mz）处理
4. **数据滤波和重采样**：使用巴特沃斯低通滤波器
5. **步态阶段检测**：自动检测足着地（Heel Strike）和离地（Toe Off）时刻
6. **数据截取**：提取站立阶段并添加前后缓冲帧
7. **OpenSim 文件输出**：生成 `.trc` 和 `.mot` 格式文件

---

## 坐标系系统

### 坐标系定义

| 坐标系 | X轴方向 | Y轴方向 | Z轴方向 | 单位 |
|--------|----------|----------|----------|------|
| **Kistler 力台本地** | 右 (Right+) | 后 (Posterior+) | 下 (Down+) | mm |
| **实验室全局** | 前 (Forward+) | 左 (Left+) | 上 (Up+) | mm |
| **OpenSim** | 前 (Forward+) | 上 (Up+) | 右 (Right+) | m |

### 坐标系转换

**1. 力台本地 → 实验室全局**

```
GRF_Lab_Fx =  type2['Fy']         (前进方向)
GRF_Lab_Fy =  type2['Fx']         (向左)
GRF_Lab_Fz =  type2['Fz']         (向上)
COP_Lab_X  = plate_cx − ay        (前进方向, mm)
COP_Lab_Y  = plate_cy − ax        (向左, mm)
Tz_lab     = −Tz_kistler          (绕Z轴向上)
```

**2. 实验室全局 → OpenSim**

等价于绕 X 轴旋转 -90°：

```
OpenSim X = Lab X      (前进方向)
OpenSim Y = Lab Z      (垂直向上)
OpenSim Z = -Lab Y     (向右)
```

---

## 文件说明

| 文件名 | 说明 |
|--------|------|
| **`c3d_to_opensim.py`** | 主转换程序，处理 8通道 Kistler 力台 C3D 文件 |
| **`c3d_to_opensim_v2.py`** | **V2** 峰值检测stance + COP斜率纠正转换程序 |
| **`c3d6_to_opensim.py`** | 6通道力数据 C3D 文件转换程序 |
| **`batch_process.py`** | 批量处理工具，支持并行处理和递归搜索 |
| **`batch_process_v2.py`** | **V2** 批量处理工具（峰值stance + COP纠正 + 绝对路径配置） |
| **`batch_c3d6_to_opensim.py`** | 6通道 C3D 文件批量处理工具 |
| **`transform_utils.py`** | 核心工具函数模块（坐标变换、滤波、文件I/O、stance检测、COP纠正等） |
| **`draw_picture_check_grf.py`** | 数据可视化工具，绘制 COP 和自由力矩曲线 |
| **`draw_check_v2.py`** | **V2** 数据可视化检查工具（stance检测图 + COP对比图 + 叠加图） |
| **`run_demo.py`** | 演示脚本，单文件转换示例 |

---

## 安装依赖

### Python 版本要求
- Python 3.7+

### 安装依赖包

```bash
pip install numpy pandas scipy ezc3d matplotlib
```

### 依赖包版本说明
- `numpy`：数组计算和数据处理
- `pandas`：CSV 文件处理
- `scipy`：信号滤波和重采样
- `ezc3d`：C3D 文件读取
- `matplotlib`：数据可视化

---

## 使用方法

### 1. 单文件转换（8通道 C3D）

使用 `c3d_to_opensim.py` 处理标准的 8通道 Kistler 力台 C3D 文件：

```bash
python c3d_to_opensim.py <c3d_file> [output_dir]
```

**示例：**
```bash
python c3d_to_opensim.py data/S15T1V11.c3d output/
```

### 2. 单文件转换（6通道 C3D）

使用 `c3d6_to_opensim.py` 处理包含 6通道力数据的 C3D 文件：

```bash
python c3d6_to_opensim.py <c3d_file> [output_dir]
```

**示例：**
```bash
python c3d6_to_opensim.py data/S15T1V11.c3d output/
```

### 3. 运行演示脚本

使用 `run_demo.py` 运行内置演示：

```bash
python run_demo.py
```

### 4. 批量处理（8通道）

使用 `batch_process.py` 批量转换多个 C3D 文件：

```bash
python batch_process.py
```

运行后会提示输入：
1. C3D 文件所在目录路径
2. 输出目录路径
3. 是否递归搜索子目录
4. 是否启用并行处理
5. 最大进程数（可选）

### 5. 批量处理（6通道）

使用 `batch_c3d6_to_opensim.py` 批量转换 6通道 C3D 文件：

```bash
python batch_c3d6_to_opensim.py <input_dir> [output_dir]
```

### 6. 数据可视化

使用 `draw_picture_check_grf.py` 绘制 COP 和自由力矩曲线：

```bash
python draw_picture_check_grf.py <input_dir> [output_dir] [pattern]
```

**参数说明：**
- `input_dir`：C3D 文件所在目录
- `output_dir`：输出图表目录（可选）
- `pattern`：C3D 文件匹配模式，如 `*.c3d` 或 `S*.c3d`（可选）

---

## 处理流程详解

### 完整处理流程（8通道）

```
Step 1: 读取 C3D 文件
  ├─ 提取 Marker 点数据 (4, n_markers, n_frames)
  ├─ 提取力台模拟数据 (n_channels, n_analog_frames)
  ├─ 提取采样频率（point_rate, analog_rate）
  └─ 提取力台参数（ORIGIN, CORNERS）

Step 2: 滤波和重采样模拟信号
  ├─ 对力台数据应用低通滤波（50 Hz）
  └─ 重采样到 Marker 频率

Step 3: Kistler Type3 → Type2 计算
  ├─ 从 8 通道原始数据计算 6 分量 + 自由力矩
  └─ 输出 Fx, Fy, Fz, ax, ay, Tz

Step 4: 力台本地坐标系 → 实验室坐标系
  └─ 应用坐标变换和 CORNERS 偏移

Step 5: 实验室坐标系 → OpenSim 坐标系
  ├─ 对 Marker 数据应用旋转矩阵
  └─ 对力台数据应用旋转矩阵

Step 6: Marker 数据滤波
  └─ 应用低通滤波（6 Hz）

Step 7: 力台数据重构和单位转换
  ├─ 重构为 OpenSim 9 列格式
  ├─ COP: mm → m
  └─ 力矩: N·mm → N·m

Step 8: 步态阶段检测和截取
  ├─ 根据垂直力检测 Heel Strike 和 Toe Off
  ├─ 添加前后缓冲帧（默认各 25 帧）
  └─ 同步截取 Marker 数据

Step 9: 写入输出文件
  ├─ 生成 .trc 文件（Marker 数据）
  ├─ 生成 .mot 文件（力台数据）
  └─ 生成 cut_records.csv（截取记录）
```

### 6通道 C3D 处理流程

6通道流程与 8通道类似，但有以下区别：
- Step 2 直接处理 6通道力数据（Fx, Fy, Fz, Mx, My, Mz）
- 跳过 Kistler Type3 → Type2 计算
- 从 6通道数据直接计算 COP 和自由力矩

---

## 参数说明

### 主要处理参数

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `marker_cutoff` | float | 6.0 | Marker 数据低通滤波截止频率 (Hz) |
| `force_cutoff` | float | 50.0 | 力台数据低通滤波截止频率 (Hz) |
| `stance_threshold` | float | 30.0 | 步态检测的垂直力阈值 (N) |
| `stance_pad_frames` | int | 25 | 着地前后添加的缓冲帧数 |

### 批量处理参数

| 参数名 | 类型 | 默认值 | 说明 |
|--------|------|--------|------|
| `recursive` | bool | True | 是否递归搜索子目录 |
| `parallel` | bool | False | 是否启用并行处理 |
| `max_workers` | int | None | 最大工作进程数（None=自动） |

---

## 输出文件

### 1. `.trc` 文件格式

```
PathFileType    4    (X/Y/Z)    {filename}
DataRate    CameraRate    NumFrames    NumMarkers    Units    OrigDataRate    OrigDataStartFrame    OrigNumFrames
{frame_rate}    {frame_rate}    {n_frames}    {n_markers}    {units}    {frame_rate}    1    {n_frames}
Frame#    Time    {marker1}            {marker2}            ...
        X1    Y1    Z1    X2    Y2    Z2    ...

{frame_num}    {time}    {x1}    {y1}    {z1}    {x2}    {y2}    {z2}    ...
```

- **数据列数**：2 + n_markers × 3
- **数据顺序**：Frame#, Time, X1, Y1, Z1, X2, Y2, Z2, ...

### 2. `.mot` 文件格式

```
{filename}
version=1
nRows={n_frames}
nColumns={n_cols}
inDegrees=yes
endheader
time    ground_force_vx    ground_force_vy    ground_force_vz    ground_force_px    ground_force_py    ground_force_pz    1_ground_force_vx    ...

{time}    {force1_vx}    {force1_vy}    {force1_vz}    {cop1_x}    {cop1_y}    {cop1_z}    ...
```

- **每个力台 9 列**：force(3) + cop(3) + torque(3)
- **列顺序**：所有力台的 force 和 cop 在前，然后是所有力台的 torque
- **特殊列**：ground_force_py 和 ground_torque_x/z 始终为 0

### 3. `cut_records.csv` 文件格式

```csv
trial,heel_strike_frame,toe_off_frame,cut_start_frame,cut_end_frame,total_frames
S15T1V11,  150,               200,        125,            225,          101
```

用于记录每次截取的步态阶段信息。

---

## 示例

### 示例 1：基本单文件转换

```python
from c3d_to_opensim import process_c3d

result = process_c3d(
    c3d_path='data/S15T1V11.c3d',
    output_dir='output/',
    marker_cutoff=6.0,
    force_cutoff=50.0,
    stance_threshold=30.0,
    stance_pad_frames=25,
)

print(f"TRC 文件: {result['trc_path']}")
print(f"MOT 文件: {result['mot_path']}")
print(f"截取信息: {result['stance_info']}")
```

### 示例 2：批量处理

```python
from batch_process import batch_process

batch_process(
    input_root='data/raw/',
    output_root='data/processed/',
    recursive=True,
    parallel=True,
    max_workers=4,
    process_params={
        'marker_cutoff': 6.0,
        'force_cutoff': 50.0,
        'stance_threshold': 30.0,
        'stance_pad_frames': 25,
    }
)
```

### 示例 3：可视化 COP 数据

```python
from draw_picture_check_grf import batch_extract_grf_data, plot_grf_data

# 批量提取数据
all_data = batch_extract_grf_data('data/raw/', pattern='*.c3d')

# 绘制图表
plot_grf_data(all_data, 'output/plots/', show_plot=False)
```

---

## 注意事项

### ⚠️ 坐标系确认

- 请确保您的实验室坐标系配置与代码中的定义一致（X=前, Y=左, Z=上）
- 如果坐标系不同，请修改 `transform_utils.py` 中的坐标转换函数

### ⚠️ 力台配置

- 确认 C3D 文件中的力台数量与预期一致
- 对于多力台设置，工具会自动检测力台配置并分别处理

### ⚠️ 步态检测

- 默认阈值 30N 可能需要根据您的实验数据调整
- 如果数据本身已只包含站立阶段，可以跳过截取步骤

### ⚠️ 滤波参数

- Marker 滤波频率（默认 6Hz）和力台滤波频率（默认 50Hz）可根据需要调整
- 滤波器使用 4 阶巴特沃斯低通滤波和零相位滤波（filtfilt）

### ⚠️ 并行处理

- 并行处理适用于大量文件的批处理
- 注意内存使用情况，避免同时处理过多大文件

---

## 参考文献

### 相关标准

- **C3D 文件格式**：[C3D.org](https://www.c3d.org/)
- **OpenSim 文件格式**：[OpenSim Documentation](https://simtk.org/projects/opensim/)
- **Kistler 力台**：[Kistler Force Plates](https://www.kistler.com/)

### 相关代码

本项目参考了以下代码实现：
- Kistler 力台本地坐标计算
- 数据处理和滤波函数
- 步态阶段截取算法
- 坐标系旋转函数

---

## V2 版本说明（Peak-based Stance Detection + COP Correction）

V2 版本对 stance 检测和 COP 处理做了重要改进，生成的文件与 V1 完全独立，方便对比。

### V1 vs V2 核心差异

| 特性 | V1 | V2 |
| ---- | --- | --- |
| Stance 检测 | 阈值法：找第一个/最后一个超阈值的帧 | 峰值法：从峰值向两侧遍历找 load/off |
| Padding 归零 | 仅归零 COP 和力矩 | 全面归零：力 + COP + 力矩均归零 |
| COP 纠正 | 无 | 斜率纠正 + 异常检测 |
| 输入方式 | 交互式输入路径 | 绝对路径配置（修改脚本顶部即可） |

### V2 Stance 检测算法

```text
1. 找到垂直力 (Fy) 的最大值 (Peak)
2. 从 Peak 向左遍历 → 第一个 < threshold 的帧 = load (heel-strike)
3. 从 Peak 向右遍历 → 第一个 < threshold 的帧 = off (toe-off)
4. 在 load/off 基础上加 padding 帧
5. padding 区域的力数据手动归零
```

### V2 COP 纠正算法

```text
1. 检测异常帧:
   - 力 < 阈值但 COP 非零 → 异常
   - 帧间 COP 跳变 > jump_threshold → 异常
2. 斜率纠正:
   - 取 stance 中间部分 (默认 30%~70%)，线性拟合得通用斜率 k
   - 从中间向两端遍历，变化率 > rate_multiplier * |k*dt| 的帧视为异常
   - 异常帧用 上一正常值 + k*dt 替换
3. 非 stance 阶段 COP 置零
```

### V2 新增函数（transform_utils.py）

| 函数名 | 说明 |
| ---- | ---- |
| `detect_stance_phase_from_peak()` | 基于峰值的 stance 检测，返回 load/off/peak 信息 |
| `detect_cop_anomalies()` | COP 异常帧检测（力阈值 + 帧间跳变双重机制） |
| `correct_cop_slope()` | 基于中间部分斜率的 COP 纠正 |

### V2 使用方法

#### 批量处理

编辑 `batch_process_v2.py` 顶部配置区域：

```python
INPUT_ROOT  = r"E:\your\data\input"
OUTPUT_ROOT = r"E:\your\data\output_v2"

PROCESS_PARAMS = {
    'stance_threshold':    30.0,   # 力阈值 (N)
    'stance_pad_frames':   25,     # 补偿帧数
    'cop_middle_ratio':    0.3,    # COP斜率纠正中间占比
    'cop_rate_multiplier': 2.0,    # 异常变化率倍数
    'cop_jump_threshold':  0.03,   # COP帧间跳变阈值 (m)
}
```

然后运行：

```bash
python batch_process_v2.py
```

#### 单文件转换

```bash
python c3d_to_opensim_v2.py <c3d_file> [output_dir]
```

#### 可视化检查

编辑 `draw_check_v2.py` 顶部配置区域，然后运行：

```bash
python draw_check_v2.py
```

生成的图表包括：

- `*_stance.png` — Stance 检测可视化（标注 Peak/Load/Off 位置 + padding 区域）
- `*_combined.png` — Fy/COPx/COPz 组合图
- `*_compare.png` — COP 原始 vs 纠正后对比图
- 叠加折线图（多文件叠加对比）

### V2 参数说明

| 参数 | 默认值 | 说明 |
| ---- | ------ | ---- |
| `stance_threshold` | 30.0 N | stance 检测力阈值 |
| `stance_pad_frames` | 25 | 前后补偿帧数 |
| `cop_middle_ratio` | 0.3 | COP 斜率纠正中间部分占比 |
| `cop_rate_multiplier` | 2.0 | COP 变化率异常判定倍数 |
| `cop_jump_threshold` | 0.03 m | COP 帧间跳变阈值 |

---

## 许可证

本项目仅用于学术和研究目的。

## 联系方式

如有问题或建议，请通过 GitHub Issues 联系。
