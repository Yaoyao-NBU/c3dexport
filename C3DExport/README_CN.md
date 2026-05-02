# C3DExport

**一个将 C3D 生物力学运动捕捉文件转换为 OpenSim 兼容格式的 Python 库。**

C3DExport 可将 C3D 文件（生物力学运动捕捉中常用格式）转换为 OpenSim 兼容的 `.trc`（标记点轨迹）和 `.mot`（地面反作用力）文件。支持 Kistler 测力台的 8 通道和 6 通道配置。

---

## 目录

- [功能特性](#功能特性)
- [项目结构](#项目结构)
- [安装](#安装)
- [快速开始](#快速开始)
- [模块：`C3DExport.core`](#模块c3dexportcore)
- [模块：`C3DExport.io`](#模块c3dexportio)
- [模块：`C3DExport.utils`](#模块c3dexportutils)
- [输出文件格式](#输出文件格式)
- [依赖项](#依赖项)
- [许可证](#许可证)

---

## 功能特性

- **三种转换管线**，适配不同测力台配置
- **自动站立相检测**（阈值法和峰值法）
- **压力中心（COP）异常检测与校正**
- **坐标系变换**：Kistler 测力台局部坐标系 -> 实验室全局坐标系 -> OpenSim Y 轴向上坐标系
- **巴特沃斯低通滤波**，用于标记点和力数据
- **多相重采样**，同步模拟信号和标记点数据的采样率
- **批处理**支持，通过 `Transform/` 目录中的辅助脚本实现

---

## 项目结构

```
c3dexport/
├── C3DExport/                 # 核心库包
│   ├── __init__.py            # 包入口，导出主转换函数
│   ├── core.py                # 主转换管线
│   ├── io.py                  # OpenSim 文件 I/O（.trc 和 .mot 写入器）
│   └── utils.py               # 工具函数（旋转、滤波、检测）
│
├── Transform/                 # 辅助脚本和批处理
│   ├── c3d_to_opensim.py      # 独立 8 通道转换器
│   ├── c3d_to_opensim_v2.py   # 独立 8 通道转换器（v2 版）
│   ├── c3d6_to_opensim.py     # 独立 6 通道转换器
│   ├── batch_process.py       # 批处理脚本
│   ├── batch_process_v2.py    # 批处理脚本（v2 版）
│   ├── batch_c3d6_to_opensim.py # 6 通道批处理
│   ├── draw_picture_check_grf.py # 地面反作用力可视化
│   ├── draw_check_v2.py       # 可视化（v2 版）
│   ├── transform_utils.py     # 附加变换工具
│   └── README.md              # Transform 模块文档
│
└── README.md                  # 英文文档
└── README_CN.md               # 本文件（中文文档）
```

---

## 安装

### 前置要求

- Python 3.8+
- 依赖包：

```bash
pip install numpy scipy ezc3d pandas
```

### 安装方式

将 `C3DExport` 文件夹复制到你的项目目录中，然后导入：

```python
from C3DExport import convert_c3d, convert_c3d_v2, convert_c3d6
```

---

## 快速开始

```python
from C3DExport import convert_c3d

# 转换 8 通道 Kistler C3D 文件
result = convert_c3d("data/trial.c3d", "output/")

print(result['trc_path'])      # -> "output/trial.trc"
print(result['mot_path'])      # -> "output/trial.mot"
print(result['csv_path'])      # -> "output/cut_records.csv"
print(result['stance_info'])   # -> 站立相元数据字典
```

### 使用 V2 版本（峰值站立检测 + COP 校正）

```python
from C3DExport import convert_c3d_v2

result = convert_c3d_v2("data/trial.c3d", "output/",
                         cop_middle_ratio=0.3,
                         cop_rate_multiplier=2.0,
                         cop_jump_threshold=0.03)
```

### 使用 6 通道转换器

```python
from C3DExport import convert_c3d6

result = convert_c3d6("data/trial.c3d", "output/")
```

---

## 模块：`C3DExport.core`

核心模块提供主转换管线。三个转换器遵循相同的通用流程：

1. **读取 C3D** - 从 C3D 文件中解析标记点和测力台数据
2. **重采样** - 将模拟（力）数据同步到标记点帧率
3. **力计算** - 将原始测力台通道转换为 Fx/Fy/Fz/COP/Tz
4. **坐标变换** - 测力台局部坐标 -> 实验室坐标 -> OpenSim Y 轴向上坐标
5. **滤波** - 对标记点数据应用巴特沃斯低通滤波
6. **站立相检测** - 识别着地（heel-strike）和离地（toe-off）事件
7. **输出** - 写入 `.trc`、`.mot` 和 CSV 记录

### 主转换函数

#### `convert_c3d(c3d_path, output_dir, ...)`

将 **8 通道 Kistler** C3D 转换为 OpenSim 文件，使用**阈值法**进行站立相检测。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `c3d_path` | `str` | 必填 | 输入 `.c3d` 文件路径 |
| `output_dir` | `str` | 必填 | 输出文件目录 |
| `marker_cutoff` | `float` | `6.0` | 标记点低通滤波截止频率（Hz） |
| `force_cutoff` | `float` | `50.0` | 力数据低通滤波截止频率（Hz） |
| `stance_threshold` | `float` | `30.0` | 站立相检测的垂直力阈值（N） |
| `stance_pad_frames` | `int` | `25` | 站立相前后的填充帧数 |
| `verbose` | `bool` | `True` | 是否打印进度信息 |

**返回值：** `dict`，包含以下键：
- `trc_path` (`str`): 输出 `.trc` 文件路径
- `mot_path` (`str`): 输出 `.mot` 文件路径
- `csv_path` (`str`): 切割记录 CSV 文件路径
- `stance_info` (`dict`): 站立相元数据（着地帧、离地帧、切割范围等）

---

#### `convert_c3d_v2(c3d_path, output_dir, ...)`

将 **8 通道 Kistler** C3D 转换为 OpenSim 文件，使用**峰值法**进行站立相检测，并带有 **COP 斜率校正**。

除 `convert_c3d` 的参数外，还有以下额外参数：

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `cop_middle_ratio` | `float` | `0.3` | COP 斜率拟合的中间段比例（30%-70%） |
| `cop_rate_multiplier` | `float` | `2.0` | 异常值检测的比率乘数 |
| `cop_jump_threshold` | `float` | `0.03` | 帧间 COP 跳变阈值（米） |

**返回值：** `dict`，与 `convert_c3d` 相同的键，`stance_info` 中额外包含：
- `peak_frame` (`int`): 峰值垂直力的帧索引
- `peak_value` (`float`): 峰值垂直力大小（N）
- `copx_slope` (`float`): 拟合的 COPx 斜率
- `copz_slope` (`float`): 拟合的 COPz 斜率
- `copx_outlier_count` (`int`): COPx 异常帧校正数量
- `copz_outlier_count` (`int`): COPz 异常帧校正数量

---

#### `convert_c3d6(c3d_path, output_dir, ...)`

将 **6 通道**测力台 C3D 转换为 OpenSim 文件。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `c3d_path` | `str` | 必填 | 含 6 通道力数据的 `.c3d` 文件路径 |
| `output_dir` | `str` | 必填 | 输出文件目录 |
| `marker_cutoff` | `float` | `6.0` | 标记点低通滤波截止频率（Hz） |
| `force_cutoff` | `float` | `50.0` | 力数据低通滤波截止频率（Hz） |
| `stance_threshold` | `float` | `30.0` | 站立相检测的垂直力阈值（N） |
| `stance_pad_frames` | `int` | `25` | 站立相前后的填充帧数 |
| `verbose` | `bool` | `True` | 是否打印进度信息 |

**返回值：** `dict`，与 `convert_c3d` 相同的键。

### 内部辅助函数

这些是转换器内部使用的私有函数：

| 函数 | 说明 |
|------|------|
| `_read_c3d(c3d_path)` | 解析 C3D 文件，提取标记点/力数据、测力台几何参数 |
| `_resample_analog(...)` | 将模拟信号重采样到标记点帧率 |
| `_compute_type2_8ch(...)` | 从 8 通道 Kistler 原始数据计算 Type 2 力数据 |
| `_compute_type2_6ch(...)` | 从 6 通道原始数据计算 Type 2 力数据 |
| `_transform_forces(...)` | 力数据坐标变换：测力台局部 -> 实验室 -> OpenSim |
| `_rotate_markers(...)` | 将标记点坐标旋转到 OpenSim 坐标系 |
| `_filter_markers(...)` | 对标记点应用巴特沃斯低通滤波 |
| `_build_structured_plates_v1(...)` | 构建带 COP 阈值过滤的结构化测力台数据 |
| `_detect_active_plate(...)` | 检测垂直力最大的测力台 |
| `_cut_and_zero_padding(...)` | 将数据切割到站立相范围，站立相外置零 |
| `_markers_to_flat(...)` | 将标记点从 (3, N, T) 重塑为 (T, N*3)，用于 .trc 输出 |
| `_build_mot_plates(...)` | 将结构化测力台数据转换为 `write_mot` 所需格式 |
| `_write_csv_record(...)` | 将站立相信息追加到 CSV 记录文件 |

---

## 模块：`C3DExport.io`

OpenSim 兼容输出文件的 I/O 函数。

### 函数

#### `write_trc(filepath, marker_data, marker_labels, frame_rate, units='mm', orig_start_frame=1)`

将标记点轨迹数据写入 OpenSim `.trc` 文件。

| 参数 | 类型 | 说明 |
|------|------|------|
| `filepath` | `str` | 输出文件路径 |
| `marker_data` | `ndarray (n_frames, n_markers*3)` | 标记点 XYZ 坐标，交错排列 |
| `marker_labels` | `list[str]` | 标记点名称（如 `['R.ASIS', 'L.ASIS', ...]`） |
| `frame_rate` | `float` | 采样率（Hz） |
| `units` | `str` | `'mm'` 或 `'m'`（默认：`'mm'`） |
| `orig_start_frame` | `int` | 原始起始帧号（默认：`1`） |

`.trc` 文件格式包含：
- 文件元数据头（数据率、帧数、标记点数、单位）
- 列标题（标记点名称和轴标签 `X1, Y1, Z1, X2, Y2, Z2, ...`）
- 数据行（帧号、时间戳、每个标记点的 XYZ 坐标）

---

#### `write_mot(filepath, force_data_per_plate, n_plates, frame_rate, filename=None)`

将地面反作用力数据写入 OpenSim `.mot` 文件。

| 参数 | 类型 | 说明 |
|------|------|------|
| `filepath` | `str` | 输出文件路径 |
| `force_data_per_plate` | `list[dict]` | 每个测力台的数据，每个 dict 包含键：`force` (N,3), `cop` (N,3), `torque` (N,3) |
| `n_plates` | `int` | 测力台数量 |
| `frame_rate` | `float` | 采样率（Hz） |
| `filename` | `str` 或 `None` | 头部文件名（默认：文件路径的基本名称） |

`.mot` 文件列布局：
```
time
ground_force_vx  ground_force_vy  ground_force_vz    (FP1)
ground_force_px  ground_force_py  ground_force_pz    (FP1 压力中心)
1_ground_force_vx ...                                (FP2，如果存在)
ground_torque_x  ground_torque_y  ground_torque_z    (FP1)
1_ground_torque_x ...                                (FP2，如果存在)
```

---

## 模块：`C3DExport.utils`

核心工具函数，分为五个类别。

### 旋转矩阵

#### `rotation_matrix(axis, angle_deg)`

生成 3x3 旋转矩阵（右手定则）。

| 参数 | 类型 | 说明 |
|------|------|------|
| `axis` | `str` | 旋转轴：`'X'`、`'Y'` 或 `'Z'` |
| `angle_deg` | `float` | 旋转角度（度） |

**返回值：** `ndarray (3, 3)` - 旋转矩阵

---

#### `chain_rotations(*rotations)`

链接多个旋转（从左到右依次应用）。

| 参数 | 类型 | 说明 |
|------|------|------|
| `*rotations` | `tuples` | `(axis, angle_deg)` 元组序列 |

**返回值：** `ndarray (3, 3)` - 组合旋转矩阵

**示例：**
```python
R = chain_rotations(('X', -90), ('Y', 90))
```

---

#### `apply_rotation(data_3xN, R)`

将 3x3 旋转矩阵应用于 (3, N) 数据。

| 参数 | 类型 | 说明 |
|------|------|------|
| `data_3xN` | `ndarray (3, N)` | 输入数据（3 行 x N 列） |
| `R` | `ndarray (3, 3)` | 旋转矩阵 |

**返回值：** `ndarray (3, N)` - 旋转后的数据

---

### Kistler 测力台计算

#### `compute_kistler_channel8(channels_8, a, b, az0)`

从 **8 通道 Kistler Type 3** 原始数据计算 Type 2 测力台数据。

| 参数 | 类型 | 说明 |
|------|------|------|
| `channels_8` | `ndarray (8, N)` | 原始通道：`[fx12, fx34, fy14, fy23, fz1, fz2, fz3, fz4]` |
| `a` | `float` | 局部 Y 方向传感器偏移量（前后方向），mm |
| `b` | `float` | 局部 X 方向传感器偏移量（内外方向），mm |
| `az0` | `float` | 顶板偏移量，mm（通常为负值） |

**返回值：** `dict`，包含以下键：
- `Fx`, `Fy`, `Fz` (`ndarray`): 地面反作用力（N）
- `ax`, `ay` (`ndarray`): 压力中心距测力台中心的坐标（mm）
- `Tz` (`ndarray`): 自由垂直力矩（N*mm）

**算法：**
1. 汇总原始通道得到 Fx_raw, Fy_raw, Fz_raw
2. 根据传感器几何参数计算力矩 Mx, My, Mz
3. 使用 `az0` 偏移量将力矩转移到板面
4. 计算压力中心：`ax = -Myp/Fz`，`ay = Mxp/Fz`（仅当 |Fz| >= 20N 时有效）
5. 计算自由垂直力矩：`Tz = Mz - Fy*ax + Fx*ay`
6. 取反以获得地面反作用力惯例

---

#### `compute_kistler_channel6(channels_6, az0)`

从 **6 通道**测力台数据计算压力中心和自由力矩。

| 参数 | 类型 | 说明 |
|------|------|------|
| `channels_6` | `ndarray (6, N)` | 原始通道：`[Fx, Fy, Fz, Mx, My, Mz]` |
| `az0` | `float` | 顶板偏移量，mm |

**返回值：** `dict`，与 `compute_kistler_channel8` 相同的键。

---

### 坐标变换

#### `plate_local_to_lab(type2, corners)`

将力数据从 Kistler 测力台局部坐标系转换为实验室全局坐标系。

| 坐标系 | X | Y | Z |
|--------|---|---|---|
| Kistler 局部 | 右+ | 后+ | 下+ |
| 实验室全局 | 前+ | 左+ | 上+ |

| 参数 | 类型 | 说明 |
|------|------|------|
| `type2` | `dict` | `compute_kistler_channel8` 或 `compute_kistler_channel6` 的输出 |
| `corners` | `ndarray (3, 4)` | 来自 C3D `FORCE_PLATFORM:CORNERS` 的测力台角点 |

**返回值：** `dict`，包含键：`Fx`, `Fy`, `Fz` (N), `COPx`, `COPy`, `COPz` (mm), `Tz` (N*mm)

---

#### `lab_to_opensim_force(data_lab)`

将力数据从实验室全局坐标系转换为 OpenSim 坐标系。

| 坐标系 | X | Y | Z |
|--------|---|---|---|
| 实验室全局 | 前 | 左 | 上 |
| OpenSim | 前 | 上 | 右 |

| 参数 | 类型 | 说明 |
|------|------|------|
| `data_lab` | `dict` | 实验室坐标系下的力数据 |

**返回值：** `dict`，相同键，值为 OpenSim 坐标系下的值。

---

### 滤波与重采样

#### `butter_lowpass_filter(data, cutoff, fs, order=4)`

4 阶巴特沃斯低通滤波器（零相位，通过 `filtfilt` 实现）。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `data` | `ndarray (N,)` | 必填 | 输入信号 |
| `cutoff` | `float` | 必填 | 截止频率（Hz） |
| `fs` | `float` | 必填 | 采样频率（Hz） |
| `order` | `int` | `4` | 滤波器阶数 |

**返回值：** `ndarray (N,)` - 滤波后的信号

---

#### `resample_to_target_rate(data_1d, src_rate, tgt_rate)`

使用多相重采样将 1D 信号从源采样率重采样到目标采样率。

| 参数 | 类型 | 说明 |
|------|------|------|
| `data_1d` | `ndarray (N,)` | 输入信号 |
| `src_rate` | `float` | 源采样率（Hz） |
| `tgt_rate` | `float` | 目标采样率（Hz） |

**返回值：** `ndarray (M,)` - 重采样后的信号

---

### 站立相检测

#### `detect_stance_phase(vertical_force, threshold=30.0, pad_frames=25)`

使用**阈值法**检测站立相。

**算法：** 找到 |vertical_force| > threshold 的第一帧和最后一帧，然后添加填充帧。

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `vertical_force` | `ndarray (N,)` | 必填 | 垂直地面反作用力（OpenSim Y 轴，向上为正） |
| `threshold` | `float` | `30.0` | 力阈值（N） |
| `pad_frames` | `int` | `25` | 站立相前后的额外帧数 |

**返回值：** `(start_idx, end_idx, hs_idx, toe_idx)` - 均为 0 索引帧号

**异常：** 若未检测到站立相，抛出 `ValueError`

---

#### `detect_stance_phase_from_peak(vertical_force, threshold=30.0, pad_frames=25)`

使用**峰值法**检测站立相。

**算法：**
1. 找到峰值（最大绝对垂直力）
2. 从峰值向左遍历直到力 < threshold -> 着地（heel-strike）
3. 从峰值向右遍历直到力 < threshold -> 离地（toe-off）
4. 添加填充帧

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `vertical_force` | `ndarray (N,)` | 必填 | 垂直地面反作用力（OpenSim Y 轴，向上为正） |
| `threshold` | `float` | `30.0` | 力阈值（N） |
| `pad_frames` | `int` | `25` | 填充帧数 |

**返回值：** `(start_idx, end_idx, hs_idx, toe_idx, peak_idx, peak_val)` - 均为 0 索引

**异常：** 若峰值 <= threshold，抛出 `ValueError`

---

### COP 异常检测与校正

#### `detect_cop_anomalies(copx, copy, time, force_vertical, threshold=20.0, jump_threshold=0.03)`

检测站立相期间的 COPx/COPy 异常帧。

**双重检测策略：**
1. 力 < threshold 但 COP != 0 -> 异常
2. 帧间 COP 跳变 > jump_threshold（米） -> 异常

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `copx`, `copy` | `ndarray (N,)` | 必填 | COP 坐标（米） |
| `time` | `ndarray (N,)` | 必填 | 时间向量（秒） |
| `force_vertical` | `ndarray (N,)` | 必填 | 垂直力（N） |
| `threshold` | `float` | `20.0` | 力阈值（N） |
| `jump_threshold` | `float` | `0.03` | COP 跳变阈值（米） |

**返回值：**
- `copx_anomaly` (`ndarray (N,) bool`): COPx 异常标记
- `copy_anomaly` (`ndarray (N,) bool`): COPy 异常标记
- `info` (`dict`): 异常统计信息

---

#### `correct_cop_slope(copx, copy, time, force_vertical, threshold=20.0, middle_ratio=0.3, rate_multiplier=2.0)`

使用中间段斜率拟合校正 COP 异常值。

**算法：**
1. 确定站立相（力 >= threshold）
2. 对站立相中间段（如 30%-70%）进行线性拟合
3. 从中间段向外遍历，若帧间差值 > rate_multiplier * |k*dt|，则用斜率预测值替换
4. 非站立相 COP 值置零

| 参数 | 类型 | 默认值 | 说明 |
|------|------|--------|------|
| `copx`, `copy` | `ndarray (N,)` | 必填 | COP 坐标（米） |
| `time` | `ndarray (N,)` | 必填 | 时间向量（秒） |
| `force_vertical` | `ndarray (N,)` | 必填 | 垂直力（N） |
| `threshold` | `float` | `20.0` | 力阈值（N） |
| `middle_ratio` | `float` | `0.3` | 中间段比例 |
| `rate_multiplier` | `float` | `2.0` | 异常值比率乘数 |

**返回值：**
- `copx_corrected` (`ndarray`): 校正后的 COPx
- `copy_corrected` (`ndarray`): 校正后的 COPy
- `info` (`dict`): 斜率信息和异常帧计数

---

## 输出文件格式

### `.trc` 文件（标记点轨迹）

制表符分隔文件，包含：
- **第 1 行：** 文件类型标识（`PathFileType 4 (X/Y/Z)`）
- **第 2 行：** 头部字段名称
- **第 3 行：** 头部值（数据率、帧数、标记点数、单位等）
- **第 4 行：** 标记点名称（每个名称跨 3 列）
- **第 5 行：** 轴子标题（`X1, Y1, Z1, X2, Y2, Z2, ...`）
- **第 6 行：** 空行
- **数据行：** `帧号  时间  X1  Y1  Z1  X2  Y2  Z2  ...`

### `.mot` 文件（地面反作用力）

制表符分隔文件，包含：
- **头部：** 文件名、版本、nRows、nColumns、inDegrees、endheader
- **列标签：** `time`, `ground_force_vx/vy/vz`, `ground_force_px/py/pz`, `ground_torque_x/y/z`
- **数据行：** 力（N）、压力中心（m）和力矩（N*m）的时间序列

### `cut_records.csv`（站立相元数据）

记录每个处理试次的站立相信息的 CSV 文件：
- `trial`: 试次名称
- `heel_strike_frame`: 着地帧号
- `toe_off_frame`: 离地帧号
- `cut_start_frame`, `cut_end_frame`: 切割范围
- `total_frames`: 切割区域内的帧数

---

## 依赖项

| 包 | 版本 | 用途 |
|----|------|------|
| `numpy` | >= 1.20 | 数组运算、线性代数 |
| `scipy` | >= 1.7 | 巴特沃斯滤波器、多相重采样 |
| `ezc3d` | >= 1.5 | C3D 文件解析 |
| `pandas` | >= 1.3 | CSV 记录写入 |

安装所有依赖：

```bash
pip install numpy scipy ezc3d pandas
```

---

## 许可证

本项目仅供研究和教育用途。
