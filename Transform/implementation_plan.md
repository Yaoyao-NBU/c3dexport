# C3D → OpenSim 文件转换器

将 C3D 实验数据转换为 OpenSim 可用的 `.trc`（标记点）和 `.mot`（地面反力）文件，并将坐标系从实验室坐标系转换为 OpenSim 的 Y 轴垂直坐标系。

## User Review Required

> [!IMPORTANT]
> **坐标系映射**：根据参考代码分析，实验室坐标系为「X=前, Y=左, Z=上」，OpenSim 坐标系为「X=前, Y=上, Z=右」。映射关系：
> - OpenSim X = Lab X（前进方向）
> - OpenSim Y = Lab Z（垂直向上）
> - OpenSim Z = -Lab Y（向右）
>
> 请确认这个映射是否正确。

> [!IMPORTANT]
> **力台数量**：从模板 `.mot` 文件中看到有两个力台（`ground_force_*` 和 `1_ground_force_*`），共 19 列（time + 2×9）。请确认 C3D 文件中是否确实有 2 个力台。

> [!WARNING]
> **Kistler 力台坐标系→实验室坐标系映射**：从 `extract_kistler_global.py` 中看到，力台本地坐标系的映射关系是：
> - Kistler Y (posterior+) → Global −X
> - Kistler X (right+) → Global −Y
> - Kistler Z (down+) → Global −Z
>
> 因此力（地面反力 = 反向）：GRF_Fx = Kistler Fy, GRF_Fy = Kistler Fx, GRF_Fz = Kistler Fz
> COP 需要用 CORNERS 偏移。请确认此映射是否适合你的实验室设置。

## Proposed Changes

项目将在 `E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Transform\` 目录下创建一个完整的 Python 程序。

---

### 核心模块

#### [NEW] [c3d_to_opensim.py](file:///E:/Python_Learn/CarbonShoes_DataProcess/C3D_Data_Process/Transform/c3d_to_opensim.py)

主程序，完整的处理流水线：

**Step 1 — 读取 C3D 文件**
- 使用 `ezc3d` 读取 C3D 文件
- 提取 Marker 点数据 `(4, n_markers, n_frames)` → 取 `(3, n_markers, n_frames)` (XYZ)
- 提取力台模拟数据 `(1, n_channels, n_analog_frames)`
- 提取采样频率：`point_rate`（标记点）和 `analog_rate`（力台）
- 提取单位信息（Marker 通常为 mm）
- 提取力台参数（ORIGIN: a, b, az0）

**Step 2 — Kistler Type3 → Type2 计算**
- 从 8 通道原始传感器数据计算出 6 分量 + 自由力矩（Type 2）
- 公式来源于 `extract_kistler_local.py` 中的 `compute_kistler_type3()` 函数
- 输出：`Fx, Fy, Fz`（力）, `ax, ay`（COP，相对板中心，单位 mm）, `Tz`（自由力矩，N·mm）
- 封装为 `type2_force` 数组，每个力台一组

**Step 3 — 力台坐标系 → 实验室坐标系旋转**
- 将 Type2 力台数据从 Kistler 本地坐标系转换到实验室全局坐标系
- Kistler 本地：X=右, Y=后, Z=下
- 实验室全局：X=前, Y=左, Z=上
- 力（地面反力）：GRF_Fx = Fy_kistler, GRF_Fy = Fx_kistler, GRF_Fz = Fz_kistler
- COP：使用 CORNERS 平均值偏移到全局坐标（从板中心到实验室原点）

**Step 4 — 实验室坐标系 → OpenSim 坐标系旋转**
- 对 Marker 数据和力台数据统一进行第二次旋转
- Lab→OpenSim 映射：X→X, Z→Y（垂直）, -Y→Z
- 等价于绕 X 轴旋转 -90°
- 使用 3×3 旋转矩阵实现

**Step 5 — 力台数组重构**
- 为每个力台构建 9 列格式：
  ```
  ground_force_vx, ground_force_vy, ground_force_vz   (力 XYZ)
  ground_force_px, ground_force_py, ground_force_pz    (COP XYZ)
  ground_torque_x, ground_torque_y, ground_torque_z    (力矩 XYZ)
  ```
- 自由力矩 X 和 Z 填充为零（OpenSim 坐标系中 Y 为垂直方向，自由力矩仅 Y 有值）
- COP 垂直方向（Y）填为零

**Step 6 — 单位转换**
- COP：mm → m（÷1000）
- 力矩：N·mm → N·m（÷1000）
- 力：保持 N

**Step 7 — 重采样与滤波**
- 力台数据从 `analog_rate` 降采样至 `point_rate`（Marker 频率）
- 滤波器：4 阶巴特沃斯低通滤波器（`scipy.signal.butter` + `filtfilt`）
- Marker 数据：截止频率 6 Hz
- 力台数据：截止频率 50 Hz（滤波后再降采样）

**Step 8 — Stance 阶段截取**
- 根据力台垂直力 (Y 方向在 OpenSim 中) 判定 Stance 阶段
- 阈值 30 N：大于 30N = 着地（Heel Strike），小于 30N = 离地（Toe Off）
- 向前补 25 帧，向后补 25 帧（用于深度学习上下文）
- 输出截取帧的起止信息到 CSV
- 同步截取 Marker 数据

**Step 9 — 输出文件封装**
- **`.trc` 文件**：基于 `Maker.trc` 模板
  - Header 6 行 + 空行 + 数据
  - 填写 DataRate, CameraRate, NumFrames, NumMarkers, Units, OrigDataRate 等
  - Marker 名称行：按 C3D 中的标签顺序排列
  - 数据行：Frame#, Time, X1, Y1, Z1, X2, Y2, Z2, ...

- **`.mot` 文件**：基于 `Grf.mot` 模板
  - Header 7 行（filename, version, nRows, nColumns, inDegrees, endheader, labels）
  - 数据行：time, force/COP/torque for each plate

---

### 工具函数模块

#### [NEW] [transform_utils.py](file:///E:/Python_Learn/CarbonShoes_DataProcess/C3D_Data_Process/Transform/transform_utils.py)

公用工具函数：

1. **`rotation_matrix(axis, angle_deg)`** — 生成 3×3 旋转矩阵
2. **`apply_rotation(data_3xN, R)`** — 对 (3, N) 数组应用旋转矩阵
3. **`compute_kistler_type2(channels_8, a, b, az0)`** — Kistler Type3→Type2 计算
4. **`plate_local_to_lab(type2_data, corners)`** — 力台本地→实验室坐标系转换
5. **`lab_to_opensim(data_3xN)`** — 实验室→OpenSim 坐标系转换
6. **`butter_lowpass_filter(data, cutoff, fs, order=4)`** — 巴特沃斯低通滤波器
7. **`resample_to_target_rate(data, src_rate, tgt_rate)`** — 重采样（使用 scipy.signal.resample）
8. **`detect_stance_phase(vertical_force, threshold, pad_frames)`** — 检测 Stance 阶段
9. **`write_trc(filepath, markers, labels, rate, units)`** — 写 .trc 文件
10. **`write_mot(filepath, force_data, labels, rate, filename)`** — 写 .mot 文件

---

### 运行入口

#### [NEW] [run_demo.py](file:///E:/Python_Learn/CarbonShoes_DataProcess/C3D_Data_Process/Transform/run_demo.py)

Demo 脚本：
- 输入：`E:\...\Data\S15T1V11.c3d`
- 输出：`E:\...\Transform\output\S15T1V11.trc` 和 `S15T1V11_grf.mot`
- 截取记录：`E:\...\Transform\output\cut_records.csv`
- 打印关键步骤日志

## Open Questions

> [!IMPORTANT]
> 1. 实验室坐标系是否确认为 X=前, Y=左, Z=上？如有不同请告知。
> 2. C3D 文件中是否确有 2 个力台？如果只有 1 个力台，`.mot` 文件中第二个力台的数据将全部填零。
> 3. Stance 阶段截取：如果文件本身就已经只包含 Stance 阶段（无空台阶段），是否跳过截取直接输出全部数据？

## Verification Plan

### Automated Tests
1. 运行 `run_demo.py`，确认无报错
2. 检查输出 `.trc` 和 `.mot` 文件是否可被 OpenSim 打开（通过检查 header 格式）
3. 验证输出文件的行数和列数与 header 声明一致

### Manual Verification
1. 检查输出的 `.trc` 文件 header 格式与模板 `Maker.trc` 一致
2. 检查输出的 `.mot` 文件 header 格式与模板 `Grf.mot` 一致
3. 检查 Marker 点名称是否与模板中一致（41 个 Marker）
4. 验证坐标系转换后，垂直方向数据在 Y 列（OpenSim 约定）
