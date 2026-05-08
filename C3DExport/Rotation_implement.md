力台局部坐标系至全局坐标系旋转矩阵 ($R_{plate \to lab}$) 计算规范

本规范定义了通过力台四个角点坐标精确推导三维旋转矩阵的标准化算法流程。此流程包含强制正交化步骤，旨在消除物理设备制造或标定误差导致的非正交畸变，满足运动生物力学数据解析工作流的严谨性要求。

1. 算法输入 (Inputs)

变量名: corners

数据结构: 形状为 $(3, 4)$ 的浮点数矩阵（例如 NumPy ndarray）。

物理定义: 包含 4 个角点在全局实验室坐标系中的三维坐标，按列分别记为 $C_1, C_2, C_3, C_4$。

注: 硬件物理顺序需对齐厂商规范。本推导基于行业通用标准：向量 $C_1 \to C_2$ 定义局部 X 轴，向量 $C_1 \to C_4$ 定义局部 Y 轴。

2. 核心推导步骤 (Operational Flow)

步骤 1：提取基础方向向量 (Extract Basis Vectors)

计算定义力台平面的两个初始方向向量 $\vec{v}_x$ 和 $\vec{v}_y$：


$$\vec{v}_x = C_2 - C_1$$

$$\vec{v}_y = C_4 - C_1$$

步骤2：归一化（变成单位向量）：将向量除以它们的模长，使其长度为 1。$$\hat{i} = \frac{V_x}{||V_x||}, \quad \hat{j} = \frac{V_y}{||V_y||}$$

步骤 3：计算正交法向量 (Orthogonal Normal Vector Calculation)

为确保坐标系的严格正交性，通过叉乘（Cross Product）求得垂直于力台表面的局部 Z 轴单位向量 $\hat{k}$：


$$\vec{v}_z = \vec{v}_x \times \vec{v}_y$$

$$\hat{k} = \frac{\vec{v}_z}{||\vec{v}_z||}$$

步骤 4：强制正交化次轴 (Secondary Axis Orthogonalization)

利用已确认绝对正交的 $\hat{k}$ 和 $\hat{i}$，根据右手定则反向通过叉乘计算出严格正交的局部 Y 轴单位向量 $\hat{j}$：


$$\hat{j} = \hat{k} \times \hat{i}$$

工程说明: 此步骤产生的 $\hat{j}$ 必然为单位向量，且绝对垂直于 X 轴与 Z 轴，有效消除了 $C_1, C_2, C_4$ 构成的非直角物理误差。

步骤 5：组装旋转矩阵 (Matrix Assembly)

将计算得到的三个正交单位向量作为列向量，按顺序拼装为 $3 \times 3$ 的方向余弦矩阵 $R_{plate \to lab}$：


$$R_{plate \to lab} = \begin{bmatrix} \hat{i}_x & \hat{j}_x & \hat{k}_x \\ \hat{i}_y & \hat{j}_y & \hat{k}_y \\ \hat{i}_z & \hat{j}_z & \hat{k}_z \end{bmatrix}$$

3. 附加输出：计算物理中心点 (Center Point Calculation)

同步计算力台的绝对几何中心 $P_{center}$，作为空间仿射变换的原点偏移量（Translation Vector）：


$$P_{center} = \frac{1}{4} \sum_{n=1}^{4} C_n$$

4. 编码实施指引 (Implementation Directives)

Agent 在生成代码时必须遵循以下约束：

依赖库: 必须使用 numpy，特别是 np.cross 进行叉乘，np.linalg.norm 进行模长计算。

容错机制 (Division by Zero Guard): 在执行步骤 2 和步骤 3 的除法归一化前，必须校验向量模长是否极小（例如 $< 1e-6$）。若低于阈值，需抛出明确的数值异常 (ValueError)。

输出格式: 函数需返回一个元组 (R_plate2lab, P_center)，分别为 (3, 3) 的旋转矩阵和 (3,) 的中心点向量。