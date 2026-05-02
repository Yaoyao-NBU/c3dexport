# implement file：   本项目的主要作用就是把C3D的实验数据转换为Opensim可以使用的文件！将C3D的坐标系转换为Opensim 的Y轴垂直的坐标系。
1."使用ezc3d的函数读取C3D文件，将Marker点和力台的信号数据读取出来，并存在np的数组中。同时提取出来对应文件中Marker点的采样频率和力台数据的采样频率。也需要提取出来整个marker的单位和力台数据的单位。"
2.对于取出来的两个数组，先对力台的数组进行力学数据的求解，将基础的ksilter的8通道的type3型的传感器数据计算成6外加自由力矩的Type2型力台数据。并重新封装一个叫做type2_force的数组.
（需要写一个旋转矩阵函数，例如rotate（（x，-90），（y，90））就是代表先绕X轴旋转-90度，再绕Y轴旋转90度这样。也可以这个函数一次只能旋转一次也是可以的）
3.对于type2_force进行一次矩阵转换（2*2旋转矩阵），将力台的坐标系旋转到实验室坐标系。
4.对于包含标记点和Type2_force的力的数组（统一坐标系）进行第二次矩阵转换（2*2旋转矩阵），转换为opensim的坐标系。
5.对于力台的数组，需要重新重构一下；生成这样的数组结构（ground_force_vx	ground_force_vy	ground_force_vz ground_force_px	ground_force_py	ground_force_pz ground_torque_x ground_torque_y ground_torque_z）自由力矩X和Y填充为零，Cop列的垂直方向Y为零。
6.对旋转好的力台文件进行单位的换算，将COp换算成（M单位）单位换力矩算成Nm。
7.将力台文件按照帧率进行重采样，采样为Marker点的频率。并写一个用于滤波的函数，滤波器类型为巴特沃斯四阶低通滤波器，滤波频率可以自选（对于Marker点的数据滤波频率为6hz，力台数据为50Hz）。
8.还需要写一个截取数据长度的函数，主要作用就是仿真整个不太的stance阶段，可以根据力台的Y轴（垂直方向的力）来进行定义，当力大于30N时为落地时刻，当力小于30N的时候为toe Off时刻！为了仿真准确，和后面做深度学习数据更好增加上下文长度，可以手动在落地时刻向前补25帧，在落地后向后手动不25帧！（可以的话导出一个记录灭个截取文件的开始和结束帧数）。在通过力学数据确定截取的帧数长度范围后，需要对marker数据进行同样的截取！
9.开始封装opensim需要的文件MArker数组的数据封装倒.trc文件中，力台的数据封装到.mot文件中。这两个文件都有一个各自的Header需要写好。我会制作好一个母文件，你只需要写入进去并填好对应的一些信息就可以。
    我简单的介绍一下.trc/.mot文件的header.
    1:.trc file:  这个是需要注意的表头，按照点输入maker数据的时候需要注意按照顺序Marker1（X1，y1，z1） Marker2（X1，y2，z2）排好顺序，里面的帧率和marker数量，单位都按照C3d中的单位和数据去填好。同时需要注意记得在这个填好之后空一排再去排入对应的marker的数据。
    PathFileType	4	(X/Y/Z)	subject01_walk1.trc							
    DataRate	CameraRate	NumFrames	NumMarkers	Units	OrigDataRate	OrigDataStartFrame	OrigNumFrames			
    60	60	151	41	mm	60	1	151			
    Frame#	Time	R.ASIS			L.ASIS			V.Sacral		
            X1	Y1	Z1	X2	Y2	Z2	X3	Y3	Z3
    
    2: .mot file
    subject01_walk1_grf.mot	 "这里指的是文件的名字"							
    version=1				"这里指的是版本1，是力台的版本"				
    nRows=1501				"这里指的是帧数（也就是数据有的帧数）"				
    nColumns=19					"这里指有多少列的数据"			
    inDegrees=yes				"默认就好"				
    endheader						"不变"		
    time	ground_force_vx	ground_force_vy	ground_force_vz	ground_force_px	ground_force_py	ground_force_pz	1_ground_force_vx	1_ground_force_vy  "这里指的时间列和对应的Lable名字"
    

数据位置：
1. .mot文件（总）：E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Data\Grf.mot
2. .trc文件（总）：E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Data\Marker.trc
3. 一个C3D数据文件：E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\Data\S15T1V11.c3d
参考代码：
1.计算力台六通道数据的代码：E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\kistler_c3d_extractor\extract_kistler_local.py  / E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\kistler_c3d_extractor\extract_kistler_global.py
2.滤波函数： E:\Python_Learn\CarbonShoes_DataProcess\OPenSIm\Data_ProcessFunction.py
3.截取代码：E:\Python_Learn\CarbonShoes_DataProcess\OPenSIm\Batch_CUT_TrcSto.py 和 E:\Python_Learn\CarbonShoes_DataProcess\OPenSIm\Data_ProcessFunction.py
4.旋转函数：E:\Python_Learn\CarbonShoes_DataProcess\C3D_Data_Process\kistler_c3d_extractor\rotate_lab_coordinates.py
