import ezc3d as ed
import numpy as np

# Create a c3d file 
c3d = ed.c3d()

# Set the parameters and data
c3d['parameters']['POINT']['UNITS']['value'] = ['mm']
c3d['parameters']['POINT']['RATE']['value'] = [100]
c3d['parameters']['POINT']['LABELS']['value'] = ['jiayin01', 'jiayin02', 'jiayin03']
c3d['data']['points'] = np.random.rand(4, 3, 1500)  # 3D points, 1 point, 1 frame
c3d['data']['points'][1, :, :] = 2  # Set the last point to a specific value
c3d['data']['points'][2, :, :] = 3  # Set the first point in the first frame to a specific value

# Set the analog data
c3d["parameters"]["ANALOG"]["RATE"]["value"] = [1000]
c3d["parameters"]["ANALOG"]["LABELS"]["value"] = ('jiayinFx','jiayinFy','jiayinFz','jiayingMx','jiayingMy','jiayingMz')
c3d["data"]["analogs"] = np.random.rand(1,6, 15000)  # 6 analog channels, 15000 frames
c3d["data"]["analogs"][0, :, :] = 4  # Set the first analog channel to a specific value
c3d["data"]["analogs"][0, 1, :] = 5  # Set the second analog channel to a specific value
c3d["data"]["analogs"][0, 2, :] = 6
c3d["data"]["analogs"][0, 3, :] = 7
c3d["data"]["analogs"][0, 4, :] = 8
c3d["data"]["analogs"][0, 5, :] = 9


# Write the c3d file
c3d.add_parameter("POINT", "newParam",[1,2,3])
c3d.add_parameter("newGroup", "newParam",["MyParam1","MyParam2"])
c3d.write("test_ezc3d_floder/jiaying_is_so_cute.c3d")    