import h5py
import numpy as np

data = {}

for field in ["Density", "Vx", "Vy", "Vz"]:
    f = h5py.File("mhd512li.%s" % field)
    data[field] = f["/mhd512li.%s" % field][:]
    f.close()

np.random.seed(0x4d3d3d3)

data["Px"] = data["Density"] * data["Vx"]
data["Py"] = data["Density"] * data["Vy"]
data["Pz"] = data["Density"] * data["Vz"]
data["Energy"] = np.random.random((512, 512, 512))

with open("UM_IC", "wb") as f:
    for field in ["Density", "Px", "Py", "Pz", "Energy"]:
        print("Writing %s" % field)
        np.float32(data[field]).tofile(f)

