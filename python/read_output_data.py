import numpy as np
import struct

fmt = "d"*100

with open("../data/numpy_output_z.bin", "rb") as f:
    np_z = np.array(struct.unpack_from(fmt, f.read())).reshape([10, 10])

with open("../data/eigen_output_z.bin", "rb") as f:
    eg_z = np.array(struct.unpack_from(fmt, f.read())).reshape([10, 10])
