import numpy as np
import struct

fmt = "d"*100

with open("../data/input_x.bin", "rb") as f:
    x = np.array(struct.unpack_from(fmt, f.read())).reshape([10, 10])

with open("../data/input_y.bin", "rb") as f:
    y = np.array(struct.unpack_from(fmt, f.read())).reshape([10, 10])


np.sum(x.T[0] * y[0])
z = np.dot(x, y)
z_bin = struct.pack(fmt, *z.ravel())

with open("../data/numpy_output_z.bin", "wb") as f:
    f.write(z_bin)
