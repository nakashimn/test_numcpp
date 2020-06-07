import numpy as np
import struct

x = np.random.random([10, 10])
y = np.random.random([10, 10])

fmt = "d"*100
x_bin = struct.pack(fmt, *x.ravel())
y_bin = struct.pack(fmt, *y.ravel())

with open("../data/input_x.bin", "wb") as f:
    f.write(x_bin)

with open("../data/input_y.bin", "wb") as f:
    f.write(y_bin)
