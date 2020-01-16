import os
import numpy as np

try:
    path = os.path.dirname(__file__)
except:
    path = "."

gprof_weights = np.loadtxt(os.path.join(path, "gprof_weights.txt"))
gprof_weights /= gprof_weights.sum()
