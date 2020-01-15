from netCDF4 import Dataset
import os
import numpy as np
import matplotlib.pyplot as plt

def create_output_file(path):
    f = Dataset(path, "w", clobber = True)
    f.createDimension("swaths", size = None)
    f.createDimension("scenes", size = None)
    f.createDimension("swath_width", size = 221)
    f.createDimension("date_string", size = 14)
    f.createDimension("x_size", size = 1900)
    f.createDimension("y_size", size = 2200)

    f.createVariable("scene_id", "i8", ("scenes"))
    f.createVariable("start_time", "c", ("swaths", "date_string"))
    f.createVariable("end_time", "c", ("swaths", "date_string"))
    f.createVariable("mask", "i4", ("y_size", "x_size"))

    # Location
    f.createVariable("lon", "f4", ("swaths", "swath_width"))
    f.createVariable("lat", "f4", ("swaths", "swath_width"))

    # GPROF
    g = f.createGroup("gprof")
    g.createVariable("tcwv_index", "i4", ("swaths", "swath_width"))
    g.createVariable("t2m_index", "i4", ("swaths", "swath_width"))
    g.createVariable("st_index", "i4", ("swaths", "swath_width"))
    g.createVariable("surface_precipitation", "f4", ("swaths", "swath_width"))


    # 1c
    g = f.createGroup("1c")
    g.createDimension("channels", size = 13)
    g.createVariable("brightness_temperatures", "f4", ("swaths", "swath_width", "channels"))

    # Opera
    g = f.createGroup("opera")
    g.createVariable("precipitation_5", "f4", ("swaths", "swath_width",))
    g.createVariable("precipitation_10", "f4", ("swaths", "swath_width",))

    # Combined
    g = f.createGroup("combined")
    g.createDimension("swaths_combined", None)
    g.createDimension("swath_width_combined", 49)
    g.createVariable("scene_id", "i4", ("swaths_combined",))
    g.createVariable("lat", "f4", ("swaths_combined", "swath_width_combined"))
    g.createVariable("lon", "f4", ("swaths_combined", "swath_width_combined"))
    g.createVariable("surface_precipitation", "f4", ("swaths_combined", "swath_width_combined"))

    return f

