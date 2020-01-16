"""
GPM-Opera colocations

Script to match up GPM and Opera ground-radar observations.
"""
try:
    from IPython import get_ipython
    ip = get_ipython()
    if not ip is None:
        ip.magic("%load_ext autoreload")
        ip.magic("%autoreload 2")
except:
    pass

from calendar import monthrange
from datetime import datetime
from cloud_colocations.colocations.products import set_cache
from cloud_colocations.colocations.formats import GPMGMI1C, GPM_2B_CMB, GPM_2A_GPROF,\
    OperaRainfall
from cloud_colocations.shapes import opera
from shapely.geometry import MultiPoint
from colocations import create_output_file

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import argparse

################################################################################
# CLI arguments
################################################################################

desc = "Download and match GPM-Opera colocations for given year and month."
parser = argparse.ArgumentParser(description=desc)
parser.add_argument('year', metavar='year', type=int, nargs=1,
                    help='Year for which to generate the colocations.')
parser.add_argument('month', metavar='month', type=int, nargs=1,
                    help='Month for which to generate the colocations.')
parser.add_argument('output_file', metavar='output_file', type=str, nargs=1,
                    help='Filename of output file.')
parser.add_argument('--cache', metavar='cache_folder', type=str, nargs=1,
                    default="../data", help='Folder to use as cache.')



args = parser.parse_args()
year = args.year[0]
month = args.month[0]
output_file = args.output_file[0]
cache = args.cache[0]

################################################################################
# Preparations
################################################################################

set_cache(cache)


output_file = create_output_file(output_file)
_, days = monthrange(year, month)

################################################################################
# Iterate over days of the month
################################################################################

scene_id = 0
for day in range(days):
    t = datetime(year, month, day + 1)
    files_gprof = GPM_2A_GPROF.get_files_by_day(t)
    files_cmb = GPM_2B_CMB.get_files_by_day(t)
    files_1c = GPMGMI1C.get_files_by_day(t)

    for f_gprof, f_cmb, f_1c in zip(files_gprof, files_cmb, files_1c):

        #
        # Check if obs over Europe.
        #

        # Get pixels within Europe
        lons = f_gprof.lon
        lats = f_gprof.lat
        coords = np.hstack([lons.reshape(-1, 1), lats.reshape(-1, 1)])
        points = MultiPoint(coords)
        inds = np.array(list(map(opera.contains, points)))
        inds = np.reshape(inds, lons.shape)

        if ~np.any(inds):
            continue

        i_start = np.where(inds)[0][0]
        i_end = np.where(inds)[0][-1]
        ts = f_gprof.get_time(i_start)
        te = f_gprof.get_time(i_end)

        #
        # Interpolate opera data.
        #

        files_opera = OperaRainfall.get_files_in_range(ts,
                                                    te,
                                                    t0_inclusive = True,
                                                    t1_inclusive = True)
        if not all(map(lambda x: getattr(x, "valid"), files_opera)):
            continue
        t = np.array([(of.mean_time - f_gprof.start_time).total_seconds() for of in files_opera])
        x = files_opera[0].x
        y = files_opera[0].y[::-1]

        x_size = files_opera[0].x_size
        y_size = files_opera[0].y_size

        # 5km smoothed
        data_5 = np.zeros((int(x_size), int(y_size), t.size))
        for i, of in enumerate(files_opera):
            data_5[:, :, i] = of.get_smoothed_data(5e3)[0].T[:, ::-1]

        int_5 = sp.interpolate.RegularGridInterpolator((x, y, t),
                                                    data_5,
                                                    bounds_error = False,
                                                    fill_value = np.nan)
        x_i, y_i = files_opera[0].projection(lons[i_start : i_end, :],
                                            lats[i_start : i_end, :])
        t_i = f_gprof.time[i_start : i_end, :]
        data_5_i = int_5((x_i, y_i, t_i))

        # 10km smoothed
        data_10 = np.zeros((int(x_size), int(y_size), t.size))
        for i, of in enumerate(files_opera):
            data_10[:, :, i] = of.get_smoothed_data(5e3)[0].T[:, ::-1]


        int_10 = sp.interpolate.RegularGridInterpolator((x, y, t),
                                                    data_10,
                                                    bounds_error = False,
                                                    fill_value = np.nan)
        x_i, y_i = files_opera[0].projection(lons[i_start : i_end, :],
                                            lats[i_start : i_end, :])
        t_i = f_gprof.time[i_start : i_end, :]
        data_10_i = int_10((x_i, y_i, t_i))

        #
        # Store data.
        #
        m = output_file.dimensions["swaths"].size
        n = m + i_end - i_start

        output_file["start_time"][scene_id] = ts.strftime("%Y%m%d%H%M%S")
        output_file["end_time"][scene_id] = te.strftime("%Y%m%d%H%M%S")

        output_file["scene_id"][m : n] = scene_id
        output_file["lon"][m : n, :] = lons[i_start : i_end, :]
        output_file["lat"][m : n, :] = lats[i_start : i_end, :]

        # gprof
        g = output_file["gprof"]
        g["tcwv_index"][m : n, :] = f_gprof.tcwv_index[i_start : i_end, :]
        g["t2m_index"][m : n, :] = f_gprof.t2m_index[i_start : i_end, :]
        g["st_index"][m : n, :] = f_gprof.st_index[i_start : i_end, :]
        g["surface_precipitation"][m : n, :] = f_gprof.precip[i_start : i_end, :]

        # 1c
        g = output_file["1c"]
        g["brightness_temperatures"][m : n, :, :9] = f_1c.y_s1[i_start : i_end, :]
        g["brightness_temperatures"][m : n, :, 9:] = f_1c.y_s2[i_start : i_end, :]

        # opera
        g = output_file["opera"]
        g["precipitation_5"][m : n, :] = data_5_i
        g["precipitation_10"][m : n, :] = data_10_i


        # combined
        # Get pixels within Europe
        lons = f_cmb.lon
        lats = f_cmb.lat
        coords = np.hstack([lons.reshape(-1, 1), lats.reshape(-1, 1)])
        points = MultiPoint(coords)
        inds = np.array(list(map(opera.contains, points)))
        inds = np.reshape(inds, lons.shape)

        if np.any(inds):
            i_start = np.where(inds)[0][0]
            i_end = np.where(inds)[0][-1]

            g = output_file["combined"]
            m = g.dimensions["swaths_combined"].size
            n = m + i_end - i_start
            g["scene_id"][m : n] = scene_id
            g["lat"][m : n, :] = f_cmb.lat[i_start : i_end, :]
            g["lon"][m : n, :] = f_cmb.lon[i_start : i_end, :]
            g["surface_precipitation"][m : n, :] = f_cmb.precip[i_start : i_end, :]

        scene_id = scene_id + 1
