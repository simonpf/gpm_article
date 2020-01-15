from matplotlib.colors import LogNorm
from matplotlib.cm import Greys
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
from netCDF4 import Dataset
import numpy as np

def plot_colocation(f, index):

    norm = LogNorm(1e-2, 5e1)
    cmap = "magma"

    indices = f["scene_id"][:] == index
    p_opera = f["opera"]["precipitation_5"][indices, :]
    p_gprof = f["gprof"]["surface_precipitation"][indices, :]
    lons = f["lon"][indices, :]
    lats = f["lat"][indices, :]

    proj_opera = ccrs.LambertAzimuthalEqualArea(central_longitude=10,
                                                central_latitude=55,
                                                false_easting=1950000,
                                                false_northing=-2100000)
    proj_pc = ccrs.PlateCarree()

    ll = np.array([-10.434, 31.746])
    ur = np.array([57.81, 67.62])
    ll_t = proj_opera.transform_point(ll[0], ll[1], proj_pc)
    ur_t = proj_opera.transform_point(ur[0], ur[1], proj_pc)


    # OPERA
    plt.figure(figsize = (12, 4), dpi= 200)
    gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 0.05])
    ax = plt.subplot(gs[0], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)

    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.2)
    ax.pcolormesh(lons, lats, p_opera, norm = norm, cmap = cmap, transform = proj_pc)
    ax.set_title("(a) Opera ground radar")

    # GPROF
    ax = plt.subplot(gs[1], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)

    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.5)
    img = ax.pcolormesh(lons, lats, p_gprof, norm = norm, cmap = cmap, transform = proj_pc)

    i = np.where(p_gprof >= 0)[1][0]
    ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.5)
    i = np.where(p_gprof >= 0)[1][-1]
    ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.5)
    ax.set_title("(b) GPROF GMI")


    ax = plt.subplot(gs[2], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)

    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.5)

    g = f["combined"]
    indices = g["scene_id"][:] == index

    if np.any(indices):
        p_cmb = g["surface_precipitation"][indices, :]
        lons = g["lon"][indices, :]
        lats = g["lat"][indices, :]

        ax.pcolormesh(lons, lats, p_cmb, norm = norm, cmap = cmap, transform = proj_pc)

        i = np.where(p_cmb >= 0)[1][0]
        ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.5)
        i = np.where(p_cmb >= 0)[1][-1]
        ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.5)

    ax.set_title("(c) GPM Combined")

    ax = plt.subplot(gs[3])
    plt.colorbar(img, cax = ax, label = r"Rainfall rate $[mm\ h^{-1}]$")
    plt.tight_layout()

