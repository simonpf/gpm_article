from matplotlib.colors import LogNorm
from matplotlib.cm import Greys
from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
import cartopy
import cartopy.crs as ccrs
import scipy as sp
from scipy.signal import convolve
from netCDF4 import Dataset
import numpy as np

from gprof_article import gprof_weights

def plot_colocation(f, index):
    """
    Plot colocation scene with given index from file.

    Args:
        f(netCDF4.Dataset): NetCDF file handle containing the
            colocations.
        index(int): The scene_id value identifying the scene
            to plot.
    """
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


    ################################################################################
    # OPERA
    ################################################################################

    plt.figure(figsize = (12, 4), dpi= 200)
    gs = GridSpec(1, 4, width_ratios = [1.0, 1.0, 1.0, 0.05])
    ax = plt.subplot(gs[0], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)

    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.2)
    ax.set_title("(a) Opera ground radar", loc = "left")

    ax.pcolormesh(lons, lats, p_opera, norm = norm, cmap = cmap, transform = proj_pc)

    #
    # Boundary
    #

    mask_opera = np.logical_and(p_gprof >= 0, np.isfinite(p_opera)).astype(np.float)
    mask_opera[0, :] = 0.0
    mask_opera[-1, :] = 0.0
    ax.contour(lons, lats, mask_opera, levels = [0.0, 1.0], colors = "k",
               linewidths = 0.8, transform = proj_pc)

    ################################################################################
    # GPROF
    ################################################################################

    ax = plt.subplot(gs[1], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)
    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.5)
    ax.set_title("(b) GPROF GMI", loc = "left")

    img = ax.pcolormesh(lons, lats, p_gprof, norm = norm, cmap = cmap, transform = proj_pc)

    #
    # Boundary
    #

    i = np.where(p_gprof >= 0)[1][0]
    ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8)
    i = np.where(p_gprof >= 0)[1][-1]
    ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8)
    ax_gprof = ax

    ################################################################################
    # Combined
    ################################################################################

    ax = plt.subplot(gs[2], projection = proj_opera)
    ax.set_extent([ll_t[0], ur_t[0], ll_t[1], ur_t[1]], crs = proj_opera)
    ax.stock_img()
    ax.coastlines(resolution="10m", linewidth=0.5)
    ax.set_title("(c) GPM Combined", loc = "left")

    g = f["combined"]
    indices = g["scene_id"][:] == index

    if np.any(indices):
        p_cmb = g["surface_precipitation"][indices, :]
        p_cmb_s = convolve(p_cmb, gprof_weights, "same")
        p_cmb_s[p_cmb_s < 1e-2] = np.nan
        lons = g["lon"][indices, :]
        lats = g["lat"][indices, :]

        ax.pcolormesh(lons, lats, p_cmb_s, norm = norm, cmap = cmap, transform = proj_pc)

        #
        # Boundary
        #

        i = np.where(p_cmb >= 0)[1][0]
        ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8)
        ax_gprof.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8, ls = "--")
        i = np.where(p_cmb >= 0)[1][-1]
        ax.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8)
        ax_gprof.plot(lons[:, i], lats[:, i], c = "k", transform = proj_pc, lw = 0.8, ls = "--")

    ################################################################################
    # Colorbar
    ################################################################################

    ax = plt.subplot(gs[3])
    plt.colorbar(img, cax = ax, label = r"Rainfall rate $[mm\ h^{-1}]$")

    plt.tight_layout()

def scatter_plot(f, index):
    """
    Plot colocation scene with given index from file.

    Args:
        f(netCDF4.Dataset): NetCDF file handle containing the
            colocations.
        index(int): The scene_id value identifying the scene
            to plot.
    """
    norm = LogNorm(1e-2, 5e1)
    cmap = "magma"

    indices = f["scene_id"][:] == index
    p_opera = f["opera"]["precipitation_5"][indices, :]
    p_gprof = f["gprof"]["surface_precipitation"][indices, :]

    mask = np.logical_and(p_gprof >= 0.0, np.isfinite(p_opera))

    p_opera = p_opera[mask]
    p_gprof = p_gprof[mask]

    f = plt.figure(figsize = (5, 4))
    gs = GridSpec(1, 2, width_ratios = [1.0, 0.05])
    ax = plt.subplot(gs[0])

    bins = np.logspace(-2, 2, 21)
    img, xx, yy = np.histogram2d(p_gprof, p_opera, bins = bins)
    img = ax.pcolormesh(xx, yy, img.T)
    ax.plot(xx, xx, c = "white", ls = "--")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_aspect(1.0)
    ax.set_xlabel("Opera rainfall [$mm\ h^{-1}$]")
    ax.set_ylabel("GPROF Surface precip. [$mm\ h^{-1}$]")

    ax = plt.subplot(gs[1])
    plt.colorbar(img, cax = ax, label = "Frequency")
