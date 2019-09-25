from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import os
import warnings
import matplotlib.cbook
import radcos_fio

warnings.filterwarnings("ignore",category=matplotlib.cbook.mplDeprecation)
def healpix_converter(filename):
    data_array = radcos_fio.data()
    # file pixels may be organized as ring or nest. Currently on nest
    ra = (np.pi * 2 / 360) * data_array[:,0]
    dec = (np.pi * 2 / 360) * data_array[:,1]
    flux = data_array[:,2]
    ra[np.where(ra > np.pi)] -= 2 * np.pi
    pixel_refs = hp.pixelfunc.ang2pix(32, ra, dec, nest=False, lonlat=False)
    print pixel_refs

    return

    if nest_or_ring is 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    if nest_or_ring is 'nest':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)

    # plot of Galactic gas with coordinate projection
    min_ra = np.min(ra)
    print min_ra
    max_ra = np.max(ra)
    print max_ra
    mean_ra = np.mean(ra)
    min_dec = np.min(dec)
    print min_dec
    max_dec = np.max(dec)
    print max_dec
    mean_dec = np.mean(dec)

    fig = plt.figure()
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # llc and urc are x & y ranges, and are specific to a location.
    # The latitude and longitude settings are part of basemap.
    m = Basemap(projection='hammer', llcrnrlon=-11, llcrnrlat=-15, urcrnrlon=13.5, urcrnrlat=-37, resolution='h', epsg=5520)
    # m = Basemap(projection='hammer', lon_0=mean_ra, lat_0=mean_dec, llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra, urcrnry=max_dec, resolution='h', epsg=5520)
    x, y = m(ra, dec)
    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    m.drawparallels(np.arange(-90., 120., 10.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0., 420., 10.), labels=[0, 0, 0, 1])
    # creates a scatter plot of the selected data on a globe.
    m.scatter(x, y, 3, marker='o', linewidths=.1, c=data, cmap=plt.cm.coolwarm)
    m.colorbar()
    plt.show(plotfile)
