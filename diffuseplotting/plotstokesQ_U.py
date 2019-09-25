# Use Stokes Q and U files(linear polarizations) to plot a histogramatic plot which shows the polarization angle
# and magnitude of polarization, where color shows density of points.
# Create a map, plotted on a spherical, RA-dec grid, which plots a third variable
# of the data (theta, K, I) on a color scale.
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
import re
import plothealpix_map
import healpix_functions
import fits_data_extraction
import stokes_math

"""
def plotstokesQ_U(ra, dec, signal_Q, signal_U, signal_I, map_var='theta',
                  histogram_plotfile_base=None, map_plotfile_base=None, plotdirectory=None,
                  save_show='show', file_extension='.png', projection='ortho'):
"""
def plotstokesQ_U(signal_Q, signal_U, signal_I, ra, dec, histogram_title, map_var='theta',
                  directory=None,
                  save_show='show', file_extension='.png', projection='ortho'):
    #signal_Q, signal_U, signal_I, ra, dec =\
    #    fits_data_extraction.fits_extractor(filename_Q, filename_U, filename_I)
    """
    # we set information from a fits file. theta is in radians.
    # x_stokes and y_stokes are the x and y components of points found using theta and magnitude.
    # Kratio[np.where(I < K)] = 1

    outliers_K = [K[np.all([K >= .1], axis=0)]]
    outliers_theta = [theta[np.all([K >= .1], axis=0)]]
    x_stokes_out = np.cos(outliers_theta) * outliers_K
    y_stokes_out = np.sin(outliers_theta) * outliers_K

    # setting file and path names
    filename_Q = os.path.split(filename_Q)[-1]
    filename_U = os.path.split(filename_U)[-1]
    obsID = re.findall('\d+',filename_Q)
    obsID = 'ObsID_' + ''.join(obsID)
    if histogram_plotfile_base is None:
        histogram_plotfile = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        histogram_plotfile_base = obsID + '_hist_angle_vs_magnitude'
    if map_plotfile_base is None:
        map_plotfile = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        map_plotfile_base = obsID
    if plotdirectory is None:
        plotdirectory = os.getcwd()
    if not os.path.isdir(plotdirectory):
        os.makedirs(plotdirectory)

    path_histogram_plotfile_base = os.path.join(plotdirectory, histogram_plotfile_base)
    path_map_plotfile_base = os.path.join(plotdirectory, map_plotfile_base)
    histogram_plotfile = path_histogram_plotfile_base + file_extension
    map_plotfile = path_map_plotfile_base + file_extension
    """
    K, theta, x_stokes, y_stokes = stokes_math.math(signal_I, signal_Q, signal_U)

    # plots selected third variable to map
    if map_var is 'theta':
        map_data_var = theta
    elif map_var is 'K':
        map_data_var = K
    elif map_var is 'x_stokes':
        map_data_var = x_stokes
    elif map_var is 'xval':
        map_data_var = xval
    elif map_var is 'y_stokes':
        map_data_var = y_stokes
    elif map_var is 'yval':
        map_data_var = yval
    elif map_var is 'I':
        map_data_var = I
    elif map_var is 'Q':
        map_data_var = Q
    elif map_var is 'U':
        map_data_var = U
    else:
        raise ValueError("unrecognized map_var")

    plothealpix_map.mapping(ra, dec, map_data_var, map_var, newplotfile_base=histogram_title, projection=projection,
                                            save_show=save_show, file_extension=file_extension)
    # manipulate histogram for number of bins, color scheme.
    # histogram is of 'angle' of polarization between Stokes Q and Stokes U.
    plt.hist2d(x_stokes, y_stokes, bins=150, norm=LogNorm(), cmap=plt.cm.plasma)
    plt.colorbar()
    plt.title(histogram_title)
    plt.title('ObsID {} polarization angle vs magnitude'.format([int(s) for s in re.findall('\d+', histogram_title)]))
    plt.title(os.path.split(histogram_title)[-1])

    # either show the histogram, or save it to a location.
    if save_show == 'show':
        plt.show(histogram_title)
    elif save_show == 'none':
        print "histogram not displayed"
    elif save_show == 'save':
        plt.savefig(histogram_title)
        print 'saved polarization histogram to ' + histogram_title
    else:
        raise ValueError('save_show needs to be equal to "save" or "show" to save or show the image.')

    # x_stokes and y_stokes are the vector components necessary to create a drapery plot.
    return x_stokes, y_stokes
