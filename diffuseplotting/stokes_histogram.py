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
import healpix_functions
import fits_data_extraction
import stokes_math


def plotstokesQU(obsID, x_stokes, y_stokes, ra, dec,
                 new_histogram_filename=None, directory=None,
                 save_show='show', projection='ortho'):

    # manipulate histogram for number of bins, color scheme.
    # histogram is of 'angle' of polarization between Stokes Q and Stokes U.
    plt.hist2d(x_stokes, y_stokes, bins=150, norm=LogNorm(), cmap=plt.cm.plasma)
    plt.colorbar()
    plt.title(('{} histogram'.format(obsID)))
    #plt.title('ObsID {} polarization angle vs magnitude'.format([int(s) for s in re.findall('\d+', new_histogram_filename)]))
    #plt.title(os.path.split(new_histogram_filename)[-1])

    # either show the histogram, or save it to a location.
    if save_show == 'show':
        plt.show(new_histogram_filename)
    elif save_show == 'none':
        print "histogram not displayed"
    elif save_show == 'save':
        plt.savefig(new_histogram_filename)
        print 'saved polarization histogram to ' + new_histogram_filename
    else:
        raise ValueError('save_show needs to be equal to "save" or "show" to save or show the image.')

    # x_stokes and y_stokes are the vector components necessary to create a drapery plot.
    return ra, dec, x_stokes, y_stokes
