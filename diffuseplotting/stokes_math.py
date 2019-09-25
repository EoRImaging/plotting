from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from mpl_toolkits.basemap import Basemap
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
import os
import re
"""
Takes Stokes Q, U, and I signals in ra/dec coordinates and computes x_stokes and y_stokes.
Theta is the polarization angle of the linearly polarized portion of the data.
x_stokes and y_stokes are the vector components of theta projected onto an x-y axis.
"""
def math(I, Q, U):
    # Finding x and y from Stokes parameters U and Q
    K = np.sqrt(U**2 + Q**2)
    Kratio = K / I
    theta = np.arctan(U/(np.sqrt(Q**2 + U**2)+Q))
    theta[np.where(theta < 0)] = theta[np.where(theta < 0)] + np.pi
    x_stokes = np.cos(theta) * K
    y_stokes = np.sin(theta) * K
    return K, theta, x_stokes, y_stokes
