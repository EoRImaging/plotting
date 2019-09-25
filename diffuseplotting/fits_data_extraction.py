from scipy.io.idl import readsav
from netCDF4 import Dataset, num2date
from astropy.io import fits
import numpy as np
import healpy as hp
import scipy as sp
import matplotlib
from matplotlib.colors import LogNorm
import os
import healpix_functions
"""
This module extracts variable data lists from .fits files.  It has an option to
convert coordinates to RA and dec, or leave as healpix pixels.
"""

# this function selects an area in a healpix file in a shape with x vertices.
def cutout_shape(nside, pixelnum, radius, ordering='ring'):
    nside = nside
    ordering = ordering
    #vertices = [[(1,2),(3,4),(5,6),(7,8)]]
    print "before the funct"
    if ordering.lower() == 'ring':
        nest=False
    elif ordering.lower() == 'nested':
        nest=True

    cutout_pixels = hp.query_disc(nside, (0,0,0), radius, inclusive=False, nest=nest)
    cutout_pixelnum = len(cutout_pixels)


    return nside, cutout_pixelnum, cutout_pixels, ordering


def cutout_square(ra, dec):
    mean_ra = np.mean(ra)
    mean_dec = np.mean(dec)
    min_ra = np.min(ra)
    max_ra = np.max(ra)
    min_dec = np.min(dec)
    max_dec = np.max(dec)

    center = (mean_ra, mean_dec)
    ra_radius = mean_ra - min_ra
    ra_a = np.sqrt((ra_radius**2) / 2)
    dec_radius = mean_dec - min_dec
    dec_a = np.sqrt((dec_radius**2) / 2)
    upper_ra_boundary = mean_ra + ra_a
    lower_ra_boundary = mean_ra - ra_a
    upper_dec_boundary = mean_dec + dec_a
    lower_dec_boundary = mean_dec - dec_a
    #upper_left = (mean_ra - ra_a, mean_dec + dec_a)
    #upper_left = center + (-ra_a, dec_a)
    #upper_right = center + (ra_a, dec_a)
    #lower_left = center + (-ra_a, -dec_a)
    #lower_right = center + (ra_a, -dec_a)
    #print "corners clockwise from upper left:"
    #print upper_left, upper_right, lower_right, lower_left
    print "upper ra, lower ra, upper dec, lower dec:"
    print upper_ra_boundary, lower_ra_boundary, upper_dec_boundary, lower_dec_boundary
    return upper_ra_boundary, lower_ra_boundary, upper_dec_boundary, lower_dec_boundary
# this function converts healpix coordinate systems to ra and dec.
# def healpix_to_ra_dec(nside, pixelnum, plotfile, data, ordering='ring', projection='ortho', save_show='show'):


def healpix_to_RA_dec(nside, pixelnum, ordering='ring'):
    # file pixels may be organized as ring or nest
    if ordering.lower() == 'ring':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=False, lonlat=True)
    elif ordering.lower() == 'nested':
        ra, dec = hp.pixelfunc.pix2ang(int(nside), pixelnum, nest=True, lonlat=True)
    return ra, dec


def fits_extractor(filename_Q, filename_U, filename_I, coords='RA_dec'):
    if not filename_U.endswith('.fits'):
        raise ValueError('Stokes U file is not a fits file. Files must be .fits files')
    if not filename_Q.endswith('.fits'):
        raise ValueError('Stokes Q file is not a fits file. Files must be .fits files')
    if not filename_I.endswith('.fits'):
        raise ValueError('Stokes I file is not a fits file. Files must be .fits files')
    # get header and data info from the Stokes fits files
    contents_Q = fits.open(filename_Q)
    pixelnum_Q = contents_Q[1].header['naxis1']
    data_Q = contents_Q[1].data
    nside_Q = contents_Q[1].header['nside']
    ordering_Q = contents_Q[1].header['ordering']

    contents_U = fits.open(filename_U)
    pixelnum_U = contents_U[1].header['naxis1']
    data_U = contents_U[1].data
    nside_U = contents_U[1].header['nside']
    ordering_U = contents_U[1].header['ordering']

    contents_I = fits.open(filename_I)
    pixelnum_I = contents_I[1].header['naxis1']
    data_I = contents_I[1].data
    nside_I = contents_I[1].header['nside']
    ordering_I = contents_I[1].header['ordering']


    if pixelnum_U != pixelnum_Q:
        raise ValueError('files do not have same number of pixels.')

    if nside_U != nside_Q:
        raise ValueError('files do not have same nside.')

    if ordering_U != ordering_Q:
        raise ValueError('files do not have same ordering.')

    # extract data from specified files
    pixels_Q = data_Q.field('PIXEL')
    signal_Q = data_Q.field('SIGNAL')
    pixels_U = data_U.field('PIXEL')
    signal_U = data_U.field('SIGNAL')
    pixels_I = data_I.field('PIXEL')
    signal_I = data_I.field('SIGNAL')

    if np.any(pixels_Q != pixels_U):
        raise ValueError('files do not have same set of pixels.')
    """
    if cutout is True:
        #nside, cutout_pixelnum, cutout_pixels, ordering
        nside_Q, pixelnum_Q, pixels_Q, ordering_Q = cutout_shape(nside_Q, pixels_Q, 5, ordering=ordering_Q)
        print "pixels selected in shape with vertices"
    elif cutout is False:
        nside_Q, pixels_Q, ordering_Q = nside_Q, pixels_Q, ordering_Q
        print "all pixels selected"
    """
    # converting to RA and dec if option is selected
    if coords is 'RA_dec':
        ra, dec = healpix_to_RA_dec(nside_Q, pixels_Q, ordering=ordering_Q)
        return signal_Q, signal_U, signal_I, ra, dec
    elif coords is 'pixels':
        return signal_Q, signal_U, signal_I, nside_Q, pixels_Q, ordering_Q
