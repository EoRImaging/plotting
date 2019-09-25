# gaussfuncs.py - contains functions for plotting Gaussian functions on
# HealPIX pixel maps.
from math import exp
import healpy as hp
import numpy as np
from havnp import haversine
import itertools
from numba import autojit

@autojit
def computeGauss(amp, dist, var):
    prod =  amp*(np.exp(-np.square(dist)/(2*(var**2))))
    #if any(prod < 0):
    #    print('computeGauss returning negative number' + str(prod) + 'with inputs:')
    #    print('amp: '+str(amp))
    #    print('dist: '+str(dist))
    #    print('var: '+str(var))
    return prod
    
def makeGauss(nside, ra, dec, amp, varparam=np.sqrt(2*np.pi)): #bounds=None):#radius=6447797.0):
    print('RA length:' + str(len(ra)))
    print('dec length:' + str(len(dec)))
    print('amp length:' + str(len(amp)))
    raMin = min(ra)-2.0
    raMax = max(ra)+2.0
    decMin = min(dec)-2.0
    decMax = max(dec)+2.0
    #m = np.arange(hp.nside2npix(nside))
    vertices = hp.ang2vec([raMin,raMin,raMax,raMax],[decMin,decMax,decMax,decMin],lonlat=True)
    m = hp.query_polygon(nside,vertices)
    m_ra, m_dec = hp.pix2ang(nside, m, lonlat=True)
    m_amp = np.zeros_like(m,dtype='float64')
    #ra = ra*1.0
    #dec = dec*1.0
    numpix = m.size
    
    #if True: #bounds != None:
    #    raMin = min(ra)-2.0
    #    raMax = max(ra)+2.0
    #    decMin = min(dec)-2.0
    #    decMax = min(dec)+2.0
        # Compute distance of center of data to each pixel
        # centDist = haversine(np.full_like(m,raCent),np.full_like(m,decCent),m_ra,m_dec)
        # Get indices of pixels to delete
    #    toDelete = [ k for (i,j,k) in zip(m_ra,m_dec,m) if i < raMin or i > raMax or j < decMin or j > decMax ]
    #    #m = np.delete(m, toDelete)
    #    m_ra = np.delete(m_ra, toDelete)
    #    m_dec = np.delete(m_dec, toDelete)
    #    m_amp = np.delete(m_amp, toDelete)
    #    numdeleted = str(len(toDelete))
    #    print(numdeleted + ' uninteresting pixels deleted of '+str(numpix))
    
    ### TO IMPLEMENT: Check ra, dec, amp are of same size
    #if isinstance(ra, np.ndarray):
    for rai, deci, ampi in itertools.izip(ra, dec, amp):
        # Used to convert from ra/dec to pixel to vec; should not be necessary
        #nearpix = hp.pix2vec(nside,hp.ang2pix(nside,rai,deci,lonlat=True))
        nearpix = hp.ang2vec(rai,deci,lonlat=True)
        # Calculate Gaussian on pixels within following radius (in degrees)
        radius = 3.0/60.0 # Hard-coded; remove eventually
        pxls = hp.query_disc(nside,nearpix,radius*np.pi/180.0)
        indices = np.searchsorted(m,pxls)
        near_ra, near_dec = hp.pix2ang(nside,pxls,lonlat=True)
        #print('Using nearest '+str(len(near_ra))+' pixels to source.')
        #near_ra = m_ra[pxls]
        #near_dec = m_dec[pxls]
        m_dist = np.zeros_like(near_ra)
        var = varparam*1.0
        m_dist = haversine(np.full_like(near_ra,rai), np.full_like(near_ra,deci), near_ra, near_dec)
        m_amp[indices] += computeGauss(ampi, m_dist, var)
    print('makeGauss output length: '+str(m_ra.size))
    return m_ra, m_dec, m_dist, m_amp
    #else:
    #    #m_dist = np.empty_like(m)
    #    var = varparam*amp*1.0
    #    print('RA: ' + str(ra))
    #    print('dec: ' + str(dec))
    #    m_dist = haversine(np.full_like(m,ra), np.full_like(m,dec), m_ra, m_dec)
    #    m_amp += computeGauss(amp, m_dist, var) #amp * (np.exp(-np.square(m_dist) / (2*(var**2))))     
        
        
    #centDist = np.zeros_like(m)
    #if bounds != None:
    #    raMin = bounds[0]*1.0
    #    raMax = bounds[1]*1.0
    #    decMin = bounds[2]*1.0
    #    decMax = bounds[3]*1.0
    #    # Compute distance of center of data to each pixel
    #    # centDist = haversine(np.full_like(m,raCent),np.full_like(m,decCent),m_ra,m_dec)
    #    # Get indices of pixels to delete
    #    toDelete = [ k for (i,j,k) in zip(m_ra,m_dec,m) if i < raMin or i > raMax or j < decMin or j > decMax ]
    #    #m = np.delete(m, toDelete)
    #    m_ra = np.delete(m_ra, toDelete)
    #    m_dec = np.delete(m_dec, toDelete)
    #    m_amp = np.delete(m_amp, toDelete)
    #    numdeleted = str(len(toDelete))
    #    print(numdeleted + ' uninteresting pixels deleted of '+str(numpix))
