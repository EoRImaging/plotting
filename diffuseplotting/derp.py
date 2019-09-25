# derp.py
# Diffuse/Extended Radio source Plotter

import sys, os
import numpy as np 
import scipy.io
#filename = '1130781304_run1_catalog.sav';


# Set to True to enable verbose output
debugging = False

def dprint(str):
    # Verbose mode method
    if debugging:
        print(str)
        
def buildSources(fname):
    sav = scipy.io.readsav(fname)
    cat = sav[sav.keys()[0]]
    sourceList = [] 
    for element in cat:
        sourceList.append(Source(element['id'], element['X'], element['Y'], element['ra'], element['dec'], element['ston'], element['freq'], element['alpha'], element['gain'], element['flag'], element['extend'], element['flux']['XX'], element['flux']['YY'], element['flux']['XY'], element['flux']['YX'], element['flux']['I'], element['flux']['Q'], element['flux']['U'], element['flux']['V']))
        print('Source ID # ' + str(element['id']) + ' loaded.')
    return sourceList
### END def buildSources

# Aitoff-Hammer projection function - deprecated w/ use of Basemap
#def toAH(RA, dec):
#   z = np.sqrt(1+np.cos(dec*np.pi/180)*np.cos(RA*np.pi/360))
#   x = np.cos(dec*np.pi/180)*np.sin(RA*np.pi/360)/z
#   y = np.sin(dec*np.pi/180)/z
#   print('X: ')
#   print(x)
#   print('Y: ')
#   print(y)
#   return (x,y)

class Source(object):
    """
    
    Document!!!
    
    """

    def __init__(self, ID, X, Y, RA, dec, ston, freq, alpha, gain, flag, extend, XX, YY, XY, YX, I, Q, U, V):
        self.ID = ID
        self.X = X
        self.Y = Y
        self.RA = RA
        self.dec = dec
        self.ston = ston
        self.freq = freq
        self.alpha = alpha
        self.gain = gain
        self.flag = flag
        self.XX = XX
        self.YY = YY
        self.XY = XY
        self.YX = YX
        self.I = I
        self.Q = Q
        self.U = U
        self.V = V
        self.extend = []
        if extend != None:
            for element in extend:
                self.extend.append({'ID' : element['ID'], 'X' : element['X'], 'Y' : element['Y'], 'RA' : element['RA'], 'dec' : element['dec'], 'ston' : element['ston'], 'freq' : element['freq'], 'alpha' : element['alpha'], 'gain' : element['gain'], 'flag' : element['flag'], 'XX' : element['flux']['XX'], 'YY' : element['flux']['YY'], 'XY' : element['flux']['XY'], 'YX' : element['flux']['YX'], 'I' : element['flux']['I'], 'Q' : element['flux']['Q'], 'U' : element['flux']['U'], 'V' : element['flux']['V']})
    # def Gaussian(self, raArray, decArray, cap):
    
    def getSourceData(self,args):
        validQuantities =['ID','X','Y','RA','dec','ston','freq','alpha','gain','flag','XX','YY','XY','YX','I','Q','U','V']
        filtArgs = []
        numArgs = 0
            # args must be a list of strings matching valid data product identifiers:
        # ID, X, Y, RA, dec, ston, freq, alpha, gain, flag, XX, YY, XY, YX, I, Q, U, V.
        # getComps returns an NxM numpy array where N is the number of sources in an
        # extended object (1 for non-extended objects) and M is the number of requested
        # data products. The "primary" data for an extended object are included as the
        # first row of the returned array.
        for quantity in args:
            if not(quantity in validQuantities):
                raise ValueError(quantity + ' is not a valid data identifier. The available data are: \nID: Source identification number \nX: Ortho-slant X-coordinate \nY: Ortho-slant y-coordinate \nRA: Right ascension \ndec: Declination \nston: Signal to noise ratio \nfreq: Observation frequency \nalpha: Characteristic number for ston variance with frequency? \ngain: Telescope relative gain \nflag: No idea! \nXX, YY, XY, YX, I, Q, U, V: Polarization and flux data.')
            numArgs = numArgs + 1
        
        n = 0
        outArray = np.empty((len(self.extend)+1,numArgs))
        for quantity in args:
                outArray[0,n] = self.__dict__[quantity]
                n = n+1
        m = 1
        n = 0
        for subSource in self.extend:
            for quantity in args:
                outArray[m,n] = subSource[quantity]
                n = n+1
            m = m+1
            n=0
        return outArray
        
    ### END def getSourceData
    
### END Source class definition
