from gaussfuncs import makeGauss
import healpy as hp
import numpy as np
import scipy.io
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import matplotlib.colors as colors
#from derp import *
import os, errno, re
import cartopy.crs as ccrs
from cartopy.mpl.gridliner import LONGITUDE_FORMATTER, LATITUDE_FORMATTER, Gridliner
import types
import warnings
from jdb import *

def crdir(directory):
    """
crdir - Create directory. Creates a directory if it doesn't already exist.
If the directory already exists, it handles the error intelligently.
If the directory cannot be created for any other reason, an error is raised.
Input: string containing path of directory to be created. 
Returns nothing.   
    """
    try: 
        os.makedirs(directory)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
### END def crdir

#class MyGridliner(Gridliner):
#    def _assert_can_draw_ticks(self):
#        return True
### END class MyGridliner

def numpize(fname,saveData=True,points=True,centroids=False,components=False,*quantities):
    """Convert IDL .sav file of catalog data with RA, dec, flux to numpy, .npz
    
    Arguments:
    fname - String containing path to IDL .sav file
    saveData - Boolean; saves .npz file in same dir as fname
    """
    ###
    # DRY helper function - generates 1xN np array from source of quantities
    def getData(source,ext=False):
        arr = np.array([[]])
        srcparams = ['id','ID','x','X','y','Y','ra','RA','dec','DEC','ston','STON','freq','FREQ','alpha','ALPHA','gain','GAIN','flag','FLAG','extend','EXTEND','flux','FLUX']
        fluxparams = ['xx','XX','yy','YY','xy','XY','yx','YX','i','I','q','Q','u','U','v','V']
        if not(ext):
            for q in quantities:
                if q in fluxparams:
                    arr = np.append(arr,source['flux'][q])
                elif q in srcparams:
                    arr = np.append(arr,source[q])
                else:
                    raise ValueError(q + ' is not a quantity in catalog files.')
            return np.reshape(arr,(1,len(quantities)))
        else:
            for subsrc in src['extend']:
                for q in quantities:
                    if q in fluxparams:
                        arr = np.append(arr, subsrc['flux'][q])
                    elif q in srcparams:
                        arr = np.append(arr, subsrc[q])
                    else:
                        raise ValueError(q + ' is not a quantity in catalog files.')
            numSources = len(src['extend'])
            return np.reshape(arr,(numSources,len(quantities)))
    ### END def getData

    # Load IDL file into dict
    cat = scipy.io.readsav(fname)
    # Shed dict wrapper, get np structured dtype
    cat = cat[cat.keys()[0]]
    # Initialize output array
    data = np.empty((0,len(quantities)))
    for src in cat:
        if src['extend']==None:
            if points:
                data = np.append(data, getData(src),axis=0)
        else:
            if centroids:
                data = np.append(data,getData(src),axis=0)
            if components:
                data = np.append(data,getData(src,ext=True),axis=0)
    i=0
    outDict = {}
    for q in quantities:
        outDict[q] = data[:,i]
        i += 1
    if saveData:
        outname = fname[:-4]
        if points:
            outname+='_points'
        if centroids:
            outname+='_centroids'
        if components:
            outname+='_components'
        np.savez_compressed(outname+'.npz',**outDict)
    return outDict
    ###    
    #for source in sourceList:
        #print('Reading source ' + str(source.ID) + ' with RA ' + str(source.RA) + ', dec ' + str(source.dec) + ', flux ' + str(source.I)) 
    #    if flag == 0:
    #        dataArray = source.getSourceData(['RA','dec','I'])
    #        flag = 1
    #    else:
    #        newArray = source.getSourceData(['RA','dec','I'])
    #        dprint('Appending array of shape '+str(np.shape(newArray))+' to dataArray of shape '+str(np.shape(dataArray)))
    #        dataArray = np.append(dataArray, newArray,axis=0)

    #print('RA, dec, I extracted')
    #raout = dataArray[:,0]
    #decout = dataArray[:,1]
    #fluxout = dataArray[:,2]
    #if saveData:
    #    np.savez_compressed(fname[:-3]+'npz',ra=raout,dec=decout,flux=fluxout)
    #return raout, decout, fluxout

def gaussify(sra, sdec, sflux, var, nside, obsID, saveData, nameMod='',radius=2.0):
    # Gaussed data output file name
    gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+nameMod+'.npz'
    
    #TODO: Get rid of following awful hack
    if nameMod == '_components':
        print('Detecting lengthy extended-object component calculation; reducing Gaussian computation radius.')
        radius = radius * 0.1
        print('Using Gaussian computation radius of ' + str(radius)+' for ' + str(len(sra)) + ' extended components.')

    # Normalize to 1.0 to avoid getting negative numbers from log10
    # TODO: Implement +=1, multiplicative normalization switch, compare
    print('Min flux: '+str(min(sflux)))
    print('Max flux: '+str(max(sflux)))
    normConst = 1.0/min(sflux)
    print('Normalizing flux to 1.0 by multiplying by '+str(normConst))
    sflux *= normConst
    vdiag('gaussify', sra=sra, sdec=sdec, log10sflux=np.log10(sflux))
    ra, dec, amp = makeGauss(nside,sra,sdec,np.log10(sflux),var,radius)
    #print('[gaussify] Values returned from makeGauss: ')
    #print('RA: ' + str(ra))
    #print('dec: ' + str(dec))
    #print('amp: ' + str(amp))
    if saveData:
        print('Saving ' + gfname)
        np.savez_compressed(gfname,ra=ra,dec=dec,amp=amp,normConst=normConst)

    return ra, dec, amp, normConst
    ### END def gaussify

def gaussplot(ra,dec,amp,cap,var,nside,obsID,normConst=1.0,savePlot=True,plotPath=None,points=True,centroids=False,components=False,cmap=['Greys'],alpha=[1.0]):
    """ DOCUMENT ME PLZ """
    #print('[gaussplot] Brightest single point has amplitude ' + str(max(amp)))
    
    # Set figure size to be nice and big
    mpl.rcParams['figure.figsize'] = [24.0, 18.0]
    
    # Configure cartopy projection
    proj = ccrs.PlateCarree() # TODO figure out why this isn't Mollweide
    # TODO: Fix ugly central_long hack on next 2 lines
    central_longitude = (max(ra[ra.keys()[0]])+min(ra[ra.keys()[0]]))*0.5
    print('[gaussplot] Min|mid|max RA:'+str(min(ra[ra.keys()[0]]))+'|'+str(central_longitude)+'|'+str(max(ra[ra.keys()[0]])))
    ax = plt.axes(projection=ccrs.Mollweide(central_longitude=central_longitude))

    # TODO: Implement singular source type case
    # TODO: Make alpha, cmap dicts w/ src keywords early on
    # TODO: Make cap a list/dict
    # Plot data set(s)
    i = 0
    for src in ra.keys():
        if src == 'points':
            srcCmap = cmap[0]
        elif src == 'centroids':
            srcCmap = cmap[1]
        elif src == 'components':
            srcCmap = cmap[2]
        print('[gaussplot] Plotting ' + src)
        ax.scatter(ra[src],dec[src],c=np.clip(amp[src],None,cap),s=0.0001,cmap=srcCmap,alpha=alpha[i],transform=proj)
    
        # Create colorbar in a funny way 'cuz Cartopy
        minflux = (10**min(amp[src]))/normConst[src]
        maxflux = (10**max(amp[src]))/normConst[src]
        smap = plt.cm.ScalarMappable(cmap=srcCmap,norm=colors.LogNorm(vmin=minflux,vmax=min(maxflux,10**cap)))
        smap._A = []
        plt.colorbar(smap,ax=ax,alpha=alpha[i], label=src+'(Jy)')
        i += 1
    
    # Invert x-axis since RA increases eastward, opposite of longitude
    plt.gca().invert_xaxis()

    # TODO: Implement axis labels

    plt.title('obsID: '+str(obsID)+'\n var = '+str(var)+'$^\circ$ cap = '+str(cap)+' nside='+str(nside))
    #ax.set_xticks([45,60,75,90,105,120])
    #ax.set_yticks([-20,-30,-40,-50,-60,-70])
    gl = ax.gridlines(crs=proj)#MyGridliner(ax, crs=proj)
    gl._assert_can_draw_ticks = types.MethodType(lambda x: True, gl)
    gl.xlabels_bottom=True
   # gl.xlocator = mticker.LinearLocator(6)
   # gl.xformatter = LONGITUDE_FORMATTER
    #gl.xlines=True
    #gl.ylines=True
    gl.ylabels_left=True
   # gl.ylocator = mticker.LinearLocator(4)
   # gl.yformatter = LATITUDE_FORMATTER
    
    # Plot location of Pictor A galaxy
    #plt.plot((5.0+19.0/60.0+49.7/3600.0)*360.0/24.0,-(45.0+46.0/60.0+44.0/3600.0),'ro',ms=1.0)

    if savePlot:
        if plotPath==None:
            plotPath = './gaussplots/'+obsID+'_nside'+str(nside)+'_var'+str(var)+'_cap' + str(cap)
            for stype in ra.keys():
                plotPath += '_' + stype
            plotPath += '.png'
        print('Saving '+plotPath)
        plt.savefig(plotPath)
        print('Saved.')
        plt.close()
        return
    else:
        print('ACK NOT IMPLEMENTED DO NOT DO THAT')
        return
        #return plt, ax, gl
#   

def makeFromScratch(fname, nside, var, cap, obsID, 
                    saveData, savePlot, plotPath, # I/O opts
                    points, centroids, components, # Source type opts
                    cmaps, alpha, radius): # Plotting opts
    srcTypes = []
    if points:
        srcTypes.append('points')
    if centroids:
        srcTypes.append('centroids')
    if components:
        srcTypes.append('components')
    ra={}
    dec={}
    flux={}
    amp={}
    normConst={}
    print('Extracting ra, dec, flux (I) from IDL .sav file')
    for src in srcTypes:
        if src=='points':
            d = numpize(fname, saveData, True, False, False, 'RA','dec','I')
        if src=='centroids':
            d = numpize(fname, saveData, False, True, False, 'RA','dec','I')
        if src=='components':
            d = numpize(fname, saveData, False, False, True, 'RA','dec','I')
        #ra[src] = d['RA']
        #dec[src] = d['dec']
        #flux[src] = d['I']
        mask = d['dec'] < -20
        ra[src] = d['RA'][mask]
        dec[src] = d['dec'][mask]
        flux[src] = d['I'][mask]
        print('Computing Gaussians of '+ src)
        # Note: At this next line, ra and dec stop being the ra and dec
        # of the catalog sources and start being the gaussian-convolved data
        # on a HealPIX grid.
        ra[src], dec[src], amp[src], normConst[src] = gaussify(ra[src], dec[src], flux[src], var, nside, obsID, saveData,nameMod='_'+src,radius=radius)
    print('Plotting Gaussians')
    gaussplot(ra, dec, amp, cap, var, nside, obsID, normConst, savePlot, plotPath, points, centroids, components, cmaps, alpha)
    print('makeGaussPlot done')
    return
    
def handleFuckery(fname, nside, var, cap, obsID, saveData, loadData, savePlot):
    print('')
    print('makeGaussPlot started using following parameters:\n')
    print('obsID (GPS seconds): ' + obsID)
    print('nside - HealPIX grid resolution: ' + str(nside))
    print('var (degrees) - Variance for each Gaussian: ' + str(var))
    print('cap (arb. units) - Upper limit on plotted values: ' + str(cap))
    #if colorbarLevels > 998:
    #    raise ValueError('Cannot have more than 998 color levels in output contour plot.')
    #    #colorbarLevels = 998
    #elif colorbarLevels < 0:
    #    raise ValueError('Cannot have negative number of colorbar levels.')
    #    #colorbarLevels = 10
    #else:
    #print('Output colorbar levels: ' + str(colorbarLevels))
    if saveData:
        print('saveData is True; numpy ra/dec/flux arrays for sources and their Gaussian convolution will be saved under ./gaussdata')
    else:
        print('saveData is False; no intermediate data will be saved to disk.')
    if loadData:
        print('loadData is True; existing ra/dec/flux numpy arrays of Gaussian-convolved or catalog data will be loaded if available.')
    else:
        print('loadData is False; the plot will be created from scratch starting with the .sav file.')
    if savePlot:
        print('savePlot is True; the plot will be saved under ./gaussplots')
    else:
        raise NotImplementedError('savePlot is False; this option is not implemented.')
        #savePlot = True
    print('')


def makeGaussPlot(fname, nside, var, cap, 
                  saveData=True, loadData=True, savePlot=True, plotPath=None,
                  points=True, centroids=False, components=False,
                  cmaps=['Greys'],alpha=[1.0],radius=5.0):
    """Make truncated-Gaussian plot of EoR catalog .sav data.
    
    Arguments:
    fname -- Path to IDL .sav file of source catalog data to plot
    nside -- Resolution of HealPIX grid to convolve & plot upon
    var -- Variance of Gaussian in decimal degrees
    cap -- Upper limit on plotted values, arbitrary units
    
    Keyword arguments: TODO
    """
    warnings.filterwarnings('ignore', '.*Python3 style.*', RuntimeWarning)
    warnings.filterwarnings('ignore','.*lementwise .= comparison failed.*',FutureWarning)
    
    # Extract obsID from fname using newbie regex
    # DEBUGGING NOTE: Any instance of '113' elsewhere in the file name followed
    # by 7 digits will cause an error!
    searchObj = re.search( r'113\d{7}', fname)
    obsID = searchObj.group()
    handleFuckery(fname, nside, var, cap, obsID, saveData, loadData, savePlot)
    # Create output directories if they don't already exist
    if saveData:
        crdir('./gaussdata')
    if savePlot:
        crdir('./gaussplots')
    
    # Make plots from scratch if loadData is False
    if not(loadData):
        makeFromScratch(fname, nside, var, cap, obsID, saveData, savePlot, plotPath, points, centroids, components, cmaps, alpha, radius)
        return
    
    srcTypes = []
    if points:
        srcTypes.append('points')
    if centroids:
        srcTypes.append('centroids')
    if components:
        srcTypes.append('components')

    ra = {}
    dec = {}
    amp = {}
    normConst = {}
    # Gauss'd file name
    for src in srcTypes:
        gfname='./gaussdata/'+obsID+'_gauss_nside'+str(nside)+'_var'+str(var)+'_'+src+'.npz'
    # If Gauss'd output file already exists, just plot it
        if os.path.isfile(gfname):
            print('Loading existing Gaussian-convolved HealPIX-gridded data')
            npdict = np.load(gfname)
            ra[src] = npdict['ra']
            dec[src] = npdict['dec']
            amp[src] = npdict['amp']
            normConst[src] = npdict['normConst']
            print('Plotting Gaussians')
            print('Plot data length: ' + str(amp[src].size))
            #gaussplot(ra, dec, amp, cap, var, nside, obsID, savePlot)
            #print('makeGaussPlot done')
            #return
        else: # If Gauss'd data doesn't exist:
            # If ra, dec, flux (I) data is already saved in np format, use it
            if os.path.isfile(fname[:-3]+'_'+src+'npz'):
                print('Loading existing .npz file of RA, dec, flux (I)')
                npdict=np.load(fname[:-3]+'_'+src+'npz')
                sra = npdict['ra']
                sdec = npdict['dec']
                sflux = npdict['flux']
                print('Computing Gaussians')
                ra[src], dec[src], amp[src], normConst[src] = gaussify(sra, sdec, sflux, var, nside, obsID, saveData,nameMod='_'+src,radius=radius)
                print('Plotting Gaussians')
                #gaussplot(ra,dec,amp,cap,var,nside,obsID,savePlot)
                #print('makeGaussPlot done')
                #return
            else: # If np format ra, dec, flux don't exist, extract from .sav
                #makeFromScratch(fname, nside, var, cap, obsID, saveData, savePlot,plotPath,points,centroids,components)
                
                if src=='points':
                    d = numpize(fname, saveData, True, False, False, 'RA','dec','I')
                if src=='centroids':
                    d = numpize(fname, saveData, False, True, False, 'RA','dec','I')
                if src=='components':
                    d = numpize(fname, saveData, False, False, True, 'RA','dec','I')
                mask = d['dec'] < -20
                ra[src] = d['RA'][mask]
                dec[src] = d['dec'][mask]
                flux = d['I'][mask]
                print('Computing Gaussians of '+ src)
                # Note: At this next line, ra and dec stop being the ra and dec
                # of the catalog sources and start being the gaussian-convolved data
                # on a HealPIX grid.
                ra[src], dec[src], amp[src], normConst[src] = gaussify(ra[src], dec[src], flux, var, nside, obsID, saveData,nameMod='_'+src,radius=radius)
    gaussplot(ra,dec,amp,cap,var,nside,obsID,normConst,savePlot,plotPath,points,centroids,components, cmaps, alpha)
    print('makeGaussPlot done.')
# A cool thing Bryna did that I may need to reference sometime
#ra[np.where(ra > 180)] -= 360
