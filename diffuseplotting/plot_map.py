from mpl_toolkits.basemap import Basemap, projection_params
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from matplotlib.colors import LinearSegmentedColormap
import warnings
import matplotlib.cbook
import haversine_function
import healpix_functions
import fits_data_extraction
import matplotlib.tri as tri
import matplotlib.colors as colors
import os
import re
import seaborn as sns


warnings.filterwarnings("ignore", category=matplotlib.cbook.mplDeprecation)
# This module can map various types of data onto a spherical projection.
# plotfile is the name of the saved plot. data_vals is the third dimension to be plotted.
# def mapping(ra, dec, data_vals, var_name, map_file_name=None, projection='ortho', save_show='show', file_extension='.png'):
def mapping(ra, dec, data_vals, var_name, obsID=None, map_file_name=None,
            projection='ortho', save_show='show', full_image=True):

    ra[np.where(ra > 180)] -= 360
    # plot of Galactic gas with coordinate projection
    if full_image is True:
        min_ra = np.min(ra)
        max_ra = np.max(ra)
        mean_ra = np.mean(ra)
        min_dec = np.min(dec)
        max_dec = np.max(dec)
        mean_dec = np.mean(dec)

    elif full_image is False:
        max_ra, min_ra, max_dec, min_dec =\
            fits_data_extraction.cutout_square(ra, dec)
        mean_ra = (max_ra + min_ra) / 2
        mean_dec = (min_dec + max_dec) / 2

    # print('RA min, max, mean: ', min_ra, max_ra, mean_ra)
    # print('Dec min, max, mean: ', min_dec, max_dec, mean_dec)

    #print('{p} projection'.format(p=projection))
    if 'no corners' in projection_params[projection]:
        if 'bounding_lat' in projection_params[projection]:
            if projection in ['nplaea', 'npstere', 'npaeqd']:
                # northern polar projection
                if max_dec < 0:
                    print('polar projection in wrong hemisphere: data is in '
                          'the southern hemisphere')
                    boundinglat = 10
                else:
                    boundinglat = abs(min_dec)
            else:
                # southern polar projection
                if max_dec > 0:
                    print('polar projection in wrong hemisphere: data is in '
                          'the northern hemisphere')
                    boundinglat = -10
                else:
                    boundinglat = max_dec

            m = Basemap(projection=projection, lon_0=mean_ra,
                        boundinglat=boundinglat, resolution=None)
        else:
            # examples are 'moll', 'hammer', 'eck4'
            print('corners are not supported in this projection')
            # these projections can only make full-globe maps, any corner information is ignored
            m = Basemap(projection=projection, lon_0=mean_ra, lat_0=mean_dec,
                        resolution=None)
    elif 'llcrnrx' in projection_params[projection]:
        # example is 'orth'
        # can't use llcrnrlon, llcrnrlat, urcrnrlon, urcrnrlat even though the docs say you can
        # you also can't use width or height even though the docs say you can
        # you have to use llcrnrx, llcrnry, urcrnrx, urcrnry but I'm not sure how to calculate them....
        # these numbers are based on some online examples and experimentation
        # ra_deg_scale_factor = 115600
        dec_deg_scale_factor = haversine_function.haversine(min_ra, mean_dec, max_ra, mean_dec) / 40
        ra_deg_scale_factor = haversine_function.haversine(min_ra, mean_dec, max_ra, mean_dec) / 50
        llcrnrx = -1 * (max_ra - min_ra) * ra_deg_scale_factor
        urcrnrx = (max_ra - min_ra) * ra_deg_scale_factor
        if projection not in ['geos']:
            llcrnry = -1 * (max_dec - min_dec) * dec_deg_scale_factor
            urcrnry = (max_dec - min_dec) * dec_deg_scale_factor
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        resolution=None, llcrnrx=llcrnrx, llcrnry=llcrnry,
                        urcrnrx=urcrnrx, urcrnry=urcrnry)
        else:
            # need a new scale factor for some reason
            deg_scale_factor = 175000
            llcrnry = min_dec * deg_scale_factor / 2
            urcrnry = max_dec * deg_scale_factor / 2
            m = Basemap(projection=projection, lon_0=np.mean(ra), resolution=None,
                        llcrnrx=llcrnrx, llcrnry=llcrnry, urcrnrx=urcrnrx,
                        urcrnry=urcrnry)
    elif 'lat_1,lat_2,lon_1,lon_2' in projection_params[projection]:
        m = Basemap(projection=projection, lon_0=mean_ra, lat_0=mean_dec,
                    lon_1=mean_ra, lat_1=min_dec, lon_2=mean_ra, lat_2=max_dec,
                    llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                    urcrnrlat=max_dec, resolution=None)
    elif 'lon_0' in projection_params[projection]:
        if projection in ['rotpole']:
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        o_lat_p=0, o_lon_p=0,
                        llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                        urcrnrlat=max_dec, resolution=None)
        else:
            # example: 'aeqd'
            m = Basemap(projection=projection, lon_0=np.mean(ra), lat_0=np.mean(dec),
                        llcrnrlon=min_ra, llcrnrlat=min_dec, urcrnrlon=max_ra,
                        urcrnrlat=max_dec, resolution=None)
    else:
        m = Basemap(projection=projection, llcrnrlon=min_ra, llcrnrlat=min_dec,
                    urcrnrlon=max_ra, urcrnrlat=max_dec, resolution=None)

    fig = plt.figure()
    x, y = m(ra, dec)
    ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
    # types of contour plots are pcolor, and contourf.


    cdict1 = {'red':   ((0.0, 0.0, 0.0),
                        (0.25, 0.5, 0.5),
                        (0.5, 1.0, 1.0),
                        (0.75, 0.5, 0.5),
                        (1.0, 0.0, 0.0)),

             'green':   ((0.0, 0.0, 0.0),
                         (0.25, 0.5, 0.5),
                         (0.5, 1.0, 1.0),
                         (0.75, 0.5, 0.5),
                         (1.0, 0.0, 0.0)),

             'blue':    ((0.0, 0.0, 0.0),
                         (0.25, 0.5, 0.5),
                         (0.5, 1.0, 1.0),
                         (0.75, 0.5, 0.5),
                         (1.0, 0.0, 0.0)),
             }

    circular = LinearSegmentedColormap('BlueRed1', cdict1)


    if var_name is 'sources':
        #plt.plot(x, y, ".")
        image = m.scatter(x, y, 3, marker='.', linewidths=.1, c=data_vals) #, cmap=plt.cm.plasma)  # norm=LogNorm())
        print "scatter"
    else:
        #m.contourf(x, y, data_vals, tri=True, cmap=plt.cm.plasma, bins='log', norm=colors.LogNorm())
        #m.contourf(x, y, data_vals, tri=True, cmap=plt.cm.plasma, norm=LogNorm())
        image = m.contourf(x, y, data_vals, tri=True, cmap=circular)
        fig.gca().invert_xaxis()
        #image = m.pcolor(x, y, data_vals, tri=True, cmap=circular)
        #image.set_edgecolor("face")
        for c in image.collections:
            c.set_edgecolor("face")
        m.colorbar()
    # draw parallels and meridians. Labels are 1/0 as [Top,bottom,right,left]
    m.drawparallels(np.arange(-90., 120., 5.), labels=[1, 0, 0, 0])
    m.drawmeridians(np.arange(0., 420., 5.), labels=[0, 0, 0, 1])
    # var_name = ''.join(var_name)
    #plt.imshow((x, y), aspect='auto', extent=(min_ra, max_ra, min_dec, max_dec))

    if obsID is not None:
        plt.title('{} {} map: {} proj.'.format(obsID, var_name, projection))
    else:
        plt.title('{} map: {} proj.'.format(var_name, projection))
    # either show the graphs, or save them to a location.
    if save_show == 'show':
        plt.show()
    elif save_show == 'none':
        print "polarization map not displayed"
    elif save_show == 'save':
        new_map_file = map_file_name
        plt.savefig(new_map_file, dpi=300)
        print 'saved mapped data to ' + new_map_file
    else:
        raise ValueError('save_show needs to be equal to "save" or "show" to save or show the image.')
