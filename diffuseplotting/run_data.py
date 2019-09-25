import lic_maps
import plotstokesQ_U
import stokes_math
import fits_data_extraction
import plot_map
import stokes_histogram
import os
import re
import numpy as np
from matplotlib import gridspec
import matplotlib.pyplot as plt


def run_data(filename_Q, filename_U, filename_I, plot_variable=None,
             save_show='show', graph_selection='all', file_extension='.png', projection='cyl', transparency=1,
             interp_theta=False, lic_length=1000, lic_width=8., polarization=1, directory=None,
             histogram_file_basename=None, map_file_basename=None, drapery_file_basename=None):
    signal_Q, signal_U, signal_I, ra, dec =\
        fits_data_extraction.fits_extractor(filename_Q, filename_U, filename_I)
    ra[np.where(ra > 180)] -= 360

    K, theta, x_stokes, y_stokes = stokes_math.math(signal_I, signal_Q, signal_U)

    plot_variables = {'theta': theta, 'K': K, 'I': signal_I, 'Q': signal_Q, 'U': signal_U}

    if plot_variable not in plot_variables.keys():
        raise ValueError('Plot variable is not in dict')

    map_data_var = plot_variables[plot_variable]

    """
    Filename and directory information
    """
    if directory is None:
        directory = os.getcwd()
    if not os.path.isdir(directory):
        os.makedirs(directory)

    filename_Q = os.path.split(filename_Q)[-1]
    filename_U = os.path.split(filename_U)[-1]
    obsID = re.findall('\d+',filename_Q)
    obsID = 'ObsID_' + ''.join(obsID)

    if histogram_file_basename is None:
        # histogram_filename = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        histogram_file_basename = obsID + '_theta_vs_magnitude_hist'
    histogram_filename = histogram_file_basename + file_extension

    if map_file_basename is None:
        # histogram_filename = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        map_file_basename = obsID + '_map'

    map_filename = map_file_basename + '_' + plot_variable + file_extension

    if drapery_file_basename is None:
        # histogram_filename = os.path.splitext(filename_Q)[0] + '_' + os.path.splitext(filename_U)[0]
        drapery_file_basename = obsID + '_drapery'

    drapery_filename = drapery_file_basename + file_extension
    # A whole bunch of filename and directory meshing
    path_map_filename = os.path.join(directory, map_filename)
    map_filename = path_map_filename
    path_new_histogram_filename = os.path.join(directory, histogram_filename)
    histogram_filename = path_new_histogram_filename
    path_new_drapery_filename = os.path.join(directory, drapery_filename)
    drapery_filename = path_new_drapery_filename

    # All three options are plotted
    if graph_selection is 'all':
        stokes_histogram.plotstokesQU(obsID, x_stokes, y_stokes, ra, dec,
                          new_histogram_filename=histogram_filename, directory=None,
                          save_show=save_show, projection=projection)
        if plot_variable is not None:
            plot_map.mapping(ra, dec, map_data_var, plot_variable, obsID=obsID,
                                    map_file_name=map_filename, projection=projection, save_show=save_show, full_image=False)

        else:
            raise ValueError('Must have mapping variable and variable name inputs to plot map.')

        lic_maps.LIC(obsID, x_stokes, y_stokes, ra, dec, 300, 1000, length=lic_length, width=lic_width, full_image=False,
                     disp_drapery=save_show, name_of_plot=drapery_filename, transparency=transparency, interp_theta=interp_theta, polarization=polarization)

    # Stokes theta/magnitude histogram WORKS
    elif graph_selection is 'stokes_histogram':
        stokes_histogram.plotstokesQU(obsID, x_stokes, y_stokes, ra, dec,
                          new_histogram_filename=histogram_filename, directory=None,
                          save_show=save_show, projection=projection)

    # A selected variable mapped onto some map projection
    elif graph_selection is 'map':
        if plot_variable is not None:
            plot_map.mapping(ra, dec, map_data_var, plot_variable, obsID=obsID,
                                    map_file_name=map_filename, projection=projection, save_show=save_show, full_image=False)

        else:
            raise ValueError('Must have mapping variable and variable name inputs to plot map.')

    elif graph_selection is 'overlay':
        plot_map.mapping(ra, dec, map_data_var, plot_variable, obsID=obsID,
                                map_file_name=map_filename, projection=projection, save_show=save_show, full_image=False)
        lic_maps.LIC(obsID, x_stokes, y_stokes, ra, dec, 300, 1000, length=lic_length, width=lic_width, full_image=False,
                     disp_drapery=save_show, name_of_plot=drapery_filename, transparency=transparency, interp_theta=interp_theta, polarization=polarization)

    # Plot of the direction of polarization created with the line integral convolution method
    elif graph_selection is 'drapery':
        lic_maps.LIC(obsID, x_stokes, y_stokes, ra, dec, 300, 1000, length=lic_length, width=lic_width, full_image=False,
                     disp_drapery=save_show, name_of_plot=drapery_filename, transparency=transparency, interp_theta=interp_theta, polarization=polarization)

    else:
        raise ValueError('Graph selection options are "all", "stokes_histogram",\
                         "map", "overlay", and "drapery".')
