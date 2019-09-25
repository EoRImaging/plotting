## Synopsis

The purpose of these plotting modules is to be able to analyze the linearly polarized components
of a set of Stokes data files, I, Q, and U.

## Functions

### Stokes Plotting and Line Integral Convolution (LIC) Plots
The run_data module is a wrapper that allows you to create multiple kinds of plots of the same
data.  It expects to be given fits files.

```
run_data(filename_Q, filename_U, filename_I, plot_variable=None,
            save_show='show', graph_selection='all', file_extension='.png', projection='cyl', transparency=1,
            interp_theta=False, lic_length=1000, lic_width=8., polarization=1, directory=None,
            histogram_file_basename=None, map_file_basename=None, drapery_file_basename=None)
```
**Parameter Options**

  **Filename_Q, filename_U, filename_I**:  3 files of a typical 4pol data set.  They are configured to
  use .fits files.

  **plot_variable**: the variable to be plotted for the plot_map module.  The other two plotting
  modules plot specific variables.  Default is None.

  **save_show**: Options are 'save' and 'show'. 'Save' option saves the plots you've chosen, as whatever type
  of file has been specified in **file_extension**.

  **graph_selection**: Options are 'map', 'stokes_histogram', 'drapery', 'overlay,'and 'all'.  These options are linked to the modules **plot_map, stokes_histogram** and **lic_maps**.  The option 'overlay' plots drapery and map, and is useful if you want to combine the images.

  **file_extension**: Type of file to be saved. Not necessary if save_show option 'show' is selected.
  Default is .png.

  **projection**: The map projection for the mapping module plot_map. The default for run_data is 'cyl', which is a flattened projection. Projection 'ortho' is good if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot the best.
  Projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.

  **transparency**: this dictates the transparency of the drapery plot.  This is useful for the 'overlay' option, if you want to use the map with the drapery of the same data on top of it.

  **interp_theta**: If plotting the drapery, if True, this will make a basic plot of the data after it's been interpolated to a regular grid.

  **lic_length**: Dictates the length of the "lines" for a drapery plot.  Large values (1000s) will make flow patterns more obvious, while small values show more detail.  This is affected by how large the area of data is.

  **lic_width**: Dictates the size of the randomly generated pixels used in the LIC.  Larger values will reduce detail but make flow lines thicker and more easily followed by eye.  Must be a float. Default is 8., which increases the area of each randomly generated 'pixel' by 64.

  **polarization**: Options are 1 or 2.  This basically allows you to switch which direction the flow lines in the LIC will flow, by 90 degrees.  This is helpful if you see a lot of lines that seem to be parallel to each other and very short.  Possibly the flow will be more visible by trying the other polarization. Default is 1.  If you choose the 2 option, the interpolated_theta plot will not work, but the color-scheme of the map is such that it won't matter.  Just use a map plotted using plot_map or lic_maps where polarization=1.

  **histogram_file_basename, map_file_basename, drapery_file_basename**: Defaults are None.
  If not given an input, run_data will extract the ObsID from filename_Q and append it with variable information.
  (The ObsID extraction finds an integer in the filename.  If the file does not contain an integer number,
  you must set a basename for all plots you're creating.)

  **directory**: If you want to put the plots in a different directory than the current working directory.
  If a directory with the name inputted doesn't exist, one will be created.


```
plot_map.mapping(ra, dec, data_vals, var_name, obsID=None, map_file_name=None,
            projection='ortho', save_show='show', full_image=True)
```

  This module creates a map projected onto a sphere, showing a third dimension in color.  All projections that work for run_data work for this.

  **ra, dec**: The coordinates of the data points.

  **data_vals**: The 3rd dimension data. This info will be plotted with color.

  **var_name**: The reference of the 3rd dimension; This will only be used to title the plot.
  If 'sources' is selected, the generated map will be a scatterplot, instead of a contourf plot.  'Sources' really only works if you don't have a huge number of points.  If there are a lot, they will overlap.

  **obsID**: The Obs ID of the area of data. It's used to title the plot.

  **map_file_name**: The location of the file the plot will go in.

  **projection**: The map projection for the mapping module plot_map. The default for run_data is 'cyl'.
  'ortho' is a good projection if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot.
  Other projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.

  **save_show**: Options are 'save' and 'show'. 'Save' option saves the plots you've chosen, as whatever type
  of file has been specified.

  **full_image**: If True, the window will attempt to show all the data.  If False, the maximum rectangular shape that fits inside the circle of data is selected. This area is selected with the module fits_data_extraction.cutout_square(ra, dec).


```
stokes_histogram.plotstokesQU((obsID, x_stokes, y_stokes, ra, dec,
                 new_histogram_filename=None, directory=None,
                 save_show='show', projection='ortho'))
```

  This module takes data lists 'x_stokes' and 'y_stokes', which are obtained using the module stokes_math.math(signal_I, signal_Q, signal_U), and creates a 2D histogram of the data.  

  **x_stokes and y_stokes** are the projections onto a single axis, of the relationship between the Q and U axes.  These axes are at 45-degree angles to each other, making this transform necessary.

  The output plot shows the density of points at some angle between 0 and pi, and some magnitude (distance from zero).  Ideally, the data should appear as a roughly semicircular shape.

  **obsID**: The Obs ID of the area of data

  **ra, dec**: The coordinates of the data points.

  **new_histogram_filename**: The name of the written-out file. If you don't use the run_data wrapper, this needs to be a complete path, including extension.

  **directory**: location of the file.  Default is current working directory.

  **save_show**: Options are 'save' and 'show'.

  **projection**: The map projection for the mapping module plot_map. The default for run_data is 'cyl'.
  'ortho' is a good projection if you want to see curved RA and Dec lines, while 'cyl' matches up with the drapery plot.
  Other projections include 'ortho', 'cyl', 'merc', 'moll' and more, detailed in the matplotlib package Basemap.


```
lic_maps.LIC(obsID, x_stokes, y_stokes, ra, dec, dpi, size, length=1000, width=8., full_image=True,
        disp_drapery='save', name_of_plot='flow-image.png', transparency=1, interp_theta=False, rotate=False,
        name_of_interp_plot='interpolated_theta.jpg', polarization=1)
```
This module creates a 'drapery' plot of stokes linear polarizations.  The direction of polarization at every point is represented visually as a texture or draping of the angles.

  **obsID**: The Obs ID of the area of data

  **x_stokes and y_stokes** are the projections onto a single axis, of the relationship between the Q and U axes.  These axes are at 45-degree angles to each other, making this transform necessary.  You get them using the stokes_math module.

  **ra, dec**: The coordinates of the data points.

  **dpi**: The "resolution" of the output image.  Recommended dpi is 200-300.

  **size**: The approximate number of dots across. LIC finds an aspect ratio from this information and uses it to scale the x and y directions.  Size will be the x direction.

  **length**: Dictates the length of the "lines" for a drapery plot.  Large values (1000s) will make flow patterns more obvious, while small values show more detail.  This is affected by how large the area of data is.

  **width**: Dictates the size of the randomly generated pixels used in the LIC.  Larger values will reduce detail but make flow lines thicker and more easily followed by eye.  Must be a float. Default is 2., which increases the length of each pixel side by 2x, and increases the area of each randomly generated 'pixel' by 4.

  **full_image**: If True, the window will attempt to show all the data.  If False, the maximum rectangular shape that fits inside the circle of data is selected. This area is selected with the module fits_data_extraction.cutout_square(ra, dec).

  **disp_drapery**: Can be 'save' or 'show'.

  **name_of_plot**: The location of the file the plot will go in. Needs to include the extension if not using run_data module.

  **transparency**: this dictates the transparency of the drapery plot.  This is useful for the 'overlay' option, if you want to use the map with the drapery of the same data on top of it.

  **interp_theta**: If True, this will make a basic plot of the data after it's been interpolated to a regular grid.

  **rotate**: Default is false.  If True, "rotates" the polarizations of the points.

  **name_of_interp_plot**: If interp_theta is True, the name of the plot to be saved out.

  Note: The option interp_theta default is to save the file.  This is because there are some artifacts that sometimes appear when plotting the data using the plot_map module, that are only visible if the plot is saved to a publication-quality file type such as jpg, pdf or eps.

  **polarization**: Options are 1 or 2.  This basically allows you to switch which direction the flow lines will flow by 90 degrees.  This is helpful if you see a lot of lines that seem to be parallel to each other and very short.  Possibly the flow will be more visible by trying the other polarization. Default is 1.  If you choose the 2 option, the interpolated_theta plot will not work, but the color-scheme of the map is such that it won't matter.  Just use a map plotted using plot_map or lic_maps where polarization=1.

### Minor Modules that make the above work

**fits_data_extraction**

```
cutout_square(ra, dec)
```

This function selects the largest rectangular shape that fits inside an area of data.  This is what is called when the full_image option is set to False.

```
healpix_to_RA_dec(nside, pixelnum, ordering='ring')
```

This function converts a healpix formatted file to RA and Dec instead of nside and pixelnumber.  If you choose the wrong ordering, the data will appear as one long streak of data instead of a rectangular or ovoid shape.  Options for **ordering** are 'ring' and 'nested'.

```
fits_extractor(filename_Q, filename_U, filename_I, coords='RA_dec')
```

This function takes fits data files and returns signal info.

**returns**:
signal_Q, signal_U, signal_I, nside_Q, pixels_Q, ordering_Q
All three fits files need to be the same nside, pixel number, and ordering style in order to use this function.

```
math(I, Q, U)
```

Inputs signals of I, Q, and U, and outputs K, theta, x_stokes, and y_stokes.

**K** is the magnitude of Q and U combined. In other words, it is the sum of all linear polarization.

**theta** is the angle in reference to an x-y coordinate system of a signal's Q and U components.

**x_stokes and y_stokes** are the projections onto the x and y axes of theta.

```
haversine(lon1, lat1, lon2, lat2)
```
This is basically only necessary in order to make basemap work in plot_map. It calculates the great-circle distance between two points on Earth.
