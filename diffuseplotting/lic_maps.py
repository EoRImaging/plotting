import numpy as np
import pylab as plt
from scipy.interpolate import griddata
import math
from numpy import interp
import lic_internal
import stokes_math
import plothealpix_map
import fits_data_extraction
import copy
import matplotlib.pyplot as plt
import matplotlib.axes as axes
from matplotlib.colors import LinearSegmentedColormap
# automatically focuses at center


def LIC(obsID, x_stokes, y_stokes, ra, dec, dpi, size, length=1000, width=8., full_image=True,
        disp_drapery='save', name_of_plot='flow-image.png', transparency=1, interp_theta=False, rotate=False,
        name_of_interp_plot='interpolated_theta.jpg', polarization=2):

    if polarization == 1:
        x_stokes = x_stokes
        y_stokes = y_stokes
    elif polarization == 2:
        x_stokes_rotated = -y_stokes
        y_stokes_rotated = x_stokes
        x_stokes = x_stokes_rotated
        y_stokes = y_stokes_rotated
    dpi = dpi
    if full_image is True:
        mean_ra = np.mean(ra)
        mean_dec = np.mean(dec)
        lower_ra = np.min(ra)
        upper_ra = np.max(ra)
        lower_dec = np.min(dec)
        upper_dec = np.max(dec)
    elif full_image is False:
        upper_ra, lower_ra, upper_dec, lower_dec =\
            fits_data_extraction.cutout_square(ra, dec)
        mean_ra = (upper_ra + lower_ra) / 2
        mean_dec = (upper_dec + lower_dec) / 2

    video = False

    x_range = upper_ra - lower_ra
    y_range = upper_dec - lower_dec
    aspect_ratio = y_range / x_range
    x_size = size
    y_size = int(round(size * aspect_ratio))

    xi = np.linspace(0, x_size, x_size)
    yi = np.linspace(0, y_size, y_size)
    # scale = divide by lenght of x, multoply by the max-min, add min ra. then check in right range
    ra_reg = (xi / max(xi))*(upper_ra - lower_ra) + lower_ra
    dec_reg = (yi / max(yi))*(upper_dec - lower_dec) + lower_dec

    # yi = ((yi / max(yi)) * (upper_range_dec - lower_range_dec)) + lower_range_dec

    indices = (np.where((ra >= lower_ra) & (ra <= upper_ra) &
                        (dec >= lower_dec) & (dec <= upper_dec)))
    #print indices

    x_stokes = x_stokes[indices]
    y_stokes = y_stokes[indices]
    #x_stokes_rotated = x_stokes_rotated[indices]
    #y_stokes_rotated = y_stokes_rotated[indices]

    #if rotate=True:
    #    x_stokes = x_stokes_rotated
    #    y_stokes = y_stokes_rotated

    ra = ra[indices]
    dec = dec[indices]
    #K2 = np.sqrt(x_stokes**2 + y_stokes**2)
    #theta2 = np.arcsin(x_stokes / K2)
    #plothealpix_map.mapping(ra, dec, theta2, 'theta', obsID=None, map_file_name=None,
    #        projection='ortho', save_show='show', full_image=True)
    #plothealpix_map.mapping(ra, dec, K, 'K', )
    x = ra_reg[:, np.newaxis]
    y = dec_reg[np.newaxis, :]
    x = np.repeat(x, y_size, axis=1).flatten()
    y = np.repeat(y, x_size, axis=0).flatten()
    #print ra_reg
    #print dec_reg
    # theta = np.arctan(y_stokes / x_stokes) + np.pi/2
    xval = griddata((ra, dec), x_stokes, (x, y), method='nearest')
    yval = griddata((ra, dec), y_stokes, (x, y), method='nearest')


    xval = np.reshape(xval, (x_size, y_size))
    yval = np.reshape(yval, (x_size, y_size))

    if interp_theta is True:
        K3 = np.sqrt(xval**2 + yval**2)
        theta3 = np.arcsin(xval / K3)
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
        print theta3.shape
        theta3 = np.rot90(theta3, 1, axes=(0,1))
        goober = plt.imshow(theta3, cmap=circular, aspect='equal', extent=(lower_ra, upper_ra, lower_dec, upper_dec))
        plt.gca().invert_xaxis()
        plt.colorbar()
        #plt.imsave(name_of_interp_plot, theta3, cmap=circular, aspect='equal', extent=(lower_ra, upper_ra, lower_dec, upper_dec))
        plt.savefig(name_of_interp_plot)
        print "saved interpolated data to " + name_of_interp_plot

    vectors = np.zeros((len(ra_reg), len(dec_reg), 2), dtype=np.float32)
    vectors[:, :, 0] = xval
    vectors[:, :, 1] = yval

    # lic_internal is configured for texture shape in the order (y, x)
    # instead of (x, y)
    original_texture = np.random.rand(y_size, x_size).astype(np.float32)
    texture = np.random.rand(int(np.ceil(y_size / width)), int(np.ceil(x_size / width))).astype(np.float32)
    text = np.repeat(np.repeat(texture, int(width), axis=0), int(width), axis=1)
    text = text[0:y_size, 0:x_size]
    fig = plt.figure()

    plt.bone()
    frame = 0

    if video:
        kernellen = length
        for t in np.linspace(0, 1, 16 * 5):
            kernel = np.sin(np.arange(kernellen)*np.pi/kernellen)*(1+np.sin(2*np.pi*5*(np.arange(kernellen)/float(kernellen)+t)))

            kernel = kernel.astype(np.float32)

            image = lic_internal.line_integral_convolution(vectors, text, kernel)

            plt.clf()
            plt.axis('off')
            plt.figimage(image)
            plt.gcf().set_size_inches((x_size / float(dpi), y_size / float(dpi)))
            plt.savefig("flow-%04d.png" % frame, dpi=dpi)
            frame += 1
    else:
        kernellen = length
        kernel = np.sin(np.arange(kernellen) * np.pi / kernellen)
        kernel = kernel.astype(np.float32)
        image = lic_internal.line_integral_convolution(vectors, text, kernel)

        #plt.clf()
        #plt.axis('equal')
        plt.imshow(image, aspect='equal', alpha=transparency, extent=(lower_ra, upper_ra, lower_dec, upper_dec))
        plt.gca().invert_xaxis()
        #plt.figimage(image, resize=True)
        #plt.gcf().set_size_inches((x_size / float(dpi), y_size / float(dpi)))

        if disp_drapery == 'save':
            plt.savefig(name_of_plot)
            print 'saved drapery to ' + name_of_plot
        elif disp_drapery == 'show':
            plt.show()
        elif disp_drapery == 'none':
            print "drapery plot not created"
        #return plt.imshow(image, aspect='auto', alpha=transparency, extent=(lower_ra, upper_ra, lower_dec, upper_dec))
