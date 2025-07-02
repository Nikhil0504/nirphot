import astropy.units as u
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from astropy.io import fits
from astropy.visualization import ImageNormalize, LinearStretch, PercentileInterval

import logging

logger = logging.getLogger(__name__)

def imshow_with_mouseover(image, ax=None, *args, **kwargs):
    """Wrapper for matplotlib imshow that displays the value under the
    cursor position

    Wrapper for pyplot.imshow that sets up a custom mouseover display
    formatter so that mouse motions over the image are labeled in the
    status bar with pixel numerical value as well as X and Y coords.
    """
    if ax is None:
        ax = plt.gca()
    ax.imshow(image, *args, **kwargs)
    aximage = ax.images[0].properties()['array']
    # need to account for half pixel offset of array coordinates for mouseover relative to pixel center,
    # so that the whole pixel from e.g. ( 1.5, 1.5) to (2.5, 2.5) is labeled with the coordinates of pixel (2,2)

    # We use the extent and implementation to map back from the data coord to pixel coord
    # There is probably an easier way to do this...
    imext = ax.images[0].get_extent()  # returns [-X, X, -Y, Y]
    imsize = ax.images[0].get_size()  # returns [sY, sX]g

    def report_pixel(x, y):
        # map data coords back to pixel coords
        # and be sure to clip appropriately to avoid array bounds errors
        img_y = np.floor((y - imext[2]) / (imext[3] - imext[2]) * imsize[0])
        img_y = int(img_y.clip(0, imsize[0] - 1))

        img_x = np.floor((x - imext[0]) / (imext[1] - imext[0]) * imsize[1])
        img_x = int(img_x.clip(0, imsize[1] - 1))

        return "(%6.3f, %6.3f)     %-12.6g" % (x, y, aximage[img_y, img_x])

    ax.format_coord = report_pixel
    return ax

def display_psf(hdulist_or_filename, ext=0, vmin=1e-7, vmax=1e-1,
                scale='log', cmap=None, title=None, imagecrop=None,
                adjust_for_oversampling=False, normalize='None',
                crosshairs=False, markcentroid=False, colorbar=True,
                colorbar_orientation='vertical', pixelscale='PIXELSCL',
                ax=None, return_ax=False, interpolation=None, cube_slice=None,
                angular_coordinate_unit=u.arcsec):
    """Display nicely a PSF from a given hdulist or filename

    This is extensively configurable. In addition to making an attractive display, for
    interactive usage this function provides a live display of the pixel value at a
    given (x,y) as you mouse around the image.

    Parameters
    ----------
    hdulist_or_filename : fits.hdulist or string
        FITS file containing image to display.
    ext : int
        FITS extension. default = 0
    vmin, vmax : float
        min and max for image display scaling
    scale : str
        'linear' or 'log', default is log
    cmap : matplotlib.cm.Colormap instance or None
        Colormap to use. If not given, taken from user's
        `poppy.conf.cmap_sequential` (Default: 'gist_heat').
    title : string, optional
        Set the plot title explicitly.
    imagecrop : float
        size of region to display (default is whole image)
    adjust_for_oversampling : bool
        rescale to conserve surface brightness for oversampled PSFs?
        (Making this True conserves surface brightness but not
        total flux.) Default is False, to conserve total flux.
    normalize : string
        set to 'peak' to normalize peak intensity =1, or to 'total' to
        normalize total flux=1. Default is no normalization.
    crosshairs : bool
        Draw a crosshairs at the image center (0, 0)? Default: False.
    markcentroid : bool
        Draw a crosshairs at the image centroid location?
        Centroiding is computed with the JWST-standard moving box
        algorithm. Default: False.
    colorbar : bool
        Draw a colorbar on the image?
    colorbar_orientation : 'vertical' (default) or 'horizontal'
        How should the colorbar be oriented? (Note: Updating a plot and
        changing the colorbar orientation is not supported. When replotting
        in the same axes, use the same colorbar orientation.)
    pixelscale : str or float
        if str, interpreted as the FITS keyword name for the pixel scale in arcsec/pixels.
        if float, used as the pixelscale directly.
    ax : matplotlib.Axes instance
        Axes to display into.
    return_ax : bool
        Return the axes to the caller for later use? (Default: False)
        When True, this function returns a matplotlib.Axes instance, or a
        tuple of (ax, cb) where the second is the colorbar Axes.
    interpolation : string
        Interpolation technique for PSF image. Default is None,
        meaning it is taken from matplotlib's `image.interpolation`
        rcParam.
    cube_slice : int or None
        if input PSF is a datacube from calc_datacube, which slice
        of the cube should be displayed?
    angular_coordinate_unit : astropy Unit
        Coordinate unit to use for axes display. Default is arcseconds.
    """
    if isinstance(hdulist_or_filename, str):
        hdulist = fits.open(hdulist_or_filename)
    elif isinstance(hdulist_or_filename, fits.HDUList):
        hdulist = hdulist_or_filename
    else:
        raise ValueError("input must be a filename or FITS HDUList object")

    # Get a handle on the input image
    if hdulist[ext].data.ndim == 2:
        im0 = hdulist[ext].data
        psf_array_shape = hdulist[ext].data.shape
    elif hdulist[ext].data.ndim == 3:
        if cube_slice is None:
            raise ValueError("To display a PSF datacube, you must set cube_slice=<#>.")
        else:
            im0 = hdulist[ext].data[cube_slice]
            psf_array_shape = hdulist[ext].data.shape[1:]
    else:
        raise RuntimeError("Unsupported image dimensionality.")

    # Normalization
    if adjust_for_oversampling:
        try:
            scalefactor = hdulist[ext].header['OVERSAMP'] ** 2
        except KeyError:
            logger.error("Could not determine oversampling scale factor; "
                       "therefore NOT rescaling fluxes.")
            scalefactor = 1
        im = im0 * scalefactor
    else:
        # don't change normalization of actual input array, work with a copy!
        im = im0.copy()

    if normalize.lower() == 'peak':
        logger.debug("Displaying image normalized to peak = 1")
        im /= im.max()
    elif normalize.lower() == 'total':
        logger.debug("Displaying image normalized to PSF total = 1")
        im /= im.sum()
    
    # make 0s NaNs to avoid autoscaling issues
    im[im == 0] = np.nan

    if scale == 'linear':
        norm = matplotlib.colors.Normalize(vmin=vmin, vmax=vmax)
    else:
        norm = matplotlib.colors.LogNorm(vmin=vmin, vmax=vmax)

    if isinstance(pixelscale, str):
        pixelscale = hdulist[ext].header[pixelscale]
    else:
        pixelscale = float(pixelscale)

    if angular_coordinate_unit != u.arcsec:
        coordinate_rescale = (1 * u.arcsec).to_value(angular_coordinate_unit)
        pixelscale *= coordinate_rescale
    else:
        coordinate_rescale = 1
    halffov_x = pixelscale * psf_array_shape[1] / 2.0
    halffov_y = pixelscale * psf_array_shape[0] / 2.0

    unit_label = str(angular_coordinate_unit)
    extent = [-halffov_x, halffov_x, -halffov_y, halffov_y]

    if cmap is None:
        import poppy
        cmap = getattr(matplotlib.cm, poppy.conf.cmap_sequential)
    # update and get (or create) image axes
    ax = imshow_with_mouseover(
        im,
        extent=extent,
        cmap=cmap,
        norm=norm,
        ax=ax,
        interpolation=interpolation,
        origin='lower'
    )
    ax.set_xlabel(unit_label)

    if imagecrop is not None:
        halffov_x = min((imagecrop / 2.0, halffov_x))
        halffov_y = min((imagecrop / 2.0, halffov_y))
    ax.set_xbound(-halffov_x, halffov_x)
    ax.set_ybound(-halffov_y, halffov_y)
    if crosshairs:
        ax.axhline(0, ls=':', color='k')
        ax.axvline(0, ls=':', color='k')
    if title is None:
        try:
            fspec = "%s, %s" % (hdulist[ext].header['INSTRUME'], hdulist[ext].header['FILTER'])
        except KeyError:
            fspec = str(hdulist_or_filename)
        title = "PSF sim for " + fspec
    ax.set_title(title)

    if colorbar:
        if ax.images[0].colorbar is not None:
            # Reuse existing colorbar axes (Issue #21)
            colorbar_axes = ax.images[0].colorbar.ax
            cb = plt.colorbar(
                ax.images[0],
                ax=ax,
                cax=colorbar_axes,
                orientation=colorbar_orientation
            )
        else:
            cb = plt.colorbar(
                ax.images[0],
                ax=ax,
                orientation=colorbar_orientation
            )
        if scale.lower() == 'log':
            ticks = np.logspace(np.log10(vmin), np.log10(vmax), int(np.round(np.log10(vmax / vmin) + 1)))
            if colorbar_orientation == 'horizontal' and vmax == 1e-1 and vmin == 1e-8:
                ticks = [1e-8, 1e-6, 1e-4, 1e-2, 1e-1]  # looks better
            cb.set_ticks(ticks)
            cb.set_ticklabels(ticks)
        if normalize.lower() == 'peak':
            cb.set_label('Intensity relative to peak pixel')
        else:
            cb.set_label('Fractional intensity per pixel')

    if return_ax:
        if colorbar:
            return ax, cb
        else:
            return ax


def save_psf_plots(hdu, path, filt, ext=None):
    from matplotlib.backends.backend_pdf import PdfPages

    with PdfPages(path) as pdf:
        if ext is None:
            ext = ['DET_DIST', 'OVERDIST', 'ROTATED_OVERDIST', 'ROTATED_DET_DIST']
        for e in ext:
            display_psf(hdu, ext=e, ax=plt.gca(), title=f'{filt} PSF {e}', return_ax=False)
            pdf.savefig()
            plt.close()
    
    logger.info(f"PSF plots saved to {path}")


def display_kernel(kernel, title=None, ax=None, return_ax=False, interpolation=None, save_path=None, save=False):
    """Display nicely a PSF kernel

    This is extensively configurable. In addition to making an attractive display, for
    interactive usage this function provides a live display of the pixel value at a
    given (x,y) as you mouse around the image.

    Parameters
    ----------
    kernel : numpy.ndarray
        2D array containing the kernel to display.
    title : string, optional
        Set the plot title explicitly.
    ax : matplotlib.Axes instance
        Axes to display into.
    return_ax : bool
        Return the axes to the caller for later use? (Default: False)
        When True, this function returns a matplotlib.Axes instance.
    interpolation : string
        Interpolation technique for PSF image. Default is None,
        meaning it is taken from matplotlib's `image.interpolation`
        rcParam.
    """
    if ax is None:
        ax = plt.gca()

    # make 0s NaNs to avoid autoscaling issues
    kernel[kernel == 0] = np.nan
    norm = ImageNormalize(
        kernel, interval=PercentileInterval(98), stretch=LinearStretch()
        )

    ax.imshow(kernel, cmap='viridis', interpolation=interpolation, origin='lower', norm=norm)
    ax.set_title(title)
    ax.set_xlabel('X [pixels]')
    ax.set_ylabel('Y [pixels]')

    if return_ax:
        return ax
    
    if save:
        try:
            plt.savefig(save_path)
            logger.info(f"Kernel plot saved to {save_path}")
        except Exception as e:
            logger.error(f"Error saving kernel plot: {e}")