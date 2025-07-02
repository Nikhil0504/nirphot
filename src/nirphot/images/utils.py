import logging
from scipy.ndimage import rotate
import numpy as np
import astropy

logger = logging.getLogger(__name__)


def rotate_image(
    image: np.ndarray,
    angle: float,
    interp_order: int = 1,
    reshape: bool = False,
    prefilter: bool = False,
):
    """
    Rotate a 2D image array by a specified angle.

    Parameters
    ----------
    image : numpy.ndarray
        The input 2D image array to rotate.
    angle : float
        The angle (in degrees) to rotate the image. Positive values rotate counter-clockwise.
    interp_order : int, optional
        The order of the spline interpolation. Default is 1 (bilinear).
    reshape : bool, optional
        If True, the output shape is adapted so that the input array is contained entirely in the output.
        If False, the output shape matches the input shape. Default is False.
    prefilter : bool, optional
        Determines if the input is prefiltered with spline_filter before interpolation. Default is False.

    Returns
    -------
    numpy.ndarray
        The rotated image array.
    """
    logger.info(f"Rotating image by {angle} degrees")
    logger.debug(f"Interpolation order: {interp_order}")
    logger.debug(f"Reshape: {reshape}")
    logger.debug(f"Prefilter: {prefilter}")
    return rotate(
        image, -1.0 * angle, order=interp_order, reshape=reshape, prefilter=prefilter
    )


def get_pixscl(hdu: astropy.io.fits.ImageHDU, ext: int | str = 0):
    """
    Calculate the pixel scale of a FITS image.

    For JWST it's sqrt(|CDELT1 * CDELT2| - 0) * 60 * 60


    Parameters
    ----------
    hdu : astropy.io.fits.ImageHDU
        The FITS image HDU.
    ext : int|str, optional
        The extension number or name of the FITS image. Default is 0.

    Returns
    -------
    float
        The pixel scale in arcseconds per pixel.
    """
    try:
        logger.info(f"Calculating pixel scale for extension {ext}")
        pixscl = (
            np.sqrt(abs(hdu[ext].header["CDELT1"] * hdu[ext].header["CDELT2"]) - abs(0))
            * 60
            * 60
        )
        logger.debug(f"Pixel scale calculated: {pixscl} arcseconds/pixel")
    except KeyError:
        logger.error(f"CDELT1 or CDELT2 not found in header for extension {ext}")
        raise ValueError("CDELT1 or CDELT2 not found in header")
    return pixscl
