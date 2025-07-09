import stpsf as webbpsf
from astropy import time
from astropy.io import fits

import logging
from nirphot.images.utils import get_pixscl, rotate_image
from nirphot.plotting.psfs import save_psf_plots

webbpsf.setup_logging()
logger = logging.getLogger(__name__)


def generate_webbpsf(filter, img_hdu, ext, pixscl, oversample, fov_pixels, fov_arcsec):
    """Generate a WebbPSF for the given filter.

    Parameters:
    -----------
    filter : str
        The filter to generate the WebbPSF for.
    img_hdu : astropy.io.fits.hdu.hdulist.HDUList
        The image HDUList object.
    ext : int
        The FITS extension to extract the pixel scale from.
    pixscl : float
        The pixel scale for the WebbPSF.
    oversample : int
        The oversampling factor for the WebbPSF.
    fov_pixels : int
        The field of view in pixels for the WebbPSF.
    fov_arcsec : float
        The field of view in arcseconds for the WebbPSF.

    Returns:
    --------
    psf : astropy.io.fits.hdu.hdulist.HDUList
        The WebbPSF HDUList object.
    """
    logger.info(f"Generating WebbPSF for filter {filter}")

    nc = webbpsf.NIRCam()
    nc.filter = filter
    nc.pixelscale = pixscl if pixscl else get_pixscl(img_hdu, ext)

    t = time.Time(img_hdu[ext].header["MJD-AVG"], format="mjd")
    nc.load_wss_opd_by_date(t.iso.replace(" ", "T"))

    # output format
    nc.options["output_mode"] = "both"

    try:
        psf = nc.calc_psf(
            oversample=oversample, fov_pixels=fov_pixels, fov_arcsec=fov_arcsec
        )
        logger.info(f"WebbPSF generated for filter {filter}")
        logger.info(
            f"Pixscale for the oversampled PSF: {psf['OVERSAMP'].header['PIXELSCL']}"
        )
        logger.info(
            f"Pixscale for the detector-sampled PSF: {psf['DET_SAMP'].header['PIXELSCL']}"
        )

    except Exception as e:
        logger.error(
            f"Failed to generate WebbPSF for filter {filter}: {e}", exc_info=True
        )
        raise RuntimeError(f"Failed to generate WebbPSF for filter {filter}: {e}")
    return psf


def rotate_webbpsf(psf, img_hdu=None, ext=None, pa=None):
    """Rotate the WebbPSF by the given position angle.

    Parameters:
    -----------
    psf : astropy.io.fits.hdu.hdulist.HDUList
        The WebbPSF HDUList object.
    img_hdu : astropy.io.fits.hdu.hdulist.HDUList
        The image HDUList object.
    ext : int
        The FITS extension to extract the pixel scale from.
    pa : float
        The position angle in degrees.

    Returns:
    --------
    psf : astropy.io.fits.hdu.hdulist.HDUList
        The rotated WebbPSF HDUList object.
    """

    try:
        angle = pa if pa else img_hdu[ext].header["PA_APER"]
        logger.info(f"Rotating WebbPSF by PA={angle}")
        rotated_psf_overdist = rotate_image(psf["OVERDIST"].data, angle, reshape=False)
        rotated_psf_detdist = rotate_image(psf["DET_DIST"].data, angle, reshape=False)

        # add these two files to the hdu list
        psf.append(
            fits.ImageHDU(
                data=rotated_psf_overdist,
                name="ROTATED_OVERDIST",
                header=psf["OVERDIST"].header.copy(),
            )
        )
        psf.append(
            fits.ImageHDU(
                data=rotated_psf_detdist,
                name="ROTATED_DET_DIST",
                header=psf["DET_DIST"].header.copy(),
            )
        )

        # copy the header from the input image HDU to the rotated PSF HDU and add a comment about the rotation
        psf["ROTATED_OVERDIST"].header["COMMENT"] = f"Rotated by PA={angle} degrees"
        psf["ROTATED_DET_DIST"].header["COMMENT"] = f"Rotated by PA={angle} degrees"

    except Exception as e:
        logger.error(f"Failed to rotate WebbPSF by PA={angle}: {e}")
        raise RuntimeError(f"Failed to rotate WebbPSF by PA={angle}: {e}")

    logger.info(f"WebbPSF rotated by PA={angle}")

    return psf
