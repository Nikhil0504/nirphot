import numpy as np
import astropy.units as u


# not detected: mag =  99; magerr = 1-sigma upper limit assuming zero flux
# not observed: mag = -99; magerr = 0
def fluxes2mags(flux, fluxerr):
    nondet = flux < 0  # Non-detection if flux is negative
    unobs = (fluxerr <= 0) + (fluxerr == np.inf)  # Unobserved if flux uncertainty is negative or infinity

    mag = flux.to(u.ABmag)
    magupperlimit = fluxerr.to(u.ABmag)  # 1-sigma upper limit if flux=0

    mag = np.where(nondet, 99 * u.ABmag, mag)
    mag = np.where(unobs, -99 * u.ABmag, mag)

    magerr = 2.5 * np.log10(1 + fluxerr / flux)
    magerr = magerr.value * u.ABmag

    magerr = np.where(nondet, magupperlimit, magerr)
    magerr = np.where(unobs, 0 * u.ABmag, magerr)

    return mag, magerr


PROPERTIES = [
    "label",
    "xcentroid",
    "ycentroid",
    "sky_centroid",
    "area",
    "semimajor_sigma",
    "semiminor_sigma",
    "fwhm",
    "ellipticity",
    "orientation",
    "gini",
    "kron_radius",
    "local_background",
    "segment_flux",
    "segment_fluxerr",
    "kron_flux",
    "kron_fluxerr",
]
