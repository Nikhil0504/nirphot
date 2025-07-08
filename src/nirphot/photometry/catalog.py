from astropy.io import fits
from astropy.stats import gaussian_fwhm_to_sigma
from astropy.convolution import Gaussian2DKernel, convolve

import astropy.units as u
import astropy.wcs as wcs  

import numpy as np

from photutils.background import Background2D
from photutils.segmentation import SourceFinder, SourceCatalog, make_2dgaussian_kernel
from photutils.utils import calc_total_error

from nirphot.paths import DATA_DIR, PLOT_DIR
from nirphot.utils import generate_filename

from nirphot.photometry.utils import PROPERTIES, fluxes2mags

import logging

logger = logging.getLogger(__name__)

JWST_UNITS = u.MJy / u.sr

class PhotometryDetectionCatalog:
    def __init__(self, detection_file, weight_file, run='run', verbose=True):
        self.verbose = verbose
        self.run = run

        self.detection_file = detection_file
        self.weight_file = weight_file
        self.detection_hdu = fits.open(self.detection_file)
        self.weight_hdu = fits.open(self.weight_file)

        self.detection_data = self.detection_hdu[0].data
        self.weight_data = self.weight_hdu[0].data
        self.convolved_image_data = self._smooth_img()

        self.imwcs = wcs.WCS(self.detection_hdu[0].header)
        
        self.pixscale = wcs.utils.proj_plane_pixel_scales(self.imwcs)[0] 
        self.pixscale *= self.imwcs.wcs.cunit[0].to('arcsec')
        self.flux_units = JWST_UNITS * (self.pixscale * u.arcsec)**2

    def _smooth_img(self, smooth_fwhm=2, kernel_size=5):
        self.smooth_kernel = make_2dgaussian_kernel(smooth_fwhm, kernel_size)
        self.smooth_kernel.normalize()
        
        return convolve(self.detection_data, self.smooth_kernel)
        

    def measure_background(self, bkg_size=50, filter_size=3):
        if self.verbose:
            logger.info(f"Measuring background with size {bkg_size} and filter size {filter_size}")

        self.background_map = Background2D(self.detection_data, (bkg_size, bkg_size), filter_size=filter_size)

    def detect_sources(self, n_sigma=3, npixels=10):
        # threshold = (sigma * background_rms) + background
        detection_threshold = (n_sigma * self.background_map.background_rms) + self.background_map.background
        finder = SourceFinder(npixels=npixels, progress_bar=True)
        
        if self.verbose:
            logger.info("Detecting sources...")

        self.segmentation_map = finder(self.convolved_image_data, detection_threshold)

        # save segmentation map
        segmentation_hdu = fits.PrimaryHDU(self.segmentation_map)
        segmentation_hdu.header.update(self.imwcs.to_header())
        
        fp = generate_filename(f'{self.run}_segmentation_map', 'fits', DATA_DIR)
        segmentation_hdu.writeto(fp, overwrite=True)
        
        if self.verbose:
            logger.info(f"Detected {self.segmentation_map.nlabels} sources")
            logger.info("Segmentation map saved.")
       
        segmentation_regions = self.segmentation_map.to_regions() 
        segmentation_regions.write(fp.replace('.fits', '.reg'), format="ds9")
        
        
    def measure_source_properties(self, exposure_time, local_background_width=24):
        # "effective_gain" = exposure time map (conversion from data rate units to counts)
        # weight = inverse variance map = 1 / sigma_background**2 (without sources)
        self.exposure_time_map = exposure_time * self.background_map.background_rms_median**2 * self.weight_data
        
        background_rms = 1 / np.sqrt(self.weight_data)
        
        # make sure to not have 0s in the exposure time map
        # If effective_gain is zero (or contains zero values in an array), 
        # then the source Poisson noise component will not be included.
        self.data_rms = calc_total_error(self.detection_data, background_rms, self.exposure_time_map+1e-8)
        
        self.catalog = SourceCatalog(self.detection_data-self.background_map.background, 
            self.segmentation_map, 
            convolved_data=self.convolved_image_data,
            error=self.data_rms,
            background=self.background_map.background,
            wcs=self.imwcs,
            localbkg_width=local_background_width,
            progress_bar=True
        )
        
        self.catalog_table = self.catalog.to_table(columns=PROPERTIES)
        
        for aperture in ['segment', 'kron']:
            flux    = self.catalog_table[aperture+'_flux']    * self.flux_units.to(u.nJy)
            fluxerr = self.catalog_table[aperture+'_fluxerr'] * self.flux_units.to(u.nJy)
            mag, magerr = fluxes2mags(flux, fluxerr)
            
            self.catalog_table[aperture+'_flux']    = flux
            self.catalog_table[aperture+'_fluxerr'] = fluxerr
            self.catalog_table[aperture+'_mag']     = mag
            self.catalog_table[aperture+'_magerr']  = magerr
            
        fp_catalog = generate_filename(f'{self.run}_catalog', 'ecsv', DATA_DIR)
        self.catalog_table.write(fp_catalog, format='ascii.ecsv')
