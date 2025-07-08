from nirphot.photometry.catalog import PhotometryDetectionCatalog
from nirphot.photometry.utils import fluxes2mags
from nirphot.utils import generate_filename
from nirphot.paths import DATA_DIR

import numpy as np

class ComputePhotometry:
    def __init__(self, detection_filter, sci_images, wht_images, exposure_times, run=''):
        self.detection_filter = detection_filter
        self.sci_images = sci_images
        self.wht_images = wht_images
        self.exposure_times = exposure_times
        self.run = run
        
        self.filters = sci_images.keys()
        
        self.detection_image = self.sci_images[self.detection_filter]
        self.wht_image = self.wht_images[self.detection_filter]
        self.exposure_time = self.exposure_times[self.detection_filter]
 
        self.detection_catalog = self._run_detection_catalog(self.detection_image, self.wht_image, self.exposure_time, detection_mode=True, save=True)
        
    def _run_detection_catalog(self, sci_img, wht_img, exposure_time, detection_mode=True, segmentation_map=None, save=False):
        catalog = PhotometryDetectionCatalog(sci_img, wht_img, run=self.run, save=save)
        catalog.measure_background()
        
        if segmentation_map is not None and not detection_mode:
            # use the provided segmentation map
            catalog.segmentation_map = segmentation_map
        elif detection_mode and segmentation_map is None:
            # creates the segmentation map
            catalog.detect_sources()
            catalog.segmentation_map.remove_masked_labels(catalog.detection_img_mask)
        else:
            raise ValueError("Provide either segmentation_map or enable detection_mode, not both.")
            
        catalog.measure_source_properties(exposure_time)
        return catalog
    
    def get_multiband_photometry(self):
        # isophotal aperture photometry
        self.source_table = self.detection_catalog.catalog_table
        
        for filter in self.filters:
            sci_img_file = self.sci_images[filter]
            wht_img_file = self.wht_images[filter]
            exposure_time = self.exposure_times[filter]
            
            filter_catalog = self._run_detection_catalog(
                sci_img_file, 
                wht_img_file, 
                exposure_time,
                detection_mode=False, 
                segmentation_map=self.detection_catalog.segmentation_map,
            )
            
            
            self.source_table[filter+'_flux']    = filter_catalog.catalog_table['segment_flux']
            self.source_table[filter+'_fluxerr'] = filter_catalog.catalog_table['segment_fluxerr']
        
            self.source_table[filter+'_mag']     = filter_catalog.catalog_table['segment_mag']
            self.source_table[filter+'_magerr']  = filter_catalog.catalog_table['segment_magerr']
            
        self._add_flux_corrections()
        
        fp = generate_filename(f'{self.run}_isophotal_catalog', 'ecsv', DATA_DIR)
        self.source_table.write(fp, format='ascii.ecsv', overwrite=True)
        
    
    def _add_flux_corrections(self):
        reference_flux_auto = self.source_table['kron_flux']    # Kron total flux estimate
        reference_flux_iso  = self.source_table['segment_flux'] # flux in isophotal aperture defined by detection segment
        
        kron_flux_corrections = reference_flux_auto / reference_flux_iso
        self.source_table['total_flux_cor'] = kron_flux_corrections
        
        for filter in self.filters:
            self.source_table[filter+'_flux']    *= kron_flux_corrections
            self.source_table[filter+'_fluxerr'] *= kron_flux_corrections
            #total_flux_table[filt+'_mag'] = total_flux_table[filt+'_flux'].to(u.ABmag)  # doesn't handle non-detections
            self.source_table[filter+'_mag'] = fluxes2mags(self.source_table[filter+'_flux'], self.source_table[filter+'_fluxerr'])[0]
            # magnitude uncertainty magerr stays the same
    