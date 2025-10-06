from nirphot.photometry.compute import ComputePhotometry

detection_filter = "f444w"

sci_images = {
    "f090w": "/Users/ng27753/Data/plckg165/convolved/f090w_convolved_to_f444w_median_filter_20250709.fits",
    "f150w": "/Users/ng27753/Data/plckg165/convolved/f150w_convolved_to_f444w_median_filter_20250709.fits",
    "f200w": "/Users/ng27753/Data/plckg165/convolved/f200w_convolved_to_f444w_median_filter_20250709.fits",
    "f277w": "/Users/ng27753/Data/plckg165/convolved/f277w_convolved_to_f444w_median_filter_20250709.fits",
    "f356w": "/Users/ng27753/Data/plckg165/convolved/f356w_convolved_to_f444w_median_filter_20250709.fits",
    "f444w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f444w_30mas_20250527_drz.fits",
}

wht_images = {
    "f090w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f090w_30mas_20250527_wht.fits",
    "f150w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f150w_30mas_20250527_wht.fits",
    "f200w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f200w_30mas_20250527_wht.fits",
    "f277w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f277w_30mas_20250527_wht.fits",
    "f356w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f356w_30mas_20250527_wht.fits",
    "f444w": "/Users/ng27753/Data/plckg165/20250527/ep5/30mas/mosaic_plckg165_ep5_jwst_nircam_f444w_30mas_20250527_wht.fits",
}

exp_times = {
    "f090w": 2490.932,
    "f150w": 1889.672,
    "f200w": 2104.408,
    "f277w": 2104.408,
    "f356w": 1889.672,
    "f444w": 2490.932,
}

cp = ComputePhotometry(
    detection_filter, sci_images, wht_images, exp_times, run="test1_ep5_psf_conv_f444w"
)
cp.get_multiband_photometry()
