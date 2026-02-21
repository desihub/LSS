# Add predicted DECam mags

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits


# Coefficients for EDR3
coeffs = dict(
    g = [-0.1125681175, 0.3506376997, 0.9082025788, -1.0078309266,
        -1.4212131445, 4.5685722177, -4.5719415419, 2.3816887292,
        -0.7162270722, 0.1247021438, -0.0114938710, 0.0003949585,
        0.0000051647],
    r = [0.1431278873, -0.2999797766, -0.0553379742, 0.1544273115,
        0.3068634689, -0.9499143903, 0.9769739362, -0.4926704528,
        0.1272539574, -0.0133178183, -0.0008153813, 0.0003094116,
        -0.0000198891],
    z = [0.5173814296, -1.0450176704, 0.1529797809, 0.1856005222,
        -0.2366580132, 0.1018331214, -0.0189673240, 0.0012988354])

bprp_min, bprp_max = -0.5, 4.7


gaia = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_dr9.fits'))

for band in ['g', 'r', 'z']:
    mag = np.copy(gaia['PHOT_G_MEAN_MAG'])
    for order, c in enumerate(coeffs[band]):
        x = gaia['PHOT_BP_MEAN_MAG']-gaia['PHOT_RP_MEAN_MAG']
        x = np.clip(x, bprp_min, bprp_max)
        mag += c * (x)**order
    gaia['decam_mag_'+band] = mag
    
mask = (gaia['PHOT_BP_MEAN_MAG']==0) | (gaia['PHOT_RP_MEAN_MAG']==0)
for band in ['g', 'r', 'z']:
    gaia['decam_mag_'+band][mask] = np.nan

gaia = gaia[['decam_mag_g', 'decam_mag_r', 'decam_mag_z']]
gaia.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_predict_decam_dr9.fits')

