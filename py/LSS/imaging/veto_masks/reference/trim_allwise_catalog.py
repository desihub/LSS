# Based on AllWISE W1MPRO<13.3 catalog
# Select stars within 1.5 degree of a brick center

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

import healpy as hp
from astropy import units as u
from astropy.coordinates import SkyCoord

allwise_path = '/global/project/projectdirs/desi/users/rongpu/useful/w1_bright-13.3.fits'
output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-13.3-dr9.fits'

search_radius = 1.5 * 3600.
search_radius_init = 4 * 3600.  # search radius for the initial healpix selection

allwise = Table(fitsio.read(allwise_path))
print(len(allwise))
allwise['idx'] = np.arange(len(allwise))

bricks_south = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/survey-bricks-dr9-south.fits.gz'))
bricks_north = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/survey-bricks-dr9-north.fits.gz'))
bricks = vstack([bricks_south, bricks_north])
_, idx = np.unique(bricks['brickid'], return_index=True)
bricks = bricks[idx]
print(len(bricks))

ra2 = np.array(bricks['ra'])
dec2 = np.array(bricks['dec'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

########################## Pre-selection using healpix ##########################
nside = 32
print('Healpix pixel size (square deg): {:.5f}'.format(hp.nside2pixarea(nside, degrees=True)))

npix = hp.nside2npix(nside)
hp_ra, hp_dec = hp.pix2ang(nside, np.arange(npix), nest=True, lonlat=True)

# Select healpix pixels within 4 degrees of a brick center
sky1 = SkyCoord(hp_ra*u.degree, hp_dec*u.degree, frame='icrs')
_, d2d, _ = sky1.match_to_catalog_sky(sky2)
d2d = d2d.to_value('arcsec')
mask = d2d < search_radius_init
hp_idx = np.where(mask)[0]
print(len(hp_idx))

wise_hp_idx = hp.ang2pix(nside, allwise['RA'], allwise['DEC'], nest=True, lonlat=True)
mask = np.in1d(wise_hp_idx, hp_idx)
allwise = allwise[mask]
print(len(allwise))

########################## Find the stars near a brick ##########################

ra1 = np.array(allwise['RA'])
dec1 = np.array(allwise['DEC'])
sky1 = SkyCoord(allwise['RA']*u.degree, allwise['DEC']*u.degree, frame='icrs')

_, d2d, _ = sky1.match_to_catalog_sky(sky2)
d2d = d2d.to_value('arcsec')
mask = d2d < search_radius
allwise = allwise[mask]
print(len(allwise))

allwise.write(output_path, overwrite=True)
