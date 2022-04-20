# Based on AllWISE+2MASS catalog provided by Aaron Meisner
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

allwise_path = '/global/cfs/cdirs/desi/users/rongpu/useful/unwise_bitmask_bright_stars/w1_bright-2mass.fits.gz'
output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-dr9.fits'

search_radius = 1.5 * 3600.

allwise = Table(fitsio.read(allwise_path))

bricks_south = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/survey-bricks-dr9-south.fits.gz'))
bricks_north = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/survey-bricks-dr9-north.fits.gz'))
bricks = vstack([bricks_south, bricks_north])
_, idx = np.unique(bricks['brickid'], return_index=True)
bricks = bricks[idx]
print(len(bricks))

ra2 = np.array(bricks['ra'])
dec2 = np.array(bricks['dec'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

ra1 = np.array(allwise['RA'])
dec1 = np.array(allwise['DEC'])
sky1 = SkyCoord(allwise['RA']*u.degree, allwise['DEC']*u.degree, frame='icrs')

idx1, _, _, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
idx1 = np.unique(idx1)
allwise = allwise[idx1]
# allwise['idx'] = idx1

allwise.write(output_path, overwrite=True)
