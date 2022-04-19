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

twomass_path = '/global/cfs/cdirs/desi/users/rongpu/useful/2mass_psc/2mass_psc_j_12.fits'
output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/2mass_psc/2mass_psc_j_12-dr9.fits'

search_radius = 1.5 * 3600.

twomass = Table(fitsio.read(twomass_path))
print(len(twomass))

twomass.rename_columns(['RAJ2000', 'DEJ2000'], ['RA', 'DEC'])

bricks_south = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/survey-bricks-dr9-south.fits.gz'))
bricks_north = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/survey-bricks-dr9-north.fits.gz'))
bricks = vstack([bricks_south, bricks_north])
_, idx = np.unique(bricks['brickid'], return_index=True)
bricks = bricks[idx]
print(len(bricks))

ra2 = np.array(bricks['ra'])
dec2 = np.array(bricks['dec'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

ra1 = np.array(twomass['RA'])
dec1 = np.array(twomass['DEC'])
sky1 = SkyCoord(twomass['RA']*u.degree, twomass['DEC']*u.degree, frame='icrs')

idx1, _, _, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)
idx1 = np.unique(idx1)
twomass = twomass[idx1]
twomass['idx'] = idx1
print(len(twomass))

twomass.write(output_path, overwrite=True)
