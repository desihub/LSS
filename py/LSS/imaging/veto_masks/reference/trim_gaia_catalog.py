# Select G<18 stars within 1.5 degree of a brick center

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


# gaia_dir = '/project/projectdirs/cosmo/data/gaia/dr2/healpix'
# output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/gaia_dr2_g_18_dr9.fits'

gaia_dir = '/project/projectdirs/cosmo/data/gaia/edr3/healpix'
output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_dr9.fits'

search_radius = 1.5 * 3600.
search_radius_init = 4 * 3600.  # search radius for the initial healpix selection

bricks_south = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/south/survey-bricks-dr9-south.fits.gz'))
bricks_north = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/north/survey-bricks-dr9-north.fits.gz'))
bricks = vstack([bricks_south, bricks_north])
_, idx = np.unique(bricks['brickid'], return_index=True)
bricks = bricks[idx]
print(len(bricks))

ra2 = np.array(bricks['ra'])
dec2 = np.array(bricks['dec'])
sky2 = SkyCoord(ra2*u.degree, dec2*u.degree, frame='icrs')

gaia_nside = 32
print('Healpix pixel size (square deg): {:.5f}'.format(hp.nside2pixarea(gaia_nside, degrees=True)))

gaia_npix = hp.nside2npix(gaia_nside)
gaia_hp_ra, gaia_hp_dec = hp.pix2ang(gaia_nside, np.arange(gaia_npix), nest=True, lonlat=True)

# Select healpix pixels within 4 degrees of a brick center
sky1 = SkyCoord(gaia_hp_ra*u.degree, gaia_hp_dec*u.degree, frame='icrs')
idx1, _, _, _ = sky2.search_around_sky(sky1, seplimit=search_radius_init*u.arcsec)
gaia_list = np.unique(idx1)
print(len(gaia_list))

gaia = []
for index, hp_index in enumerate(gaia_list):
    gaia_fn = str(hp_index).zfill(5)
    if index%100==0:
        print('{}/{}, {}'.format(index, len(gaia_list), gaia_fn))
    tmp = Table(fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), columns=['RA', 'DEC', 'PHOT_G_MEAN_MAG']))
    mask_mag = tmp['PHOT_G_MEAN_MAG']<18.0
    sky1 = SkyCoord(tmp['RA'][mask_mag]*u.degree, tmp['DEC'][mask_mag]*u.degree, frame='icrs')
    idx1, _, _, _ = sky2.search_around_sky(sky1, seplimit=search_radius*u.arcsec)

    if len(idx1)>0:
        idx = np.where(mask_mag)[0][idx1]
        tmp = Table(fitsio.read(os.path.join(gaia_dir, 'healpix-{}.fits'.format(gaia_fn)), rows=idx, columns=['SOURCE_ID', 'RA', 'DEC', 'PHOT_G_MEAN_MAG', 'PHOT_BP_MEAN_MAG', 'PHOT_RP_MEAN_MAG', 'PHOT_G_MEAN_FLUX_OVER_ERROR', 'ASTROMETRIC_EXCESS_NOISE']))
        gaia.append(tmp)
gaia = vstack(gaia)
print(len(gaia))

mask = gaia['PHOT_G_MEAN_FLUX_OVER_ERROR']>0
gaia = gaia[mask]
print(len(gaia))

gaia.write(output_path, overwrite=True)
