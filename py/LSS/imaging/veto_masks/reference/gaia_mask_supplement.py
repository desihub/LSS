# Catalog of DR9 reference stars that have mask_mag<8 and are brighter than the new reference catalog by more than 0.05

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
from match_coord import search_around, match_coord

# DR9 bright star mask catalog
dr9ref = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/gaia-mask-dr9.fits', columns=['mask_mag']))
mask = dr9ref['mask_mag']<8
idx = np.where(mask)[0]
dr9ref = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/tmp/gaia-mask-dr9.fits', rows=idx))
print(len(dr9ref))

# New mask catalog
cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits', columns=['mask_mag']))
mask = cat['mask_mag']<10
idx = np.where(mask)[0]
cat = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits', rows=idx))


idx1, idx2, d2d, d_ra, d_dec = search_around(cat['RA'], cat['DEC'], dr9ref['ra'], dr9ref['dec'], search_radius=5.)
print(len(dr9ref)-len(np.unique(idx2)))

# Find the brightest mag from new catalog for each for the old catalog star
idx2_unique = np.unique(idx2)
mask_mag = np.zeros(len(idx2_unique))
for ii, index2 in enumerate(idx2_unique):
    mask = idx2==index2
    mask_mag[ii] = np.min(cat['mask_mag'][idx1[mask]])

mask = mask_mag-dr9ref['mask_mag'][idx2_unique] > 0.05
dr9ref_add = dr9ref[idx2_unique[mask]].copy()
print(len(dr9ref_add))

dr9ref_add.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_suppl_dr9.fits')

