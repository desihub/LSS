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
from match_coord import search_around


##################################### Tycho-2 stars missing from GAIA EDR3 #####################################

tycho2 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/tycho2-reference-dr9.fits'))
print(len(tycho2))

mask = tycho2['MAG_VT']<10.
mask &= (tycho2['MAG_VT']!=0) | (tycho2['MAG_HP']!=0)
print(np.sum(mask)/len(mask))
tycho2 = tycho2[mask]

gaia_pm_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_14_pm_dr9.fits'
gaia_pm = Table(fitsio.read(gaia_pm_path))

mask = (tycho2['EPOCH_RA']!=0) & (tycho2['EPOCH_DEC']!=0)
tycho_ra_j1991 = tycho2['RA'].copy()
tycho_dec_j1991 = tycho2['DEC'].copy()
tycho_ra_j1991[mask] = (tycho2['RA'] - (tycho2['EPOCH_RA']-1991.) * tycho2['PM_RA'] * 1/3600 / np.cos(np.radians(tycho2['DEC'])))[mask]
tycho_dec_j1991[mask] = (tycho2['DEC'] - (tycho2['EPOCH_DEC']-1991.) * tycho2['PM_DEC'] * 1/3600)[mask]

gaia_pm['RA_J1991'] = gaia_pm['RA'] - 16 * gaia_pm['PMRA'] * 1e-3/3600 / np.cos(np.radians(gaia_pm['DEC']))
gaia_pm['DEC_J1991'] = gaia_pm['DEC'] - 16 * gaia_pm['PMDEC'] * 1e-3/3600

idx1, idx2, d2d, d_ra, d_dec = search_around(gaia_pm['RA_J1991'], gaia_pm['DEC_J1991'], tycho_ra_j1991, tycho_dec_j1991, search_radius=1.)

mask_missing = np.full(len(tycho2), True)
mask_missing[idx2] = False
print(np.sum(mask_missing), np.sum(mask_missing)/len(tycho2))

tycho2 = tycho2[mask_missing]

tycho2['mask_mag'] = -99
tycho2['mask_mag'] = np.min([tycho2['ggguess'], tycho2['zguess']+1], axis=0)

mask = np.isnan(tycho2['zguess'])
mask_vt = mask & (tycho2['MAG_VT']!=0)
mask_hp = mask & (tycho2['MAG_VT']==0) & (tycho2['MAG_HP']!=0)
tycho2['mask_mag'][mask_vt] = tycho2['MAG_VT'][mask_vt]
tycho2['mask_mag'][mask_hp] = tycho2['MAG_HP'][mask_hp]
print(np.sum(tycho2['mask_mag']==-99), np.sum(np.isnan(tycho2['mask_mag'])))

##################################### GAIA EDR3 #####################################

gaia_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_dr9.fits'
gaia_decam_mags_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_edr3_g_18_predict_decam_dr9.fits'

gaia_columns = ['RA', 'DEC', 'PHOT_G_MEAN_MAG', 'PHOT_G_MEAN_FLUX_OVER_ERROR']

gaia = Table(fitsio.read(gaia_path, columns=gaia_columns))
gaia_decam_mags = Table(fitsio.read(gaia_decam_mags_path))
gaia = hstack([gaia, gaia_decam_mags], join_type='exact')
print(len(gaia))

gaia['mask_mag'] = np.min([gaia['PHOT_G_MEAN_MAG'], gaia['decam_mag_z']+1], axis=0)
mask = np.isnan(gaia['decam_mag_z'])
gaia['mask_mag'][mask] = gaia['PHOT_G_MEAN_MAG'][mask]
print(np.sum(np.isnan(gaia['mask_mag'])))

############################################################################################################

tycho2 = tycho2[['RA', 'DEC', 'mask_mag', 'ggguess', 'zguess']]
tycho2.rename_column('ggguess', 'gaia_g')
tycho2['is_tycho2'] = True

gaia = gaia[['RA', 'DEC', 'mask_mag', 'PHOT_G_MEAN_MAG', 'decam_mag_z']]
gaia.rename_columns(['PHOT_G_MEAN_MAG', 'decam_mag_z'], ['gaia_g', 'zguess'])
gaia['is_tycho2'] = False

gaia = vstack([tycho2, gaia], join_type='exact')

gaia.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits')

