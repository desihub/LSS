# Predict the DECam z and GAIA G magnitudes using Tycho-2 and 2MASS photometry

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

sys.path.append(os.path.expanduser('~/git/Python/user_modules/'))
import match_coord


tycho2_path = '/global/project/projectdirs/cosmo/staging/tycho2/tycho2.kd.fits'
twomass_path = '/global/cfs/cdirs/desi/users/rongpu/useful/2mass_psc/2mass_psc_j_12.fits'

output_path = '/global/cfs/cdirs/desi/users/rongpu/useful/tycho2-reference.fits'

tycho2 = Table(fitsio.read(tycho2_path))
twomass = Table(fitsio.read(twomass_path))

mask_bad = tycho2['MAG_VT']==0

mask = (tycho2['EPOCH_RA']!=0) & (tycho2['EPOCH_DEC']!=0)
tycho_ra_j2000 = tycho2['RA'].copy()
tycho_dec_j2000 = tycho2['DEC'].copy()
tycho_ra_j2000[mask] = (tycho2['RA'] - (tycho2['EPOCH_RA']-2000.) * tycho2['PM_RA'] * 1/3600 / np.cos(np.radians(tycho2['DEC'])))[mask]
tycho_dec_j2000[mask] = (tycho2['DEC'] - (tycho2['EPOCH_DEC']-2000.) * tycho2['PM_DEC'] * 1/3600)[mask]

# Give invalid mags lowest priority
mag_vt = tycho2['MAG_VT'].copy()
mag_vt[mask_bad] = 99.
mask_hp = mask_bad & (tycho2['MAG_HP']!=0)  # use MAG_HP (Hipparcos mag) when VT is unavailable
mag_vt[mask_hp] = tycho2['MAG_HP'][mask_hp]

idx1, idx2, d2d, d_ra, d_dec = match_coord.match_coord(twomass['RAJ2000'], twomass['DEJ2000'], tycho_ra_j2000, tycho_dec_j2000, priority2=-mag_vt, search_radius=5., plot_q=False)
print(len(idx1)/len(twomass))
print(len(idx1)/len(tycho2))

tycho2['Jmag'] = np.nan
tycho2['Hmag'] = np.nan
tycho2['Kmag'] = np.nan
tycho2['zguess'] = np.nan
tycho2['ggguess'] = np.nan

twomass = twomass[idx1]
tycho2['Jmag'][idx2] = twomass['Jmag']
tycho2['Hmag'][idx2] = twomass['Hmag']
tycho2['Kmag'][idx2] = twomass['Kmag']

coeffs_z = [-0.01835938, -0.68084937, 0.49222576]
coeffs_gg = [0.00445346, -0.07819228, -0.07145574, 0.00278177]
xmin, xmax = -1, 8

x = tycho2['MAG_VT'][idx2]-twomass['Jmag']

pz = np.poly1d(coeffs_z)
tycho2['zguess'][idx2] = pz(np.clip(x, xmin, xmax)) + tycho2['MAG_VT'][idx2]
tycho2['zguess'][mask_bad] = np.nan

pgg = np.poly1d(coeffs_gg)
tycho2['ggguess'][idx2] = pgg(np.clip(x, xmin, xmax)) + tycho2['MAG_VT'][idx2]
tycho2['ggguess'][mask_bad] = np.nan

tycho2.write(output_path, overwrite=True)
