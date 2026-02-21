# Create Gaia reference catalog

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits


output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_lrg_mask_v1.fits'

from scipy.interpolate import interp1d

max_mag = 18.

# South
mags = np.array([4.0, 9.0, 10.0, 10.5, 11.5, 12.0, 12.5, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 17.0, 18.0])
radii = np.array([429.18637985, 80.95037032, 57.98737129, 36.80882682,
        26.36735446, 25.29190318, 21.40616169, 15.33392671,
        13.74150366, 13.56870306, 12.03092488, 11.10823009,
         9.79334208, 7.01528803, 5.02527796])
log_radii = np.log10(radii)
f_radius_log_south = interp1d(mags, log_radii, bounds_error=False, fill_value='extrapolate')
f_radius_south = lambda mags: 10**f_radius_log_south(mags)

# North
mags = np.array([4.0, 9.0, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0, 14.5, 15.0, 15.5, 16.0, 17.0, 18.0])
radii = np.array([429.18637985, 80.95037032, 60., 60.,
        60., 47.46123803, 38.68173428, 32.73883553,
        27.70897871, 23.45188791, 19.84883862, 16.79934664,
        13.67150555, 11.57107301, 7.83467367, 5.61223042,
         4.02022236])
log_radii = np.log10(radii)
f_radius_log_north = interp1d(mags, log_radii, bounds_error=False, fill_value='extrapolate')
f_radius_north = lambda mags: 10**f_radius_log_north(mags)

gaia_columns = ['RA', 'DEC', 'mask_mag']
gaia_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_dr9.fits'
gaia_suppl_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_reference_suppl_dr9.fits'

hdu = fits.open(gaia_path)
gaia = Table()
for col in gaia_columns:
    gaia[col] = np.copy(hdu[1].data[col])

suppl = Table(fitsio.read(gaia_suppl_path, columns=['ra', 'dec', 'mask_mag', 'istycho', 'zguess', 'gaia_phot_g_mean_mag']))
suppl.rename_columns(['ra', 'dec', 'istycho', 'gaia_phot_g_mean_mag'], ['RA', 'DEC', 'is_tycho2', 'gaia_g'])
suppl = suppl[gaia_columns]
gaia = vstack([gaia, suppl], join_type='exact')

mask = gaia['mask_mag']<max_mag

gaia['radius_south'] = 0.
gaia['radius_north'] = 0.

gaia['radius_south'][mask] = f_radius_south(gaia['mask_mag'][mask])
gaia['radius_north'][mask] = f_radius_north(gaia['mask_mag'][mask])

gaia.write(output_path)
