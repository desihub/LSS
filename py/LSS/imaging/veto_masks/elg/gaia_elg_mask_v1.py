# Create Gaia reference catalog

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits


output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/gaia_elg_mask_v1.fits'

from scipy.interpolate import interp1d

max_mag = 18.

mags = np.array([4., 5., 6., 7.] + np.arange(7.5, 18.5, 0.5).tolist())

# South
# radius vs magnitude relation from Mehdi
radii_south = [700., 500., 350., 260., 200.,
               180., 170., 120., 100., 75., 65.,
               55., 50., 40., 30., 25., 22.,
               19., 18., 14., 11., 10.0, 7.5,
               6.0, 5.5, 4.2]
log_radii = np.log10(radii_south)
f_radius_log_south = interp1d(mags, log_radii, bounds_error=False, fill_value='extrapolate')
f_radius_south = lambda mags: np.maximum(10**f_radius_log_south(mags), 1630. * 1.396**(-mags))

# North
# radius vs magnitude relation from Mehdi
radii_north = np.array([690., 590., 550., 510., 310.,
                        270., 260., 250., 150., 120., 95.,
                        85., 70., 65., 40., 35.,
                        30., 25., 20., 18., 15.,
                        12., 10., 10., 8.0, 7.0])
log_radii = np.log10(radii_north)
f_radius_log_north = interp1d(mags, log_radii, bounds_error=False, fill_value='extrapolate')
f_radius_north = lambda mags: np.maximum(10**f_radius_log_north(mags), 1630. * 1.396**(-mags))

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
