from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
# import matplotlib.pyplot as plt
import numpy as np
from astropy.table import Table, vstack, hstack
import fitsio
from astropy.io import fits

from scipy.interpolate import interp1d


output_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-lrg_mask_v1.fits'

# WISE mask
w1_mags = [0, 0.5, 1, 1.5, 2, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0]
w1_radii = [600, 600, 550, 500, 475, 425, 400, 400, 390, 392.5, 395, 370, 360, 330, 275, 240, 210, 165, 100, 75, 60]
w1_max_mag = 10.0

f_radius = interp1d(w1_mags, w1_radii, bounds_error=False, fill_value='extrapolate')

wise_path = '/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-13.3-dr9.fits'
wise = Table(fitsio.read(wise_path))
# print(len(wise))

wise['w1ab'] = np.array(wise['W1MPRO']) + 2.699
mask = wise['w1ab']<w1_max_mag

wise['radius'] = 0.
wise['radius'][mask] = f_radius(wise['w1ab'][mask])
wise.write(output_path)
