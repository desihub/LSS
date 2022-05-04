# Add fainter stars (to W1MPRO=13.3) to the W1+2MASS catalog (which only goes to W1MPRO=9.5)

from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack
import fitsio
# from astropy.io import fits

t1 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-dr9.fits', columns=['RA', 'DEC', 'W1MPRO']))
t2 = Table(fitsio.read('/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-13.3-dr9.fits', columns=['RA', 'DEC', 'W1MPRO']))

mask = t2['W1MPRO']>=9.5
t2 = t2[mask]

cat = vstack([t1, t2])
cat.write('/global/cfs/cdirs/desi/users/rongpu/desi_mask/w1_bright-2mass-13.3-dr9.fits')
