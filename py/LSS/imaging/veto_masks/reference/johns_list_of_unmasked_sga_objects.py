from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
# import matplotlib
# matplotlib.use("Agg")
# import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio
# from astropy.io import fits


columns = ['SGA_ID', 'REF_CAT', 'GALAXY', 'RA', 'DEC', 'DIAM']

ee = Table(fitsio.read('/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/masking/SGA-ellipse-v3.0.kd.fits', columns=columns))
nn = Table(fitsio.read('/global/cfs/cdirs/cosmo/staging/largegalaxies/v3.0/SGA-ellipse-v3.0.fits', columns=columns))
old = ee[ee['SGA_ID'] > -1]
new = nn[nn['SGA_ID'] > -1]
old = old[np.argsort(old['SGA_ID'])]
new = new[np.argsort(new['SGA_ID'])]
ww = (old['REF_CAT'] == '') * (new['REF_CAT'] == 'L3')

cat = new[ww]['SGA_ID', 'GALAXY', 'RA', 'DEC', 'DIAM']

t = cat[['RA', 'DEC']]
t['radius'] = cat['DIAM']/2. * 60.
# t.write('/global/u2/r/rongpu/temp/missing_sga.fits', overwrite=True)

for index in range(len(t)):
    print('{}, {}, {}'.format(t['RA'][index], t['DEC'][index], t['radius'][index]))
