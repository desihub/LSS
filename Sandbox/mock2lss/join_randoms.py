import astropy.io.fits as pf
from astropy.table import Table, hstack, vstack, Column
import numpy as np
files = ['mockRandom_1100_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1300_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1700_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2000_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1900_FirstGen_CutSky_alltracers_sv3bits.fits']
fits_tables = []
maxT = 0
for f in files:
    print('------------', f)
    temp = Table.read(f,hdu=1)
    print('First ', temp['TARGETID'])
    temp['TARGETID'] += maxT
    maxT = np.max(temp['TARGETID'])
    print('After ', temp['TARGETID'])

    fits_tables.append(temp)
new = vstack(fits_tables)
new.write('mockRandom_5X_FirstGen_CutSky_alltracers_sv3bits.fits', overwrite = True)
pf.setval('mockRandom_5X_FirstGen_CutSky_alltracers_sv3bits.fits', 'EXTNAME', value='TARGETS', ext=1)

