import astropy.io.fits as pf
from astropy.table import Table, hstack, vstack, Column
import numpy as np
# first: files = ['mockRandom_1100_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1300_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1700_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2000_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1900_FirstGen_CutSky_alltracers_sv3bits.fits']
# second: files = ['mockRandom_1200_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1400_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1500_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1600_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_1800_FirstGen_CutSky_alltracers_sv3bits.fits']
#files = ['mockRandom_2100_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2200_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2300_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2400_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2500_FirstGen_CutSky_alltracers_sv3bits.fits']
#files = ['mockRandom_2600_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2700_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2800_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_2900_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3000_FirstGen_CutSky_alltracers_sv3bits.fits']
#files = ['mockRandom_3100_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3200_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3300_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3400_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3500_FirstGen_CutSky_alltracers_sv3bits.fits']
#files = ['mockRandom_3600_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3700_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3800_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_3900_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_4000_FirstGen_CutSky_alltracers_sv3bits.fits']
files = ['mockRandom_4600_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_4700_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_4800_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_4900_FirstGen_CutSky_alltracers_sv3bits.fits', 'mockRandom_5000_FirstGen_CutSky_alltracers_sv3bits.fits']

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
new.write('mockRandom_5X_8_FirstGen_CutSky_alltracers_sv3bits.fits', overwrite = True)
pf.setval('mockRandom_5X_8_FirstGen_CutSky_alltracers_sv3bits.fits', 'EXTNAME', value='TARGETS', ext=1)

