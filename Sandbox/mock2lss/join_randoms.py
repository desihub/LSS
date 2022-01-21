import astropy.io.fits as pf
from astropy.table import Table, hstack, vstack, Column
import numpy as np
import random

def select_random_Ns(lst, n):
    random.shuffle(lst)
    result = []
    for i in range(0, len(lst), n):
        result.append(lst[i:i + n])
    return result


SNUMS = np.linspace(100, 5000, num=50, dtype=np.int)

sets = select_random_Ns(SNUMS, 5)

file_name = 'mockRandom_{SNUM}_FirstGen_CutSky_alltracers_sv3bits.fits'


info_randoms = open('info_join_randoms','w')

for i,s in enumerate(sets):
    fits_tables = []
    maxT = 0
    info_randoms.write('%d '%i)
    for f in s:
        info_randoms.write('%d '%f)
        print('------------', f)
        temp = Table.read(file_name.format(SNUM=f), hdu=1)
        print('First ', temp['TARGETID'])
        temp['TARGETID'] += maxT
        maxT = np.max(temp['TARGETID'])
        print('After ', temp['TARGETID'])

        fits_tables.append(temp)
    info_randoms.write('\n')
    new = vstack(fits_tables)
    new.write('mockRandom_5X_{ID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(ID=i), overwrite = True)
    pf.setval('mockRandom_5X_{ID}_FirstGen_CutSky_alltracers_sv3bits.fits'.format(ID=i), 'EXTNAME', value='TARGETS', ext=1)
info_randoms.close()
