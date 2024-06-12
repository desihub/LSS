import numpy as np
from astropy.table import Table,vstack
import os

program = 'dark'

rmin = 0
rmax = 64

#path = '/pscratch/sd/a/acarnero/test_main/altmtl{MOCKNUM}/Univ000' 
path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/altmtl21_R64/Univ{MOCKNUM}' 

extratiles = Table.read('/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4/aux_data/extra_{PROGRAM}.ecsv'.format(PROGRAM = program), format='ascii.ecsv')

tileref = extratiles['TILEID'][-1]
print(tileref)

for i in range(rmin, rmax):
    input_track = os.path.join(path, 'mainsurvey-{PRG}obscon-TileTracker.ecsv').format(MOCKNUM = "%03d" % i, PRG=program.upper())
    tiles = Table.read(input_track, format='ascii.ecsv')
    tiles.meta['amtldir'] = path.format(MOCKNUM = i)
    if tiles['TILEID'][-1] != tileref:
        print('merging for mock', i)
        newtable = vstack([tiles, extratiles])
        newtable.meta = tiles.meta
        newtable.write(input_track, overwrite=True)
    else:
        print(i, 'it has already been merged')
