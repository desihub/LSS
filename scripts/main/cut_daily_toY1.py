#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import fitsio
import glob
import argparse
import gc
#gc.enable()
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
import healpy as hp

#import tracemalloc

#tracemalloc.start()

from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 



mainp = main(type)

mt = mainp.mtld
tiles = mainp.tiles

datemax = 20220620
wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == 'dark'
wd &= mt['LASTNIGHT'] <= datemax
    
mtld = mt[wd]

wb = mt['SURVEY'] == 'main'
wb &= mt['ZDONE'] == 'true'
wb &= mt['FAPRGRM'] == 'bright'
wb &= mt['LASTNIGHT'] <= datemax

mtlb = mt[wb]

tps_dark = ['QSO','LRG','ELG','ELG_LOP','ELG_LOPnotqso']

for tp in tps_dark:
    dd = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/datcomb_'+tp+'_tarspecwdup_zdone.fits')
    sel = np.isin(dd['TILEID'],mtld['TILEID'])
    print(tp,len(dd),len(dd[sel]))
    fno = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/daily/datcomb_'+tp+'_tarspecwdup_zdone.fits'
    common.write_LSS(dd[sel],fno)
    
tps_bright = ['BGS_ANY','BGS_BRIGHT']
    
for tp in tps_bright:
    dd = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/datcomb_'+tp+'_tarspecwdup_zdone.fits')
    sel = np.isin(dd['TILEID'],mtlb['TILEID'])
    print(tp,len(dd),len(dd[sel]))
    fno = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/daily/datcomb_'+tp+'_tarspecwdup_zdone.fits'
    common.write_LSS(dd[sel],fno)

#print('found 