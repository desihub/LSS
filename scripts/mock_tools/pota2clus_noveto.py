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
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
#import LSS.mkCat_singletile.fa4lsscat as fa
#from LSS.globals import main

if os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--tracer", help="tracer type to be selected")
parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_output", help="base directory for output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/test/')
parser.add_argument("--random_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/')

parser.add_argument("--minr", help="minimum number for random files",default=1,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=2,type=int) 

args = parser.parse_args()
print(args)

rm = int(args.minr)
rx = int(args.maxr)

notqso = ''
#if args.notqso == 'y':
#    notqso = 'notqso'

tracer = args.tracer

if tracer == 'LRG':
    zmin = 0.4
    zmax = 1.1

in_data_fn = args.mockdir+'pota_'+tracer+'.fits'
out_data_fn = args.base_output+tracer+'_complete_noveto_clustering.dat.fits'
mock_data = fitsio.read(in_data_fn)
selcoll = mock_data['COLLISION'] == False
mock_data = mock_data[selcoll]
selz = mock_data['RSDZ'] > zmin
selz &= mock_data['RSDZ'] < zmax
mock_data = mock_data[selz]
mock_data = Table(mock_data)
mock_data.rename_column('RSDZ', 'Z')
common.write_LSS(mock_data,out_data_fn)

ran_samp_cols = ['Z','WEIGHT']

def ran_col_assign(randoms,data,sample_columns)
    inds = np.random.choice(len(data),len(randoms))
    dshuf = data[inds]
    for col in sample_columns:
        randoms[col] =  dshuf[col]
    return randoms

for rann in range(rm,rx):
    in_ran_fn = args.random_dir+tracer+'_'+str(rann)+'_full_noveto.ran.fits'
    out_ran_fn = args.base_output+tracer++'_complete_noveto_'+str(rann)+'clustering.ran.fits'
    ran = Table(fitsio.read(in_ran_fn,columns=['RA','DEC']))
    ran = ran_col_assign(ran,mock_data,ran_samp_cols)
    common.write_LSS(ran,out_ran_fn)




