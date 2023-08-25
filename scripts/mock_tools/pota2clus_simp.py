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
parser.add_argument("--realization",type=int)
parser.add_argument("--prog", default="DARK")
parser.add_argument("--veto",default='_noveto')
#parser.add_argument("--mockdir", help="directory when pota mock data is",default='/global/cfs/cdirs/desi/users/acarnero/y1mock/SecondGen/clustering/')
parser.add_argument("--base_dir", help="base directory for input/output",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit/')
parser.add_argument("--random_dir",help="where to find the data randoms",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/LSScats/v0.6/')
parser.add_argument("--specdata_dir",help="where to find the spec data ",default='/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/iron/')
parser.add_argument("--minr", help="minimum number for random files",default=0,type=int)
parser.add_argument("--maxr", help="maximum for random files, default is 1, but 40 are available (use parallel script for all)",default=1,type=int) 

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

elif tracer == 'ELG_LOP':
    zmin = 0.8
    zmax = 1.6

elif tracer == 'QSO':
    zmin = 0.8
    zmax = 2.1

else:
    sys.exit('tracer type '+args.tracer+' not supported (yet)')


mockdir = args.base_dir+'mock'+str(args.realization)+'/'
in_data_fn = mockdir+'pota-'+args.prog+'.fits'
print(in_data_fn)
out_data_fn = mockdir+tracer+'_complete'+args.veto+'_clustering.dat.fits'
out_data_froot = mockdir+tracer+'_complete'+args.veto+'_'
cols = ['LOCATION',
 'FIBER',
 'TARGETID',
 'RA',
 'DEC','RSDZ',
 'PRIORITY_INIT',
 'PRIORITY',
 'DESI_TARGET','BRICKID','NOBS_G',
 'NOBS_R',
 'NOBS_Z',
 'MASKBITS','ZWARN',
 'COLLISION',
 'TILEID']
mock_data = fitsio.read(in_data_fn,columns=cols)
selcoll = mock_data['COLLISION'] == False
mock_data = mock_data[selcoll]

if args.prog == 'DARK':
    bit = targetmask.desi_mask[args.tracer]
    desitarg='DESI_TARGET'

ndattot = len(mock_data)
seltar = mock_data[desitarg] & bit > 0
mock_data = mock_data[seltar]
lmockdat_noveto = len(mock_data)
print('length before/after cut to target type '+args.tracer)
print(ndattot,lmockdata_noveto)

'''
PUT IN SOMETHING HERE TO MASK TO GOODHARDLOC AS AN OPTION
'''
if 'gtl' in args.veto:
    tilelocid = 10000*mock_data['TILEID']+mock_data['LOCATION']
    specfo = args.specdata_dir+'datcomb_'+args.program.lower()+'_spec_zdone.fits'
    print('loading specf file '+specfo)
    specf = Table(fitsio.read(specfo))
    print(len(np.unique(specf['TILEID'])))
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    print('loaded specf file '+specfo)
    specfc = common.cut_specdat(specf,badfib=mainp.badfib)
    gtl = np.unique(specfc['TILELOCID'])
    goodtl = np.isin(tilelocid,gtl)
    mock_data = mock_data[goodtl]
    print(lmock_data_noveto,len(mock_data))
    

selz = mock_data['RSDZ'] > zmin
selz &= mock_data['RSDZ'] < zmax
mock_data = mock_data[selz]
mock_data = Table(mock_data)
mock_data = unique(mock_data,keys=['TARGETID'])
print('length after cutting to redshift and unique targetid',len(mock_data))
mock_data.rename_column('RSDZ', 'Z')
mock_data['WEIGHT'] = 1
common.write_LSS(mock_data,out_data_fn)

def splitGC(flroot,datran='.dat',rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+'.fits'
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+'.fits'

    fn = Table(fitsio.read(flroot+app))
    c = SkyCoord(fn['RA']* u.deg,fn['DEC']* u.deg,frame='icrs')
    gc = c.transform_to('galactic')
    sel_ngc = gc.b > 0
    outf_ngc = flroot+'NGC_'+app
    common.write_LSS(fn[sel_ngc],outf_ngc)
    outf_sgc = flroot+'SGC_'+app
    common.write_LSS(fn[~sel_ngc],outf_sgc)

splitGC(out_data_froot,'.dat')

ran_samp_cols = ['Z','WEIGHT']

def ran_col_assign(randoms,data,sample_columns):
    inds = np.random.choice(len(data),len(randoms))
    dshuf = data[inds]
    for col in sample_columns:
        randoms[col] =  dshuf[col]
    return randoms

for rann in range(rm,rx):
    in_ran_fn = args.random_dir+args.tracer'_'+str(rann)+'_full_noveto.ran.fits' #all noveto have same ra,dec, tracer becomes important for LRG imaging veto
    out_ran_fn = out_data_froot+str(rann)+'_clustering.ran.fits'
    ran = Table(fitsio.read(in_ran_fn,columns=['RA','DEC','TILELOCID']))
    if 'gtl' in args.veto:
        goodtl = np.isin(ran['TILELOCID'],gtl)
        ran = ran[goodtl]
    ran = ran_col_assign(ran,mock_data,ran_samp_cols)
    common.write_LSS(ran,out_ran_fn)
    splitGC(out_data_froot,'.ran',rann)




