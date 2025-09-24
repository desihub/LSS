#standard python
import sys
import os
import shutil
import unittest
from datetime import datetime
import json
import numpy as np
import healpy as hp
import fitsio
import glob
import argparse
from astropy.table import Table,join,unique,vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot
from desitarget import targetmask

import logging
# create logger
logname = 'add_night_info'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)

logger.info('script is starting')


#import logging
#logging.getLogger().setLevel(logging.ERROR)


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
from LSS.globals import main
from LSS.qso_cat_utils import qso_catalog_maker,build_qso_catalog_from_healpix,build_qso_catalog_from_tiles

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="the full path of the file you are adding the columns to")
parser.add_argument("--outfile", help="the full path to where output gets saved")
parser.add_argument("--verspec",help="version for redshifts",default='jura')
parser.add_argument("--verspecrel",help="version for redshifts",default='v1')
parser.add_argument("--survey",help="e.g., main, sv3",default='main')




args = parser.parse_args()
print(args)


specrel = args.verspec

#required columns for importing from zcatalogs, add any as needed
columns = ['TARGETID','ZWARN','ZERR','SPECTYPE','TSNR2_QSO', 'TSNR2_LYA','TSNR2_ELG','TSNR2_LRG']


surpipe = args.survey

reldir = '/global/cfs/cdirs/desi/spectro/redux/'+specrel


def add_lastnight(qf,prog='dark'):
    qf['LASTNIGHT'] = np.zeros(len(qf),dtype=int)
    targetid2index = {targetid:index for index,targetid in enumerate(qf["TARGETID"])}
    tilefn = reldir+'/zcatalog/ztile-'+surpipe+'-'+prog+'-cumulative.fits'
    t=Table(fitsio.read(tilefn,columns=['TARGETID','LASTNIGHT']))
    selection=np.in1d(t["TARGETID"],qf["TARGETID"])
    if np.sum(selection)==0 :
        print("no intersection")
    ii=[targetid2index[tid] for tid in t["TARGETID"][selection]]
    qf["LASTNIGHT"][ii] = np.maximum(qf["LASTNIGHT"][ii],t["LASTNIGHT"][selection])
    print(np.sum(qf["LASTNIGHT"]==0),"entries without LASTNIGHT info")
    return qf

def add_fminfo(qf,expinfo):
    # sort expinfo by NIGHT in descending order
    print('getting FIRST info')
    expinfo.sort('NIGHT')
    expinfo_first = unique(expinfo,keys=['TARGETID'])
    expinfo_first['NIGHT'].name = 'COADD_FIRSTNIGHT'
    expinfo_first['MJD'].name = 'COADD_FIRSTMJD'
    expinfo_first.keep_columns(['TARGETID','COADD_FIRSTNIGHT','COADD_FIRSTMJD'])
    qf = join(qf,expinfo_first,keys=['TARGETID'],join_type='left')
    del expinfo_first
    print('getting LAST info')
    expinfo_last = unique(expinfo,keys=['TARGETID'],keep='last')
    expinfo_last['NIGHT'].name = 'COADD_LASTNIGHT'
    expinfo_last['MJD'].name = 'COADD_LASTMJD'
    expinfo_last.keep_columns(['TARGETID','COADD_LASTNIGHT','COADD_LASTMJD'])
    qf = join(qf,expinfo_last,keys=['TARGETID'],join_type='left')
    del expinfo_last
    print('getting mean info')
    expinfo.sort('TARGETID')
    tids = np.unique(expinfo['TARGETID'])
    meanmjd = np.zeros(len(tids))
    ti = 0
    i = 0
    while i < len(expinfo):
        mjds = 0
        mjdw = 0
        while expinfo[i]['TARGETID'] == tids[ti]:
            mjds += expinfo[i]['MJD']*expinfo[i]['EXPTIME']
            mjdw += expinfo[i]['EXPTIME']
            i += 1
            if i == len(expinfo):
                break
        meanmjd[ti] = mjds/mjdw
        ti += 1
    meantab = Table()
    meantab['TARGETID'] = tids
    meantab['COADD_MEANMJD'] = meanmjd
    qf = join(qf,meantab,keys=['TARGETID'],join_type='left')
    return qf

    

#load the dark time healpix zcatalog, to be used for getting extra columns
logger.info('loading zcat')
zcat = Table(fitsio.read(reldir.replace('global','dvs_ro')+'/zcatalog/'+args.verspecrel+'/zpix-'+surpipe+'-dark.fits',columns=columns))
logger.info('loading exp info')
expinfo = Table(fitsio.read(reldir.replace('global','dvs_ro')+'/zcatalog/'+args.verspecrel+'/zpix-'+surpipe+'-dark.fits', 'EXP_FIBERMAP', columns=['TARGETID', 'NIGHT', 'MJD','EXPTIME']))
logger.info('loading '+args.infile+' to add columns to')
qf = fitsio.read(args.infile.replace('global','dvs_ro'))
qcols = list(qf.dtype.names)
kc = ['TARGETID']
for col in columns:
    if col not in qcols:
        kc.append(col)
if len(kc) > 1:
    zcat.keep_columns(kc)
    qf = join(qf,zcat,keys=['TARGETID'])
#get night/tile info from tiles zcat
#add_lastnight(qf,prog='dark')
qf = add_fminfo(qf,expinfo)
common.write_LSS_scratchcp(qf,args.outfile,logger=logger)
    
