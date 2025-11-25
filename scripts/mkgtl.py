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

import logging

# create logger
logname = 'LSSran'
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


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
import LSS.common_tools as common
#import LSS.mkCat_singletile.fa4lsscat as fa
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    logger.info('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--prog", help="program to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')

args = parser.parse_args()
logger.info(args)




basedir = args.basedir
version = args.version
specrel = args.verspec


if 'bright' in args.prog:
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

if '1b' in args.prog:
    pr += '1B'
    pdir += '1b'


ldirspec = basedir+'/LSS/'+specrel+'/'

dirout = ldirspec+'LSScats/'+version+'/'

if 'dark' in args.prog:
    globtype = 'LRG'
if 'bright' in args.prog:
    globtype = 'BGS'
mainp = main(globtype,args.verspec)

mt = mainp.mtld
tiles = mainp.tiles
imbits = mainp.imbits #mask bits applied to targeting
ebits = mainp.ebits #extra mask bits we think should be applied


tsnrcut = mainp.tsnrcut
dchi2 = mainp.dchi2
tnsrcol = mainp.tsnrcol        
zmin = mainp.zmin
zmax = mainp.zmax
badfib = mainp.badfib


wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
if args.survey == 'Y1':
    wd &=mt['ZDATE'] < 20220900

if args.survey == 'DA2':
    wd &=mt['ZDATE'] < 20240410


mtld = mt[wd]
#logger.info('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
logger.info('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

selt = np.isin(tiles['TILEID'],mtld['TILEID'])
ta = Table()
ta['TILEID'] = tiles[selt]['TILEID']
ta['RA'] = tiles[selt]['RA']
ta['DEC'] =tiles[selt]['DEC']



#if mkfullr or combr:
specfo = ldirspec+'datcomb_'+pdir+'_spec_zdone.fits'
logger.info('loading specf file '+specfo)
specf = Table(fitsio.read(specfo.replace('global','dvs_ro')))
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    
logger.info('loaded specf file '+specfo)
#specfc = common.cut_specdat(specf,badfib=mainp.badfib,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status)
if specrel == 'daily':
    specfc = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=False,mask_petal_nights=False,logger=logger)
else:
    specfc = common.cut_specdat(specf,badfib=mainp.badfib_td,tsnr_min=tsnrcut,tsnr_col=tnsrcol,fibstatusbits=mainp.badfib_status,remove_badfiber_spike_nz=True,mask_petal_nights=True,logger=logger)
gtl = np.unique(specfc['TILELOCID'])
filena = dirout+'/'+pdir+'_unique_good_TILELOCID.txt'
np.savetxt(filena, np.array([gtl]).astype(np.int64).T, fmt='%d')
