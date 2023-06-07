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
    scratch = '$CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = '$PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 


parser = argparse.ArgumentParser()
parser.add_argument("--prog", help="tracer type to be selected",choices=['dark','bright'])
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--survey", help="e.g., Y1 or DA02",default='Y1')
parser.add_argument("--verspec",help="version for redshifts",default='iron')

args = parser.parse_args()
print(args)

basedir = args.basedir
specrel = args.verspec


pdir = args.prog

pd = pdir

mainp = main('BBB') #tp argument should be irrelevant

mt = mainp.mtld
print(len(mt))
tiles = mainp.tiles

if args.survey == 'Y1':
    datemax = 20220620

wd = mt['SURVEY'] == 'main'
print(np.sum(wd))
wd &= mt['ZDONE'] == 'true'
print(np.sum(wd))
wd &= mt['FAPRGRM'] == pdir
print(np.sum(wd))
wd &= mt['LASTNIGHT'] <= datemax
print(np.sum(wd))
    
mtld = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
print('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

selt = np.isin(tiles['TILEID'],mtld['TILEID'])
ta = Table()
ta['TILEID'] = tiles[selt]['TILEID']
ta['RA'] = tiles[selt]['RA']
ta['DEC'] =tiles[selt]['DEC']

#tiles4comb = Table()
#tiles4comb['TILEID'] = mtld['TILEID']
#tiles4comb['ZDATE'] = mtld['LASTNIGHT']

sdir = basedir +'/'+args.survey+'/LSS/'





ldirspec = sdir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)



print(len(ta))

print(specrel)

specfo = ldirspec+'datcomb_'+pdir+'_spec_zdone.fits'
specf = Table.read(specfo)
sel = np.isin(specf['TILEID'],mtld['TILEID'])
specf = specf[sel]
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
	
print('loaded specf file '+specfo)
atl = np.unique(specf['TILELOCID'])
fout = ldirspec + '/alltileloc_zcat_'+pdir+'.txt'
np.savetxt(fout,atl)
print('wrote to '+fout)

specfc = common.cut_specdat(specf)
gtl = np.unique(specfc['TILELOCID'])

fout = ldirspec + '/goodtileloc_zmtl_'+pdir+'.txt'
np.savetxt(fout,gtl)
print('wrote to '+fout)