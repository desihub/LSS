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

from desitarget.internal import sharedmem

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


parser = argparse.ArgumentParser()
parser.add_argument("--type", help="tracer type to be selected")
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--ranmtl", help="make a random mtl file for the tile",default='n')
parser.add_argument("--rfa", help="run randoms through fiberassign",default='n')
parser.add_argument("--par", help="run different random number in parallel?",default='y')
parser.add_argument("--rann", help="number for input random file, 0 through 19",default=0,type=int)
parser.add_argument("--ran_ind",help='index for the input randoms, just 1 by default',default=1,type=int)

args = parser.parse_args()
print(args)

type = args.type
basedir = args.basedir
version = args.version
specrel = args.verspec
par = True
if args.par == 'n':
    par = False




if type[:3] == 'BGS' or type == 'bright' or type == 'MWS_ANY':
    pr = 'BRIGHT'
    pdir = 'bright'
else:
    pr = 'DARK'
    pdir = 'dark'

pd = pdir

globtype = args.type
if args.type == 'dark':
    globtype = 'LRG'
if args.type == 'bright':
    globtype == 'BGS'
mainp = main(globtype,args.verspec)

mt = mainp.mtld
tiles = mainp.tiles


wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == pdir
if specrel != 'daily':
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+pdir+'-cumulative.fits')
        wd &= np.isin(mt['TILEID'],np.unique(specf['TILEID']))
mtld = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
print('found '+str(len(mtld))+' '+pdir+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

selt = np.isin(tiles['TILEID'],mtld['TILEID'])
ta = Table()
ta['TILEID'] = tiles[selt]['TILEID']
ta['RA'] = tiles[selt]['RA']
ta['DEC'] =tiles[selt]['DEC']

maindir = basedir +'/main/LSS/'

ran_out = (args.ran_ind-1)*20+args.rann
common.printlog('ran_out is '+str(ran_out),logger)
randir = maindir+'random'+str(ran_out)

if not os.path.exists(randir):
        os.mkdir(randir)
        print('made '+randir+' random directory')
        

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)


dirout = ldirspec+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)



print(len(ta))

print(specrel)




    

dirrt = '/global/cfs/cdirs/desi/target/catalogs/dr9/2.4.0/randoms/resolve/'

common.printlog('making random target files for tiles',logger)
trad = desimodel.focalplane.get_tile_radius_deg()*1.1 #make 10% greater just in case

nd = 0
sel_tile = np.zeros(len(ta),dtype=bool)
for i in range(0,len(ta)):
	fname = dirout+str(ran_out)+'/tilenofa-'+str(ta['TILEID'][i])+'.fits'
	if os.path.isfile(fname):
		#print(fname +' already exists')
		pass
	else:
		sel_tile[i] = True
tiles = ta[sel_tile]
if len(tiles) == 0:
	sys.exit('no tiles to process ')

common.printlog('creating files for '+str(len(tiles))+' tiles',logger)
	#for i in range(0,len(tiles)):
def _create_rantile(ind):
	fname = dirout+str(ran_out)+'/tilenofa-'+str(tiles['TILEID'][ind])+'.fits'
	rtw = read_targets_in_tiles(dirrt+'randoms-'+str(args.ran_ind)+'-'+str(ii),tiles[ind])
	#print('creating '+fname)
	#tdec = tiles['DEC'][ind]
	#decmin = tdec - trad
	#decmax = tdec + trad
	#wdec = (rtall['DEC'] > decmin) & (rtall['DEC'] < decmax)
	#print(len(rt[wdec]))
	#inds = desimodel.footprint.find_points_radec(tiles['RA'][ind], tdec,rtall[wdec]['RA'], rtall[wdec]['DEC'])
	#print('got indexes')
	#rtw = rtall[wdec][inds]
	rmtl = Table(rtw)
	#print('made table for '+fname)
	del rtw
	#rmtl['TARGETID'] = np.arange(len(rmtl))
	#print(len(rmtl['TARGETID'])) #checking this column is there
	rmtl['DESI_TARGET'] = np.ones(len(rmtl),dtype=int)*2
	rmtl['NUMOBS_INIT'] = np.zeros(len(rmtl),dtype=int)
	rmtl['NUMOBS_MORE'] = np.ones(len(rmtl),dtype=int)
	rmtl['PRIORITY'] = np.ones(len(rmtl),dtype=int)*3400
	rmtl['OBSCONDITIONS'] = np.ones(len(rmtl),dtype=int)*516#tiles['OBSCONDITIONS'][i]
	rmtl['SUBPRIORITY'] = np.random.random(len(rmtl))
	#print('added columns for '+fname)
	rmtl.write(fname,format='fits', overwrite=True)
	del rmtl
	common.printlog('added columns, wrote to '+fname,logger)
	#nd += 1
	#print(str(nd),len(tiles))
inds = np.arange(len(tiles))
#for ind in inds:
#   _create_rantile(ind)
from multiprocessing import Pool
with Pool() as pool:
	res = pool.map(_create_rantile, inds)

    
