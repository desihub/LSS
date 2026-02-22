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
from desimodel.footprint import is_point_in_desi

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='everest')


args = parser.parse_args()
print(args)

basedir = args.basedir
specrel = args.verspec
prog = args.prog
progu = prog.upper()

mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wd = mt['SURVEY'] == 'main'
#wd &= mt['EFFTIME_SPEC']/mt['GOALTIME'] > 0.85
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == prog
if specrel != 'daily':
    #wd &= mt['LASTNIGHT'] < 20210801
    if specrel == 'everest':
        specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/everest/zcatalog/ztile-main-'+prog+'-cumulative.fits')
        wd &= np.isin(mt['TILEID'],np.unique(specf['TILEID']))
mtd = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
print('found '+str(len(mtd))+' '+prog+' time main survey tiles with zdone true for '+specrel+' version of reduced spectra')

tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['LASTNIGHT']

if len(tiles4comb) > 0:
    ral = []
    decl = []
    mtlt = []
    fal = []
    obsl = []
    pl = []
    #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
    for tile in tiles4comb['TILEID']:
        ts = str(tile).zfill(6)
        fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
        ral.append(fht['TILERA'])
        decl.append(fht['TILEDEC'])
        mtlt.append(fht['MTLTIME'])
        fal.append(fht['FA_RUN'])
        obsl.append(fht['OBSCON'])
    tiles4comb['RA'] = ral
    tiles4comb['DEC'] = decl
    tiles4comb['MTLTIME'] = mtlt
    tiles4comb['FA_RUN'] = fal
    tiles4comb['OBSCON'] = obsl

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/DA02/LSS/'

ldirspec = maindir+specrel+'/'

outf = ldirspec +'DA02tiles-'+progu+'.fits'   
tiles4comb.write(outf,overwrite=True,format='fits')    