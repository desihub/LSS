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

#import logging
#logging.getLogger().setLevel(logging.ERROR)


#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=os.environ[scratch])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='DA02')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='newQSOtemp')
parser.add_argument("--verspec_comp",help="version for redshifts to compare to",default='guadalupe')
parser.add_argument("--doqso",help="whether or not to combine qso data",default='y')
parser.add_argument("--mkemlin",help="whether or not to make emission line files",default='y')
parser.add_argument("--dospec",help="whether or not to combine spec data",default='y')
parser.add_argument("--redospec",help="whether or not to combine spec data from beginning",default='n')
#parser.add_argument("--counts_only",help="skip to just counting overlaps",default='n')
#parser.add_argument("--combpix",help="if n, just skip to next stage",default='y')
#parser.add_argument("--redotarspec",help="re-join target and spec data even if no updates",default='n')
#parser.add_argument("--fixspecf",help="search for problem tiles and fix them in spec comb file",default='n')
#parser.add_argument("--subguad",help="replace daily data with guadalupe tiles with gauadlupe info",default='n')




args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec
prog = args.prog
progu = prog.upper()



mainp = main(prog,specver=args.verspec_comp)

if args.verspec == 'newQSOtemp':
    specdir = '/global/cfs/cdirs/desi/users/rongpu/redux/guadalupe/cumulative_new_qso_templates/'

if args.verspec == 'newQSOtemp_tagged':
    specdir = '/global/cfs/cdirs/desi/users/rongpu/redux/guadalupe/cumulative_new_qso_templates_tagged/'


mt = mainp.mtld
tiles = mainp.tiles
badfib = mainp.badfib

wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
print('number of tiles with zdone true '+str(len(mt[wd])))
wd &= mt['ARCHIVEDATE'] > 0
print('and with archivedate > 0 '+str(len(mt[wd])))
wd &= mt['FAPRGRM'] == prog
print('and in '+prog+' '+str(len(mt[wd])))

specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+args.verspec_comp+'/zcatalog/ztile-main-'+prog+'-cumulative.fits')
wd &= np.isin(mt['TILEID'],np.unique(specf['TILEID']))
mtd = mt[wd]
print('found '+str(len(mtd))+' '+prog+' time main survey tiles with zdone true for '+args.verspec_comp+' version of reduced spectra')


tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
#tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']

tiles.keep_columns(['TILEID','RA','DEC'])
#print(tiles.dtype.names)

tiles4comb = join(tiles4comb,tiles,keys=['TILEID'])

print('check that length of tiles4comb matches '+str(len(tiles4comb)))

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

if not os.path.exists(basedir +'/'+args.survey):
    os.mkdir(basedir +'/'+args.survey)
    print('made '+basedir +'/'+args.survey)

if not os.path.exists(maindir):
    os.mkdir(maindir)
    print('made '+maindir+'/logs')


if not os.path.exists(maindir):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')


ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)

print('specrel is '+specrel)

if args.doqso == 'y':
    outf = ldirspec+'QSO_catalog.fits'
    ct.combtile_qso_alt(tiles4comb,outf,coaddir=specdir)

if args.mkemlin == 'y':
    outf = ldirspec+'emlin_catalog.fits'
    ct.combtile_em_alt(tiles4comb,outf,coaddir=specdir)
    
    


if args.dospec == 'y':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'

    newspec = ct.combtile_spec_alt(tiles4comb,specfo,redo=args.redospec,prog=prog,coaddir=specdir)
    #specf = Table.read(specfo)
    if newspec:
        print('new tiles were found for spec dataso there were updates to '+specfo)
    else:
        print('no new tiles were found for spec data, so no updates to '+specfo)
