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
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--zmtl",help="if yes, only concatenate zmtl file",default='y')

args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
prog = args.prog
progu = prog.upper()

mt = Table.read('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/ops/tiles-specstatus.ecsv')
wd = mt['SURVEY'] == 'main'
#wd &= mt['EFFTIME_SPEC']/mt['GOALTIME'] > 0.85
wd &= mt['ZDONE'] == 'true'
wd &= mt['FAPRGRM'] == prog
mtd = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
print('found '+str(len(mtd))+' '+prog+' time main survey tiles with zdone true')

tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['LASTNIGHT']

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/main/LSS/'




if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    print('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    print('made '+maindir+'/LSScats')

dirout = maindir+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    print('made '+dirout)


#outf = maindir+'datcomb_'+prog+'_spec_premtlup.fits'
if args.zmtl == 'n':
    outf = maindir+'datcomb_'+prog+'_spec_zdone.fits'
    md = ''
if args.zmtl == 'y':
    outf = maindir+'datcomb_'+prog+'_zmtl_zdone.fits'   
    md = 'zmtl'
ct.combtile_spec(tiles4comb,outf,md=md)
