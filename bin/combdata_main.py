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

sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')

args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
prog = args.prog
progu = prog.upper()

mt = Table.read('/global/cfs/cdirs/desi/spectro/redux/daily/tiles.csv')
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
tarfo = maindir+'datcomb_'+prog+'_tarwdup_zdone.fits'
ct.combtiles_wdup(ta,mdir,tarf)
specfo = maindir+'datcomb_'+prog+'_spec_zdone.fits'
ct.combtile_spec(tiles4comb,specfo)
tarf = Table.read(tarfo)
tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
remcol = ['PRIORITY','Z','ZWARN','FIBER']
for col in remcol:
	try:
		tarf.remove_columns([col] )#we get this where relevant from spec file
	except:
		print('column '+col +' was not in tarwdup file')    
specf = Table.read(specfo)
specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
,'FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','FIBERSTATUS','PRIORITY'\
,'DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','NIGHT','EXPID','MJD','TILEID','INTEG_COADD_FLUX_B',\
'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','Z_QN','Z_QN_CONF','IS_QSO_QN'])
specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
tj.write(maindir+'datcomb_'+prog+'_tarspecwdup_zdone.fits',format='fits', overwrite=True)
tc = ct.count_tiles_better('dat',prog,specrel=specrel)
tc.write(maindir+'Alltiles_'+prog+'_tilelocs.dat.fits',format='fits', overwrite=True)

