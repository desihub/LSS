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
import desimodel.footprint as foot
from desitarget import targetmask

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='everest')


args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
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
#tiles4comb['ZDATE'] = mtd['LASTNIGHT']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']

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

ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    print('made '+ldirspec)
if not os.path.exists(ldirspec+'healpix'):
    os.mkdir(ldirspec+'healpix')
    print('made '+ldirspec+'healpix')

#outf = maindir+'datcomb_'+prog+'_spec_premtlup.fits'
#tarfo = ldirspec+'datcomb_'+prog+'_tarwdup_zdone.fits'
#ct.combtiles_wdup(tiles4comb,tarfo)
hpxs = foot.tiles2pix(8, tiles=tiles4comb)
npx = 0
for px in hpxs:
    print('combining target data for pixel '+str(px)+' '+str(npx)+' out of '+str(len(hpxs)))
    tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
    ct.combtiles_wdup_hp(px,tiles4comb,tarfo)
    npx += 1

if specrel == 'daily':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    ct.combtile_spec(tiles4comb,specfo)
    specf = Table.read(specfo)
    sel = specf['COADD_FIBERSTATUS'] == 999999
    specf['COADD_FIBERSTATUS'][sel] = specf['FIBERSTATUS'][sel]
    specf.write(specfo,overwrite=True,format='fits')
    specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','ZWARN_MTL','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
    ,'FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','COADD_FIBERSTATUS','PRIORITY'\
    ,'DELTA_X','DELTA_Y','PSF_TO_FIBER_SPECFLUX','EXPTIME','OBJTYPE','NIGHT','EXPID','MJD','TILEID','INTEG_COADD_FLUX_B',\
    'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
    'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
    'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
    'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','Z_QN','Z_QN_CONF','IS_QSO_QN'])
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    #tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    s = 0
    if prog == 'dark':
        tps = ['LRG','ELG','QSO']
    if prog == 'bright':
        tps = ['BGS_ANY','MWS_ANY']   
    for tp in tps:
        for px in hpxs:                
            tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
            if os.path.isfile(tarfo):
				tarf = Table.read(tarfo)
				tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
				remcol = ['PRIORITY','Z','ZWARN','FIBER','ZWARN_MTL']
				for col in remcol:
					try:
						tarf.remove_columns([col] )#we get this where relevant from spec file
					except:
						print('column '+col +' was not in tarwdup file')    
				sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
				if s == 0:
					tarfn = tarf[sel]
					s = 1
				else:
					tarfn = vstack([tarfn,tarf[sel]],metadata_conflicts='silent')
				print(len(tarfn),tp)
			else:
			    print('file '+tarfo+' not found')
        tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left') 
        tj.write(ldirspec+'datcomb_'+tp+'_tarspecwdup_zdone.fits',format='fits', overwrite=True)
        tc = ct.count_tiles_better('dat',tp,specrel=specrel) 
        tc.write(ldirspec+tp+'_tilelocs.dat.fits',format='fits', overwrite=True)


if specrel == 'everest':
    specf.keep_columns(['TARGETID','CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
    ,'LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
    ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
    ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
    'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
    'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG'])
    specfo = ldirspec+'datcomb_'+prog+'_zmtl_zdone.fits'
    ct.combtile_spec(tiles4comb,specfo,md='zmtl')
    fzmtl = fitsio.read(specfo)
    specf = join(specf,fzmtl,keys=['TARGETID','TILEID'])


    tarf = Table.read(tarfo)
    tarf.remove_columns(['ZWARN_MTL'])
    tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']

    #tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','FIBER'],join_type='left')
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']


#tj.write(ldirspec+'datcomb_'+prog+'_tarspecwdup_zdone.fits',format='fits', overwrite=True)
#tc = ct.count_tiles_better('dat',prog,specrel=specrel)
#tc.write(ldirspec+'Alltiles_'+prog+'_tilelocs.dat.fits',format='fits', overwrite=True)

