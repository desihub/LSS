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

#sys.path.append('../py') #this requires running from LSS/bin, *something* must allow linking without this but is not present in code yet

#from this package
#try:
import LSS.main.cattools as ct
from LSS.globals import main

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is CSCRATCH",default=os.environ['CSCRATCH'])
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='everest')
parser.add_argument("--counts_only",help="skip to just counting overlaps",default='n')
parser.add_argument("--combpix",help="if n, just skip to next stage",default='y')



args = parser.parse_args()
print(args)

basedir = args.basedir
version = args.version
specrel = args.verspec
prog = args.prog
progu = prog.upper()

combpix = True
if args.combpix == 'n':
    combpix = False

mainp = main(prog)

mt = mainp.mtld
tiles = mainp.tiles

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
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']

tiles.keep_columns(['TILEID','RA','DEC'])
#print(tiles.dtype.names)

tiles4comb = join(tiles4comb,tiles,keys=['TILEID'])

print('check that length of tiles4comb matches '+str(len(tiles4comb)))

speccols = ['TARGETID','CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE', 'NCOEFF',\
'DELTACHI2', 'PETAL_LOC','DEVICE_LOC','LOCATION','FIBER','TARGET_RA','TARGET_DEC','PMRA','PMDEC',\
'REF_EPOCH','LAMBDA_REF','FA_TARGET','FA_TYPE','OBJTYPE','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY',\
'SUBPRIORITY','OBSCONDITIONS','RELEASE','BRICKID','BRICK_OBJID','MORPHTYPE','FLUX_G','FLUX_R',\
'FLUX_Z','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','MASKBITS','REF_ID','REF_CAT',\
'GAIA_PHOT_G_MEAN_MAG','GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG','PARALLAX','BRICKNAME','EBV',\
'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',\
'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'SERSIC', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2',\
'PHOTSYS', 'PRIORITY_INIT', 'NUMOBS_INIT', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',\
'PLATE_RA', 'PLATE_DEC', 'TILEID',\
'INTEG_COADD_FLUX_B', 'MEDIAN_COADD_FLUX_B', 'MEDIAN_COADD_SNR_B', 'INTEG_COADD_FLUX_R',\
'MEDIAN_COADD_FLUX_R', 'MEDIAN_COADD_SNR_R', 'INTEG_COADD_FLUX_Z', 'MEDIAN_COADD_FLUX_Z',\
'MEDIAN_COADD_SNR_Z', 'TSNR2_ELG_B', 'TSNR2_LYA_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B',\
'TSNR2_ELG_R', 'TSNR2_LYA_R', 'TSNR2_BGS_R', 'TSNR2_QSO_R', 'TSNR2_LRG_R', 'TSNR2_ELG_Z',\
'TSNR2_LYA_Z', 'TSNR2_BGS_Z', 'TSNR2_QSO_Z', 'TSNR2_LRG_Z', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS',\
'TSNR2_QSO', 'TSNR2_LRG', 'ZWARN_MTL', 'Z_QN', 'Z_QN_CONF', 'IS_QSO_QN', 'TSNR2_GPBDARK_B',\
'TSNR2_GPBBRIGHT_B', 'TSNR2_GPBBACKUP_B', 'TSNR2_GPBDARK_R', 'TSNR2_GPBBRIGHT_R',\
'TSNR2_GPBBACKUP_R', 'TSNR2_GPBDARK_Z', 'TSNR2_GPBBRIGHT_Z', 'TSNR2_GPBBACKUP_Z',\
'TSNR2_GPBDARK', 'TSNR2_GPBBRIGHT', 'TSNR2_GPBBACKUP', 'COADD_FIBERSTATUS', 'COADD_NUMEXP',\
'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y',\
'RMS_DELTA_Y', 'MEAN_FIBER_RA', 'STD_FIBER_RA', 'MEAN_FIBER_DEC', 'STD_FIBER_DEC',\
'MEAN_PSF_TO_FIBER_SPECFLUX', 'MEAN_FIBER_X', 'MEAN_FIBER_Y']



# if len(tiles4comb) > 0:
#     ral = []
#     decl = []
#     mtlt = []
#     fal = []
#     obsl = []
#     pl = []
#     #for tile,pro in zip(mtld['TILEID'],mtld['PROGRAM']):
#     for tile in tiles4comb['TILEID']:
#         ts = str(tile).zfill(6)
#         fht = fitsio.read_header('/global/cfs/cdirs/desi/target/fiberassign/tiles/trunk/'+ts[:3]+'/fiberassign-'+ts+'.fits.gz')
#         ral.append(fht['TILERA'])
#         decl.append(fht['TILEDEC'])
#         mtlt.append(fht['MTLTIME'])
#         fal.append(fht['FA_RUN'])
#         obsl.append(fht['OBSCON'])
#     tiles4comb['RA'] = ral
#     tiles4comb['DEC'] = decl
#     tiles4comb['MTLTIME'] = mtlt
#     tiles4comb['FA_RUN'] = fal
#     tiles4comb['OBSCON'] = obsl
    


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
if args.counts_only != 'y' and combpix:
	for px in hpxs:
		print('combining target data for pixel '+str(px)+' '+str(npx)+' out of '+str(len(hpxs)))
		tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
		ct.combtiles_wdup_hp(px,tiles4comb,tarfo)
		npx += 1

if specrel == 'daily':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    if os.path.isfile(specfo):
        specf = Table.read(specfo)
        dt = specf.dtype.names
        wo = 0
        if np.isin('FIBERSTATUS',dt):
            sel = specf['COADD_FIBERSTATUS'] == 999999
            specf['COADD_FIBERSTATUS'][sel] = specf['FIBERSTATUS'][sel]
            wo = 1
        if np.isin('NUMEXP',dt):
            sel = specf['COADD_NUMEXP'] == 999999
            sel |= specf['COADD_NUMEXP'] == 16959
            specf['COADD_NUMEXP'][sel] = specf['NUMEXP'][sel]
            wo = 1
        
        specf.keep_columns(speccols)
        dtn = specf.dtype.names
        if len(dtn) != len(dt):
            wo = 1
        if wo == 1:
            specf.write(specfo,overwrite=True,format='fits')

    newspec = ct.combtile_spec(tiles4comb,specfo)
    specf = Table.read(specfo)
    if newspec:
        print('new tiles were found for spec dataso there were updates to '+specfo)
        specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','ZWARN_MTL','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
        ,'FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','COADD_FIBERSTATUS','PRIORITY'\
        ,'OBJTYPE','TILEID','INTEG_COADD_FLUX_B','MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y',\
        'RMS_DELTA_Y', 'MEAN_FIBER_RA', 'STD_FIBER_RA',\
        'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
        'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
        'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
        'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','Z_QN','Z_QN_CONF','IS_QSO_QN'])
    else:
        print('no new tiles were found for spec data, so no updates to '+specfo)
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    #tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    
    if prog == 'dark':
        tps = ['LRG','ELG','QSO','ELG_LOP','ELG_LOP']
        notqsos = ['','','','','notqso']
    if prog == 'bright':
        tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']  
        notqsos = ['',''] 
    for tp,notqso in zip(tps,notqsos):
        #first test to see if we need to update any
        print('now doing '+tp+notqso)
        print(len(tiles4comb['TILEID']))
        outf = ldirspec+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'
        update = True
        uptileloc = True
        hpxsn = hpxs
        s = 0
        if os.path.isfile(outf):
            fo = fitsio.read(outf,columns=['TARGETID','TILEID','ZWARN','ZWARN_MTL'])
            nstid = len(tiles4comb['TILEID'])
            notid = len(np.unique(fo['TILEID']))
            print('there are '+str(nstid-notid)+ ' tiles that need to be added to '+outf)
            if nstid == notid:
                update = False
                print('we will not update '+outf+' because there are no new tiles')
            else:
                
                tidc = ~np.isin(tiles4comb['TILEID'],np.unique(fo['TILEID']))
                #print('the new tileids are '+str(tiles4comb['TILEID'][tidc]))
                print(len(tiles4comb[tidc]))
                hpxsn = foot.tiles2pix(8, tiles=tiles4comb[tidc])
                          
            if os.path.isfile(outtc) and update == False:
                ftc = fitsio.read(outtc,columns=['TARGETID'])
                fc = ct.cut_specdat(fo)
                ctid = np.isin(fc['TARGETID'],ftc['TARGETID'])
                if len(ctid) == sum(ctid):
                    print('all targetids are in '+outtc+' and all tileids are in '+outf+' so '+outtc+' will not be updated')
                    uptileloc = False
                del ftc
                del fc
            del fo
        if args.counts_only != 'y' and update:
            print('updating '+outf)
            if os.path.isfile(outf):
                tarfn = fitsio.read(outf)
                theta, phi = np.radians(90-tarfn['DEC']), np.radians(tarfn['RA'])
                tpix = hp.ang2pix(8,theta,phi,nest=True)
                pin = np.isin(tpix,hpxsn)
                tarfn = tarfn[~pin] #remove the rows for the healpix that will updated
                s = 1
            
            npx =0 
            for px in hpxsn:                
                tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
                if os.path.isfile(tarfo):
                    tarf = fitsio.read(tarfo)
                    #tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
                    remcol = ['PRIORITY','Z','ZWARN','FIBER','ZWARN_MTL']
                    for col in remcol:
                        try:
                            tarf.remove_columns([col] )#we get this where relevant from spec file
                        except:
                            print('column '+col +' was not in tarwdup file')    
                    if tp == 'BGS_BRIGHT':
                        sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
                    else:
                        sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
                    if notqso == 'notqso':
                        sel &= (tarf['DESI_TARGET'] & 4) == 0
                    if s == 0:
                        tarfn = tarf[sel]
                        s = 1
                    else:
                        #tarfn = vstack([tarfn,tarf[sel]],metadata_conflicts='silent')
                        tarfn = np.hstack((tarfn,tarf[sel]))
                    print(len(tarfn),tp+notqso,npx,len(hpxsn))
                else:
                    print('file '+tarfo+' not found')
                np += 1    
            tarfn = Table(tarfn)
            tarfn['TILELOCID'] = 10000*tarfn['TILEID'] +tarfn['LOCATION']
            tarfn.write(outf,format='fits', overwrite=True)
            tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left') 
            tj.write(outfs,format='fits', overwrite=True)
        if uptileloc:
            tc = ct.count_tiles_better('dat',tp+notqso,specrel=specrel) 
            tc.write(outtc,format='fits', overwrite=True)


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

