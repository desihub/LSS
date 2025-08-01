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
import LSS.common_tools as common
from LSS.globals import main

if os.environ['NERSC_HOST'] == 'cori':
    scratch = 'CSCRATCH'
elif os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    print('NERSC_HOST is not cori or permutter but is '+os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding') 

parser = argparse.ArgumentParser()
parser.add_argument("--basedir", help="base directory for output, default is SCRATCH",default=scratch)
parser.add_argument("--version", help="catalog version; use 'test' unless you know what you are doing!",default='test')
parser.add_argument("--survey", help="e.g., main (for all), DA02, any future DA",default='main')
parser.add_argument("--prog", help="dark or bright is supported",default='dark')
parser.add_argument("--verspec",help="version for redshifts",default='daily')
parser.add_argument("--check_date_only",help="whether or not to stop after maximum night is found",default='n')
parser.add_argument("--make_tile_file",help="whether or not to make a tile file",default='n')
parser.add_argument("--doqso",help="whether or not to combine qso data",default='n')
parser.add_argument("--redoqso",help="whether or not to combine qso data, starting over",default='n')
parser.add_argument("--mkemlin",help="whether or not to make emission line files",default='n')
parser.add_argument("--dospec",help="whether or not to combine spec data",default='y')
parser.add_argument("--dotarspec",help="whether or not to combine spec and tar data per type, for non-daily data",default='y')
parser.add_argument("--redospec",help="whether or not to combine spec data from beginning",default='n')
parser.add_argument("--redo_zmtl",help="whether or not to remake zmtl file",default='n')
parser.add_argument("--counts_only",help="skip to just counting overlaps",default='n')
parser.add_argument("--combpix",help="if n, just skip to next stage",default='y')
parser.add_argument("--get_petalsky",help="if y, combine info across tiles to get dispersion in sky fibers",default='n')
parser.add_argument("--comb_petalqa",help="if y, combine petal qa info across tiles ",default='n')
parser.add_argument("--par",help="if y, using multiprocessing ",default='n')

parser.add_argument("--redotardup",help="re-run the potential assignments concatenation",default='n')
parser.add_argument("--redotarspec",help="re-join target and spec data even if no updates",default='n')
parser.add_argument("--fixspecf",help="search for problem tiles and fix them in spec comb file",default='n')
parser.add_argument("--subguad",help="replace daily data with guadalupe tiles with gauadlupe info",default='n')
parser.add_argument("--tracer", help="tracer type",default='all')
parser.add_argument("--notqso", help="tracer type",default='')



args = parser.parse_args()
print(args)

import logging
# create logger
logname = 'comb_inputs'
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


basedir = args.basedir
version = args.version
specrel = args.verspec
prog = args.prog
progu = prog.upper()
redoqso = False
if args.redoqso == 'y':
    redoqso = True

combpix = True
if args.combpix == 'n':
    combpix = False

redotarspec = False
if args.redotarspec == 'y':
    redotarspec = True

mainp = main(prog,specver=specrel)

mt = mainp.mtld
tiles = mainp.tiles
badfib = mainp.badfib

wd = mt['SURVEY'] == 'main'
#wd &= mt['EFFTIME_SPEC']/mt['GOALTIME'] > 0.85
wd &= mt['ZDONE'] == 'true'
logger.info('number of tiles with zdone true '+str(len(mt[wd])))
wd &= mt['ARCHIVEDATE'] > 0
logger.info('and with archivedate > 0 '+str(len(mt[wd])))
wd &= mt['FAPRGRM'] == prog
logger.info('and in '+prog+' '+str(len(mt[wd])))
if specrel != 'daily':
    #wd &= mt['LASTNIGHT'] < 20210801
    #if specrel == 'everest':
    
    specrell = specrel.split('-')
    coaddir = '/global/cfs/cdirs/desi/spectro/redux/'+specrell[0]+'/tiles/cumulative/'
    specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/'+specrell[0]+'/zcatalog/'+specrell[1]+'/ztile-main-'+prog+'-cumulative.fits')
    wd &= np.isin(mt['TILEID'],np.unique(specf['TILEID']))
else:
    coaddir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/'

if args.survey == 'Y1':
    wd &= mt['ZDATE'] < 20220900
if args.survey == 'DA2':
    wd &= mt['ZDATE'] < 20240410

mtd = mt[wd]
#print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
logger.info('found '+str(len(mtd))+' '+prog+' time '+args.survey+' survey tiles with zdone true for '+specrel+' version of reduced spectra')


tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
if args.verspec == 'daily':
    tiles4comb['THRUDATE'] = mtd['ZDATE']#mtd['LASTNIGHT']
else:
    tiles4comb['THRUDATE'] = mtd['LASTNIGHT']#mtd['LASTNIGHT']


logger.info('The last night of data that will be processed is for '+args.prog+' is '+str(np.max(tiles4comb['THRUDATE'] )))
logger.info('Is that what was expected based on MTL updates?')
if args.check_date_only == 'y':
    sys.exit()

tiles.keep_columns(['TILEID','RA','DEC'])
#print(tiles.dtype.names)

tiles4comb = join(tiles4comb,tiles,keys=['TILEID'])

logger.info('check that length of tiles4comb matches '+str(len(tiles4comb)))

#share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir +'/'+args.survey+'/LSS/'

if not os.path.exists(maindir+'/logs'):
    os.mkdir(maindir+'/logs')
    logger.info('made '+maindir+'/logs')

if not os.path.exists(maindir+'/LSScats'):
    os.mkdir(maindir+'/LSScats')
    logger.info('made '+maindir+'/LSScats')

dirout = maindir+'LSScats/'+version+'/'
if not os.path.exists(dirout):
    os.mkdir(dirout)
    logger.info('made '+dirout)

dailydir = maindir+'daily/'
ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    logger.info('made '+ldirspec)
if not os.path.exists(ldirspec+'healpix'):
    os.mkdir(ldirspec+'healpix')
    logger.info('made '+ldirspec+'healpix')

if args.make_tile_file == 'y':
    tiles4comb.write(maindir+'tiles-'+prog.upper()+'.fits',overwrite=True,format='fits')

logger.info('specrel is '+specrel)
if specrel == 'daily':
    #specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    specfo = basedir+'main/LSS/daily/datcomb_'+prog.replace('1b','')+'_spec_zdone.fits'
    #if not os.path.isfile(specfo) and args.subguad != 'y':
    if os.path.isfile(specfo)  and args.redospec == 'n':# and args.subguad == 'n':
        specf = fitsio.read(specfo)    
    #else:
    #    specf = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/datcomb_'+prog+'_spec_zdone.fits')

    speccols = list(specf.dtype.names)
    logger.info(str(speccols))
    #spec_cols_4tar = ['TARGETID','Z','ZERR','ZWARN','ZWARN_MTL','SPECTYPE','DELTACHI2'\
    #,'LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','TILELOCID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
    #,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
    #spec_cols_4tar = ['TARGETID','ZWARN','ZWARN_MTL','LOCATION','FIBER','TILEID','TILELOCID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
    spec_cols_4tar = ['TARGETID','ZWARN','ZWARN_MTL','LOCATION','FIBER','TILEID','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY']
    logger.info(str(spec_cols_4tar))
    if args.subguad == 'y':
        dz = Table(fitsio.read(specfo))
        dz.keep_columns(speccols)
        specf = Table(specf)
        gr = np.isin(dz['TILEID'],specf['TILEID'])
        dng = dz[~gr]
        print('length before/after removing guad tiles from daily')
        print(len(dz),len(dng))
        ng = vstack([specf,dng])
        print('length after adding guad data back in '+str(len(ng)))
        print(ng.dtype.names)
        ng.write(specfo,format='fits',overwrite=True)
    del specf

regl = ['N','S']

# speccols = ['TARGETID','CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE', 'NCOEFF',\
# 'DELTACHI2', 'PETAL_LOC','DEVICE_LOC','LOCATION','FIBER','TARGET_RA','TARGET_DEC','PMRA','PMDEC',\
# 'REF_EPOCH','LAMBDA_REF','FA_TARGET','FA_TYPE','OBJTYPE','FIBERASSIGN_X','FIBERASSIGN_Y','PRIORITY',\
# 'SUBPRIORITY','OBSCONDITIONS','RELEASE','BRICKID','BRICK_OBJID','MORPHTYPE','FLUX_G','FLUX_R',\
# 'FLUX_Z','FLUX_IVAR_G','FLUX_IVAR_R','FLUX_IVAR_Z','MASKBITS','REF_ID','REF_CAT',\
# 'GAIA_PHOT_G_MEAN_MAG','GAIA_PHOT_BP_MEAN_MAG','GAIA_PHOT_RP_MEAN_MAG','PARALLAX','BRICKNAME','EBV',\
# 'FLUX_W1', 'FLUX_W2', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'FIBERFLUX_G', 'FIBERFLUX_R', 'FIBERFLUX_Z',\
# 'FIBERTOTFLUX_G', 'FIBERTOTFLUX_R', 'FIBERTOTFLUX_Z', 'SERSIC', 'SHAPE_R', 'SHAPE_E1', 'SHAPE_E2',\
# 'PHOTSYS', 'PRIORITY_INIT', 'NUMOBS_INIT', 'DESI_TARGET', 'BGS_TARGET', 'MWS_TARGET', 'SCND_TARGET',\
# 'PLATE_RA', 'PLATE_DEC', 'TILEID',\
# 'INTEG_COADD_FLUX_B', 'MEDIAN_COADD_FLUX_B', 'MEDIAN_COADD_SNR_B', 'INTEG_COADD_FLUX_R',\
# 'MEDIAN_COADD_FLUX_R', 'MEDIAN_COADD_SNR_R', 'INTEG_COADD_FLUX_Z', 'MEDIAN_COADD_FLUX_Z',\
# 'MEDIAN_COADD_SNR_Z', 'TSNR2_ELG_B', 'TSNR2_LYA_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B',\
# 'TSNR2_ELG_R', 'TSNR2_LYA_R', 'TSNR2_BGS_R', 'TSNR2_QSO_R', 'TSNR2_LRG_R', 'TSNR2_ELG_Z',\
# 'TSNR2_LYA_Z', 'TSNR2_BGS_Z', 'TSNR2_QSO_Z', 'TSNR2_LRG_Z', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS',\
# 'TSNR2_QSO', 'TSNR2_LRG', 'ZWARN_MTL', 'Z_QN', 'Z_QN_CONF', 'IS_QSO_QN', 'TSNR2_GPBDARK_B',\
# 'TSNR2_GPBBRIGHT_B', 'TSNR2_GPBBACKUP_B', 'TSNR2_GPBDARK_R', 'TSNR2_GPBBRIGHT_R',\
# 'TSNR2_GPBBACKUP_R', 'TSNR2_GPBDARK_Z', 'TSNR2_GPBBRIGHT_Z', 'TSNR2_GPBBACKUP_Z',\
# 'TSNR2_GPBDARK', 'TSNR2_GPBBRIGHT', 'TSNR2_GPBBACKUP', 'COADD_FIBERSTATUS', 'COADD_NUMEXP',\
# 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'COADD_NUMTILE', 'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y',\
# 'RMS_DELTA_Y', 'MEAN_FIBER_RA', 'STD_FIBER_RA', 'MEAN_FIBER_DEC', 'STD_FIBER_DEC',\
# 'MEAN_PSF_TO_FIBER_SPECFLUX', 'MEAN_FIBER_X', 'MEAN_FIBER_Y']



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
    



#outf = maindir+'datcomb_'+prog+'_spec_premtlup.fits'
#tarfo = ldirspec+'datcomb_'+prog+'_tarwdup_zdone.fits'
#ct.combtiles_wdup(tiles4comb,tarfo)
if specrel == 'daily' and args.survey == 'main':
    hpxs = foot.tiles2pix(8, tiles=tiles4comb)
    npx = 0
    if args.counts_only != 'y' and combpix:
        processed_tiles_file = ldirspec+'processed_tiles_'+prog+'.fits'
        if os.path.isfile(processed_tiles_file):
            tiles_proc = Table.read(processed_tiles_file)
            tidsp = np.isin(tiles4comb['TILEID'],tiles_proc['TILEID'])
            tidsnp = ~tidsp
            tiles4hp = tiles4comb[tidsnp]
        else :
            common.printlog('didnt load processed tiles file '+processed_tiles_file,logger)
            tiles4hp = tiles4comb
    
        common.printlog('will combine pixels for '+str(len(tiles4hp))+' new tiles',logger)
        if len(tiles4hp) > 0:
            for px in hpxs:
                common.printlog('combining target data for pixel '+str(px)+' '+str(npx)+' out of '+str(len(hpxs)),logger)
                tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
                ct.combtiles_wdup_hp(px,tiles4hp,tarfo)
                npx += 1
            tiles4comb.write(processed_tiles_file,format='fits',overwrite=True)

if specrel == 'daily' and args.survey == 'DA2':
    tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'
    tids = tiles4comb['TILEID']#list()
    def _tab2list(tid):
        sel = tiles4comb['TILEID'] == tid
        logger.info('at TILEID '+str(tid))
        tl_tab = tiles4comb[sel]
        tab = ct.get_tiletab(tl_tab)
        return tab

    if os.path.isfile(tarfo) == False or args.redotardup == 'y':
        logger.info('creating '+tarfo)
        tile_list = []
        tids_todo = tids
    else:
        logger.info('not remaking '+tarfo)
        tile_list = [fitsio.read(tarfo)]
        tids_c = np.unique(tile_list[0]['TILEID'])
        tids_todo_inds = ~np.isin(tids,tids_c)
        tids_todo = tids[tids_todo_inds]
    if len(tids_todo) > 0:    
        if args.par == 'y':
            
            #test of what goes in parallel
            #tid = tiles4comb['TILEID'][0]
            #sel = tiles4comb['TILEID'] == tid       
            #tl_tab = tiles4comb[sel]
            #print(tl_tab)
            #tab = ct.get_tiletab(tl_tab)
            #logger.info(str(tab.dtype.names))
            
            #from multiprocessing import Process, Manager
            #manager = Manager()
            #tile_list = manager.list()
            #def _tab2list(tlist,tid):
            #    tlist.append(tab)
            #inds = np.arange(len(tiles4comb))
            #tiles4comb['TILEID'][i])) for i in inds
            #job = [Process(target=_tab2list, args=(tile_list, tiles4comb['TILEID'][i])) for i in inds]
            #_ = [p.start() for p in job]
            #_ = [p.join() for p in job]
            from concurrent.futures import ProcessPoolExecutor
            
            with ProcessPoolExecutor() as executor:
                for tab in executor.map(_tab2list, list(tids_todo)):
                    tile_list.append(np.array(tab))
            logger.info('tiles in list of length '+str(len(tile_list)))
            logger.info('concatenating')
            logger.info(str(tile_list[0].dtype.names))
            tarsn = np.concatenate(tile_list)#vstack(tile_list,metadata_conflicts='silent')
            logger.info(str(tarsn.dtype.names))
            del tile_list
            logger.info('doing TARGETID sort')
            tarsn = Table(tarsn)
            tarsn.sort('TARGETID')
            
            logger.info('sort done')
            common.write_LSS(tarsn,tarfo)
        
        else:
            ct.combtiles_wdup(tiles4comb,fout=tarfo)
    else:
        logger.info('no new tiles to combine')
        
        
            
if  args.doqso == 'y':
    outf = ldirspec+'QSO_catalog.fits'
    if specrel == 'daily':# and args.survey == 'main':
        ct.combtile_qso(tiles4comb,outf,restart=redoqso)
    else:
        ct.combtile_qso_alt(tiles4comb,outf,coaddir=coaddir)

if  args.mkemlin == 'y':
    outf = ldirspec+'emlin_catalog.fits'
    if specrel == 'daily' and args.survey == 'main':
        outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/emtiles/'
        guadtiles = fitsio.read('/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/datcomb_'+prog+'_spec_zdone.fits',columns=['TILEID'])
        guadtiles = np.unique(guadtiles['TILEID'])
        gtids = np.isin(tiles4comb['TILEID'],guadtiles)
        tiles4em = tiles4comb[~gtids]
        ndone = 0
        for tile,zdate,tdate in zip(tiles4em['TILEID'],tiles4em['ZDATE'],tiles4em['THRUDATE']):
            outft = outdir+'emline-'+str(tile)+'.fits'
            if not os.path.isfile(outf):
                tdate = str(tdate)
                ct.combEMdata_daily_old(tile,zdate,tdate,outf=outft)
                print('wrote '+outf)
                ndone += 1
                print('completed '+str(ndone)+' tiles')
        ct.combtile_em(tiles4comb,outf)
    elif specrel != 'daily':
        if args.par == 'y':
            tl = []
            if os.path.isfile(outf):
                specd = fitsio.read(outf)
                tl.append(specd)
                tdone = np.unique(specd['TILEID'])
                tmask = ~np.isin(tiles4comb['TILEID'],tdone)
            else:
                tmask = np.ones(len(tiles4comb)).astype('bool')
            logger.info('there are '+str(np.sum(tmask))+ ' tiles to add to emline data')
            if np.sum(tmask) == 0:
                newem = False
            else:
                newem = True
                tiles_2comb = tiles4comb[tmask]
                def _get_tile(ind):
                
                    trow = tiles_2comb[ind]
                    tile,zdate,tdate = trow['TILEID'],trow['ZDATE'],trow['THRUDATE']
                    logger.info('combining emline data for TILEID '+str(tile))
                    #tspec = ct.combEMdata_daily(str(tile),str(zdate),str(tdate))
                    tspec = ct.combEMdata_rel(str(tile),str(tdate),coaddir)
                    if tspec is not None:
                        #tspec['TILEID'] = tile
                        tspec = np.array(tspec)
                        logger.info('tile '+str(tile)+' '+str(len(tspec.dtype.names))+' '+str(len(tspec)))
                        return tspec
                        #new = np.empty(len(tspec),dtype=specd.dtype)
                        #cols = specd.dtype.names
                        #for colname in cols:
                        #    new[colname][...] = tspec[colname][...]
                    else:
                        logger.info('tile '+str(tile)+' failed for emline')
                        return None
                        #new = None
                    #return new
                inds = np.arange(len(tiles_2comb))
                from concurrent.futures import ProcessPoolExecutor
            
                with ProcessPoolExecutor() as executor:
                    for specd in executor.map(_get_tile, inds):
                        if specd is not None:
                            tl.append(np.array(specd))
            #nms = list(tl[0].dtype.names)
            #for t in tl:
            #    nmsi = list(t.dtype.names)
            #    if nmsi != nms:
            #        logger.info(str(t[0]['TILEID'])+' has mismatched name list')
            #    try:
            #        temp = np.hstack([tl[0],t])
            #    except:
            #        logger.info(str(t[0]['TILEID'])+' and '+str(tl[0][0]['TILEID'])+' failed hstack')
            #        logger.info(str(tl[0]))
            #        logger.info(str(t))
            #        break
                
            specd = np.hstack(tl)
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]
            common.write_LSS(specd,outf)
            del specd
            del tl

        else:
            ct.combtile_em_alt(tiles4comb,outf,prog='dark',coaddir=coaddir)
    else:
        if args.par == 'y':
            tl = []
            if os.path.isfile(outf):
                specd = fitsio.read(outf)
                tl.append(specd)
                tdone = np.unique(specd['TILEID'])
                tmask = ~np.isin(tiles4comb['TILEID'],tdone)
            else:
                tmask = np.ones(len(tiles4comb)).astype('bool')

            logger.info('there are '+str(np.sum(tmask))+ ' tiles to add to emline data')
            if np.sum(tmask) == 0:
                newem = False
            else:
                newem = True
                tiles_2comb = tiles4comb[tmask]
                def _get_tile(ind):
                
                    trow = tiles_2comb[ind]
                    tile,zdate,tdate = trow['TILEID'],trow['ZDATE'],trow['THRUDATE']
                    logger.info('combining emline data for TILEID '+str(tile))
                    tspec = ct.combEMdata_daily(str(tile),str(zdate),str(tdate))
                    if tspec:
                        tspec['TILEID'] = tile
                        tspec = np.array(tspec)
                        new = np.empty(len(tspec),dtype=specd.dtype)
                        cols = specd.dtype.names
                        for colname in cols:
                            new[colname][...] = tspec[colname][...]
                    else:
                        new = None
                    return new
                inds = np.arange(len(tiles_2comb))
                from concurrent.futures import ProcessPoolExecutor
            
                with ProcessPoolExecutor() as executor:
                    for specd in executor.map(_get_tile, inds):
                        if specd is not None:
                            tl.append(np.array(specd))
            specd = np.hstack(tl)
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]
            common.write_LSS(specd,outf)
            del specd
            del tl
        
        else:
            ct.combtile_em_daily(tiles4comb,outf)

if args.survey == 'Y1' and args.counts_only == 'y':    
    if prog == 'dark':
        tps = ['LRG','ELG','QSO','ELG_LOP','ELG_LOP']
        notqsos = ['','','','','notqso']
    if prog == 'bright':
        tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']  
        notqsos = ['',''] 
    for tp,notqso in zip(tps,notqsos):

        tc = ct.count_tiles_better('dat',tp+notqso,specrel=specrel,survey=args.survey,badfib=badfib) 
        outtc = ldirspec+tp+notqso+'_tilelocs.dat.fits'
        tc.write(outtc,format='fits', overwrite=True)


if prog == 'dark':
    if args.tracer == 'all':
        #tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
        #notqsos = ['','','notqso','','']
        tps = ['QSO','LRG','ELG'] #only do base types because other are subset that can be cut later and this saves i/o
        notqsos = ['','','']
    else:
        tps = [args.tracer]
        notqsos = [args.notqso]    
if prog == 'dark1b':
    if args.tracer == 'all':
        #tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
        #notqsos = ['','','notqso','','']
        tps = ['LGE','QSO','LRG','ELG'] #only do base types because other are subset that can be cut later and this saves i/o
        notqsos = ['','','','']
    else:
        tps = [args.tracer]
        notqsos = [args.notqso]    

if 'bright' in prog:
    if args.tracer == 'all':
        #tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']  
        #notqsos = ['',''] 
        tps = ['BGS_ANY']#,'MWS_ANY']  
        notqsos = [''] 

    else:
        tps = [args.tracer]
        notqsos = [args.notqso]    


if specrel == 'daily' and args.dospec == 'y' and args.survey != 'main':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    dotarspec = False
    if args.redotarspec == 'y':
        dotarspec = True
    if os.path.isfile(specfo) and args.redospec == 'n':
        specf = Table.read(specfo)
        if list(specf.dtype.names) != speccols:
            logger.info('writing original file back out with subselection of columns')
            specf.keep_columns(speccols)
            common.write_LSS(specf,specfo)
        del specf
    if args.par == 'y':
        tl = []
        if os.path.isfile(specfo) :
            specd = fitsio.read(specfo)
            tl.append(specd)
            tdone = np.unique(specd['TILEID'])
            tmask = ~np.isin(tiles4comb['TILEID'],tdone)

        else:
            tmask = np.ones(len(tiles4comb)).astype('bool')
        logger.info('there are '+str(np.sum(tmask))+ ' tiles to add to spec data')
        if np.sum(tmask) == 0:
            newspec = False
        else:
            newspec = True
            tiles_2comb = tiles4comb[tmask]
            def _get_tile(ind):
                
                trow = tiles_2comb[ind]
                tile,zdate,tdate = trow['TILEID'],trow['ZDATE'],trow['THRUDATE']
                logger.info('combining spec data for TILEID '+str(tile))
                tspec = ct.combspecdata(str(tile),str(zdate),str(tdate))
                if tspec:
                    tspec['TILEID'] = tile
                    tspec = np.array(tspec)
                    new = np.empty(len(tspec),dtype=specd.dtype)
                    cols = specd.dtype.names
                    for colname in cols:
                        new[colname][...] = tspec[colname][...]
                else:
                    new = None
                return new
            inds = np.arange(len(tiles_2comb))
            from concurrent.futures import ProcessPoolExecutor
            
            with ProcessPoolExecutor() as executor:
                for specd in executor.map(_get_tile, inds):
                    if specd is not None:
                        tl.append(np.array(specd))
        if newspec:
            specd = np.hstack(tl)
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]
            common.write_LSS(specd,specfo)
            del specd
        del tl
    else:
        newspec = ct.combtile_spec(tiles4comb,specfo,redo=args.redospec,prog=prog,par=args.par)
    specf = Table.read(specfo)
    if newspec:
        logger.info('new tiles were found for spec dataso there were updates to '+specfo)
        dotarspec = True
    else:
        logger.info('no new tiles were found for spec data, so no updates to '+specfo)
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    specf.keep_columns(spec_cols_4tar)
    #tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    
    for tp,notqso in zip(tps,notqsos):
        logger.info('now doing '+tp+notqso)
        logger.info(len(tiles4comb['TILEID']))
        #outf = ldirspec+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        if os.path.isfile(outfs) == False or dotarspec:
            #update = True
            #dotarspec = True
    
            tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'
            tarf = fitsio.read(tarfo.replace('global','dvs_ro'))#,columns=cols)
            logger.info('loaded tarspecwdup file')
            #tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            if tp == 'BGS_BRIGHT':
                sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
            else:
                sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
            if notqso == 'notqso':
                sel &= (tarf['DESI_TARGET'] & 4) == 0
            tarf = Table(tarf[sel])
            logger.info('cut to target type')
            
            tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            logger.info('added TILELOCID, about to do joins')
            #tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    
            #seems to run out of memory on join
            tjl = []
            logger.info(tarf.dtype.names)
            selreg = tarf['DEC'] > 0
            logger.info(len(tarf[selreg]))
            remcol = ['LOCATION','TILEID','FIBER','PRIORITY']
            for col in remcol:
                try:
                    specf.remove_columns([col])
                except:
                    logger.info('column '+col +' was not in stacked spec table') 
            tjl.append(join(tarf[selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
            tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
            logger.info('1st join done')
            tjl.append(join(tarf[~selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
            tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
            logger.info('2nd join done')
            del tarf
            tj = vstack(tjl)
            logger.info('stacked now writing out')
            common.write_LSS_scratchcp(tj,outfs,logger=logger)
            logger.info('joined to spec data and wrote out to '+outfs)
        else:
            logger.info(outfs +' exists already, not making again')


        
if specrel == 'daily' and args.dospec == 'y' and args.survey == 'main':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    if os.path.isfile(specfo) and args.redospec == 'n':
        specf = Table.read(specfo)
        if args.fixspecf == 'y':
            ii = 0
            for tid,adate,zdate in zip(tiles4comb['TILEID'],tiles4comb['ZDATE'],tiles4comb['THRUDATE']):
                ii += 1
                if int(zdate) >  20210730:
                    td = ct.combspecdata(tid,str(adate),str(zdate))
                    kp = (td['TARGETID'] > 0)
                    td = td[kp]
                    sel = specf['TILEID'] == tid
                    fst = specf[sel]
                    if np.array_equal(fst['ZWARN_MTL'],td['ZWARN_MTL']):
                        print('tile '+str(tid)+' passed')
                    else:
                        print('tile '+str(tid)+' is mismatched')
                        specf = specf[~sel]
                        specf = vstack([specf,td])
                        print(ii,len(tiles4comb))
            specf.write(specfo,format='fits',overwrite=True)
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
            common.printlog('rewriting spec file',logger)
            specf.write(specfo,overwrite=True,format='fits')

    newspec = ct.combtile_spec(tiles4comb,specfo,redo=args.redospec,prog=prog)
    
    if newspec:
        common.printlog('new tiles were found for spec dataso there were updates to '+specfo,logger)
    else:
        common.printlog('no new tiles were found for spec data, so no updates to '+specfo,logger)
    specf = fitsio.read(specfo,columns=spec_cols_4tar)
    common.printlog('spec file '+specfo+' has '+str(len(specf))+' rows',logger)
    if '1b' in prog:
        common.printlog('adding '+specfo.replace('1b','') +' to spec info',logger)
        specfnb = fitsio.read(specfo.replace('1b',''),columns=spec_cols_4tar)
        specf = np.concatenate([specf,specfnb])
        del specfnb
    specf = Table(specf)
    
#     specf.keep_columns(['CHI2','COEFF','Z','ZERR','ZWARN','ZWARN_MTL','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
#     ,'FIBERASSIGN_X','FIBERASSIGN_Y','TARGETID','LOCATION','FIBER','COADD_FIBERSTATUS'\
#     ,'OBJTYPE','TILEID','INTEG_COADD_FLUX_B','MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y',\
#     'RMS_DELTA_Y', 'MEAN_FIBER_RA', 'STD_FIBER_RA',\
#     'MEDIAN_COADD_FLUX_B','MEDIAN_COADD_SNR_B','INTEG_COADD_FLUX_R','MEDIAN_COADD_FLUX_R','MEDIAN_COADD_SNR_R','INTEG_COADD_FLUX_Z',\
#     'MEDIAN_COADD_FLUX_Z','MEDIAN_COADD_SNR_Z','TSNR2_ELG_B','TSNR2_LYA_B','TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
#     'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
#     'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','Z_QN','Z_QN_CONF','IS_QSO_QN'])

    
    specf['TILELOCID'] = 10000*specf['TILEID'] +specf['LOCATION']
    #specf.keep_columns(spec_cols_4tar)
    #tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    
#     if prog == 'dark':
#         if args.tracer == 'all':
#             #tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
#             #notqsos = ['','','notqso','','']
#             tps = ['QSO','LRG','ELG'] #only do base types because other are subset that can be cut later and this saves i/o
#             notqsos = ['','','']
# 
#         else:
#             tps = [args.tracer]
#             notqsos = [args.notqso]    
#     if prog == 'bright':
#         if args.tracer == 'all':
#             #tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']  
#             #notqsos = ['',''] 
#             tps = ['BGS_ANY']#,'MWS_ANY']  
#             notqsos = [''] 
# 
#         else:
#             tps = [args.tracer]
#             notqsos = [args.notqso]    

    for tp,notqso in zip(tps,notqsos):
        #first test to see if we need to update any
        common.printlog('now doing '+tp+notqso,logger)
        print(len(tiles4comb['TILEID']))
        outf = ldirspec+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'
        update = True
        uptileloc = False
        dotarspec = True
        hpxsn = hpxs
        s = 0
        if os.path.isfile(outf):
            fo = fitsio.read(outf,columns=['TARGETID','TILEID'])
            common.printlog(outf +' has '+str(len(fo))+ ' rows',logger)
            nstid = len(tiles4comb['TILEID'])
            notid = len(np.unique(fo['TILEID']))
            test_tid = np.isin(tiles4comb['TILEID'],np.unique(fo['TILEID']))
            common.printlog('there are '+str(nstid-np.sum(test_tid))+ ' tiles that need to be added to '+outf,logger)
            if nstid == np.sum(test_tid):
                update = False
                print('we will not update '+outf+' because there are no new tiles')
            else:
                
                tidc = ~np.isin(tiles4comb['TILEID'],np.unique(fo['TILEID']))
                #print('the new tileids are '+str(tiles4comb['TILEID'][tidc]))
                print(len(tiles4comb[tidc]))
                hpxsn = foot.tiles2pix(8, tiles=tiles4comb[tidc])
            del fo
        if os.path.isfile(outfs):
            fo = fitsio.read(outfs,columns=['TARGETID','TILEID','ZWARN','ZWARN_MTL'])
            stids = np.unique(fo['TILEID'])
            if len(stids) == notid:      
                dotarspec = False   
            if os.path.isfile(outtc) and update == False and redotarspec == False and dotarspec == False:
                ftc = fitsio.read(outtc,columns=['TARGETID'])
                fc = ct.cut_specdat(fo)
                ctid = np.isin(fc['TARGETID'],ftc['TARGETID'])
                if len(ctid) == sum(ctid):
                    common.printlog('all targetids are in '+outtc+' and all tileids are in '+outf+' so '+outtc+' will not be updated',logger)
                    uptileloc = False
                del ftc
                del fc
            del fo
        if args.counts_only != 'y' and update:
            common.printlog('updating '+outf,logger)
            cols = None
            if os.path.isfile(outf):
                tarfn = fitsio.read(outf)
                cols = tarfn.dtype.names
                if np.isin('TILELOCID',tarfn.dtype.names):
                    print('reloading '+outf+' without reading TILELOCID column')
                    #sel = cols != 'TILELOCID'
                    #cols = cols[sel]
                    cols = []
                    for col in tarfn.dtype.names:
                        if col != 'TILELOCID':
                            cols.append(col)
                    tarfn = fitsio.read(outf,columns=cols)
                    print(tarfn.dtype.names)
                theta, phi = np.radians(90-tarfn['DEC']), np.radians(tarfn['RA'])
                tpix = hp.ang2pix(8,theta,phi,nest=True)
                pin = np.isin(tpix,hpxsn)
                tarfn = tarfn[~pin] #remove the rows for the healpix that will updated
                s = 1
            
            npx =0 
            for px in hpxsn:                
                tarfo = ldirspec+'healpix/datcomb_'+prog+'_'+str(px)+'_tarwdup_zdone.fits'
                if os.path.isfile(tarfo):
                    if cols is not None:
                        tarf = fitsio.read(tarfo,columns=cols)
                    else:
                        tarf = fitsio.read(tarfo)
                    #tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
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
                        common.printlog(tarfn.dtype.names,logger)
                        common.printlog(tarf.dtype.names,logger)
                        tarfn = np.hstack((tarfn,tarf[sel]))
                    common.printlog(str(len(tarfn))+','+tp+notqso+','+str(npx)+','+str(len(hpxsn)),logger)
                else:
                    common.printlog('file '+tarfo+' not found',logger)
                npx += 1    
            tarfn = Table(tarfn)           
            remcol = ['Z','ZWARN','FIBER','ZWARN_MTL']
            for col in remcol:
                try:
                    tarfn.remove_columns([col] )#we get this where relevant from spec file
                except:
                    print('column '+col +' was not in stacked tarwdup table')    

            common.write_LSS_scratchcp(tarfn,outf,logger=logger)
            #tarfn.write(outf,format='fits', overwrite=True)
            common.printlog('wrote out '+outf+' '+str(len(tarfn))+' rows',logger)
            
            #try:
            #    specf.remove_columns(['PRIORITY'])
            #except:
            #    print('column PRIORITY was not in spec table')  
            
            #join to spec info; now only do so after updating 1b tiles
            if '1b' in prog:
                tarfn['TILELOCID'] = 10000*tarfn['TILEID'] +tarfn['LOCATION']
                print('added TILELOCID, about to do joins')
                #tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
    
                #seems to run out of memory on join
                tjl = []
                print(tarfn.dtype.names)
                selreg = tarfn['DEC'] > 0
                print(len(tarfn[selreg]))
                remcol = ['LOCATION','TILEID']
                for col in remcol:
                    try:
                        specf.remove_columns([col])
                    except:
                        print('column '+col +' was not in stacked spec table') 
                if np.sum(selreg) > 0:
                    tjl.append(join(tarfn[selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
                    tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
                    common.printlog('1st join done',logger)
                if np.sum(~selreg) > 0:
                    tjl.append(join(tarfn[~selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
                    tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
                    common.printlog('2nd join done',logger)
                del tarfn
                if len(tjl) > 1:
                    tj = vstack(tjl)
                else:
                    tj = tjl[0]
                del tjl
                common.printlog('stacked now writing out',logger)
                #for reg in regl:                
                #    sel = tarfn['PHOTSYS'] == reg
                #    tjr = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left') 
                #tj.write(outfs,format='fits', overwrite=True)
                common.write_LSS_scratchcp(tj,outfs,logger=logger)
                common.printlog('joined to spec data and wrote out to '+outfs,logger)
        elif redotarspec or dotarspec:
            common.printlog('joining spec info to target info',logger)
            tarfn = fitsio.read(outf)
            tarfn = Table(tarfn)
            tarfn['TILELOCID'] = 10000*tarfn['TILEID'] +tarfn['LOCATION']
            remcol = ['LOCATION','TILEID']
            for col in remcol:
                try:
                    specf.remove_columns([col])
                except:
                    common.printlog('column '+col +' was not in stacked spec table',logger) 
            common.printlog('added TILELOCID, about to do joins',logger)
            #tj = join(tarfn,specf,keys=['TARGETID','TILELOCID'],join_type='left')
            tjl = []
            selreg = tarfn['DEC'] > 0
            if np.sum(selreg) > 0:
                tjl.append(join(tarfn[selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
                tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
                common.printlog('1st join done',logger)
            if np.sum(~selreg) > 0:
                tjl.append(join(tarfn[~selreg],specf,keys=['TARGETID','TILELOCID'],join_type='left'))
                tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
                common.printlog('2nd join done',logger)
            if len(tjl) > 1:
                tj = vstack(tjl)
            else:
                tj = tjl[0]
            del tjl
            del tarfn
            #tj = np.concatenate(tjl)
            common.printlog('stacked now writing out',logger)
            #tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left') 
            #print(np.unique(tj['ZWARN'],return_counts=True))
            common.write_LSS_scratchcp(tj,outfs,logger=logger)
            #tj.write(outfs,format='fits', overwrite=True)
            common.printlog('joined to spec data and wrote out to '+outfs,logger)

        if uptileloc:
            common.printlog('counting tiles',logger)
            tc = ct.count_tiles_better('dat',tp+notqso,specrel=specrel) 
            common.printlog('writing tile counts',logger)
            common.write_LSS_scratchcp(tj,outtc,logger=logger)
            #tc.write(outtc,format='fits', overwrite=True)


if args.get_petalsky == 'y':
    petalsky_fn = ldirspec+'tile_petal_skydisp_'+prog+'.fits'
    ct.combtile_skystd(tiles4comb,petalsky_fn,specver=specrel,clip=3)

if args.comb_petalqa == 'y':
    petalqa_fn = ldirspec+'tile_petal_qa_'+prog+'.fits'
    ct.combtile_petalqa(tiles4comb,petalqa_fn,specver=specrel)

if specrel != 'daily' and args.dospec == 'y':
    specf.keep_columns(['TARGETID','CHI2','COEFF','Z','ZERR','ZWARN','NPIXELS','SPECTYPE','SUBTYPE','NCOEFF','DELTACHI2'\
    ,'LOCATION','FIBER','COADD_FIBERSTATUS','TILEID','FIBERASSIGN_X','FIBERASSIGN_Y','COADD_NUMEXP','COADD_EXPTIME','COADD_NUMNIGHT'\
    ,'MEAN_DELTA_X','MEAN_DELTA_Y','RMS_DELTA_X','RMS_DELTA_Y','MEAN_PSF_TO_FIBER_SPECFLUX','TSNR2_ELG_B','TSNR2_LYA_B'\
    ,'TSNR2_BGS_B','TSNR2_QSO_B','TSNR2_LRG_B',\
    'TSNR2_ELG_R','TSNR2_LYA_R','TSNR2_BGS_R','TSNR2_QSO_R','TSNR2_LRG_R','TSNR2_ELG_Z','TSNR2_LYA_Z','TSNR2_BGS_Z',\
    'TSNR2_QSO_Z','TSNR2_LRG_Z','TSNR2_ELG','TSNR2_LYA','TSNR2_BGS','TSNR2_QSO','TSNR2_LRG','PRIORITY','DESI_TARGET','BGS_TARGET','TARGET_RA','TARGET_DEC','LASTNIGHT'])
    specfo = ldirspec+'datcomb_'+prog+'_zmtl_zdone.fits'
    outfs = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    if args.redo_zmtl == 'y':
        ct.combtile_spec(tiles4comb,specfo,md='zmtl',specver=specrell[0])   
        fzmtl = fitsio.read(specfo)
        specf = join(specf,fzmtl,keys=['TARGETID','TILEID'])
        specf.write(outfs,format='fits', overwrite=True)

    #if os.path.isfile(outfs):
    specf = Table(fitsio.read(outfs.replace('global','dvs_ro')))
    #else:
    specf.remove_columns(['DESI_TARGET','BGS_TARGET','TARGET_RA','TARGET_DEC','PRIORITY']) #remove these columns because they are in the targets already
    if specrel == 'everest' or specrel =='guadalupe':
        #tarfo = ldirspec+'datcomb_'+prog+'_tarwdup_zdone.fits'
        tps = [prog]
        notqsos = ['']
    else:
        #tar
        if prog == 'dark':
            if args.tracer == 'all':
                tps = ['LRG','ELG','QSO','ELG_LOP','ELG_LOP']
                notqsos = ['','','','','notqso']
            else:
                tps = [args.tracer.strip('notqso')]
                notqsos = ['']
                if 'notqso' in args.tracer:
                    notqsos = ['notqso']
        if prog == 'bright':
            tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']  
            notqsos = ['',''] 
    if args.dotarspec == 'y':
        for tp,notqso in zip(tps,notqsos):
            #first test to see if we need to update any
            logger.info('now doing '+tp+notqso)
            #logger.info(str(len(tiles4comb['TILEID'])))
            #tarfo = dailydir+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
            outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
            #outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'

            #tarf = Table.read(tarfo)
            tarfo = dailydir+'/datcomb_'+prog+'_tarwdup_zdone.fits'
            tarf = fitsio.read(tarfo.replace('global','dvs_ro'))#,columns=cols)
            logger.info('loaded tarwdup file')
            #tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            if tp == 'BGS_BRIGHT':
                sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
            else:
                sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
            if notqso == 'notqso':
                sel &= (tarf['DESI_TARGET'] & 4) == 0
            tarf = Table(tarf[sel])
            logger.info('cut to target type')
            
            tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            logger.info('added TILELOCID, about to do joins')

            remcol = ['Z','ZWARN','FIBER','ZWARN_MTL']
            for col in remcol:
                try:
                    tarf.remove_columns([col] )#we get this where relevant from spec file
                except:
                    logger.info('column '+col +' was not in stacked tarwdup table')    

            #tarf.remove_columns(['ZWARN_MTL'])
            tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            #specf.remove_columns(['PRIORITY'])
            tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID'],join_type='left')
            del tarf
            cols_fromspec = list(specf.dtype.names)
            for col in cols_fromspec:
                if np.ma.is_masked(tj[col]):
                    tj[col] = tj[col].filled(999999)
                    #logger.info(str(np.unique(tj[col],return_counts=True)))
            #del specf
            logger.info('joined tar and spec, now writing')
            #tj.write(outfs,format='fits', overwrite=True)
            common.write_LSS_scratchcp(tj,outfs,logger=logger)
            del tj
            #print('wrote, now counting tiles')
            #tc = ct.count_tiles_better('dat',tp+notqso,specrel=specrel,survey=args.survey) 
            #outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'
            #tc.write(outtc,format='fits', overwrite=True)


#tj.write(ldirspec+'datcomb_'+prog+'_tarspecwdup_zdone.fits',format='fits', overwrite=True)
#tc = ct.count_tiles_better('dat',prog,specrel=specrel)
#tc.write(ldirspec+'Alltiles_'+prog+'_tilelocs.dat.fits',format='fits', overwrite=True)

