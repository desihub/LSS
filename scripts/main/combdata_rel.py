# standard python
import logging
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
from astropy.table import Table, join, unique, vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desimodel.footprint import is_point_in_desi
import desimodel.footprint as foot
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
from LSS.globals import main

scratch = os.getenv('SCRATCH')
parser = argparse.ArgumentParser()
parser.add_argument(
    "--basedir", help="base directory for output, default is SCRATCH", default=scratch)
parser.add_argument(
    "--version", help="catalog version; use 'test' unless you know what you are doing!", default='test')
parser.add_argument(
    "--survey", help="e.g., DA3, DA02, any future DA", default='DA3')
parser.add_argument(
    "--prog", help="dark or bright is supported", default='dark')
parser.add_argument("--verspec", help="version for redshifts", default='daily')
parser.add_argument("--check_date_only",
                    help="whether or not to stop after maximum night is found", default='n')
parser.add_argument("--make_tile_file",
                    help="whether or not to make a tile file", action='store_true')
parser.add_argument(
    "--doqso", help="whether or not to combine qso data", action='store_true')
parser.add_argument(
    "--redoqso", help="whether or not to combine qso data, starting over", action='store_true')
parser.add_argument(
    "--mkemlin", help="whether or not to make emission line files", action='store_true')
parser.add_argument(
    "--dotarg", help="whether or not to combine spec and tar data per type, for non-daily data", action='store_true')
parser.add_argument(
    "--dotarspec", help="whether or not to combine spec and tar data per type, for non-daily data", action='store_true')
parser.add_argument(
    "--dospec", help="whether or not to combine spec data from beginning", action='store_true')

parser.add_argument(
    "--redospec", help="whether or not to combine spec data from beginning", action='store_true')
parser.add_argument("--par", help="use multiprocessing ", action='store_true')
parser.add_argument("--tracer", help="tracer type", default='all')
parser.add_argument("--notqso", help="tracer type", default='')


args = parser.parse_args()
print(args)

# create logger
logname = 'comb_inputs'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter(
    '%(asctime)s - %(name)s - %(levelname)s - %(message)s')

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

mainp = main(prog, specver=specrel)

mt = mainp.mtld
tiles = mainp.tiles
badfib = mainp.badfib

wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
logger.info('number of tiles with zdone true '+str(len(mt[wd])))
wd &= mt['ARCHIVEDATE'] > 0
logger.info('and with archivedate > 0 '+str(len(mt[wd])))
if args.doqso == 'y' or args.mkemlin == 'y':
    wd &= ((mt['FAPRGRM'] == prog) | (mt['FAPRGRM'] == prog+'1b'))
    common.printlog('tiles being considered from fa programs ' +
                    str(np.unique(mt[wd]['FAPRGRM'])), logger)
else:
    wd &= mt['FAPRGRM'] == prog
logger.info('and in '+prog+' '+str(len(mt[wd])))
if specrel != 'daily':
    specrell = specrel.split('-')
    coaddir = '/global/cfs/cdirs/desi/spectro/redux/' + \
        specrell[0]+'/tiles/cumulative/'
    specf = Table.read('/global/cfs/cdirs/desi/spectro/redux/' +
                       specrell[0]+'/zcatalog/'+specrell[1]+'/ztile-main-'+prog+'-cumulative.fits')
    wd &= np.isin(mt['TILEID'], np.unique(specf['TILEID']))
else:
    coaddir = '/global/cfs/cdirs/desi/spectro/redux/daily/tiles/archive/'

if args.survey == 'Y1':
    wd &= mt['ZDATE'] < 20220900
if args.survey == 'DA2':
    wd &= mt['ZDATE'] < 20240410
if args.survey == 'DA3':
    wd &= mt['ZDATE'] < 20260417

mtd = mt[wd]
# print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
logger.info('found '+str(len(mtd))+' '+prog+' time '+args.survey +
            ' survey tiles with zdone true for '+specrel+' version of reduced spectra')


tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
if args.verspec == 'daily':
    tiles4comb['THRUDATE'] = mtd['ZDATE']  # mtd['LASTNIGHT']
else:
    tiles4comb['THRUDATE'] = mtd['LASTNIGHT']  # mtd['LASTNIGHT']


logger.info('The last night of data that will be processed is for ' +
            args.prog+' is '+str(np.max(tiles4comb['THRUDATE'])))
logger.info('Is that what was expected based on MTL updates?')
if args.check_date_only:
    sys.exit()

tiles.keep_columns(['TILEID', 'RA', 'DEC'])
# print(tiles.dtype.names)

tiles4comb = join(tiles4comb, tiles, keys=['TILEID'])

logger.info('check that length of tiles4comb matches '+str(len(tiles4comb)))

# share basedir location '/global/cfs/cdirs/desi/survey/catalogs'
maindir = basedir + '/'+args.survey+'/LSS/'

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

if args.make_tile_file:
    tiles4comb.write(maindir+'tiles-'+prog.upper() +
                     '.fits', overwrite=True, format='fits')
else:
    print(args.make_tile_file, ' is False, not making tile file')

logger.info('specrel is '+specrel)
if specrel == 'daily':
    specfo = basedir+'main/LSS/daily/datcomb_' + \
        prog.replace('1b', '')+'_spec_zdone.fits'

    if not os.path.isfile(specfo):
        sys.exit('daily spec file '+specfo+' not found, exiting')

    spec_cols_4tar = ['TARGETID', 'ZWARN', 'ZWARN_MTL', 'LOCATION', 'FIBER', 'TILEID',
                      'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'PRIORITY']
    logger.info(str(spec_cols_4tar))

regl = ['N', 'S']


if args.dotarg:
    tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'
    tids = tiles4comb['TILEID']  # list()

    def _tab2list(tid):
        sel = tiles4comb['TILEID'] == tid
        logger.info('at TILEID '+str(tid))
        tl_tab = tiles4comb[sel]
        tab = ct.get_tiletab(tl_tab)
        return tab

    if os.path.isfile(tarfo) == False:  # or args.redotardup::
        logger.info('creating '+tarfo)
        tile_list = []
        tids_todo = tids
    else:
        logger.info('not remaking '+tarfo)
        tile_list = [fitsio.read(tarfo)]
        tids_c = np.unique(tile_list[0]['TILEID'])
        tids_todo_inds = ~np.isin(tids, tids_c)
        tids_todo = tids[tids_todo_inds]
    if len(tids_todo) > 0:
        if args.par:
            from concurrent.futures import ProcessPoolExecutor

            with ProcessPoolExecutor() as executor:
                for tab in executor.map(_tab2list, list(tids_todo)):
                    tile_list.append(np.array(tab))
            logger.info('tiles in list of length '+str(len(tile_list)))
            logger.info('concatenating')
            logger.info(str(tile_list[0].dtype.names))
            # vstack(tile_list,metadata_conflicts='silent')
            tarsn = np.concatenate(tile_list)
            logger.info(str(tarsn.dtype.names))
            del tile_list
            logger.info('doing TARGETID sort')
            tarsn = Table(tarsn)
            tarsn.sort('TARGETID')

            logger.info('sort done')
            common.write_LSS_scratchcp(tarsn, tarfo, logger=logger)

        else:
            ct.combtiles_wdup(tiles4comb, fout=tarfo)
    else:
        logger.info('no new tiles to combine')


if args.doqso:
    outf = ldirspec+'QSO_catalog.fits'
    print('making/adding to QSO catalog '+outf)
    print(tiles4comb.dtype)
    if specrel == 'daily':  # and args.survey == 'main':
        ct.combtile_qso(tiles4comb, outf, restart=redoqso)
    else:
        ct.combtile_qso_alt(tiles4comb, outf, coaddir=coaddir)

if args.mkemlin:
    outf = ldirspec+'emlin_catalog.fits'
    if specrel == 'daily' and args.survey == 'main':
        outdir = '/global/cfs/cdirs/desi/survey/catalogs/main/LSS/daily/emtiles/'
        guadtiles = fitsio.read(
            '/global/cfs/cdirs/desi/survey/catalogs/DA02/LSS/guadalupe/datcomb_'+prog+'_spec_zdone.fits', columns=['TILEID'])
        guadtiles = np.unique(guadtiles['TILEID'])
        gtids = np.isin(tiles4comb['TILEID'], guadtiles)
        tiles4em = tiles4comb[~gtids]
        ndone = 0
        for tile, zdate, tdate in zip(tiles4em['TILEID'], tiles4em['ZDATE'], tiles4em['THRUDATE']):
            outft = outdir+'emline-'+str(tile)+'.fits'
            if not os.path.isfile(outf):
                tdate = str(tdate)
                ct.combEMdata_daily_old(tile, zdate, tdate, outf=outft)
                print('wrote '+outf)
                ndone += 1
                print('completed '+str(ndone)+' tiles')
        ct.combtile_em(tiles4comb, outf)
    elif specrel != 'daily':
        if args.par == 'y':
            tl = []
            if os.path.isfile(outf):
                specd = fitsio.read(outf)
                tl.append(specd)
                tdone = np.unique(specd['TILEID'])
                tmask = ~np.isin(tiles4comb['TILEID'], tdone)
            else:
                tmask = np.ones(len(tiles4comb)).astype('bool')
            logger.info('there are '+str(np.sum(tmask)) +
                        ' tiles to add to emline data')
            if np.sum(tmask) == 0:
                newem = False
            else:
                newem = True
                tiles_2comb = tiles4comb[tmask]

                def _get_tile(ind):

                    trow = tiles_2comb[ind]
                    tile, zdate, tdate = trow['TILEID'], trow['ZDATE'], trow['THRUDATE']
                    logger.info('combining emline data for TILEID '+str(tile))
                    # tspec = ct.combEMdata_daily(str(tile),str(zdate),str(tdate))
                    tspec = ct.combEMdata_rel(str(tile), str(tdate), coaddir)
                    if tspec is not None:
                        # tspec['TILEID'] = tile
                        tspec = np.array(tspec)
                        logger.info(
                            'tile '+str(tile)+' '+str(len(tspec.dtype.names))+' '+str(len(tspec)))
                        return tspec
                        # new = np.empty(len(tspec),dtype=specd.dtype)
                        # cols = specd.dtype.names
                        # for colname in cols:
                        #    new[colname][...] = tspec[colname][...]
                    else:
                        logger.info('tile '+str(tile)+' failed for emline')
                        return None
                        # new = None
                    # return new
                inds = np.arange(len(tiles_2comb))
                from concurrent.futures import ProcessPoolExecutor

                with ProcessPoolExecutor() as executor:
                    for specd in executor.map(_get_tile, inds):
                        if specd is not None:
                            tl.append(np.array(specd))
            # nms = list(tl[0].dtype.names)
            # for t in tl:
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
            common.write_LSS_scratchcp(specd, outf, logger=logger)
            del specd
            del tl

        else:
            ct.combtile_em_alt(tiles4comb, outf, prog='dark', coaddir=coaddir)
    else:
        if args.par == 'y':
            tl = []
            if os.path.isfile(outf):
                specd = fitsio.read(outf)
                tl.append(specd)
                tdone = np.unique(specd['TILEID'])
                tmask = ~np.isin(tiles4comb['TILEID'], tdone)
            else:
                tmask = np.ones(len(tiles4comb)).astype('bool')

            logger.info('there are '+str(np.sum(tmask)) +
                        ' tiles to add to emline data')
            if np.sum(tmask) == 0:
                newem = False
            else:
                newem = True
                tiles_2comb = tiles4comb[tmask]

                def _get_tile(ind):

                    trow = tiles_2comb[ind]
                    tile, zdate, tdate = trow['TILEID'], trow['ZDATE'], trow['THRUDATE']
                    logger.info('combining emline data for TILEID '+str(tile))
                    tspec = ct.combEMdata_daily(
                        str(tile), str(zdate), str(tdate))
                    if tspec:
                        tspec['TILEID'] = tile
                        tspec = np.array(tspec)
                        new = np.empty(len(tspec), dtype=specd.dtype)
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
            common.write_LSS_scratchcp(specd, outf, logger=logger)
            del specd
            del tl

        else:
            ct.combtile_em_daily(tiles4comb, outf)


if prog == 'dark':
    if args.tracer == 'all':
        # tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
        # notqsos = ['','','notqso','','']
        # only do base types because other are subset that can be cut later and this saves i/o
        tps = ['QSO', 'LRG', 'ELG']
        notqsos = ['', '', '']
    else:
        tps = [args.tracer]
        notqsos = [args.notqso]
if prog == 'dark1b':
    if args.tracer == 'all':
        # tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
        # notqsos = ['','','notqso','','']
        # only do base types because other are subset that can be cut later and this saves i/o
        tps = ['LGE', 'QSO', 'LRG', 'ELG']
        notqsos = ['', '', '', '']
    else:
        tps = [args.tracer]
        notqsos = [args.notqso]

if 'bright' in prog:
    if args.tracer == 'all':
        # tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']
        # notqsos = ['','']
        tps = ['BGS_ANY']  # ,'MWS_ANY']
        notqsos = ['']

    else:
        tps = [args.tracer]
        notqsos = [args.notqso]


if args.dotarspec and specrel == 'daily':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'

    specf = Table(fitsio.read(specfo, columns=spec_cols_4tar))

    specf['TILELOCID'] = 10000*specf['TILEID'] + specf['LOCATION']

    for tp, notqso in zip(tps, notqsos):
        logger.info('now doing '+tp+notqso)
        logger.info(len(tiles4comb['TILEID']))

        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        if os.path.isfile(outfs) == False or args.dotarspec:

            tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'

            tarf = fitsio.read(tarfo.replace('global', 'dvs_ro'))
            logger.info('loaded tarspecwdup file')

            if tp == 'BGS_BRIGHT':
                sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
            else:
                sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
            if notqso == 'notqso':
                sel &= (tarf['DESI_TARGET'] & 4) == 0
            tarf = Table(tarf[sel])
            logger.info('cut to target type')

            tarf['TILELOCID'] = 10000*tarf['TILEID'] + tarf['LOCATION']
            logger.info('added TILELOCID, about to do joins')
            # tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')

            # seems to run out of memory on join
            tjl = []
            logger.info(tarf.dtype.names)
            selreg = tarf['DEC'] > 0
            logger.info(len(tarf[selreg]))
            remcol = ['LOCATION', 'TILEID', 'FIBER', 'PRIORITY']
            for col in remcol:
                try:
                    specf.remove_columns([col])
                except:
                    logger.info('column '+col +
                                ' was not in stacked spec table')
            tjl.append(join(tarf[selreg], specf, keys=[
                       'TARGETID', 'TILELOCID'], join_type='left'))
            tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
            logger.info('1st join done')
            tjl.append(join(tarf[~selreg], specf, keys=[
                       'TARGETID', 'TILELOCID'], join_type='left'))
            tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
            logger.info('2nd join done')
            del tarf
            tj = vstack(tjl)
            logger.info('stacked now writing out')
            common.write_LSS_scratchcp(tj, outfs, logger=logger)
            logger.info('joined to spec data and wrote out to '+outfs)
        else:
            logger.info(outfs + ' exists already, not making again')


if specrel != 'daily' and args.dospec:
    specf.keep_columns(['TARGETID', 'CHI2', 'COEFF', 'Z', 'ZERR', 'ZWARN', 'NPIXELS', 'SPECTYPE', 'SUBTYPE', 'NCOEFF', 'DELTACHI2', 'LOCATION', 'FIBER', 'COADD_FIBERSTATUS', 'TILEID', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'MEAN_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_X', 'RMS_DELTA_Y', 'MEAN_PSF_TO_FIBER_SPECFLUX', 'TSNR2_ELG_B', 'TSNR2_LYA_B', 'TSNR2_BGS_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B',
                        'TSNR2_ELG_R', 'TSNR2_LYA_R', 'TSNR2_BGS_R', 'TSNR2_QSO_R', 'TSNR2_LRG_R', 'TSNR2_ELG_Z', 'TSNR2_LYA_Z', 'TSNR2_BGS_Z',
                        'TSNR2_QSO_Z', 'TSNR2_LRG_Z', 'TSNR2_ELG', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_QSO', 'TSNR2_LRG', 'PRIORITY', 'DESI_TARGET', 'BGS_TARGET', 'TARGET_RA', 'TARGET_DEC', 'LASTNIGHT'])
    specfo = ldirspec+'datcomb_'+prog+'_zmtl_zdone.fits'
    outfs = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    if args.redo_zmtl == 'y':
        ct.combtile_spec(tiles4comb, specfo, md='zmtl', specver=specrell[0])
        fzmtl = fitsio.read(specfo)
        specf = join(specf, fzmtl, keys=['TARGETID', 'TILEID'])
        specf.write(outfs, format='fits', overwrite=True)

    # if os.path.isfile(outfs):
    specf = Table(fitsio.read(outfs.replace('global', 'dvs_ro')))
    # else:
    # remove these columns because they are in the targets already
    specf.remove_columns(['DESI_TARGET', 'BGS_TARGET',
                         'TARGET_RA', 'TARGET_DEC', 'PRIORITY'])
    if specrel == 'everest' or specrel == 'guadalupe':
        # tarfo = ldirspec+'datcomb_'+prog+'_tarwdup_zdone.fits'
        tps = [prog]
        notqsos = ['']
    else:
        # tar
        if prog == 'dark':
            if args.tracer == 'all':
                tps = ['LRG', 'ELG', 'QSO', 'ELG_LOP', 'ELG_LOP']
                notqsos = ['', '', '', '', 'notqso']
            else:
                tps = [args.tracer.strip('notqso')]
                notqsos = ['']
                if 'notqso' in args.tracer:
                    notqsos = ['notqso']
        if prog == 'bright':
            tps = ['BGS_ANY', 'BGS_BRIGHT']  # ,'MWS_ANY']
            notqsos = ['', '']
    if args.dotarspec:
        for tp, notqso in zip(tps, notqsos):
            # first test to see if we need to update any
            logger.info('now doing '+tp+notqso)
            # logger.info(str(len(tiles4comb['TILEID'])))
            # tarfo = dailydir+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
            outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
            # outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'

            # tarf = Table.read(tarfo)
            tarfo = dailydir+'/datcomb_'+prog+'_tarwdup_zdone.fits'
            # ,columns=cols)
            tarf = fitsio.read(tarfo.replace('global', 'dvs_ro'))
            logger.info('loaded tarwdup file')
            # tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
            if tp == 'BGS_BRIGHT':
                sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
            else:
                sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
            if notqso == 'notqso':
                sel &= (tarf['DESI_TARGET'] & 4) == 0
            tarf = Table(tarf[sel])
            logger.info('cut to target type')

            tarf['TILELOCID'] = 10000*tarf['TILEID'] + tarf['LOCATION']
            logger.info('added TILELOCID, about to do joins')

            remcol = ['Z', 'ZWARN', 'FIBER', 'ZWARN_MTL']
            for col in remcol:
                try:
                    # we get this where relevant from spec file
                    tarf.remove_columns([col])
                except:
                    logger.info('column '+col +
                                ' was not in stacked tarwdup table')

            # tarf.remove_columns(['ZWARN_MTL'])
            tarf['TILELOCID'] = 10000*tarf['TILEID'] + tarf['LOCATION']
            # specf.remove_columns(['PRIORITY'])
            tj = join(tarf, specf, keys=[
                      'TARGETID', 'LOCATION', 'TILEID'], join_type='left')
            del tarf
            cols_fromspec = list(specf.dtype.names)
            for col in cols_fromspec:
                if np.ma.is_masked(tj[col]):
                    tj[col] = tj[col].filled(999999)
                    # logger.info(str(np.unique(tj[col],return_counts=True)))
            # del specf
            logger.info('joined tar and spec, now writing')
            # tj.write(outfs,format='fits', overwrite=True)
            common.write_LSS_scratchcp(tj, outfs, logger=logger)
            del tj
            # print('wrote, now counting tiles')
            # tc = ct.count_tiles_better('dat',tp+notqso,specrel=specrel,survey=args.survey)
            # outtc =  ldirspec+tp+notqso+'_tilelocs.dat.fits'
            # tc.write(outtc,format='fits', overwrite=True)


# tj.write(ldirspec+'datcomb_'+prog+'_tarspecwdup_zdone.fits',format='fits', overwrite=True)
# tc = ct.count_tiles_better('dat',prog,specrel=specrel)
# tc.write(ldirspec+'Alltiles_'+prog+'_tilelocs.dat.fits',format='fits', overwrite=True)
