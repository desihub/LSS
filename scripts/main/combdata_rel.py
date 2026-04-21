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
    "--redoqso", help="whether or not to combine qso data, starting over", default='n')
parser.add_argument(
    "--mkemlin", help="whether or not to make emission line files", action='store_true')
parser.add_argument(
    "--dotarspec", help="whether or not to combine spec and tar data per type, for non-daily data", default='y')
parser.add_argument(
    "--redospec", help="whether or not to combine spec data from beginning", default='n')
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
redoqso = False
if args.redoqso == 'y':
    redoqso = True

combpix = True
if args.combpix == 'n':
    combpix = False

redotarspec = False
if args.redotarspec == 'y':
    redotarspec = True

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
if args.check_date_only == 'y':
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


if specrel == 'daily' and args.survey == 'DA2':
    tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'
    tids = tiles4comb['TILEID']  # list()

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
            common.write_LSS(tarsn, tarfo)

        else:
            ct.combtiles_wdup(tiles4comb, fout=tarfo)
    else:
        logger.info('no new tiles to combine')


if args.doqso == 'y':
    outf = ldirspec+'QSO_catalog.fits'
    print('making/adding to QSO catalog '+outf)
    print(tiles4comb.dtype)
    if specrel == 'daily':  # and args.survey == 'main':
        ct.combtile_qso(tiles4comb, outf, restart=redoqso)
    else:
        ct.combtile_qso_alt(tiles4comb, outf, coaddir=coaddir)

if args.mkemlin == 'y':
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
            common.write_LSS(specd, outf)
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
            common.write_LSS(specd, outf)
            del specd
            del tl

        else:
            ct.combtile_em_daily(tiles4comb, outf)

if args.survey == 'Y1' and args.counts_only == 'y':
    if prog == 'dark':
        tps = ['LRG', 'ELG', 'QSO', 'ELG_LOP', 'ELG_LOP']
        notqsos = ['', '', '', '', 'notqso']
    if prog == 'bright':
        tps = ['BGS_ANY', 'BGS_BRIGHT']  # ,'MWS_ANY']
        notqsos = ['', '']
    for tp, notqso in zip(tps, notqsos):

        tc = ct.count_tiles_better(
            'dat', tp+notqso, specrel=specrel, survey=args.survey, badfib=badfib)
        outtc = ldirspec+tp+notqso+'_tilelocs.dat.fits'
        tc.write(outtc, format='fits', overwrite=True)


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


if specrel == 'daily' and args.dospec == 'y' and args.survey != 'main':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    dotarspec = False
    if args.redotarspec == 'y':
        dotarspec = True
    if os.path.isfile(specfo) and args.redospec == 'n':
        specf = Table.read(specfo)
        if list(specf.dtype.names) != speccols:
            logger.info(
                'writing original file back out with subselection of columns')
            specf.keep_columns(speccols)
            common.write_LSS(specf, specfo)
        del specf
    if args.par == 'y':
        tl = []
        if os.path.isfile(specfo):
            specd = fitsio.read(specfo)
            tl.append(specd)
            tdone = np.unique(specd['TILEID'])
            tmask = ~np.isin(tiles4comb['TILEID'], tdone)

        else:
            tmask = np.ones(len(tiles4comb)).astype('bool')
        logger.info('there are '+str(np.sum(tmask)) +
                    ' tiles to add to spec data')
        if np.sum(tmask) == 0:
            newspec = False
        else:
            newspec = True
            tiles_2comb = tiles4comb[tmask]

            def _get_tile(ind):

                trow = tiles_2comb[ind]
                tile, zdate, tdate = trow['TILEID'], trow['ZDATE'], trow['THRUDATE']
                logger.info('combining spec data for TILEID '+str(tile))
                tspec = ct.combspecdata(str(tile), str(zdate), str(tdate))
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
        if newspec:
            specd = np.hstack(tl)
            kp = (specd['TARGETID'] > 0)
            specd = specd[kp]
            common.write_LSS(specd, specfo)
            del specd
        del tl
    else:
        newspec = ct.combtile_spec(
            tiles4comb, specfo, redo=args.redospec, prog=prog, par=args.par)
    specf = Table.read(specfo)
    if newspec:
        logger.info(
            'new tiles were found for spec dataso there were updates to '+specfo)
        dotarspec = True
    else:
        logger.info(
            'no new tiles were found for spec data, so no updates to '+specfo)
    specf['TILELOCID'] = 10000*specf['TILEID'] + specf['LOCATION']
    specf.keep_columns(spec_cols_4tar)
    # tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')

    for tp, notqso in zip(tps, notqsos):
        logger.info('now doing '+tp+notqso)
        logger.info(len(tiles4comb['TILEID']))
        # outf = ldirspec+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        if os.path.isfile(outfs) == False or dotarspec:
            # update = True
            # dotarspec = True

            tarfo = ldirspec+'/datcomb_'+prog+'_tarwdup_zdone.fits'
            # ,columns=cols)
            tarf = fitsio.read(tarfo.replace('global', 'dvs_ro'))
            logger.info('loaded tarspecwdup file')
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


if specrel == 'daily' and args.dospec == 'y' and args.survey == 'main':
    specfo = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'
    if os.path.isfile(specfo) and args.redospec == 'n':
        specf = Table.read(specfo)
        if args.fixspecf == 'y':
            ii = 0
            for tid, adate, zdate in zip(tiles4comb['TILEID'], tiles4comb['ZDATE'], tiles4comb['THRUDATE']):
                ii += 1
                if int(zdate) > 20210730:
                    td = ct.combspecdata(tid, str(adate), str(zdate))
                    kp = (td['TARGETID'] > 0)
                    td = td[kp]
                    sel = specf['TILEID'] == tid
                    fst = specf[sel]
                    if np.array_equal(fst['ZWARN_MTL'], td['ZWARN_MTL']):
                        print('tile '+str(tid)+' passed')
                    else:
                        print('tile '+str(tid)+' is mismatched')
                        specf = specf[~sel]
                        specf = vstack([specf, td])
                        print(ii, len(tiles4comb))
            specf.write(specfo, format='fits', overwrite=True)
        dt = specf.dtype.names
        wo = 0
        if np.isin('FIBERSTATUS', dt):
            sel = specf['COADD_FIBERSTATUS'] == 999999
            specf['COADD_FIBERSTATUS'][sel] = specf['FIBERSTATUS'][sel]
            wo = 1
        if np.isin('NUMEXP', dt):
            sel = specf['COADD_NUMEXP'] == 999999
            sel |= specf['COADD_NUMEXP'] == 16959
            specf['COADD_NUMEXP'][sel] = specf['NUMEXP'][sel]
            wo = 1

        specf.keep_columns(speccols)
        dtn = specf.dtype.names
        if len(dtn) != len(dt):
            wo = 1
        if wo == 1:
            common.printlog('rewriting spec file', logger)
            specf.write(specfo, overwrite=True, format='fits')

    newspec = ct.combtile_spec(
        tiles4comb, specfo, redo=args.redospec, prog=prog)

    if newspec:
        common.printlog(
            'new tiles were found for spec data so there were updates to '+specfo, logger)
    else:
        common.printlog(
            'no new tiles were found for spec data, so no updates to '+specfo, logger)
    specf = fitsio.read(specfo.replace('global', 'dvs_ro'),
                        columns=spec_cols_4tar)
    common.printlog('spec file '+specfo+' has ' +
                    str(len(specf))+' rows', logger)
    if '1b' in prog:
        common.printlog('adding '+specfo.replace('1b', '') +
                        ' to spec info', logger)
        specfnb = fitsio.read(specfo.replace('1b', '').replace(
            'global', 'dvs_ro'), columns=spec_cols_4tar)
        specf = np.concatenate([specf, specfnb])
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

    specf['TILELOCID'] = 10000*specf['TILEID'] + specf['LOCATION']
    # specf.keep_columns(spec_cols_4tar)
    # tj = join(tarf,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')

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

    for tp, notqso in zip(tps, notqsos):
        # first test to see if we need to update any
        common.printlog('now doing '+tp+notqso, logger)
        # print(len(tiles4comb['TILEID']))
        outf = ldirspec+'datcomb_'+tp+notqso+'_tarwdup_zdone.fits'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup_zdone.fits'
        outtc = ldirspec+tp+notqso+'_tilelocs.dat.fits'
        update = True
        uptileloc = False
        dotarspec = True
        hpxsn = hpxs
        s = 0
        if os.path.isfile(outf):
            fo = fitsio.read(outf.replace('global', 'dvs_ro'),
                             columns=['TARGETID', 'TILEID'])
            common.printlog(outf + ' has '+str(len(fo)) + ' rows', logger)
            nstid = len(tiles4comb['TILEID'])
            notid = len(np.unique(fo['TILEID']))
            test_tid = np.isin(tiles4comb['TILEID'], np.unique(fo['TILEID']))
            common.printlog('there are '+str(nstid-np.sum(test_tid)) +
                            ' tiles that need to be added to '+outf, logger)
            if nstid == np.sum(test_tid):
                update = False
                common.printlog('we will not update '+outf +
                                ' because there are no new tiles', logger)
            else:

                tidc = ~np.isin(tiles4comb['TILEID'], np.unique(fo['TILEID']))
                # print('the new tileids are '+str(tiles4comb['TILEID'][tidc]))
                # print(len(tiles4comb[tidc]))
                hpxsn = foot.tiles2pix(8, tiles=tiles4comb[tidc])
            del fo
        if os.path.isfile(outfs):
            fo = fitsio.read(outfs.replace('global', 'dvs_ro'), columns=[
                             'TARGETID', 'TILEID', 'ZWARN', 'ZWARN_MTL'])
            stids = np.unique(fo['TILEID'])
            if len(stids) == notid:
                dotarspec = False
            if os.path.isfile(outtc) and update == False and redotarspec == False and dotarspec == False:
                ftc = fitsio.read(outtc.replace(
                    'global', 'dvs_ro'), columns=['TARGETID'])
                fc = ct.cut_specdat(fo)
                ctid = np.isin(fc['TARGETID'], ftc['TARGETID'])
                if len(ctid) == sum(ctid):
                    common.printlog('all targetids are in '+outtc+' and all tileids are in ' +
                                    outf+' so '+outtc+' will not be updated', logger)
                    uptileloc = False
                del ftc
                del fc
            del fo
        if args.counts_only != 'y' and update:
            common.printlog('updating '+outf, logger)
            cols = None
            if os.path.isfile(outf):
                tarfn = fitsio.read(outf.replace('global', 'dvs_ro'), rows=1)
                cols = tarfn.dtype.names
                # if np.isin('TILELOCID',tarfn.dtype.names):
                # common.printlog('reloading '+outf+' without reading TILELOCID column',logger)
                # sel = cols != 'TILELOCID'
                # cols = cols[sel]
                cols = []
                for col in tarfn.dtype.names:
                    if col != 'TILELOCID':
                        cols.append(col)
                tarfn = fitsio.read(outf.replace(
                    'global', 'dvs_ro'), columns=cols)
                # print(tarfn.dtype.names)
                theta, phi = np.radians(
                    90-tarfn['DEC']), np.radians(tarfn['RA'])
                tpix = hp.ang2pix(8, theta, phi, nest=True)
                pin = np.isin(tpix, hpxsn)
                # remove the rows for the healpix that will updated
                tarfn = tarfn[~pin]
                s = 1
            common.printlog(
                'after removing healpix to update, has '+str(len(tarfn)), logger)
            npx = 0
            tarfl = []
            if s == 1:
                tarfl = [tarfn]

            for px in hpxsn:
                tarfo = ldirspec+'healpix/datcomb_' + \
                    prog+'_'+str(px)+'_tarwdup_zdone.fits'
                if '1b' in prog:
                    tarfonb = tarfo.replace('1b', '')

                if os.path.isfile(tarfo):
                    if cols is not None:
                        tarf = fitsio.read(tarfo.replace(
                            'global', 'dvs_ro'), columns=cols)
                        if '1b' in prog and os.path.isfile(tarfonb):
                            tarfnb = fitsio.read(tarfonb.replace(
                                'global', 'dvs_ro'), columns=cols)
                    else:
                        tarf = fitsio.read(tarfo.replace('global', 'dvs_ro'))
                        if '1b' in prog and os.path.isfile(tarfonb):
                            tarfnb = fitsio.read(
                                tarfonb.replace('global', 'dvs_ro'))
                    if '1b' in prog and os.path.isfile(tarfonb):
                        tarf = np.hstack((tarf, tarfnb))
                    # tarf['TILELOCID'] = 10000*tarf['TILEID'] +tarf['LOCATION']
                    if tp == 'BGS_BRIGHT':
                        sel = tarf['BGS_TARGET'] & targetmask.bgs_mask[tp] > 0
                    else:
                        sel = tarf['DESI_TARGET'] & targetmask.desi_mask[tp] > 0
                    if notqso == 'notqso':
                        sel &= (tarf['DESI_TARGET'] & 4) == 0
                    tarfl.append(tarf[sel])
                    # if s == 0:
                    #    tarfn = tarf[sel]
                    #    s = 1
                    # else:
                    # tarfn = vstack([tarfn,tarf[sel]],metadata_conflicts='silent')
                    # common.printlog(tarfn.dtype.names,logger)
                    # common.printlog(tarf.dtype.names,logger)
                    #    tarfn = np.hstack((tarfn,tarf[sel]))
                    common.printlog(
                        str(len(tarf[sel]))+','+tp+notqso+','+str(npx)+','+str(len(hpxsn)), logger)
                else:
                    common.printlog('file '+tarfo+' not found', logger)
                npx += 1
            common.printlog('concatenating', logger)
            tarfn = np.hstack(tarfl)
            common.printlog('now length '+str(len(tarfn)), logger)
            tarfn = Table(tarfn)
            remcol = ['Z', 'ZWARN', 'FIBER', 'ZWARN_MTL', 'PRIORITY']
            for col in remcol:
                try:
                    # we get this where relevant from spec file
                    tarfn.remove_columns([col])
                except:
                    common.printlog(
                        'column '+col + ' was not in stacked tarwdup table', logger)

            common.write_LSS_scratchcp(tarfn, outf, logger=logger)
            # tarfn.write(outf,format='fits', overwrite=True)
            common.printlog('wrote out '+outf+' ' +
                            str(len(tarfn))+' rows', logger)

            # try:
            #    specf.remove_columns(['PRIORITY'])
            # except:
            #    print('column PRIORITY was not in spec table')

            # join to spec info; now only do so after updating 1b tiles
            if '1b' in prog:
                tarfn['TILELOCID'] = 10000*tarfn['TILEID'] + tarfn['LOCATION']
                print('added TILELOCID, about to do joins')
                # tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')

                # seems to run out of memory on join
                tjl = []
                print(tarfn.dtype.names)
                selreg = tarfn['DEC'] > 0
                print(len(tarfn[selreg]))
                remcol = ['LOCATION', 'TILEID']
                for col in remcol:
                    try:
                        specf.remove_columns([col])
                    except:
                        print('column '+col + ' was not in stacked spec table')
                if np.sum(selreg) > 0:
                    tjl.append(join(tarfn[selreg], specf, keys=[
                               'TARGETID', 'TILELOCID'], join_type='left'))
                    tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
                    common.printlog('1st join done', logger)
                if np.sum(~selreg) > 0:
                    tjl.append(join(tarfn[~selreg], specf, keys=[
                               'TARGETID', 'TILELOCID'], join_type='left'))
                    tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
                    common.printlog('2nd join done', logger)
                del tarfn
                if len(tjl) > 1:
                    tj = vstack(tjl)
                else:
                    tj = tjl[0]
                del tjl
                common.printlog('stacked now writing out has ' +
                                str(len(tj)) + ' rows', logger)
                # for reg in regl:
                #    sel = tarfn['PHOTSYS'] == reg
                #    tjr = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
                # tj.write(outfs,format='fits', overwrite=True)
                common.write_LSS_scratchcp(tj, outfs, logger=logger)
                common.printlog(
                    'joined to spec data and wrote out to '+outfs, logger)
        elif (redotarspec or dotarspec) and '1b' in prog:
            common.printlog('joining spec info to target info', logger)
            tarfn = fitsio.read(outf.replace('global', 'dvs_ro'))
            tarfn = Table(tarfn)
            remcol = ['Z', 'ZWARN', 'FIBER', 'ZWARN_MTL', 'PRIORITY']
            for col in remcol:
                try:
                    # we get this where relevant from spec file
                    tarfn.remove_columns([col])
                except:
                    common.printlog(
                        'column '+col + ' was not in stacked tarwdup table', logger)

            tarfn['TILELOCID'] = 10000*tarfn['TILEID'] + tarfn['LOCATION']
            remcol = ['LOCATION', 'TILEID']
            for col in remcol:
                try:
                    specf.remove_columns([col])
                except:
                    common.printlog(
                        'column '+col + ' was not in stacked spec table', logger)
            common.printlog('added TILELOCID, about to do joins', logger)
            # tj = join(tarfn,specf,keys=['TARGETID','TILELOCID'],join_type='left')
            tjl = []
            selreg = tarfn['DEC'] > 0
            if np.sum(selreg) > 0:
                tjl.append(join(tarfn[selreg], specf, keys=[
                           'TARGETID', 'TILELOCID'], join_type='left'))
                tjl[0]['ZWARN'] = tjl[0]['ZWARN'].filled(999999)
                common.printlog('1st join done', logger)
            if np.sum(~selreg) > 0:
                tjl.append(join(tarfn[~selreg], specf, keys=[
                           'TARGETID', 'TILELOCID'], join_type='left'))
                tjl[1]['ZWARN'] = tjl[1]['ZWARN'].filled(999999)
                common.printlog('2nd join done', logger)
            if len(tjl) > 1:
                tj = vstack(tjl)
            else:
                tj = tjl[0]
            del tjl
            del tarfn
            # tj = np.concatenate(tjl)
            common.printlog('stacked now writing out', logger)
            # tj = join(tarfn,specf,keys=['TARGETID','LOCATION','TILEID','TILELOCID'],join_type='left')
            # print(np.unique(tj['ZWARN'],return_counts=True))
            common.write_LSS_scratchcp(tj, outfs, logger=logger)
            # tj.write(outfs,format='fits', overwrite=True)
            common.printlog(
                'joined to spec data and wrote out to '+outfs, logger)

        if uptileloc:
            common.printlog('counting tiles', logger)
            tc = ct.count_tiles_better('dat', tp+notqso, specrel=specrel)
            common.printlog('writing tile counts', logger)
            common.write_LSS_scratchcp(tj, outtc, logger=logger)
            # tc.write(outtc,format='fits', overwrite=True)


if args.get_petalsky == 'y':
    petalsky_fn = ldirspec+'tile_petal_skydisp_'+prog+'.fits'
    ct.combtile_skystd(tiles4comb, petalsky_fn, specver=specrel, clip=3)

if args.comb_petalqa == 'y':
    petalqa_fn = ldirspec+'tile_petal_qa_'+prog+'.fits'
    ct.combtile_petalqa(tiles4comb, petalqa_fn, specver=specrel)

if specrel != 'daily' and args.dospec == 'y':
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
    if args.dotarspec == 'y':
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
