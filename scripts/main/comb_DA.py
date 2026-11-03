"""Build combined LSS data products for a survey release.

This script is a command-line driver that assembles several downstream catalog
products used by DESI LSS analyses. Depending on the flags provided, it can:

- identify the set of main-survey tiles that should be included for a given
  survey/release combination,
- write a tile summary file,
- combine repeated target observations into a target-level catalog,
- combine target and spectra information for per-tracer catalogs, and

It intentionally only works for data assemblies, Nevis and beyond:

- the ``--verspec`` value is treated as a release label of the form
  ``<redux>-<zcatalog>`` and uses cumulative tiles and zcatalog files from that
  release.

The script relies on helpers in ``LSS.main.cattools`` and
``LSS.common_tools`` for the heavy lifting. This file mainly orchestrates which
input tiles are eligible, which directories to use, and which combination steps
to run.

Typical usage patterns include:

- ``--make_tile_file``: write the tile list used for the requested combination.
- ``--dotarg``: build the target-with-duplicates combined catalog.
- ``--dotarspec``: join target and spectra products for specific tracers.
- ``--dospec`` / ``--dozmtl``: build spectra and zMTL-based outputs for
  non-daily releases.

The script is designed for batch execution on DESI systems and expects the
usual DESI filesystem layout under ``--basedir``.
"""

# standard python
# works for daily and DA3 onward
from desitarget.targetmask import zwarn_mask
from desispec.maskbits import fibermask
from desispec.exposure_qa import get_qa_params
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
parser = argparse.ArgumentParser(
    description=(
        'Combine DESI LSS spectra '
        'for a chosen survey/data assembly.'
    )
)
parser.add_argument(
    "--basedir", help="base directory for output, default is SCRATCH", default=scratch)
parser.add_argument(
    "--version", help="catalog version; use 'test' unless you know what you are doing!", default='test')
parser.add_argument(
    "--survey", help="e.g., DA3, DA02, any future DA", default='DA3')
parser.add_argument(
    "--prog", help="dark or bright is supported", default='dark')
parser.add_argument("--verspec", help="version for redshifts", default='nevis')
# options about what to run
parser.add_argument("--make_tile_file",
                    help="whether or not to make a tile file", action='store_true')
parser.add_argument(
    "--dotarspec", help="whether or not to combine spec and tar data per type, for non-daily data", action='store_true')
parser.add_argument(
    "--dospec", help="whether or not to combine spec data from beginning", action='store_true')


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

dailydir = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/'+args.survey+'/LSS/daily/'
ldirspec = maindir+specrel+'/'
if not os.path.exists(ldirspec):
    os.mkdir(ldirspec)
    logger.info('made '+ldirspec)
if not os.path.exists(ldirspec+'healpix'):
    os.mkdir(ldirspec+'healpix')
    logger.info('made '+ldirspec+'healpix')


mainp = main(prog, specver=specrel)

mt = mainp.mtld
tiles = mainp.tiles
badfib = mainp.badfib

wd = mt['SURVEY'] == 'main'
wd &= mt['ZDONE'] == 'true'
logger.info('number of tiles with zdone true '+str(len(mt[wd])))
wd &= mt['ARCHIVEDATE'] > 0
logger.info('and with archivedate > 0 '+str(len(mt[wd])))
wd &= mt['FAPRGRM'] == prog
logger.info('and in '+prog+' '+str(len(mt[wd])))
specrell = specrel.split('-')
coaddir = '/global/cfs/cdirs/desi/spectro/redux/' + \
    specrell[0]+'/tiles/cumulative/'

# wd &= np.isin(mt['TILEID'], np.unique(specf['TILEID']))

if args.survey == 'Y1':
    wd &= mt['ZDATE'] < 20220900
if args.survey == 'DA2':
    wd &= mt['ZDATE'] < 20240410
if args.survey == 'DA3':
    wd &= mt['ZDATE'] < 20260415

mtd = mt[wd]
# print('found '+str(len(mtd))+' '+prog+' time main survey tiles that are greater than 85% of goaltime')
logger.info('found '+str(len(mtd))+' '+prog+' time '+args.survey +
            ' survey tiles with zdone true for '+specrel+' version of reduced spectra')


tiles4comb = Table()
tiles4comb['TILEID'] = mtd['TILEID']
tiles4comb['ZDATE'] = mtd['ARCHIVEDATE']
tiles4comb['THRUDATE'] = mtd['LASTNIGHT']  # mtd['LASTNIGHT']


tiles.keep_columns(['TILEID', 'RA', 'DEC'])
# print(tiles.dtype.names)

tiles4comb = join(tiles4comb, tiles, keys=['TILEID'])

logger.info('check that length of tiles4comb matches '+str(len(tiles4comb)))


if args.make_tile_file:
    tiles4comb.write(maindir+'tiles-'+prog.upper() +
                     '.fits', overwrite=True, format='fits')
else:
    print(args.make_tile_file, ' is False, not making tile file')


def qafiberstatus_to_zwarn(qafiberstatus):
    """Compress QAFIBERSTATUS into BAD_SPECQA/BAD_PETALQA ZWARN bits."""
    status = np.asarray(qafiberstatus)
    params = get_qa_params()["exposure_qa"]

    petal_mask = fibermask.mask(params["bad_petal_mask"])
    spec_mask = fibermask.mask(params["bad_qafstatus_mask"]) & ~petal_mask

    zwarn = np.zeros(status.shape, dtype=np.int64)
    zwarn[(status & spec_mask) != 0] |= zwarn_mask.BAD_SPECQA
    zwarn[(status & petal_mask) != 0] |= zwarn_mask.BAD_PETALQA

    return zwarn


outfs = ldirspec+'datcomb_'+prog+'_spec_zdone.fits'

if args.dospec:
    primary_cols = ['TARGETID', 'TILEID', 'LASTNIGHT', 'FIBER', 'COADD_FIBERSTATUS', 'TARGET_RA',
                    'TARGET_DEC', 'FIBERASSIGN_X', 'FIBERASSIGN_Y', 'PRIORITY', 'DESI_TARGET', 'BGS_TARGET', 'COADD_NUMEXP', 'COADD_EXPTIME', 'COADD_NUMNIGHT', 'Z_BEST', 'Z_CONF', 'GOOD_SPEC']

    specf = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/spectro/redux/' +
                              specrell[0]+'/zcatalog/'+specrell[1]+'/main/ztile-main-'+prog.strip('1b')+'-cumulative.fits', columns=primary_cols))  # this will get used again below
    sel_1b = np.isin(specf['TILEID'], np.unique(mtd['TILEID']))
    test_allin = np.isin(mtd['TILEID'], np.unique(specf['TILEID']))
    if np.sum(test_allin) != len(mtd):
        sys.exit('Not all tiles in mtd are in specf')
    specf = specf[sel1b]
    logger.info('length of specf is '+str(len(specf)))

    extra_cols = [
        'TARGETID', 'TILEID', 'Z', 'ZERR', 'ZWARN', 'CHI2', 'COEFF', 'NPIXELS',
        'SPECTYPE', 'SUBTYPE', 'NCOEFF', 'DELTACHI2', 'LOCATION',
        'MEAN_DELTA_X', 'RMS_DELTA_X', 'MEAN_DELTA_Y', 'RMS_DELTA_Y',
        'MEAN_PSF_TO_FIBER_SPECFLUX',
        'TSNR2_BGS', 'TSNR2_ELG', 'TSNR2_LRG', 'TSNR2_LYA', 'TSNR2_QSO',
        'GOOD_Z_BGS', 'GOOD_Z_LRG', 'GOOD_Z_ELG', 'GOOD_Z_QSO', 'GOOD_Z_LYA', 'Z_QSO',
        'QAFIBERSTATUS', 'OII_FLUX', 'OII_FLUX_IVAR', 'OIII_FLUX', 'OIII_FLUX_IVAR'
    ]
    logger.info('reading extra file')
    specfe = Table(fitsio.read('/dvs_ro/cfs/cdirs/desi/spectro/redux/' +
                               specrell[0]+'/zcatalog/'+specrell[1]+'/main/ztile-main-'+prog.strip('1b')+'-cumulative-extra.fits', columns=extra_cols))
    # To construct zmtl ZWARN from Redrock ZWARN:
    zwarn_mtl = specfe['ZWARN'] | qafiberstatus_to_zwarn(
        specfe['QAFIBERSTATUS'])
    specfe['ZWARN_MTL'] = zwarn_mtl
    logger.info('about to join extra columns')
    specf = join(specf, specfe, keys=['TARGETID', 'TILEID'])
    del specfe
    logger.info('joined extra columns, now writing')
    common.write_LSS_scratchcp(specf, outfs, logger=logger)

else:
    specf = Table(fitsio.read(outfs.replace('global', 'dvs_ro')))
    logger.info('loaded specf from '+outfs.replace('global', 'dvs_ro'))


if prog == 'dark':
    # if args.tracer == 'all':
    # tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
    # notqsos = ['','','notqso','','']
    # only do base types because other are subset that can be cut later and this saves i/o
    tps = ['ELG', 'QSO', 'LRG']
    notqsos = ['', '', '']
    # else:
    #    tps = [args.tracer]
    #    notqsos = [args.notqso]
if prog == 'dark1b':
    # if args.tracer == 'all':
    # tps = ['QSO','LRG','ELG_LOP','ELG_LOP','ELG'] #order is not least to most memory intensive
    # notqsos = ['','','notqso','','']
    # only do base types because other are subset that can be cut later and this saves i/o
    tps = ['LGE', 'QSO', 'LRG', 'ELG']
    notqsos = ['', '', '', '']
    # else:
    #    tps = [args.tracer]
    #    notqsos = [args.notqso]

if 'bright' in prog:
    # if args.tracer == 'all':
    # tps = ['BGS_ANY','BGS_BRIGHT']#,'MWS_ANY']
    # notqsos = ['','']
    tps = ['BGS_ANY']  # ,'MWS_ANY']
    notqsos = ['']

    # else:
    #    tps = [args.tracer]
    #    notqsos = [args.notqso]


if args.dotarspec:
    specf.remove_columns(['DESI_TARGET', 'BGS_TARGET',
                         'TARGET_RA', 'TARGET_DEC', 'PRIORITY'])

    for tp, notqso in zip(tps, notqsos):

        fnb = ''
        if '1b' in prog:
            fnb = '_1b'
        outfs = ldirspec+'datcomb_'+tp+notqso+'_tarspecwdup'+fnb+'_zdone.fits'

        logger.info('now doing '+tp+notqso)
        tarfo = dailydir+'/datcomb_'+prog+'_tarwdup_zdone.fits'

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

        tj = join(tarf, specf, keys=[
            'TARGETID', 'LOCATION', 'TILEID'], join_type='left')
        del tarf
        cols_fromspec = list(specf.dtype.names)
        for col in cols_fromspec:
            if np.ma.is_masked(tj[col]):
                dtype = tj[col].dtype
                if np.issubdtype(dtype, np.floating):
                    fillval = 999999.0
                elif np.issubdtype(dtype, np.integer):
                    info = np.iinfo(dtype)
                    # use 999999 where it fits, otherwise the type's max value
                    fillval = 999999 if info.max >= 999999 else info.max
                elif np.issubdtype(dtype, np.bool_):
                    fillval = True
                else:
                    # strings/bytes and anything else: fall back to a string sentinel
                    fillval = '999999'
                tj[col] = tj[col].filled(fillval)

        logger.info('joined tar and spec, now writing')
        common.write_LSS_scratchcp(tj, outfs, logger=logger)
        del tj
