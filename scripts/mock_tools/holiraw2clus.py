# standard python
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
from astropy.table import Table, join, unique, vstack
from matplotlib import pyplot as plt
from desitarget.io import read_targets_in_tiles
from desitarget.mtl import inflate_ledger
from desitarget import targetmask
from desitarget.internal import sharedmem
from desimodel.footprint import is_point_in_desi
from desitarget import targetmask

import LSS.main.cattools as ct
import LSS.common_tools as common
import LSS.mocktools as mocktools
from LSS.globals import main

import logging
logger = logging.getLogger('mkCat')
logger.setLevel(level=logging.INFO)

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

logger.info('run started')

# import LSS.mkCat_singletile.fa4lsscat as fa
# from LSS.globals import main

if os.environ['NERSC_HOST'] == 'perlmutter':
    scratch = 'PSCRATCH'
else:
    logger.info('NERSC_HOST is not cori or permutter but is ' +
                os.environ['NERSC_HOST'])
    sys.exit('NERSC_HOST not known (code only works on NERSC), not proceeding')


parser = argparse.ArgumentParser()
parser.add_argument("--realization", default=1,
                    help="which realization to use, default is 1", type=int)
parser.add_argument("--base_dir", help="base directory for input/output",
                    default='/global/cfs/cdirs/desi/mocks/cai/holi/webjax_v4.80/')
parser.add_argument("--data_dir", help="where to find the data randoms",
                    default='/global/cfs/cdirs/desi/survey/catalogs/DA2/LSS/loa-v1/')
parser.add_argument(
    "--minr", help="minimum number for random files", default=0, type=int)
parser.add_argument(
    "--maxr", help="maximum for random files, default is all 18)", default=18, type=int)

parser.add_argument("--tracer", default='all')
parser.add_argument("--outloc", default=None)
parser.add_argument("--outmd", help='write out in h5 or fits',
                    choices=['.h5', '.fits'], default='.h5')
parser.add_argument("--par", default='y',
                    help='whether to run random steps in parallel or not')
parser.add_argument("--mkdat", default='y')
parser.add_argument("--mkran", default='y')
parser.add_argument("--nz", default='y')
parser.add_argument("--splitGC", default='y')


args = parser.parse_args()
logger.info(args)

rm = int(args.minr)
rx = int(args.maxr)


if args.tracer == 'all':
    tracers = ['QSO', 'LRG', 'ELG']
else:
    tracers = [args.tracer]

logger.info(tracers)


def read_file(fn, columns=None):
    if '.fits' in fn:
        data = Table(fitsio.read(fn.replace('global', 'dvs_ro')))
        if columns is not None:
            data.keep_columns(columns)
    if '.h5' in fn:
        data = common.read_hdf5_blosc(fn.replace(
            'global', 'dvs_ro'), columns=columns)
    return data


def splitGC(flroot, datran='.dat', rann=0):
    import LSS.common_tools as common
    from astropy.coordinates import SkyCoord
    import astropy.units as u
    app = 'clustering'+datran+args.outmd
    if datran == '.ran':
        app = str(rann)+'_clustering'+datran+args.outmd

    # Table(fitsio.read(flroot.replace('global','dvs_ro') +app))
    fn = read_file(flroot.replace('global', 'dvs_ro') + app)
    sel_ngc = common.splitGC(fn)  # gc.b > 0
    outf_ngc = flroot+'NGC_'+app

    outf_sgc = flroot+'SGC_'+app
    if args.outmd == '.fits':
        common.write_LSS_scratchcp(fn[sel_ngc], outf_ngc, logger=logger)
        common.write_LSS_scratchcp(fn[~sel_ngc], outf_sgc, logger=logger)

    if args.outmd == '.h5':
        common.write_LSShdf5_scratchcp(fn[sel_ngc], outf_ngc, logger=logger)
        common.write_LSShdf5_scratchcp(fn[~sel_ngc], outf_sgc, logger=logger)


def ran_col_assign(randoms, data, sample_columns, tracer, seed=0):

    rng = np.random.default_rng(seed=seed)

    def _resamp(selregr, selregd):
        for col in sample_columns:
            randoms[col] = np.zeros_like(data[col], shape=len(randoms))
        rand_sel = [selregr, ~selregr]
        dat_sel = [selregd, ~selregd]
        for dsel, rsel in zip(dat_sel, rand_sel):
            inds = rng.choice(len(data[dsel]), len(randoms[rsel]))
            # logger.info(str(len(data[dsel]),len(inds),np.max(inds))
            dshuf = data[dsel][inds]
            for col in sample_columns:
                randoms[col][rsel] = dshuf[col]

        rdl = []
        for dsel, rsel in zip(dat_sel, rand_sel):
            rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
            rdl.append(rd)
        rdr = rdl[0]/rdl[1]
        logger.info('norm factor is '+str(rdr))
        randoms['WEIGHT'][rand_sel[1]] *= rdr

    des_resamp = False
    if 'QSO' in tracer:
        des_resamp = True
    selregr = randoms['PHOTSYS'] == 'N'
    selregd = data['PHOTSYS'] == 'N'
    _resamp(selregr, selregd)
    rand_sel = [selregr, ~selregr]
    dat_sel = [selregd, ~selregd]

    for dsel, rsel in zip(dat_sel, rand_sel):
        rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
        logger.info('data/random weighted ratio after resampling:'+str(rd))

    if des_resamp:
        logger.info('resampling in DES region')
        from regressis import footprint
        import healpy as hp
        foot = footprint.DR9Footprint(
            256, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
        north, south, des = foot.get_imaging_surveys()
        th_ran, phi_ran = (-randoms['DEC']+90.) * \
            np.pi/180., randoms['RA']*np.pi/180.
        th_dat, phi_dat = (-data['DEC']+90.)*np.pi/180., data['RA']*np.pi/180.
        pixr = hp.ang2pix(256, th_ran, phi_ran, nest=True)
        selregr = des[pixr]
        pixd = hp.ang2pix(256, th_dat, phi_dat, nest=True)
        selregd = des[pixd]
        _resamp(selregr, selregd)
        rand_sel = [selregr, ~selregr]
        dat_sel = [selregd, ~selregd]

        for dsel, rsel in zip(dat_sel, rand_sel):
            rd = np.sum(randoms[rsel]['WEIGHT'])/np.sum(data[dsel]['WEIGHT'])
            logger.info('data/random weighted ratio after resampling:'+str(rd))

    return randoms


nproc = 18

mockdir = args.base_dir
if args.outloc == None:
    outdir = os.getenv(scratch)+'/holi-webjax/mock'+str(args.realization)+'/'
else:
    outdir = args.outloc


if not os.path.exists(outdir):
    os.makedirs(outdir)


logger.info('input directory is '+mockdir)
logger.info('output directory is '+outdir)


for tracer in tracers:
    tracerd = tracer

    out_data_fn = outdir+tracerd+'_input_clustering.dat'+args.outmd
    out_data_froot = outdir+tracerd+'_input_'
    ss = str(args.realization).zfill(4)
    mock_data_tr = common.read_hdf5_blosc(
        mockdir+'seed'+ss+'/holi_'+tracerd+'_v4.80_GCcomb_clustering.dat.h5', extname=None)

    if tracer == 'LRG':
        zmin = 0.4
        zmax = 1.1

    elif (tracer == 'ELG_LOP') or (tracer == 'ELG'):
        zmin = 0.8
        zmax = 1.6

    elif tracer == 'QSO':
        zmin = 0.8
        zmax = 3.5
    elif tracer == 'BGS_BRIGHT-21.5':
        zmin = 0.1
        zmax = 0.4
    if args.mkdat == 'y':
        if 'TARGETID' not in mock_data_tr.colnames:
            mock_data_tr['TARGETID'] = np.arange(len(mock_data_tr))
        selz = mock_data_tr['Z'] > zmin
        selz &= mock_data_tr['Z'] < zmax
        mock_data_tr = mock_data_tr[selz]
        logger.info('length after cutting to redshift range:' +
                    str(len(mock_data_tr)))
        mock_data_tr['WEIGHT_SYS'] = np.ones(len(mock_data_tr))
        mock_data_tr['WEIGHT_COMP'] = np.ones(len(mock_data_tr))
        mock_data_tr['WEIGHT_ZFAIL'] = np.ones(len(mock_data_tr))
        '''
        place to add imaging systematic weights and redshift failure weights would be here
        '''
        mock_data_tr['WEIGHT'] = mock_data_tr['WEIGHT_SYS'] * \
            mock_data_tr['WEIGHT_COMP']*mock_data_tr['WEIGHT_ZFAIL']
        if args.outmd == '.fits':
            common.write_LSS_scratchcp(
                mock_data_tr, out_data_fn, logger=logger)
        if args.outmd == '.h5':
            common.write_LSShdf5_scratchcp(
                mock_data_tr, out_data_fn, logger=logger)

        # splitGC(out_data_froot,'.dat')

    ran_samp_cols = ['Z', 'WEIGHT', 'WEIGHT_COMP',
                     'WEIGHT_SYS', 'WEIGHT_ZFAIL', 'TARGETID_DATA']

    nran = rx-rm
    tracerr = tracer
    if tracer[:3] == 'BGS':
        tracerr = 'BGS_BRIGHT'
    # ran_fname_base = args.base_dir.replace('global','dvs_ro') +tracerr+'_ffa_imaging_HPmapcut'
    ran_fname_base = args.data_dir.replace(
        'global', 'dvs_ro') + 'rands_intiles_DARK_nomask_'

    if args.mkran == 'y':
        if args.mkdat == 'n':
            # Table(fitsio.read(out_data_fn))
            mock_data_tr = read_file(out_data_fn)
        mock_data_tr.rename_column('TARGETID', 'TARGETID_DATA')

        def _mkran(rann):

            tracerr = tracer
            in_ran_fn = ran_fname_base+str(rann)+'.fits'
            out_ran_fn = out_data_froot+str(rann)+'_clustering.ran'+args.outmd
            rcols = ['RA', 'DEC']
            ran = Table(fitsio.read(in_ran_fn, columns=rcols))
            ran['TARGETID'] = np.arange(len(ran)) + (rann*10**9)

            ran = ran_col_assign(
                ran, mock_data_tr, ran_samp_cols, tracer, seed=rann)
            if args.outmd == '.fits':
                common.write_LSS_scratchcp(ran, out_ran_fn, logger=logger)
            if args.outmd == '.h5':
                common.write_LSShdf5_scratchcp(ran, out_ran_fn, logger=logger)
            del ran
            return True
            # splitGC(out_data_froot,'.ran',rann)

        inds = np.arange(nran)
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool(processes=nproc) as pool:
                res = pool.map(_mkran, inds)
        else:
            for rn in inds:  # range(rm,rx):
                _mkran(rn)

    if tracer == 'QSO':
        dz = 0.02
        P0 = 6000

    else:
        dz = 0.01

    if tracer == 'LRG':
        P0 = 10000
    if tracer[:3] == 'ELG':
        P0 = 4000
    if tracer[:3] == 'BGS':
        P0 = 7000

    regions = ['NGC', 'SGC']

    if args.nz == 'y':
        # this calculates the n(z) and then adds nbar(completeness) and FKP weights to the catalogs
        # for reg in allreg:
        fb = out_data_froot[:-1]
        fcr = fb+'_0_clustering.ran'+args.outmd
        fcd = fb+'_clustering.dat'+args.outmd
        fout = fb+'_nz.txt'
        common.mknz(fcd, fcr, fout, bs=dz, zmin=zmin, zmax=zmax, compmd='')
        common.addnbar(fb, bs=dz, zmin=zmin, zmax=zmax, P0=P0, nran=nran,
                       compmd='', par=args.par, nproc=nproc, logger=logger, exttp=args.outmd)

    if args.splitGC == 'y':
        splitGC(out_data_froot, '.dat')

        def _spran(rann):
            splitGC(out_data_froot, '.ran', rann)
        inds = np.arange(nran)
        if args.par == 'y':
            from multiprocessing import Pool
            with Pool(processes=nproc) as pool:
                res = pool.map(_spran, inds)
        else:
            for rn in inds:  # range(rm,rx):
                _spran(rn)
