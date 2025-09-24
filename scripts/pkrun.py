#!/usr/bin/env python
# coding: utf-8

"""
Note
----
The script can (should) be called with multiple processes as (e.g. on 2 nodes, 64 threads for each):
```
srun -n 128 python pkrun.py ...
```
Using -n 64 per node is recommended when applying rp-cut.
"""

# To run: srun -n 64 python pkrun.py --tracer ELG...

import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

from pypower import CatalogFFTPower, PowerSpectrumStatistics, CatalogSmoothWindow, PowerSpectrumSmoothWindow, PowerSpectrumOddWideAngleMatrix, PowerSpectrumSmoothWindowMatrix, utils, setup_logging
from LSS.tabulated_cosmo import TabulatedDESI

from xirunpc import compute_angular_weights
from LSS.cosmodesi_io_tools import read_clustering_positions_weights, concatenate_data_randoms, catalog_dir, get_regions, get_zlims, get_scratch_dir


os.environ['OMP_NUM_THREADS'] = os.environ['NUMEXPR_MAX_THREADS'] = '1'
logger = logging.getLogger('pkrun')


def barrier_idle(mpicomm, tag=0, sleep=0.01):
    """
    MPI barrier fonction that solves the problem that idle processes occupy 100% CPU.
    See: https://goo.gl/NofOO9.
    """
    import time
    size = mpicomm.size
    if size == 1: return
    rank = mpicomm.rank
    mask = 1
    while mask < size:
        dst = (rank + mask) % size
        src = (rank - mask + size) % size
        req = mpicomm.isend(None, dst, tag)
        while not mpicomm.Iprobe(src, tag):
            time.sleep(sleep)
        mpicomm.recv(None, src, tag)
        req.Wait()
        mask <<= 1


def compute_power_spectrum(edges, distance, dtype='f8', wang=None, weight_type='default', tracer='ELG', tracer2=None, recon_dir=None, rec_type=None, ells=(0, 2, 4), boxsize=5000., cellsize=6.,nmesh=None, dowin=False, mpicomm=None, mpiroot=0,P0=None, thetacut=None, **kwargs):

    autocorr = tracer2 is None
    catalog_kwargs = kwargs.copy()
    catalog_kwargs['weight_type'] = weight_type
    catalog_kwargs['concatenate'] = True
    with_shifted = rec_type is not None

    if 'angular' in weight_type and wang is None:
        #wang = compute_angular_weights(nthreads=1, dtype=dtype, weight_type=weight_type, tracer=tracer, tracer2=tracer2, mpicomm=mpicomm, mpiroot=mpiroot, **kwargs)
        # Does not run faster, why?
        # Because the number of cores is ncores // mpicomm.size
        nthreads = 64
        color = mpicomm.rank % nthreads == 0
        subcomm = mpicomm.Split(color, 0)
        if color:
            wang = compute_angular_weights(nthreads=nthreads, dtype=dtype, weight_type=weight_type, tracer=tracer, tracer2=tracer2, mpicomm=subcomm, mpiroot=0, **kwargs)
        barrier_idle(mpicomm)
        wang = mpicomm.bcast(wang, root=0)

    data_positions1, data_weights1, data_positions2, data_weights2 = None, None, None, None
    randoms_positions1, randoms_weights1, randoms_positions2, randoms_weights2 = None, None, None, None
    shifted_positions1, shifted_weights1, shifted_positions2, shifted_weights2 = None, None, None, None

    if mpicomm.rank == mpiroot:
        logger.info('compute power P0 is '+str(P0))
        data, randoms = read_clustering_positions_weights(distance, name=['data', 'randoms'], recon_dir=recon_dir, rec_type=rec_type, tracer=tracer, option=option,P0=P0, **catalog_kwargs)
        if with_shifted:
            shifted = randoms  # above returned shifted randoms
            randoms = read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer, option=option,P0=P0, **catalog_kwargs)
        (data_positions1, data_weights1), (randoms_positions1, randoms_weights1) = concatenate_data_randoms(data, randoms, **catalog_kwargs)
        if with_shifted:
            shifted_positions1, shifted_weights1 = concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]

        if not autocorr:
            data, randoms = read_clustering_positions_weights(distance, name=['data', 'randoms'], rec_type=rec_type, tracer=tracer2, option=option,P0=P0, **catalog_kwargs)
            if with_shifted:
                shifted = randoms
                randoms = read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer2, option=option,P0=P0, **catalog_kwargs)
            (data_positions2, data_weights2), (randoms_positions2, randoms_weights2) = concatenate_data_randoms(data, randoms, **catalog_kwargs)
            if with_shifted:
                shifted_positions2, shifted_weights2 = concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]

    kwargs = {}
    kwargs.update(wang or {})
    direct_edges = win_direct_selection_attrs = direct_selection_attrs = None
    if wang or mpicomm.bcast(len(data_weights1) > 1 if mpicomm.rank == mpiroot else None):
        direct_selection_attrs = {'theta': (0., 1.)}
    elif thetacut is not None:
        win_direct_selection_attrs = direct_selection_attrs = {'theta': (0., thetacut)}
        #direct_edges = {'min': 0., 'step': 0.1, 'max': 200.}  # use this to reduce the computing time for direct pair counts to a few seconds
        direct_edges = {'min': 0., 'step': 0.1}

    result = CatalogFFTPower(data_positions1=data_positions1, data_weights1=data_weights1,
                             data_positions2=data_positions2, data_weights2=data_weights2,
                             randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,
                             randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2,
                             shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1,
                             shifted_positions2=shifted_positions2, shifted_weights2=shifted_weights2,
                             edges=edges, ells=ells, boxsize=boxsize, cellsize=cellsize,nmesh=nmesh, resampler='tsc', interlacing=3,
                             position_type='rdd', dtype=dtype, direct_selection_attrs=direct_selection_attrs, direct_edges=direct_edges,
                             direct_attrs={"nthreads": 64}, **kwargs, mpicomm=mpicomm, mpiroot=mpiroot).poles
    window = wmatrix = None
    if dowin:
        windows = []
        boxsizes = [scale * boxsize for scale in [20., 5., 1.]]
        edges = {'step': 2. * np.pi / boxsizes[0]}
        for iboxsize, boxsize in enumerate(boxsizes):
            windows.append(CatalogSmoothWindow(randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,
                                               randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2,
                                               power_ref=result, edges=edges, boxsize=boxsize, position_type='rdd',
                                               direct_selection_attrs=win_direct_selection_attrs if iboxsize == 0 else None,
                                               direct_edges=direct_edges if iboxsize == 0 else None,
                                               mpicomm=mpicomm, mpiroot=mpiroot).poles)

        if mpicomm.rank == mpiroot:
            window = PowerSpectrumSmoothWindow.concatenate_x(*windows, frac_nyq=0.9)
            # Let us compute the wide-angle and window function matrix
            ellsin = [0, 2, 4] # input (theory) multipoles
            wa_orders = 1 # wide-angle order
            sep = np.geomspace(1e-4, 2e4, 1024*16) # configuration space separation for FFTlog, 2e4 > sqrt(3) * 8000
            kin_rebin = 4 # rebin input theory to save memory
            kin_lim = (0, 2e1) # pre-cut input (theory) ks to save some memory
            # Input projections for window function matrix:
            # theory multipoles at wa_order = 0, and wide-angle terms at wa_order = 1
            projsin = ellsin + PowerSpectrumOddWideAngleMatrix.propose_out(ellsin, wa_orders=wa_orders)
            # Window matrix
            wmatrix = PowerSpectrumSmoothWindowMatrix(result, projsin=projsin, window=window, sep=sep, kin_rebin=kin_rebin, kin_lim=kin_lim)
            # We resum over theory odd-wide angle
            wmatrix.resum_input_odd_wide_angle()

    return result, wang, window, wmatrix


def get_edges():
    return {'min': 0., 'step': 0.001}


def power_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, recon_dir='n', rec_type=False, weight_type='default', bin_type='lin',option=None, thetacut=None, out_dir='.',P0=None):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    if recon_dir != 'n':
        out_dir = out_dir[:-2] + recon_dir+'/pk/'
    if option:
        zmax = str(zmax) + option

    root = '{}_{}_{}_{}_{}'.format(tracer, zmin, zmax, weight_type, bin_type)
    if thetacut is not None:
        root += '_thetacut{}'.format(thetacut)
    if P0 is not None:
        root += '_P0{}'.format(P0)

    if file_type == 'npy':
        return os.path.join(out_dir, 'pkpoles_{}.npy'.format(root))
    return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


def window_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, recon_dir='n', rec_type=False, weight_type='default', bin_type='lin',option=None, thetacut=None, out_dir='.',P0=None):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    if recon_dir != 'n':
        out_dir = out_dir[:-2] + recon_dir+'/pk/'
    if option:
        zmax = str(zmax) + option

    root = '{}_{}_{}_{}_{}'.format(tracer, zmin, zmax, weight_type, bin_type)
    if thetacut is not None:
        root += '_thetacut{}'.format(thetacut)
    if P0 is not None:
        root += '_P0{}'.format(P0)

    if file_type == 'npy':
        return os.path.join(out_dir, 'window_smooth_{}.npy'.format(root))
    return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


def wmatrix_fn(region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, recon_dir='n', rec_type=False, weight_type='default', bin_type='lin',option=None, thetacut=None, out_dir='.',P0=None):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    if recon_dir != 'n':
        out_dir = out_dir[:-2] + recon_dir+'/pk/'
    if option:
        zmax = str(zmax) + option

    root = '{}_{}_{}_{}_{}'.format(tracer, zmin, zmax, weight_type, bin_type)
    if thetacut is not None:
        root += '_thetacut{}'.format(thetacut)
    if P0 is not None:
        root += '_P0{}'.format(P0)
    return os.path.join(out_dir, 'wmatrix_smooth_{}.npy'.format(root))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer(s) to be selected - 2 for cross-correlation', type=str, nargs='+', default=['ELG'])
    parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
    parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main', 'Y1','DA2'], default='Y1')
    parser.add_argument('--verspec', help='version for redshifts', type=str, default='iron')
    parser.add_argument('--version', help='catalog version', type=str, default='test')
    parser.add_argument('--ran_sw', help='extra string in random name', type=str, default='')
    parser.add_argument('--region', help='regions; by default, run on N, S; pass NS to run on concatenated N + S', type=str, nargs='*', choices=['N', 'S', 'NS', 'NGC', 'SGC', 'ALL'], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz", "fullonly"', type=str, nargs='*', default=None)
    parser.add_argument('--option', help='place to put extra options for cutting catalogs', default=None)
    parser.add_argument('--weight_type', help='types of weights to use; use "default_angular_bitwise" for PIP with angular upweighting; "default" just uses WEIGHT column', type=str, default='default')
    parser.add_argument('--boxsize', help='box size', type=float, default=8000.)
    parser.add_argument('--nmesh', help='mesh size', default=None)
    parser.add_argument('--cellsize', help='cell size', default=6.)
    parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=4)
    parser.add_argument('--outdir', help='base directory for output (default: SCRATCH)', type=str, default=None)
    parser.add_argument('--calc_win', help='also calculate window?; use "y" for yes', default='n')
    parser.add_argument('--vis', help='show plot of each pk?', action='store_true', default=False)
    parser.add_argument('--rebinning', help='whether to rebin the pk or just keep the original .npy file', default='n')
    parser.add_argument('--thetacut', help='apply this theta-cut, standard 0.05', type=float, default=None)
    parser.add_argument('--P0', help='value of P0 to use in FKP weights (None defaults to WEIGHT_FKP)', type=float, default=None)
    parser.add_argument('--recon_dir', help='if recon catalogs are in a subdirectory, put that here', type=str, default='n')
    #only relevant for reconstruction
    parser.add_argument('--rec_type', help='reconstruction algorithm + reconstruction convention', choices=['IFTrecsym', 'IFTreciso', 'MGrecsym', 'MGreciso'], type=str, default=None)

    setup_logging()
    args = parser.parse_args()
    if args.nmesh != None:
        args.nmesh = int(args.nmesh)
    if args.cellsize == 'None':
        args.cellsize = None    
    if args.cellsize != None:
        args.cellsize = float(args.cellsize)
    if args.calc_win == 'n':
        args.calc_win = False
    if args.calc_win == 'y':
        args.calc_win = True
        
    if args.rebinning == 'n':
        args.rebinning = False
    if args.rebinning == 'y':
        args.rebinning = True
    P0 = args.P0
    if P0 is not None:
        P0 = float(P0)
    from pypower import mpi
    mpicomm = mpi.COMM_WORLD
    mpiroot = 0

    if os.path.normpath(args.basedir) == os.path.normpath('/global/cfs/cdirs/desi/survey/catalogs/'):
        cat_dir = catalog_dir(base_dir=args.basedir, survey=args.survey, verspec=args.verspec, version=args.version)
    elif os.path.normpath(args.basedir) == os.path.normpath('/global/project/projectdirs/desi/users/acarnero/mtl_mock000_univ1/'):
        cat_dir = args.basedir
        args.region = ['']
    else:
        cat_dir = args.basedir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Catalog directory is {}.'.format(cat_dir))

    if args.outdir is None:
        out_dir = os.path.join(get_scratch_dir(), args.survey,args.version)
    else:
        out_dir = args.outdir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Output directory is {}.'.format(out_dir))
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('using P0 '+str(P0))

    tracer, tracer2 = args.tracer[0], None
    if len(args.tracer) > 1:
        tracer2 = args.tracer[1]
        if len(args.tracer) > 2:
            raise ValueError('Provide <= 2 tracers!')
    if tracer2 == tracer:
        tracer2 = None # otherwise counting of self-pairs
    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, survey=args.survey, cat_dir=cat_dir, recon_dir=args.recon_dir, rec_type=args.rec_type, ran_sw=args.ran_sw) # survey required for zdone
    distance = TabulatedDESI().comoving_radial_distance

    regions = args.region
    if regions is None:
        regions = get_regions(args.survey, rec=bool(args.rec_type))


    option = args.option#
    if args.zlim is None:
        zlims = get_zlims(tracer, tracer2=tracer2)
    elif not args.zlim[0].replace('.', '').isdigit():
        if option is not None:
            sys.exit('conflicting options, need to fix code if both are needed')
        option = args.zlim[0]
        zlims = get_zlims(tracer, tracer2=tracer2, option=option)
    else:
        zlims = [float(zlim) for zlim in args.zlim]
    zlims = list(zip(zlims[:-1], zlims[1:])) 
    #zlims = list(zip(zlims[:-1], zlims[1:])) + ([(zlims[0], zlims[-1])] if len(zlims) > 2 else []) # len(zlims) == 2 == single redshift range

    bin_type = 'lin'
    rebinning_factors = [1, 5, 10]
    if mpicomm.rank == mpiroot:
        logger.info('Computing power spectrum multipoles in regions {} in redshift ranges {}.'.format(regions, zlims))

    for zmin, zmax in zlims:
        base_file_kwargs = dict(tracer=tracer, tracer2=tracer2, zmin=zmin, zmax=zmax, recon_dir=args.recon_dir, rec_type=args.rec_type, weight_type=args.weight_type, bin_type=bin_type, out_dir=os.path.join(out_dir, 'pk'), thetacut=args.thetacut,option=option,P0=P0)
        for region in regions:
            if mpicomm.rank == mpiroot:
                logger.info('Computing power spectrum in region {} in redshift range {}.'.format(region, (zmin, zmax)))
            edges = get_edges()
            wang = None
            result, wang, window, wmatrix = compute_power_spectrum(edges=edges, distance=distance, nrandoms=args.nran, region=region, zlim=(zmin, zmax), weight_type=args.weight_type, boxsize=args.boxsize, cellsize=args.cellsize,nmesh=args.nmesh, wang=wang, dowin=args.calc_win, thetacut=args.thetacut, P0=P0, mpicomm=mpicomm, mpiroot=mpiroot, **catalog_kwargs)
            if mpicomm.rank == mpiroot:
                fn = power_fn(file_type='npy', region=region, **base_file_kwargs)
                result.save(fn)
                if window is not None:
                    fn = window_fn(file_type='npy', region=region, **base_file_kwargs)
                    window.save(fn)
                    fn = wmatrix_fn(region=region, **base_file_kwargs)
                    wmatrix.save(fn)

        all_regions = regions.copy()
        if mpicomm.rank == mpiroot:
            if 'N' in regions and 'S' in regions:  # let's combine
                result = sum([PowerSpectrumStatistics.load(power_fn(file_type='npy', region=region, **base_file_kwargs)) for region in ['N', 'S']])
                result.save(power_fn(file_type='npy', region='NScomb', **base_file_kwargs))
                if args.calc_win:
                    result = sum([PowerSpectrumSmoothWindow.load(window_fn(file_type='npy', region=region, **base_file_kwargs)) for region in ['N', 'S']])
                    result.save(window_fn(file_type='npy', region='NScomb', **base_file_kwargs))
                    result = sum([PowerSpectrumSmoothWindowMatrix.load(wmatrix_fn(region=region, **base_file_kwargs)) for region in ['N', 'S']])
                    result.save(wmatrix_fn(region='NScomb', **base_file_kwargs))
                all_regions.append('NScomb')
            if 'NGC' in regions and 'SGC' in regions:  # let's combine
                result = sum([PowerSpectrumStatistics.load(power_fn(file_type='npy', region=region, **base_file_kwargs)) for region in ['NGC', 'SGC']])
                result.save(power_fn(file_type='npy', region='GCcomb', **base_file_kwargs))
                if args.calc_win:
                    result = sum([PowerSpectrumSmoothWindow.load(window_fn(file_type='npy', region=region, **base_file_kwargs)) for region in ['NGC', 'SGC']])
                    result.save(window_fn(file_type='npy', region='GCcomb', **base_file_kwargs))
                    result = sum([PowerSpectrumSmoothWindowMatrix.load(wmatrix_fn(region=region, **base_file_kwargs)) for region in ['NGC', 'SGC']])
                    result.save(wmatrix_fn(region='GCcomb', **base_file_kwargs))
                all_regions.append('GCcomb')
            if args.rebinning:
                for region in all_regions:
                    txt_kwargs = base_file_kwargs.copy()
                    txt_kwargs.update(region=region)
                    result = PowerSpectrumStatistics.load(power_fn(file_type='npy', **txt_kwargs))
                    for factor in rebinning_factors:
                        #result = PowerSpectrumStatistics.load(fn)
                        rebinned = result[:(result.shape[0]//factor)*factor:factor]
                        txt_kwargs.update(bin_type=bin_type+str(factor))
                        fn_txt = power_fn(file_type='pkpoles', **txt_kwargs)
                        rebinned.save_txt(fn_txt)

                        if args.vis:
                            k, poles = rebinned(return_k=True, complex=False)
                            for pole in poles: plt.plot(k, k*pole)
                            tracers = tracer
                            if tracer2 is not None: tracers += ' x ' + tracer2
                            plt.title('{} {:.2f} < z {:.2f} in {}'.format(tracers, zmin, zmax, region))
                            plt.show()
