#!/usr/bin/env python
# coding: utf-8

# To run: srun -n 64 python pkrun.py --tracer ELG...

import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

from pypower import CatalogFFTPower, PowerSpectrumStatistics, CatalogSmoothWindow, utils, setup_logging
from LSS.tabulated_cosmo import TabulatedDESI

from xirunpc import read_clustering_positions_weights, concatenate_data_randoms, compute_angular_weights, catalog_dir, get_regions, get_zlims, get_scratch_dir


os.environ['OMP_NUM_THREADS'] = os.environ['NUMEXPR_MAX_THREADS'] = '1'
logger = logging.getLogger('pkrun')


def compute_power_spectrum(edges, distance, dtype='f8', wang=None, weight_type='default', tracer='ELG', tracer2=None, rec_type=None, ells=(0, 2, 4), boxsize=5000., nmesh=1024, dowin=False, option=None, mpicomm=None, mpiroot=None, **kwargs):

    autocorr = tracer2 is None
    catalog_kwargs = kwargs.copy()
    catalog_kwargs['weight_type'] = weight_type
    catalog_kwargs['concatenate'] = True
    with_shifted = rec_type is not None

    if 'angular' in weight_type and wang is None:

        wang = compute_angular_weights(nthreads=1, dtype=dtype, weight_type=weight_type, tracer=tracer, tracer2=tracer2, mpicomm=mpicomm, mpiroot=mpiroot, **kwargs)

    data_positions1, data_weights1, data_positions2, data_weights2 = None, None, None, None
    randoms_positions1, randoms_weights1, randoms_positions2, randoms_weights2 = None, None, None, None
    shifted_positions1, shifted_weights1, shifted_positions2, shifted_weights2 = None, None, None, None

    if mpicomm is None or mpicomm.rank == mpiroot:

        data, randoms = read_clustering_positions_weights(distance, name=['data', 'randoms'], rec_type=rec_type, tracer=tracer, option=option, **catalog_kwargs)
        if with_shifted:
            shifted = randoms  # above returned shifted randoms
            randoms = read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer, option=option, **catalog_kwargs)
        (data_positions1, data_weights1), (randoms_positions1, randoms_weights1) = concatenate_data_randoms(data, randoms, **catalog_kwargs)
        if with_shifted:
            shifted_positions1, shifted_weights1 = concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]

        if not autocorr:
            data, randoms = read_clustering_positions_weights(distance, name=['data', 'randoms'], rec_type=rec_type, tracer=tracer2, option=option, **catalog_kwargs)
            if with_shifted:
                shifted = randoms
                randoms = read_clustering_positions_weights(distance, name='randoms', rec_type=False, tracer=tracer2, option=option, **catalog_kwargs)
            (data_positions2, data_weights2), (randoms_positions2, randoms_weights2) = concatenate_data_randoms(data, randoms, **catalog_kwargs)
            if with_shifted:
                shifted_positions2, shifted_weights2 = concatenate_data_randoms(data, shifted, **catalog_kwargs)[1]

    kwargs = {}
    kwargs.update(wang or {})

    result = CatalogFFTPower(data_positions1=data_positions1, data_weights1=data_weights1,
                             data_positions2=data_positions2, data_weights2=data_weights2,
                             randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,
                             randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2,
                             shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1,
                             shifted_positions2=shifted_positions2, shifted_weights2=shifted_weights2,
                             edges=edges, ells=ells, boxsize=boxsize, nmesh=nmesh, resampler='tsc', interlacing=2,
                             position_type='rdd', dtype=dtype, direct_limits=(0., 1.), direct_limit_type='degree', # direct_limits, (0, 1) degree
                             **kwargs, mpicomm=mpicomm, mpiroot=mpiroot).poles
    window = None
    if dowin:
        window = CatalogSmoothWindow(randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,
                                     power_ref=result, edges=edges, boxsize=boxsize, position_type='rdd',
                                     **kwargs, mpicomm=mpicomm, mpiroot=mpiroot).poles
    return result, wang, window


def get_edges():
    return {'min':0., 'step':0.001}


def power_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, rec_type=False, weight_type='default', bin_type='lin', out_dir='.'):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    root = '{}_{}_{}_{}_{}'.format(tracer, zmin, zmax, weight_type, bin_type)
    if file_type == 'npy':
        return os.path.join(out_dir, 'pkpoles_{}.npy'.format(root))
    return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


def window_fn(file_type='npy', region='', tracer='ELG', tracer2=None, zmin=0, zmax=np.inf, rec_type=False, weight_type='default', bin_type='lin', out_dir='.'):
    if tracer2: tracer += '_' + tracer2
    if rec_type: tracer += '_' + rec_type
    if region: tracer += '_' + region
    root = '{}_{}_{}_{}_{}'.format(tracer, zmin, zmax, weight_type, bin_type)
    if file_type == 'npy':
        return os.path.join(out_dir, 'window_smooth_{}.npy'.format(root))
    return os.path.join(out_dir, '{}_{}.txt'.format(file_type, root))


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer(s) to be selected - 2 for cross-correlation', type=str, nargs='+', default=['ELG'])
    parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
    parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main'], default='SV3')
    parser.add_argument('--verspec', help='version for redshifts', type=str, default='guadalupe')
    parser.add_argument('--version', help='catalog version', type=str, default='test')
    parser.add_argument('--region', help='regions; by default, run on N, S; pass NS to run on concatenated N + S', type=str, nargs='*', choices=['N', 'S', 'NS'], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz", "fullonly"', type=str, nargs='*', default=None)
    parser.add_argument('--weight_type', help='types of weights to use; use "default_angular_bitwise" for PIP with angular upweighting; "default" just uses WEIGHT column', type=str, default='default')
    parser.add_argument('--boxsize', help='box size', type=float, default=8000.)
    parser.add_argument('--nmesh', help='mesh size', type=int, default=1024)
    parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=4)
    parser.add_argument('--outdir', help='base directory for output (default: SCRATCH)', type=str, default=None)
    parser.add_argument('--calc_win', help='also calculate window?; use "y" for yes', action='store_true', default='n')
    parser.add_argument('--vis', help='show plot of each pk?', action='store_true', default=False)

    #only relevant for reconstruction
    parser.add_argument('--rec_type', help='reconstruction algorithm + reconstruction convention', choices=['IFTrecsym', 'IFTreciso', 'MGrecsym', 'MGreciso'], type=str, default=None)

    setup_logging()
    args = parser.parse_args()
    if args.calc_win == 'n':
        args.calc_win = False
    if args.calc_win == 'y':
        args.calc_win = True

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
        out_dir = os.path.join(get_scratch_dir(), args.survey)
    else:
        out_dir = args.outdir
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Output directory is {}.'.format(out_dir))

    tracer, tracer2 = args.tracer[0], None
    if len(args.tracer) > 1:
        tracer2 = args.tracer[1]
        if len(args.tracer) > 2:
            raise ValueError('Provide <= 2 tracers!')
    if tracer2 == tracer:
        tracer2 = None # otherwise counting of self-pairs
    catalog_kwargs = dict(tracer=tracer, tracer2=tracer2, survey=args.survey, cat_dir=cat_dir, rec_type=args.rec_type) # survey required for zdone
    distance = TabulatedDESI().comoving_radial_distance

    regions = args.region
    if regions is None:
        regions = get_regions(args.survey, rec=bool(args.rec_type))

    if args.zlim is None:
        zlims = get_zlims(tracer, tracer2=tracer2)
    elif not args.zlim[0].replace('.', '').isdigit():
        option = args.zlim[0]
        zlims = get_zlims(tracer, tracer2=tracer2, option=option)
    else:
        zlims = [float(zlim) for zlim in args.zlim]
    zlims = list(zip(zlims[:-1], zlims[1:])) + ([(zlims[0], zlims[-1])] if len(zlims) > 2 else []) # len(zlims) == 2 == single redshift range

    bin_type = 'lin'
    rebinning_factors = [1, 5, 10]
    if mpicomm.rank == mpiroot:
        logger.info('Computing power spectrum multipoles in regions {} in redshift ranges {}.'.format(regions, zlims))

    for zmin, zmax in zlims:
        base_file_kwargs = dict(tracer=tracer, tracer2=tracer2, zmin=zmin, zmax=zmax, rec_type=args.rec_type, weight_type=args.weight_type, bin_type=bin_type, out_dir=os.path.join(out_dir, 'pk'))
        for region in regions:
            if mpicomm.rank == mpiroot:
                logger.info('Computing power spectrum in region {} in redshift range {}.'.format(region, (zmin, zmax)))
            edges = get_edges()
            wang = None
            result, wang, window = compute_power_spectrum(edges=edges, distance=distance, nrandoms=args.nran, region=region, zlim=(zmin, zmax), weight_type=args.weight_type, boxsize=args.boxsize, nmesh=args.nmesh, wang=wang, dowin=args.calc_win,mpicomm=mpicomm, mpiroot=mpiroot, **catalog_kwargs)
            fn = power_fn(file_type='npy', region=region, **base_file_kwargs)
            result.save(fn)
            if window is not None:
                fn = window_fn(file_type='npy', region=region, **base_file_kwargs)
                window.save(fn)

        all_regions = regions.copy()
        if mpicomm.rank == mpiroot:
            if 'N' in regions and 'S' in regions:  # let's combine
                result = sum([PowerSpectrumStatistics.load(power_fn(file_type='npy', region=region, **base_file_kwargs)) for region in ['N', 'S']])
                result.save(power_fn(file_type='npy', region='NScomb', **base_file_kwargs))
                all_regions.append('NScomb')
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
