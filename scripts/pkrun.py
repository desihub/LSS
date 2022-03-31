#!/usr/bin/env python
# coding: utf-8

# To run: srun -n 64 python pkrun.py --tracer ELG...

import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
from matplotlib import pyplot as plt

from pypower import CatalogFFTPower, PowerSpectrumStatistics, utils, setup_logging
from LSS.tabulated_cosmo import TabulatedDESI

from xirunpc import read_clustering_positions_weights, compute_angular_weights, catalog_dir, get_regions, get_zlims


os.environ['OMP_NUM_THREADS'] = os.environ['NUMEXPR_MAX_THREADS'] = '1'
logger = logging.getLogger('pkrun')


def compute_power_spectrum(edges, distance, dtype='f8', wang=None, weight_type='default', tracer='ELG', tracer2=None, rec_type=None, ells=(0, 2, 4), boxsize=5000., nmesh=1024, mpicomm=None, mpiroot=None, **kwargs):

    autocorr = tracer2 is None
    catalog_kwargs = kwargs.copy()
    catalog_kwargs['weight_type'] = weight_type
    with_shifted = rec_type is not None

    if 'angular' in weight_type and wang is None:

        wang = compute_angular_weights(nthreads=1, dtype=dtype, weight_type=weight_type, tracer=tracer, tracer2=tracer2, mpicomm=mpicomm, mpiroot=mpiroot, **kwargs)

    data_positions1, data_weights1, data_positions2, data_weights2 = None, None, None, None
    randoms_positions1, randoms_weights1, randoms_positions2, randoms_weights2 = None, None, None, None
    shifted_positions1, shifted_weights1, shifted_positions2, shifted_weights2 = None, None, None, None

    if mpicomm is None or mpicomm.rank == mpiroot:

        if with_shifted:
            data_positions1, data_weights1 = read_clustering_positions_weights(distance, name='data', rec_type=rec_type, tracer=tracer, **catalog_kwargs)
            shifted_positions1, shifted_weights1 = read_clustering_positions_weights(distance, name='randoms', rec_type=rec_type, tracer=tracer, **catalog_kwargs)
        else:
            data_positions1, data_weights1 = read_clustering_positions_weights(distance, name='data', rec_type=rec_type, tracer=tracer, **catalog_kwargs)
        randoms_positions1, randoms_weights1 = read_clustering_positions_weights(distance, name='randoms', rec_type=rec_type, tracer=tracer, **catalog_kwargs)

        if not autocorr:
            if with_shifted:
                data_positions2, data_weights2 = read_clustering_positions_weights(distance, name='data', rec_type=rec_type, tracer=tracer2, **catalog_kwargs)
                shifted_positions2, shifted_weights2 = read_clustering_positions_weights(distance, name='randoms', rec_type=rec_type, tracer=tracer2, **catalog_kwargs)
            else:
                data_positions2, data_weights2 = read_clustering_positions_weights(distance, name='data', rec_type=rec_type, tracer=tracer2, **catalog_kwargs)
            randoms_positions2, randoms_weights2 = read_clustering_positions_weights(distance, name='randoms', rec_type=rec_type, tracer=tracer2, **catalog_kwargs)

    kwargs = {}
    kwargs.update(wang or {})

    result = CatalogFFTPower(data_positions1=data_positions1, data_weights1=data_weights1,
                             data_positions2=data_positions2, data_weights2=data_weights2,
                             randoms_positions1=randoms_positions1, randoms_weights1=randoms_weights1,
                             randoms_positions2=randoms_positions2, randoms_weights2=randoms_weights2,
                             shifted_positions1=shifted_positions1, shifted_weights1=shifted_weights1,
                             shifted_positions2=shifted_positions2, shifted_weights2=shifted_weights2,
                             edges=edges, ells=ells, boxsize=boxsize, nmesh=nmesh, resampler='tsc', interlacing=2,
                             position_type='rdd', dtype=dtype, direct_limits=(0., 1.), direct_limit_type='degree',
                             **kwargs, mpicomm=mpicomm, mpiroot=mpiroot).poles

    return result, wang


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


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer(s) to be selected - 2 for cross-correlation', type=str, nargs='+', default=['ELG'])
    parser.add_argument('--basedir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs')
    parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main'], default='SV3')
    parser.add_argument('--verspec', help='version for redshifts', type=str, default='everest')
    parser.add_argument('--version', help='catalog version', type=str, default='test')
    parser.add_argument('--region', help='regions; by default, run on all regions', type=str, nargs='*', choices=['N', 'S', 'DN', 'DS', ''], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz", "fullonly"', type=str, nargs='*', default=None)
    parser.add_argument('--weight_type', help='types of weights to use; use "default_angular_bitwise" for PIP with angular upweighting; "default" just uses WEIGHT column', type=str, default='default')
    parser.add_argument('--boxsize', help='box size', type=float, default=5000.)
    parser.add_argument('--nmesh', help='mesh size', type=int, default=1024)
    parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=4)
    parser.add_argument('--outdir', help='base directory for output', type=str, default=None)
    parser.add_argument('--vis', help='show plot of each xi?', action='store_true', default=False)

    #only relevant for reconstruction
    parser.add_argument('--rec_type', help='reconstruction algorithm + reconstruction convention', choices=['IFTrecsym', 'IFTreciso', 'MGrecsym', 'MGreciso'], type=str, default=None)

    setup_logging()
    args = parser.parse_args()

    from pypower import mpi
    mpicomm = mpi.COMM_WORLD
    mpiroot = 0

    cat_dir = catalog_dir(base_dir=args.basedir, survey=args.survey, verspec=args.verspec, version=args.version)
    out_dir = os.path.join(os.environ['CSCRATCH'], args.survey)
    if args.outdir is not None: out_dir = args.outdir
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
        zlims = get_zlims(tracer, tracer2=tracer2, option=args.zlim[0])
    else:
        zlims = [float(zlim) for zlim in args.zlim]
    zlims = list(zip(zlims[:-1], zlims[1:])) + [(zlims[0], zlims[-1])]

    bin_type = 'lin'
    rebinning_factors = [1, 5, 10]
    if mpicomm is None or mpicomm.rank == mpiroot:
        logger.info('Computing power spectrum multipoles in regions {} in redshift ranges {}.'.format(regions, zlims))

    for zmin, zmax in zlims:
        for region in regions:
            if mpicomm is None or mpicomm.rank == mpiroot:
                logger.info('Computing power spectrum in region {} in redshift range {}.'.format(region, (zmin, zmax)))
            edges = get_edges()
            wang = None
            result, wang = compute_power_spectrum(edges=edges, distance=distance, nrandoms=args.nran, region=region, zlim=(zmin, zmax), weight_type=args.weight_type, boxsize=args.boxsize, nmesh=args.nmesh, wang=wang, mpicomm=mpicomm, mpiroot=mpiroot, **catalog_kwargs)
            file_kwargs = dict(region=region, tracer=tracer, tracer2=tracer2, zmin=zmin, zmax=zmax, rec_type=args.rec_type, weight_type=args.weight_type, bin_type=bin_type, out_dir=os.path.join(out_dir, 'pk'))
            fn = power_fn(file_type='npy', **file_kwargs)
            result.save(fn)
            txt_kwargs = file_kwargs.copy()
            for factor in rebinning_factors:
                #result = PowerSpectrumStatistics.load(fn)
                rebinned = result[:(result.shape[0]//factor)*factor:factor]
                txt_kwargs.update(bin_type=bin_type+str(factor))
                fn_txt = power_fn(file_type='pkpoles', **txt_kwargs)
                rebinned.save_txt(fn_txt)

                if args.vis and (mpicomm is None or mpicomm.rank == mpiroot):
                    k, poles = rebinned(return_k=True, complex=False)
                    for pole in poles: plt.plot(k, k*pole)
                    tracers = tracer
                    if tracer2 is not None: tracers += ' x ' + tracer2
                    plt.title('{} {:.2f} < z {:.2f} in {}'.format(tracers, zmin, zmax, region))
                    plt.show()
