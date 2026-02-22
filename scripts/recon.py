#!/usr/bin/env python
# coding: utf-8

"""
Note
----
The script can be used either with pyrecon/main (non-MPI, but OpenMP-thread version) or pyrecon/mpi (MPI version).
To use the MPI version:
```
module swap pyrecon/main pyrecon/mpi
```
In this case, the script can (should) be called with multiple processes as:
```
srun -n 32 python recon.py ...
```
As other scripts, xirunpc.py, pkrun.py, the non-MPI astropy's Table() is used for catalogs...
But in some future, it would be much better to use everywhere the MPI version of mockfactory (mpytools):
```
from mockfactory import Catalog
catalog = Catalog.read(data_fn)
catalog.write(data_fn)
```
"""

import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack

import pyrecon
from pyrecon import MultiGridReconstruction, IterativeFFTReconstruction, IterativeFFTParticleReconstruction, utils, setup_logging
from LSS.tabulated_cosmo import TabulatedDESI
from LSS.cosmodesi_io_tools import get_clustering_positions_weights, catalog_dir, catalog_fn, get_regions, get_zlims, get_scratch_dir

logger = logging.getLogger('recon')


def run_reconstruction(Reconstruction, distance, data_fn, randoms_fn, data_rec_fn, randoms_rec_fn, f=0.8, bias=1.2, boxsize=None, nmesh=None, cellsize=7, smoothing_radius=15, nthreads=64, convention='reciso', dtype='f8', mpicomm=None, **kwargs):

    root = mpicomm is None or mpicomm.rank == 0

    if np.ndim(randoms_fn) == 0: randoms_fn = [randoms_fn]
    if np.ndim(randoms_rec_fn) == 0: randoms_rec_fn = [randoms_rec_fn]
    
    data_positions, data_weights = None, None
    randoms_positions, randoms_weights = None, None

    if root:
        logger.info('Loading {}.'.format(data_fn))
        data = Table.read(data_fn)
        (ra, dec, dist), data_weights, mask = get_clustering_positions_weights(data, distance, name='data', return_mask=True, **kwargs)
        data = data[mask]
        data_positions = utils.sky_to_cartesian(dist, ra, dec, dtype=dtype)

    if mpicomm is not None:
        rec_kwargs = {'mpicomm': mpicomm, 'mpiroot': 0}
    else:
        rec_kwargs = {'fft_engine': 'fftw', 'nthreads': nthreads}
    recon = Reconstruction(f=f, bias=bias, boxsize=boxsize, nmesh=nmesh, cellsize=cellsize, los='local', positions=data_positions, boxpad=1.2, dtype=dtype, **rec_kwargs)

    recon.assign_data(data_positions, data_weights)
    for fn in randoms_fn:
        if root:
            logger.info('Loading {}.'.format(fn))
            (ra, dec, dist), randoms_weights = get_clustering_positions_weights(Table.read(fn), distance, name='randoms', **kwargs)
            randoms_positions = utils.sky_to_cartesian(dist, ra, dec, dtype=dtype)
        recon.assign_randoms(randoms_positions, randoms_weights)

    recon.set_density_contrast(smoothing_radius=smoothing_radius)
    recon.run()

    field = 'rsd' if convention == 'rsd' else 'disp+rsd'
    if type(recon) is IterativeFFTParticleReconstruction:
        data_positions_rec = recon.read_shifted_positions('data', field=field)
    else:
        data_positions_rec = recon.read_shifted_positions(data_positions, field=field)

    distance_to_redshift = utils.DistanceToRedshift(distance)
    if root:
        catalog = Table(data)
        dist, ra, dec = utils.cartesian_to_sky(data_positions_rec)
        catalog['RA'], catalog['DEC'], catalog['Z'] = ra, dec, distance_to_redshift(dist)
        logger.info('Saving {}.'.format(data_rec_fn))
        utils.mkdir(os.path.dirname(data_rec_fn))
        catalog.write(data_rec_fn, format='fits', overwrite=True)

    if convention != 'rsd':
        field = 'disp+rsd' if convention == 'recsym' else 'disp'
        for fn, rec_fn in zip(randoms_fn, randoms_rec_fn):
            if root:
                catalog = Table.read(fn)
                (ra, dec, dist), randoms_weights, mask = get_clustering_positions_weights(catalog, distance, name='randoms', return_mask=True, **kwargs)
                catalog = catalog[mask]
                randoms_positions = utils.sky_to_cartesian(dist, ra, dec, dtype=dtype)
            randoms_positions_rec = recon.read_shifted_positions(randoms_positions, field=field)
            if root:
                dist, ra, dec = utils.cartesian_to_sky(randoms_positions_rec)
                catalog['RA'], catalog['DEC'], catalog['Z'] = ra, dec, distance_to_redshift(dist)
                logger.info('Saving {}.'.format(rec_fn))
                utils.mkdir(os.path.dirname(rec_fn))
                catalog.write(rec_fn, format='fits', overwrite=True)
        
        
def get_f_bias(tracer='ELG'):
    if tracer.startswith('ELG'):
        return 0.9, 1.2
    if tracer.startswith('QSO'):
        return 0.928, 2.07
    if tracer.startswith('LRG'):
        return 0.834, 1.99
    if tracer.startswith('BGS'):
        return 0.682, 1.5

    return 0.9, 1.2


if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('--tracer', help='tracer to be selected', type=str, default='ELG')
    parser.add_argument('--indir', help='where to find catalogs', type=str, default='/global/cfs/cdirs/desi/survey/catalogs/')
    parser.add_argument('--survey', help='e.g., SV3 or main', type=str, choices=['SV3', 'DA02', 'main'], default='DA02')
    parser.add_argument('--verspec', help='version for redshifts', type=str, default='guadalupe')
    parser.add_argument('--version', help='catalog version', type=str, default='test')
    parser.add_argument('--region', help='regions; by default, run on all regions', type=str, nargs='*', choices=['NGC','SGC','N', 'S', 'DN', 'DS', ''], default=None)
    parser.add_argument('--zlim', help='z-limits, or options for z-limits, e.g. "highz", "lowz"', type=str, nargs='*', default=None)
    parser.add_argument('--weight_type', help='types of weights to use; "default" just uses WEIGHT column', type=str, default='default')
    parser.add_argument('--nran', help='number of random files to combine together (1-18 available)', type=int, default=5)
    parser.add_argument('--nthreads', help='number of threads, not used in case pyrecon/mpi is used (do srun -n 64 python recon.py ... instead)', type=int, default=64)
    parser.add_argument('--outdir',  help='base directory for output (default: SCRATCH)', type=str, default=None)
    parser.add_argument('--algorithm', help='reconstruction algorithm', type=str, choices=['MG', 'IFT', 'IFTP'], default='MG')
    parser.add_argument('--convention', help='reconstruction convention', type=str, choices=['reciso', 'recsym'], default='reciso')
    parser.add_argument('--f', help='growth rate', type=float, default=None)
    parser.add_argument('--bias', help='bias', type=float, default=None)
    parser.add_argument('--boxsize', help='box size', type=float, default=None)
    parser.add_argument('--nmesh', help='mesh size', type=int, default=None)
    parser.add_argument('--cellsize', help='cell size', type=float, default=7)
    parser.add_argument('--smoothing_radius', help='smoothing radius', type=float, default=15)
    parser.add_argument('--prepare_blinding', help='Use this flag to create a realspace catalog, that can be used as innput for RSD blinding', type=bool, default=False)

    args = parser.parse_args()

    try:
        mpicomm = pyrecon.mpi.COMM_WORLD  # MPI version
    except AttributeError:
        mpicomm = None  # non-MPI version
    root = mpicomm is None or mpicomm.rank == 0
    setup_logging()

    Reconstruction = {'MG': MultiGridReconstruction, 'IFT': IterativeFFTReconstruction, 'IFTP': IterativeFFTParticleReconstruction}[args.algorithm]

    if os.path.normpath(args.indir) == os.path.normpath('/global/cfs/cdirs/desi/survey/catalogs/'):
        cat_dir = catalog_dir(base_dir=args.indir, survey=args.survey, verspec=args.verspec, version=args.version)
    elif os.path.normpath(args.indir) == os.path.normpath('/global/project/projectdirs/desi/users/acarnero/mtl_mock000_univ1/'):
        cat_dir = args.indir
        args.region = ['']
    else:
        cat_dir = args.indir
    if root: logger.info('Input directory is {}.'.format(cat_dir))

    if args.outdir is None:
        out_dir = os.path.join(get_scratch_dir(), args.survey)
    else:
        out_dir = args.outdir
    if root: logger.info('Output directory is {}.'.format(out_dir))

    distance = TabulatedDESI().comoving_radial_distance

    f, bias = get_f_bias(args.tracer)
    if args.f is not None: f = args.f
    if args.bias is not None: bias = args.bias

    regions = args.region
    if regions is None:
        regions = get_regions(args.survey, rec=True)

    if args.zlim is None:
        zlims = get_zlims(args.tracer)
    elif not args.zlim[0].replace('.', '').isdigit():
        zlims = get_zlims(args.tracer, option=args.zlim[0])
    else:
        zlims = [float(zlim) for zlim in args.zlim]
    zlims = [(zlims[0], zlims[-1])]

    for zmin, zmax in zlims:
        for region in regions:
            if root: logger.info('Running reconstruction in region {} in redshift range {} with f, bias = {}.'.format(region, (zmin, zmax), (f, bias)))
            catalog_kwargs = dict(tracer=args.tracer, region=region, ctype='clustering', nrandoms=args.nran, survey=args.survey)
            data_fn = catalog_fn(**catalog_kwargs, cat_dir=cat_dir, name='data')
            randoms_fn = catalog_fn(**catalog_kwargs, cat_dir=cat_dir, name='randoms')
            data_rec_fn = catalog_fn(**catalog_kwargs, cat_dir=out_dir, rec_type=args.algorithm+args.convention, name='data')
            randoms_rec_fn = catalog_fn(**catalog_kwargs, cat_dir=out_dir, rec_type=args.algorithm+args.convention, name='randoms')
            if args.prepare_blinding:
                data_rec_fn = catalog_fn(**catalog_kwargs, cat_dir=out_dir, rec_type=args.algorithm+'rsd', name='data')
            run_reconstruction(Reconstruction, distance, data_fn, randoms_fn, data_rec_fn, randoms_rec_fn, f=f, bias=bias, boxsize=args.boxsize, nmesh=args.nmesh, cellsize=args.cellsize, smoothing_radius=args.smoothing_radius, nthreads=args.nthreads, convention='rsd' if args.prepare_blinding else args.convention, dtype='f8', zlim=(zmin, zmax), weight_type=args.weight_type, mpicomm=mpicomm)
