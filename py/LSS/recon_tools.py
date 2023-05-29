#!/usr/bin/env python
# coding: utf-8


import os
import argparse
import logging

import numpy as np
from astropy.table import Table, vstack
import fitsio

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
    recon = Reconstruction(f=f, bias=bias, boxsize=boxsize, nmesh=nmesh, cellsize=cellsize, los='local', positions=data_positions, dtype=dtype, **rec_kwargs)

    recon.assign_data(data_positions, data_weights)
    #if root:
    #    logger.info('random files are',str(randoms_fn))

    #for fn in randoms_fn:
    if root:
        logger.info('Loading {}.'.format(randoms_fn))
        randoms = vstack([Table(fitsio.read(fn)) for fn in randoms_fn])
        (ra, dec, dist), randoms_weights = get_clustering_positions_weights(randoms, distance, name='randoms', **kwargs)
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
    if tracer.startswith('ELG') or tracer.startswith('QSO'):
        return 0.9, 1.3
    if tracer.startswith('LRG'):
        return 0.8, 2.
    if tracer.startswith('BGS'):
        return 0.67, 1.5

    return 0.8, 1.2


