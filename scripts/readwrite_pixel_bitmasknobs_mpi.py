# Get bitmask values from pixel-level per-brick masks for a catalog
import os
import logging
import argparse

import numpy as np

from mockfactory import Catalog, RandomCutskyCatalog, setup_logging
from mockfactory.desi import get_brick_pixel_quantities

from mpi4py import MPI
start = MPI.Wtime()
mpicomm = MPI.COMM_WORLD


logger = logging.getLogger('LSS masknobs')
setup_logging(level=(logging.INFO if mpicomm.rank == 0 else logging.ERROR))


parser = argparse.ArgumentParser()
parser.add_argument( '--cat_type', default='obielg', required=False)#, choices=['obielg', 'abacus'])
parser.add_argument( '--reg', default='north', choices=['north','south'], required=False)
parser.add_argument('--abacus_version', type=str, default='v4_1', required=False)
parser.add_argument('--abacus_fa_num_min', type=int, default=0, required=False)
parser.add_argument('--abacus_fa_num_max', type=int, default=1, required=False)
parser.add_argument('--do_randoms', default='n', choices = ['n', 'y'], required=False)
parser.add_argument('--random_tracer', default='LRG', required=False)
parser.add_argument('--outdir', default='', required=False)
parser.add_argument('--overwrite', default='n', required=False)
parser.add_argument('--test', default='n', required=False)
args = parser.parse_args()
logger.info(args)


if args.cat_type == 'genran':
    #cutsky = RandomCutskyCatalog(rarange=(0., 180.), decrange=(0, 90.), csize=int(3e7), seed=44, mpicomm=mpicomm)
    cat_list = [RandomCutskyCatalog(rarange=(28., 30.), decrange=(1., 2.), csize=10000, seed=44, mpicomm=mpicomm)]
    cat_size_list = [cat.size for cat in cat_list]
    ra, dec = np.concatenate([cat['RA'] for cat in cat_list]), np.concatenate([cat['DEC'] for cat in cat_list])


if args.cat_type == 'obielg':
    input_path = f'/dvs_ro/cfs/cdirs/desi/survey/catalogs/image_simulations/ELG/dr9/Y1/{args.reg}/file0_rs0_skip0/merged/matched_input_full.fits'
    output_path = f'/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/elg_obiwan_{args.reg}_matched_input_full_masknobs.fits'
    logger.info("Output to " + output_path)

    logger.info('Read '+input_path)
    cat_list = [Catalog.read(input_path, mpicomm=mpicomm)]	
    cat_size_list = [cat.size for cat in cat_list]
    ra, dec = np.concatenate([cat['input_ra'] for cat in cat_list]), np.concatenate([cat['input_dec'] for cat in cat_list])
    #cat['RA'], cat['DEC'] = cat['input_ra'], cat['input_dec']
    logger.info('Nbr objects = ' + str(cat.csize))


if args.cat_type == 'abacus':
    if args.do_randoms == 'n':
        input_path_list = [f'/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_{args.abacus_version}/forFA{i}.fits' for i in range(args.abacus_fa_num_min, args.abacus_fa_num_max)]
        outdir = args.outdir if (args.outdir != '') else f"/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_{args.abacus_version}/"
        output_path_list = [f"{args.outdir}/forFA{i}_matched_input_full_masknobs.fits" for i in range(args.abacus_fa_num_min, args.abacus_fa_num_max)]
                
    elif args.do_randoms == 'y':
        input_path_list = [f"/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock{i}/LSScats/{args.random_tracer}_1_full_noveto.ran.fits" for i in range(args.abacus_fa_num_min, abacus_fa_num_max)]
        outdir = args.outdir if (args.outdir != '') else f"/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/"
        output_path_list = [f"{args.outdir}/mock{i}/LSScats/{args.random_tracer}_1_full_matched_input_full_masknobs.ran.fits" for i in range(args.abacus_fa_num_min, args.abacus_fa_num_max)]

    logger.info(f'Read {input_path_list}')
    cat_list = [Catalog.read(input_path, mpicomm=mpicomm) for input_path in input_path_list]
    cat_size_list = [cat.size for cat in cat_list]
    nbr_objects = np.sum([cat.csize for cat in cat_list])
    logger.info(f'Nbr objects = {nbr_objects}')
    ra, dec = np.concatenate([cat['RA'] for cat in cat_list]), np.concatenate([cat['DEC'] for cat in cat_list])


columns = {}
columns['maskbits'] = {'fn': '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{region}/coadd/{brickname:.3s}/{brickname}/legacysurvey-{brickname}-maskbits.fits.fz', 'dtype': 'i2', 'default': 1}
bl = ['g','r','z']
for band in bl:
    columns['nobs_'+band] = {'fn': '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{region}/coadd/{brickname:.3s}/{brickname}/legacysurvey-{brickname}-nexp-'+band+'.fits.fz', 'dtype': 'i2', 'default': 0}
columns['brickname'] = None
columns['photsys'] = None

logger.info(f'getting brick pixel quantities: {columns.keys()}')
catalog = get_brick_pixel_quantities(ra, dec, columns, mpicomm=mpicomm)
logger.info('Output columns are {}.'.format(list(catalog.keys())))
logger.info('Pixel-level quantities read in {:2.2f} s.'.format(MPI.Wtime() - start))

for i in range(len(cat_list)):
    cat, output_path = cat_list[i], output_path_list[i]

    idx_begin, idx_end = int(np.sum(cat_size_list[:i])), int(np.sum(cat_size_list[:i+1])) if i < len(cat_size_list) else int(np.sum(cat_size_list))
    cat['MASKBITS'] = catalog['maskbits'][idx_begin:idx_end]
    for band in bl:
         cat['NOBS_{band}'] = catalog['nobs_{band}'][idx_begin:idx_end]

    if args.overwrite == 'n' and args.test == 'n':
        raise ValueError(output_path + ' already exists!')
    if args.overwrite == 'n' and args.test == 'y':
        logger.info('No output will be written (because path exists)')
    if args.overwrite == 'y':
        logger.info('Will overwrite ' + output_path)

    ## SAVE FILE:
    if (not os.path.isfile(output_path)) or (args.overwrite == 'y'):
        logger.info('Write: {output_path}')
         cat.write(output_path, mpicomm=mpicomm)

