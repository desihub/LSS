# Get bitmask values from pixel-level per-brick masks for a catalog
import os
import logging
import argparse

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
parser.add_argument('--abacus_fa_num', default=0, required=False)
parser.add_argument('--do_randoms', default='n', choices = ['n', 'y'], required=False)
parser.add_argument('--random_tracer', default='LRG', required=False)
parser.add_argument('--mock_number', default=0, required=False)
parser.add_argument('--outdir', default='', required=False)
parser.add_argument('--overwrite', default='n', required=False)
parser.add_argument('--test', default='n', required=False)
args = parser.parse_args()
logger.info(args)


if args.cat_type == 'genran':
    #cutsky = RandomCutskyCatalog(rarange=(0., 180.), decrange=(0, 90.), csize=int(3e7), seed=44, mpicomm=mpicomm)
    cat = RandomCutskyCatalog(rarange=(28., 30.), decrange=(1., 2.), csize=10000, seed=44, mpicomm=mpicomm)
    ra, dec = cat['RA'], cat['DEC']


if args.cat_type == 'obielg':
    input_path = f'/dvs_ro/cfs/cdirs/desi/survey/catalogs/image_simulations/ELG/dr9/Y1/{args.reg}/file0_rs0_skip0/merged/matched_input_full.fits'
    output_path = f'/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/elg_obiwan_{args.reg}_matched_input_full_masknobs.fits'
    logger.info("Output to " + output_path)

    logger.info('Read '+input_path)
    cat = Catalog.read(input_path, mpicomm=mpicomm)	
    ra, dec = cat['input_ra'], cat['input_dec']
    #cat['RA'] = cat['input_ra']
    #cat['DEC'] = cat['input_dec']
    logger.info('Nbr objects = ' + str(cat.csize))


if args.cat_type == 'abacus':
    if args.do_randoms == 'n':
        fa_num = args.abacus_fa_num
        str_fa_num = str(fa_num)
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/"
        input_path = '/dvs_ro/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/forFA' + str_fa_num + '.fits'
        if args.outdir != '':
            output_path = args.outdir + "/" + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"
        elif args.outdir == '':
            output_path = input_dir + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"

    elif args.do_randoms == 'y':
        ran_tr = args.random_tracer
        mockno = args.mock_number
        logger.info("Running for Mock %s on Tracer %s"%(mockno, ran_tr))
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/"%(mockno)
        input_path = "/dvs_ro/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/%s_1_full_noveto.ran.fits"%(mockno, ran_tr)
        if args.outdir != '':
            output_path = args.outdir + "/" + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"
        elif args.outdir == '':
            output_path = input_dir + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"

    logger.info("Output to " + output_path)

    logger.info('Read '+input_path)
    cat = Catalog.read(input_path, mpicomm=mpicomm)	
    logger.info('Nbr objects = ' + str(cat.csize))
    
    # for name in cat.columns():
    #     cat[name.upper()] = cat[name]
    if 'TARGET_RA' in cat.columns():
        cat['RA'] = cat['TARGET_RA']
        cat['DEC'] = cat['TARGET_DEC']
    ra, dec = cat['RA'], cat['DEC']


## Check is the file already exist:
if os.path.isfile(output_path):
    fe = True
if args.overwrite == 'n' and args.test == 'n':
    raise ValueError(output_path + ' already exists!')
if args.overwrite == 'n' and args.test == 'y':
    logger.info('Will run and get timing, but no output will be written (because path exists)')
if args.overwrite == 'y':
    logger.info('Will overwrite ' + output_path)


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

cat['MASKBITS'] = catalog['maskbits']
for band in bl:
    cat['NOBS_'+band] = catalog['nobs_'+band]


## SAVE FILE:
if (not fe) or (args.overwrite == 'y'):
    car.write(output_path, mpicomm=mpicomm)

