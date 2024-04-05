# Get bitmask values from pixel-level per-brick masks for a catalog
from __future__ import division, print_function
import sys, os, glob, time, warnings, gc
import numpy as np
import matplotlib.pyplot as plt
from astropy.table import Table, vstack, hstack, join
import fitsio

from astropy.io import fits
from astropy import wcs

from multiprocessing import Pool
import argparse

from mockfactory.desi import get_brick_pixel_quantities


import os
import logging

import fitsio
import numpy as np
from mpi4py import MPI


# create logger
logname = 'masknobs'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)

# create console handler and set level to debug
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)

# create formatter
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')

# add formatter to ch
ch.setFormatter(formatter)

# add ch to logger
logger.addHandler(ch)



#time_start = time.time()



parser = argparse.ArgumentParser()
parser.add_argument( '--cat_type', default='obielg',required=False)#, choices=['obielg', 'abacus'])
parser.add_argument( '--reg', default='north', choices=['north','south'],required=False)
parser.add_argument('--abacus_fa_num', default = 0, required = False)
parser.add_argument('--do_randoms', default = 'n', choices = ['n','y'], required = False)
parser.add_argument('--random_tracer', default = 'LRG', required = False)
parser.add_argument('--mock_number', default = 0, required = False)
parser.add_argument('--outdir', default = '', required=False )
parser.add_argument('--overwrite', default = 'n', required=False )
parser.add_argument('--test', default = 'n', required=False )

args = parser.parse_args()

output_path = None
if args.cat_type == 'obielg':
    input_path = '/global/cfs/cdirs/desi/survey/catalogs/image_simulations/ELG/dr9/Y1/'+args.reg+'/file0_rs0_skip0/merged/matched_input_full.fits'
    output_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/LSS/elg_obiwan_'+args.reg+'_matched_input_full_masknobs.fits'
    
if args.cat_type == 'abacus':
    if args.do_randoms == 'n':
        fa_num = args.abacus_fa_num
        str_fa_num = str(fa_num)
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/"
        input_path = '/global/cfs/cdirs/desi/survey/catalogs/Y1/mocks/SecondGenMocks/AbacusSummit_v4_1/forFA' + str_fa_num + '.fits'
        if args.outdir != '':
            output_path = args.outdir + "/" + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"
        elif args.outdir == '':
            output_path = input_dir + "forFA" + str_fa_num + "_matched_input_full_masknobs.fits"
            print(output_path)
            
    elif args.do_randoms == 'y':
        ran_tr = args.random_tracer
        mockno = args.mock_number
        print("Running for Mock %s on Tracer %s"%(mockno, ran_tr))
        input_dir = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/"%(mockno)
        input_path = "/global/cfs/cdirs/desi/survey/catalogs/main/mocks/FirstGenMocks/AbacusSummit/Y1/mock%s/LSScats/%s_1_full_noveto.ran.fits"%(mockno, ran_tr)
        if args.outdir != '':
            output_path = args.outdir + "/" + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"
            print("Output to " + output_path)
        elif args.outdir == '':
            output_path = input_dir + ran_tr + "_1_full_matched_input_full_masknobs.ran.fits"
            print("Output to " + output_path)

if args.cat_type == 'genran':
    from mockfactory import RandomCutskyCatalog
    cutsky = RandomCutskyCatalog(rarange=(0., 180.), decrange=(0, 90.), csize=3e7, seed=44, mpicomm=mpicomm)
    ra, dec = cutsky['RA'], cutsky['DEC']

bitmask_dir = '/global/cfs/cdirs/cosmo/data/legacysurvey/dr9/'

# input_path = '/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits'
# output_path = '/global/cscratch1/sd/rongpu/temp/randoms-1-0-lrgmask_v1.fits'
fe = False

mpicomm = MPI.COMM_WORLD
mpiroot = 0

ra,dec = None,None

if mpicomm.rank == mpiroot:
	if os.path.isfile(output_path):
		fe = True
		if args.overwrite == 'n' and args.test == 'n':
			raise ValueError(output_path+' already exists!')
		if args.overwrite == 'n' and args.test == 'y':
			logger.info('will run and get timing, but no output will be written (because path exists)')
		if args.overwrite == 'y':
			logger.info('will overwrite '+output_path)



	if args.cat_type == 'obielg':
		cat = Table(fitsio.read(input_path,columns=['input_ra','input_dec']))
		cat.rename_columns(['input_ra', 'input_dec'], ['RA', 'DEC'])
	else:
		cat = Table(fitsio.read(input_path))
	
	logger.info('loaded catalog length '+str(len(cat)))
	
	for col in cat.colnames:
		cat.rename_column(col, col.upper())

	if 'TARGETID' not in cat.colnames:
		cat['TARGETID'] = np.arange(len(cat))
	
	if 'TARGET_RA' in cat.colnames:
		cat.rename_columns(['TARGET_RA', 'TARGET_DEC'], ['RA', 'DEC'])
	
	if 'BRICKID' not in cat.colnames:
		from desiutil import brick
		tmp = brick.Bricks(bricksize=0.25)
		cat['BRICKID'] = tmp.brickid(cat['RA'], cat['DEC'])
	
	ra, dec = cat['RA'], cat['DEC']

start = MPI.Wtime()
columns = {}
columns['maskbits'] = {'fn': '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{region}/coadd/{brickname:.3s}/{brickname}/legacysurvey-{brickname}-maskbits.fits.fz', 'dtype': 'i2', 'default': 1}
bl = ['g','r','z']
for band in bl:
    columns['nobs_'+band] = {'fn': '/dvs_ro/cfs/cdirs/cosmo/data/legacysurvey/dr9/{region}/coadd/{brickname:.3s}/{brickname}/legacysurvey-{brickname}-nexp-'+band+'.fits.fz', 'dtype': 'i2', 'default': 0}
columns['brickname'] = None
columns['photsys'] = None
catalog = get_brick_pixel_quantities(ra, dec, columns, mpicomm=mpicomm)
if mpicomm.rank == 0:
    logger.info('Output columns are {}.'.format(list(catalog.keys())))
    logger.info('Pixel-level quantities read in {:2.2f} s.'.format(MPI.Wtime() - start))


logger.info('Done!')
