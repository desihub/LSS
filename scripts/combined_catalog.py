import LSS.combined_tracer_utils as comb
import LSS.common_tools as common
from astropy.table import vstack
import numpy as np
import argparse
import os
import logging
from cosmoprimo.fiducial import DESI
cosmo = DESI()

from multiprocessing import Pool

logname = 'combined_catalog'
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

bias_dict = {'LRG': 2.0, 'ELG_LOPnotqso': 1.2, 'QSO': 2.1}

parser = argparse.ArgumentParser()
parser.add_argument('--base_dir', help='directory to load from')
parser.add_argument('--save_dir', help='directory to write combined catalogs to')
parser.add_argument('--in_tracers', help='input tracers (eg. --in_tracers "LRG" "ELG_LOPnotqso")', nargs='+', default=['LRG','ELG_LOPnotqso'], choices=bias_dict.keys())
parser.add_argument('--out_tracer', help='name of the tracer in output files', default='LRG+ELG_LOPnotqso')
parser.add_argument('--cap', help='NGC or SGC', choices=['NGC', 'SGC'])
parser.add_argument('--nrands', help='number of random files to process',default=18,type=int)
parser.add_argument('--verbose', help='True of False, prints out progress steps', type=bool, default=False)
args = parser.parse_args()

base_dir = args.base_dir
save_dir = args.save_dir
cap = args.cap
verbose = args.verbose
out_tracer = args.out_tracer
if verbose:
    logger.info(f'Loading from {base_dir}')
    logger.info(f'Saving to {save_dir}')
    logger.info(f'cap = {cap}')
    
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
nrands = args.nrands  #Number of random catalogs

# Binning & Tracer Settings
dz = 0.01
tracers = args.in_tracers
ntracers = len(tracers)
bias_list = [bias_dict[tracer] for tracer in tracers]

# Setup z-binning
zmin, zmax = 0.4, 2.1
nbins = int((zmax - zmin) / dz)
zmin_comb = np.linspace(zmin, zmax, nbins, endpoint=False)
zmax_comb = zmin_comb + dz
z_comb = (zmin_comb + zmax_comb) / 2

# Setup needed lists
comp_ntl = [None] * ntracers
nz = [None] * ntracers
dcat = [None] * ntracers
N_d = [None] * ntracers

# Read completeness and n(z)
for i, tracer in enumerate(tracers):
    fb = base_dir + f'{tracer}_{cap}'
    comp_ntl[i] = comb.get_comp(fb, logger=logger)
    nz[i] = np.loadtxt(base_dir + f'{tracer}_{cap}_nz.txt')

# get neff
neff, nz_comb_all, beff = comb.calc_neff(nz, bias_list, zmin, zmax, dz, verbose, logger=logger)
f = cosmo.growth_rate(z_comb)
P0 = cosmo.pk_kz(0.14, z_comb) * (beff**2 + 2/3*f*beff + f**2/5)

# Read data and compute weights
for i, tracer in enumerate(tracers):
    d_fn = base_dir + f'{tracer}_{cap}_clustering.dat.fits'
    dcat[i], nxfacd_i = comb.read_catalog(d_fn, comp_ntl[i], zmin, zmax, verbose, logger=logger, kind='data')
    dcat[i]['WEIGHT_FKP'] = comb.calc_fkp(nxfacd_i, dcat[i]['Z'], neff, P0, zmin, zmax, dz, tracer)
    del nxfacd_i # no longer needed, free memory
    N_d[i] = np.sum(dcat[i]['WEIGHT'] * dcat[i]['WEIGHT_FKP']) # default x FKP-weighted count of galaxies in the data catalog
    dcat[i]['WEIGHT'] *= bias_list[i] # upweight each tracer by its bias
    dcat[i]['TRACER_TYPE'] = i # mark the tracer type to e.g. easily divide the combined catalog into pieces later if needed

# Concatenate catalogs
dcat_concat = vstack(dcat)
del dcat # no longer needed, free memory

save_data_fn = save_dir + f'{out_tracer}_{cap}_clustering.dat.fits'
common.write_LSS_scratchcp(dcat_concat,save_data_fn,logger=logger)
del dcat_concat # no longer needed, free memory

def _make_rancat(rdmnb):
    '''
    
    rdmnb = int : index of random file
    '''
    # Setup needed lists
    rcat = [None] * ntracers
    N_r = [None] * ntracers

    # Read data and compute weights
    for i, tracer in enumerate(tracers):
        r_fn = base_dir + f'{tracer}_{cap}_{rdmnb}_clustering.ran.fits'
        rcat[i], nxfacr_i = comb.read_catalog(r_fn.replace('global','dvs_ro'), comp_ntl[i], zmin, zmax, verbose, logger=logger, kind='random')
        rcat[i]['WEIGHT_FKP'] = comb.calc_fkp(nxfacr_i, rcat[i]['Z'], neff, P0, zmin, zmax, dz, tracer)
        del nxfacr_i # no longer needed, free memory
        N_r[i] = np.sum(rcat[i]['WEIGHT'] * rcat[i]['WEIGHT_FKP']) # default x FKP-weighted count of randoms in the current catalog
        rcat[i]['WEIGHT'] *= bias_list[i] * N_d[i] / N_r[i] # upweight each tracer by its bias, and match the default x FKP-weighted count of randoms to the same count of data galaxies for each tracer
        rcat[i]['TRACER_TYPE'] = i # mark the tracer type to e.g. easily divide the combined catalogs into pieces later if needed

    # Concatenate catalogs
    rcat_concat = vstack(rcat)
    del rcat # no longer needed, free memory
    
    save_ran_fn = save_dir + f'{out_tracer}_{cap}_{rdmnb}_clustering.ran.fits'
    common.write_LSS_scratchcp(rcat_concat,save_ran_fn,logger=logger)


#run main func
rand_idxs = np.arange(nrands)
with Pool(processes=nrands) as pool:
    res = pool.map(_make_rancat, rand_idxs)