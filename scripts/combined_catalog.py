import LSS.combined_tracer_utils as comb
import LSS.common_tools as common
from astropy.table import Table, vstack
import numpy as np
import argparse
import os
import logging
from cosmoprimo.fiducial import DESI
cosmo = DESI()

from multiprocessing import Pool

logname = 'mkCat'
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


parser = argparse.ArgumentParser()
parser.add_argument('--base_dir', help='directory to load from')
parser.add_argument('--save_dir', help='directory to write combined catalogs to')
parser.add_argument('--in_tracers', help='input tracers (eg. --in_tracers "LRG" "ELG_LOPnotqso")', nargs='+', default=['LRG','ELG_LOPnotqso'])
parser.add_argument('--out_tracer', help='name of the tracer in output files', default='LRG+ELG')
parser.add_argument('--cap', help='NGC or SGC')
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
bias_list = [2.0, 1.2, 2.1] #LRG, ELG, QSO

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
nxfacd = [None] * ntracers
fkp = [None] * ntracers
N_d = [None] * ntracers
weights = [None] * ntracers

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
    dcat[i], nxfacd[i] = comb.read_data(d_fn, comp_ntl[i], zmin, zmax, verbose, logger=logger)
    fkp[i] = comb.calc_fkp(nxfacd[i], dcat[i]['Z'], neff, P0, zmin, zmax, dz, tracer)
    N_d[i] = np.sum(dcat[i]['WEIGHT'] * fkp[i])
    weights[i] = dcat[i]['WEIGHT'] * fkp[i] * bias_list[i]
    dcat[i]['TRACER_TYPE'] = i

# Concatenate catalogs
weight_concat = np.concatenate(weights)
fkp_concat = np.concatenate(fkp)
dcat_concat = vstack([Table(dcat[i]) for i in range(ntracers)])
dcat_concat['WEIGHT'] = weight_concat / fkp_concat
dcat_concat['WEIGHT_FKP'] = fkp_concat

save_data_fn = save_dir + f'{out_tracer}_{cap}_clustering.dat.fits'
common.write_LSS_scratchcp(dcat_concat,save_data_fn,logger=logger)

def _make_rancat(rdmnb):
    '''
    
    rdmnb = int : index of random file
    '''
    # Setup needed lists
    rcat = [None] * ntracers
    nxfacr = [None] * ntracers
    fkp_r = [None] * ntracers
    N_r = [None] * ntracers
    weights_r = [None] * ntracers

    # Read data and compute weights
    for i, tracer in enumerate(tracers):
        r_fn = base_dir + f'{tracer}_{cap}_{rdmnb}_clustering.ran.fits'
        rcat[i], nxfacr[i] = comb.read_rand(r_fn.replace('global','dvs_ro'), comp_ntl[i], zmin, zmax, verbose,logger=logger)
        fkp_r[i] = comb.calc_fkp(nxfacr[i], rcat[i]['Z'], neff, P0, zmin, zmax, dz, tracer)
        N_r[i] = np.sum(rcat[i]['WEIGHT'] * fkp_r[i])
        weights_r[i] = rcat[i]['WEIGHT'] * fkp_r[i] * bias_list[i] * N_d[i] / N_r[i]
        rcat[i]['TRACER_TYPE'] = i

    # Concatenate catalogs
    weight_concat_r = np.concatenate(weights_r)
    fkp_concat_r = np.concatenate(fkp_r)
    rcat_concat = vstack([Table(rcat[i]) for i in range(ntracers)])
    rcat_concat['WEIGHT'] = weight_concat_r / fkp_concat_r
    rcat_concat['WEIGHT_FKP'] = fkp_concat_r
    
    save_ran_fn = save_dir + f'{out_tracer}_{cap}_{rdmnb}_clustering.ran.fits'
    common.write_LSS_scratchcp(rcat_concat,save_ran_fn,logger=logger)


#run main func
rand_idxs = np.arange(nrands)
with Pool(processes=nrands) as pool:
    res = pool.map(_make_rancat, rand_idxs)