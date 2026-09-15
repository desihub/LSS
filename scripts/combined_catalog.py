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

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--base_dir', help='directory to load from')
parser.add_argument('--save_dir', help='directory to write combined catalogs to')
parser.add_argument('--in_tracers', help='input tracers (eg. --in_tracers "LRG" "ELG_LOPnotqso")', nargs='+', default=['LRG','ELG_LOPnotqso'], choices=bias_dict.keys())
parser.add_argument('--out_tracer', help='name of the tracer in output files', default='LRG+ELG_LOPnotqso')
parser.add_argument('--cap', help='NGC or SGC', choices=['NGC', 'SGC'])
parser.add_argument('--nrands', help='number of random files to process', default=18, type=int)
parser.add_argument('--zmin', help='minimum redshift to use for cuts', default=0.4, type=float)
parser.add_argument('--zmax', help='maximum redshift to use for cuts', default=2.1, type=float)
parser.add_argument('--verbose', help='True of False, prints out progress steps', type=bool, default=False)
parser.add_argument('--rand_unique_ids', help='True of False, additional step to ensure unique TARGETIDs in the final random catalogs (and keeping the 2500/deg2 density if used without redshift bounds)', type=bool, default=False)
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
zmin, zmax = args.zmin, args.zmax
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

if args.rand_unique_ids: N_d_raw = np.array([len(dcat[i]) for i in range(ntracers)]) # store the raw counts of each tracer's data catalog for later use in the random catalog subselection/weighting
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
        rcat[i]['WEIGHT'] *= bias_list[i] # upweight each tracer by its bias
        rcat[i]['WEIGHT'] *= N_d[i] / N_r[i] / (N_d[0] / N_r[0]) # additional step for randoms: match the random-to-data ratio (default x FKP weighted) for each tracer to the first tracer. could just make the ratio 1 for simplicity, but it wouldn't match the description in Equation 4.14 of https://arxiv.org/pdf/2508.05467v2, and that also seems to create an issue with thecov Fourier-space covariances (https://github.com/cosmodesi/thecov) from our experience with BGS_BRIGHT+FAINT
        rcat[i]['TRACER_TYPE'] = i # mark the tracer type to e.g. easily divide the combined catalogs into pieces later if needed

    if args.rand_unique_ids:
        np.random.seed(rdmnb) # for reproducibility, but different for each random catalog
        rcat_unique = [] # catalog pieces with unique TARGETIDs will be collected here to then concatenate into the final random catalog
        for tracer_bit_encoding in range(1, 2**ntracers): # go over all possible tracer intersections, bits encoding which tracers are included/excluded. skip the 0 case (no tracers) since that would be empty for sure
            # get the indices of the tracers included in this intersection; at least one is included since we skip the tracer_bit_encoding=0 case
            included_tracers = np.array([i for i in range(ntracers) if (tracer_bit_encoding & (1 << i)) != 0])
            # get the indices of the tracers excluded from this intersection
            excluded_tracers = np.array([i for i in range(ntracers) if (tracer_bit_encoding & (1 << i)) == 0])

            # find the TARGETIDs strictly in this intersection
            current_targetids = rcat[included_tracers[0]]['TARGETID']
            for i in included_tracers[1:]: current_targetids = np.intersect1d(current_targetids, rcat[i]['TARGETID'], assume_unique=True) # TARGETIDs should be unique within each tracer catalog, so assume_unique=True should be safe and may be faster
            for i in excluded_tracers: current_targetids = np.setdiff1d(current_targetids, rcat[i]['TARGETID'], assume_unique=True) # TARGETIDs should be unique within each tracer catalog, so assume_unique=True should be safe and may be faster

            if len(current_targetids) == 0: continue # nothing to be done if this strict intersection is empty

            if len(included_tracers) == 1:
                # if only one tracer is included, just keep the rows with the unique TARGETIDs in that tracer's catalog
                mask = np.isin(rcat[included_tracers[0]]['TARGETID'], current_targetids)
                rcat_unique.append(rcat[included_tracers[0]][mask])
                del mask # no longer needed, free memory
                continue # done for this intersection, move on to the next one

            # case of multiple tracers included in the intersection, need to select which tracer catalog to draw random from for each of current_targetids
            n_random_goals = len(current_targetids) * N_d_raw[included_tracers] / N_d_raw[included_tracers].sum() # goal number of randoms to draw is proportional to the number of data in each tracer catalog (before consequent rounding to integer). this is a simple and reasonable choice. an alternative could be to use the number of data in the sky area corresponding to the intersection, but we don't seem to have a good way to compute that, and it may not be worth the effort anyway
            n_random_split = np.rint(np.cumsum(n_random_goals)[:-1]).astype(int) # get the indices to split the shuffled randoms. round each to the nearest integer
            np.random.shuffle(current_targetids) # randomly shuffle the TARGETIDs in-place to then split and select for each tracer catalog
            targetids_sel_all = np.split(current_targetids, n_random_split) # split the shuffled indices according to the number of randoms to draw for each sample
            for i, targetids_sel in zip(included_tracers, targetids_sel_all):
                if len(targetids_sel) == 0: continue # check just in case some of the splits are empty for very small intersections, though that should be rare. this presents a bit of a problem for the tracer sky density, but hopefully only could happen for very small-area intersections
                mask = np.isin(rcat[i]['TARGETID'], targetids_sel)
                rcat_unique.append(rcat[i][mask])
                rcat_unique[-1]['WEIGHT'] *= len(current_targetids) / len(targetids_sel) # upweight the selected randoms to account for the fact that we are keeping only targetids_sel out of current_targetids for this tracer. in-place multiplication is fine, as this set of random should not be encountered again. NB: this simple number-based upscaling could cause additional fluctuations in weighted random density in redshift and/or on sky; using the weight ratio may be better in that respect, but may have an issue of overly fine tuning for small intersections (and we also may need to be more careful about multiplying the weights in place)
                del mask # no longer needed, free memory
            del current_targetids, targetids_sel_all # no longer needed, free memory
        del rcat # no longer needed, free memory
        rcat_concat = vstack(rcat_unique) # concatenate the catalog pieces with unique TARGETIDs. the order may be a bit strange, but that should not matter. for unscrambling the catalogs, we have the TRACER_TYPE column
        del rcat_unique # no longer needed, free memory
    else:
        rcat_concat = vstack(rcat) # simply concatenate catalogs
        del rcat # no longer needed, free memory
    
    save_ran_fn = save_dir + f'{out_tracer}_{cap}_{rdmnb}_clustering.ran.fits'
    common.write_LSS_scratchcp(rcat_concat,save_ran_fn,logger=logger)


#run main func
rand_idxs = np.arange(nrands)
with Pool(processes=nrands) as pool:
    res = pool.map(_make_rancat, rand_idxs)