import LSS.combined_tracer_utils as comb
import LSS.common_tools as common
from astropy.table import Table, vstack
import numpy as np
import argparse
import os
import logging

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
parser.add_argument('--cap', help='NGC or SGC')
parser.add_argument('--nrands', help='number of random files to process',default=18,type=int)
parser.add_argument('--verbose', help='True of False, prints out progress steps', type=bool, default=False)
args = parser.parse_args()

base_dir = args.base_dir
save_dir = args.save_dir
cap = args.cap
verbose = args.verbose
if verbose:
    print(f'Loading from {base_dir}')
    print(f'Saving to {save_dir}')
    print(f'cap = {cap}')
    
if not os.path.exists(save_dir):
    os.makedirs(save_dir)
nrands = args.nrands  #Number of random catalogs
P0 = 6000
dz = 0.01
bias_list = [2.0, 1.2] #LRG, ELG

#get comp_ntl
#if mpicomm.rank == 0:
fb1 = base_dir + f'LRG_{cap}'
fb2 = base_dir + f'ELG_LOPnotqso_{cap}'
comp_ntl1 = comb.get_comp(fb1)
comp_ntl2 = comb.get_comp(fb2)
#else:
#    comp_ntl1, comp_ntl2, comp_ntl3, comp_ntl4 = None, None, None, None
#comp_ntl1 = mpicomm.bcast(comp_ntl1, root=0)
#comp_ntl2 = mpicomm.bcast(comp_ntl2, root=0)

#get nz
nz_file1 = np.loadtxt(base_dir + f'LRG_{cap}_nz.txt')
nz_file2 = np.loadtxt(base_dir + f'ELG_LOPnotqso_{cap}_nz.txt')
nz_list = [nz_file1, nz_file2]

#zmin, zmax = comb.setup_binning(nz_list, verbose) #setup_binning currently not working as intended
zmin, zmax = 0.8, 1.1

#get neff
neff, nz_comb_all = comb.calc_neff(nz_list, bias_list, zmin, zmax, dz, verbose)

#Load Data
d1_fn = base_dir + f'LRG_{cap}_clustering.dat.fits'
d2_fn = base_dir + f'ELG_LOPnotqso_{cap}_clustering.dat.fits'

dcat1, nxfacd1 = comb.read_data(d1_fn, comp_ntl1, zmin, zmax, verbose)
dcat2, nxfacd2 = comb.read_data(d2_fn, comp_ntl2, zmin, zmax, verbose)

#Calculate data fkp weights
fkp_d1 = comb.calc_fkp(nxfacd1, dcat1['Z'], neff, P0, zmin, zmax, dz)
fkp_d2 = comb.calc_fkp(nxfacd2, dcat2['Z'], neff, P0, zmin, zmax, dz)

#Total number of data
N_d1 = np.sum(dcat1['WEIGHT'] * fkp_d1)
N_d2 = np.sum(dcat2['WEIGHT'] * fkp_d2)

#New data weights
weight_d1 = dcat1['WEIGHT'] * fkp_d1 * bias_list[0]
weight_d2 = dcat2['WEIGHT'] * fkp_d2 * bias_list[1]

#Concatenate data tracers into one catalog
weight_d = np.concatenate((weight_d1, weight_d2))
fkp_d = np.concatenate((fkp_d1, fkp_d2))

#Data Catalog
dcat = vstack([Table(dcat1), Table(dcat2)])
dcat['WEIGHT'] = weight_d / fkp_d
dcat['WEIGHT_FKP'] = fkp_d

save_data_fn = save_dir + f'LRG+ELG_LOPnotqso_{cap}_clustering.dat.fits'

common.write_LSS_scratchcp(dcat,save_data_fn,logger=logger)

def _make_rancat(rdmnb):
    '''
    
    rdmnb = int : index of random file
    '''
    #Load Randoms
    r1_fn = base_dir + f'LRG_{cap}_{rdmnb}_clustering.ran.fits'
    r2_fn = base_dir + f'ELG_LOPnotqso_{cap}_{rdmnb}_clustering.ran.fits'
    
    rcat1, nxfacr1 = comb.read_rand(r1_fn.replace('global','dvs_ro'), comp_ntl1, zmin, zmax, verbose)
    rcat2, nxfacr2 = comb.read_rand(r2_fn.replace('global','dvs_ro'), comp_ntl2, zmin, zmax, verbose)
    
    
    #Calculate random fkp weights
    fkp_r1 = comb.calc_fkp(nxfacr1, rcat1['Z'], neff, P0, zmin, zmax, dz)
    fkp_r2 = comb.calc_fkp(nxfacr2, rcat2['Z'], neff, P0, zmin, zmax, dz)
    
    
    #Total number of randoms
    N_r1 = np.sum(rcat1['WEIGHT'] * fkp_r1)
    N_r2 = np.sum(rcat2['WEIGHT'] * fkp_r2)
    
    
    #New random weights
    weight_r1 = rcat1['WEIGHT'] * fkp_r1 * bias_list[0]
    weight_r2 = rcat2['WEIGHT'] * fkp_r2 * bias_list[1] * N_d2 * N_r1 / (N_d1 * N_r2)
    
    
    #Concatenate random tracers into one catalog
    weight_r = np.concatenate((weight_r1, weight_r2))
    fkp_r = np.concatenate((fkp_r1, fkp_r2))
    
    
    #Random Catalog
    rcat = vstack([Table(rcat1), Table(rcat2)])
    rcat['WEIGHT'] = weight_r / fkp_r
    rcat['WEIGHT_FKP'] = fkp_r
    
    #Save Catalogs
    save_rand_fn = save_dir + f'LRG+ELG_LOPnotqso_{cap}_{rdmnb}_clustering.ran.fits'
    
    common.write_LSS_scratchcp(rcat,save_rand_fn,logger=logger)


#run main func
rand_idxs = np.arange(nrands)
with Pool(processes=nrands) as pool:
    res = pool.map(_make_rancat, rand_idxs)