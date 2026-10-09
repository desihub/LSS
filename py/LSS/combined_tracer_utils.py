import numpy as np
from astropy.table import Table
import fitsio

def read_catalog(fn, comp_ntl, zmin, zmax, verbose=False, logger=None, kind='data'):
    '''
    Reads data or random clustering catalog
    '''
    dcat = Table(fitsio.read(fn))
    mask = (dcat['Z'] > zmin) & (dcat['Z'] < zmax)
    dcat = dcat[mask]
    nxfacd = comp_ntl[dcat['NTILE']-1]
    if verbose: logger.info(f'Loaded {kind} file at ' + fn)
    return dcat, nxfacd

def get_comp(fb, verbose=False, logger=None):
    '''
    Calculates completeness per number of tiles
    '''
    fn = fb + '_clustering.dat.fits'
    fd = Table(fitsio.read(fn))
    mean_comp = 1/np.mean(fd['WEIGHT_COMP'])
    if verbose: logger.info(f'mean completeness = {mean_comp}')
    weight_ntl = np.bincount(fd['NTILE']-1, weights=fd['WEIGHT_COMP']) / np.bincount(fd['NTILE']-1) # mean of WEIGHT_COMP for each (positive integer) NTILE in the data. Note that the NTILE values are shifted down by 1 to avoid guaranteed division by zero for NTILE=0
    comp_ntl = 1 / weight_ntl
    fran = fitsio.read(fb+'_0_clustering.ran.fits',columns=['NTILE','FRAC_TLOBS_TILES'])
    fttl = np.bincount(fran['NTILE']-1, weights=fran['FRAC_TLOBS_TILES']) / np.bincount(fran['NTILE']-1) # mean of FRAC_TLOBS_TILES for each (positive integer) NTILE in randoms. Note that the NTILE values are shifted down by 1 to avoid guaranteed division by zero for NTILE=0
    comp_ntl = comp_ntl*fttl
    if verbose: logger.info(f'completeness per ntile: \n{np.arange(len(comp_ntl))+1} \n{comp_ntl}')
    return comp_ntl

def setup_binning(nz_list, verbose=False, logger=None):
    '''
    Sets up zbins
    
    returns min and max of z combined zbins
    '''
    ntracers = len(nz_list)
    zmin_ti = []
    zmax_ti = []
    for i in range(ntracers):
        nz_file = nz_list[i]
        zmin_ti.append(min(nz_file[:,1]))
        zmax_ti.append(max(nz_file[:,2]))
    zmin = min(zmin_ti)
    zmax = max(zmax_ti)
    if verbose: logger.info(f'zmin = {zmin}, zmax = {zmax}')
    return zmin, zmax

def rebin_nz(z_comb, nz_t, zmin_arr_t, dz_t):
    '''
    Rebins n(z) for a given tracer to the combined zbins
    
    z_comb = centers of combined zbins
    nz_t = n(z) for given tracer
    zmin_arr_t = lower edges of tracer zbins
    dz_t = size of tracer zbins
    '''
    zmin_t = zmin_arr_t[0]
    zmax_t = zmin_arr_t[-1] + dz_t
    nbins = len(z_comb)
    nz_t_comb = np.zeros(nbins, dtype=float)
    for i in range(nbins):
        z = z_comb[i]
        if z > zmin_t and z < zmax_t:
            i_t = int((z - zmin_t) / dz_t)
            nz_t_comb[i] = nz_t[i_t]
    return nz_t_comb

def calc_neff(nz_list, bias_list, zmin, zmax, dz, verbose=False, logger=None):
    '''
    Combines n(z) files of each tracer into one neff(z) and rebinned n_t(z)
    
    returns neff[zbin], nz_comb_all[tracer, zbin], and beff[zbin]
    '''
    ntracers = len(nz_list)
    nbins = int((zmax - zmin) / dz)
    zmin_comb = np.linspace(zmin, zmax, nbins, endpoint=False)
    zmax_comb = zmin_comb + dz
    z_comb = (zmin_comb + zmax_comb) / 2
    nz_comb_all = np.zeros((ntracers,nbins), dtype=float)
    bnz_comb_all = np.zeros((ntracers,nbins), dtype=float)
    b2nz_comb_all = np.zeros((ntracers,nbins), dtype=float)
    
    for i in range(ntracers):
        nz_file = nz_list[i]
        bias = bias_list[i]
        zmin_t = nz_file[:,1]
        nz_t = nz_file[:,3]
        dz_t = zmin_t[1] - zmin_t[0]
        nz_t_comb = rebin_nz(z_comb, nz_t, zmin_t, dz_t)
        nz_comb_all[i,:] = nz_t_comb
        bnz_comb_all[i,:] = bias * nz_t_comb
        b2nz_comb_all[i,:] = bias**2 * nz_t_comb
    
    beff = np.sum(b2nz_comb_all, axis=0) / np.sum(bnz_comb_all, axis=0)
    neff = np.sum(bnz_comb_all, axis=0) / beff
    if verbose: logger.info(f'neff(z) = {neff}\n\n beff(z) = {beff}')
    return neff, nz_comb_all, beff

def calc_fkp(nxfac, zall, neff, P0, zmin, zmax, dz, tracer):
    '''
    Calculates new fkp weights for combined catalog
    '''
    fkp = np.zeros(len(zall))
    for i in range(len(zall)):
        z = zall[i]
        zind = int((z - zmin) / dz)
        fkp[i] = 1. / (1. + neff[zind]*nxfac[i]*P0[zind])
    return fkp