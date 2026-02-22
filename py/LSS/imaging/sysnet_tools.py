# Copied from Alberto Rosario's tools for preparing data for sysnet
# first provide data, randoms, hpixmaps with systematics
# clean data by zcuts (only data)
# project data to HEALPix
# project randoms to HEALPix
# read file with imaging properties maps
# assign each imaging attribute to data by converting ra,dec to hpix
# 
import logging
import pandas as pd
import numpy as np
import healpy as hp
from astropy.table import Table
import fitsio as ft

import LSS.common_tools as common

bands = ['R','G','Z']
maps_dr9 = ['EBV','STARDENS'] + [f'GALDEPTH_{b}' for b in bands] + [f'PSFSIZE_{b}' for b in bands]

def prep4sysnet(data, rands, sys, allsky_rands=None, zcolumn='Z_not4clus', zmin=0.6, zmax=1.6, nran_exp=None,
                nside=256, nest=True, use_obiwan=False, columns=maps_dr9,wtmd='fracz',tp='ELG'):
    logger = logging.getLogger('prep4sysnet')
    #if zcolumn == 'Z_not4clus':
    data = do_zcut(data, zmin, zmax, zcolumn,tp=tp)
    cols = list(data.dtype.names)
    weights = np.ones_like(data[zcolumn])
    weights_ran = np.ones(len(rands))


    if wtmd == 'fracz':
        if 'FRACZ_TILELOCID' in cols:
            print('using 1/FRACZ_TILELOCID based completeness weights')
            wts = 1/data['FRACZ_TILELOCID']
        else:
            wts = data['WEIGHT_COMP']
        #if 'FRAC_TLOBS_TILES' in cols:
        #    print('using FRAC_TLOBS_TILES')
        wts *= 1/data['FRAC_TLOBS_TILES']
    if wtmd == 'probobs':
        wts = 129/(data['PROB_OBS']*128+1)
    if wtmd == 'wt_iip':
        wts = data['WEIGHT_IIP']

    if wtmd == 'wt':
        wts = data['WEIGHT']
        weights_ran = rands['WEIGHT']
    if wtmd == 'wt_comp':
        wts = data['WEIGHT_COMP']

    if 'WEIGHT_ZFAIL' in cols:
        wts *= data['WEIGHT_ZFAIL']
    # only do true for data ???

    weights *= wts
    if use_obiwan:
        weights *= data['OBI_WEIGHT']#ut.get_nn_weights(data, run, zmin, zmax, nside, hpix=None, version=version)
    
    hpmaps = create_sysmaps(sys, nest=nest, columns=columns)
    
    if allsky_rands is None:
        data_hpmap, rands_hpmap = hpixelize(nside, data, rands, weights=weights, weights_ran=weights_ran,nest=False,
                                            return_mask=False, nest2ring=False) 
        prep_table = hpdataset(data_hpmap, rands_hpmap, hpmaps, columns, fracmd='old', nran_exp=nran_exp)
    else: 
        print("using all sky randoms in fracgood")
        # now this assumes that allsky_rands is a healpix map
        data_hpmap  = hpixsum(nside,data['RA'],data['DEC'],weights=weights,nest=False,nest2ring=False)
        rands_hpmap = hpixsum(nside,rands['RA'],rands['DEC'],nest=False,nest2ring=False)
        #all_rands_hpmap = hpixsum(nside,allsky_rands['RA'],allsky_rands['DEC'],nest=False,nest2ring=False)
        all_rands_hpmap = allsky_rands
        frac_area_hpmap = rands_hpmap / all_rands_hpmap
        
        prep_table = hpdataset(data_hpmap, frac_area_hpmap, hpmaps, columns, fracmd='new',nran_exp=None)
    
    return prep_table
    
def hpdataset(data_hpmap, rands_hpmap, hpmaps, columns, nran_exp=None, frac_min=0.0,fracmd='old'):
    logger = logging.getLogger('hpdataset')

    features = hpmaps[columns].values
    is_seen = features!=hp.UNSEEN
    seen_mask = np.ones(features.shape[0], '?')
    for i,col in enumerate(columns):
        seen_mask &= is_seen[:,i]
        
    if fracmd == 'old':
        if nran_exp is None:
            nran_exp = np.mean(rands_hpmap[rands_hpmap>0])
        print(f'nran_exp: {nran_exp}')
        logger.info(f'nran_exp: {nran_exp}')
        frac = rands_hpmap / nran_exp
    if fracmd == 'new':
        # this method assumes that rands_hpmap is already a healpix map of 
        # the fractional area of each pixel
        frac = rands_hpmap
    frac_mask = frac > frac_min
    mask = frac_mask & seen_mask
    hpix = np.argwhere(mask).flatten()
        
    dtype = [('features', ('f8', features.shape[1])), 
             ('label', 'f8'),
             ('fracgood', 'f8'),
             ('hpix', 'i8')]    
        
    dataset = np.zeros(data_hpmap[mask].size, dtype=dtype)
    dataset['label'] = data_hpmap[mask]
    dataset['fracgood'] = frac[mask]
    dataset['features'] = features[mask,:]
    dataset['hpix'] = hpix
    
    return dataset

def do_zcut(data, zmin, zmax, zcolumn,tp='ELG'):
    zgood = (data[zcolumn] > zmin) & (data[zcolumn] < zmax)
    if zcolumn == 'Z_not4clus':
        zgood &= common.goodz_infull(tp,data,zcolumn)
        zgood &= data['ZWARN'] != 999999
    print(f"# removed from quality and zcut {zmin}<{zmax}: {data[zcolumn].size - zgood.sum()}, {100 * (data[zcolumn].size - zgood.sum()) / data[zcolumn].size:.2f}%")
    return data[zgood]

def create_sysmaps(hpmaps, nest=True, columns=maps_dr9):
    logger = logging.getLogger('create_sysmaps')
    sysmaps = {}
    for prop in columns:
        if nest:
            logger.info(f"changing nest to ring ordering {prop}")
            sysmaps[prop] = hp.reorder(hpmaps[prop], inp='nest', out='ring')
        else:
            sysmaps[prop] = hpmaps[prop]
    return pd.DataFrame(sysmaps)
    
def hpixelize(nside, data, randoms, weights=None,weights_ran=None,return_mask=False, nest=False, nest2ring=False):
    if weights is None:
        data_hpmap = hpixsum(nside, data['RA'], data['DEC'], nest=nest, nest2ring=nest2ring)
    else:
        data_hpmap = hpixsum(nside, data['RA'], data['DEC'], weights=weights, nest=nest, nest2ring=nest2ring)
    
    rands_hpmap = hpixsum(nside, randoms['RA'], randoms['DEC'], weights=weights_ran,nest=nest, nest2ring=nest2ring)
    rands_mask = rands_hpmap > 0.0
    mask = rands_mask
    if return_mask:
        return data_hpmap, rands_hpmap, mask
    else:
        return data_hpmap, rands_hpmap

def make_hp(value, hpix, nside, fill_with=np.nan):
    """ A Function to create a HEALPix map
    """
    m_ = np.zeros(12*nside*nside)
    m_[:] = fill_with
    m_[hpix] = value
    
    return m_

def value2hpix(nside, ra, dec, value, statistic='mean', nest=False, nest2ring=False):
    """
    Aggregates a quantity (value) onto HEALPix with nside and ring ordering
    using `statistic`
    
    parameters
    ----------
    nside : int
    
    ra : array_like
    
    dec : array_like
    
    value : array_like
    
    statistic : str
        (optional), default is 'mean', but can work with 'min', 'max', etc
        
    nest : bool, optional
        if True, assume NESTED pixel ordering, otherwise, RING pixel ordering

    returns
    -------
    value_hp : array_like

    """
    hpix = radec2hpix(nside, ra, dec, nest=nest, nest2ring=nest2ring)
    npix = hp.nside2npix(nside)
    value_hp = binned_statistic(hpix, value, statistic=statistic,
                                bins=np.arange(0, npix+1, 1))[0]
    return value_hp

def hpixsum(nside, ra, dec, weights=None, nest=False, nest2ring=False):
    """
    Aggregates ra and dec onto HEALPix with nside and ring ordering.

    credit: Yu Feng, Ellie Kitanidis, ImagingLSS, UC Berkeley


    parameters
    ----------
    nside: int
    
    ra: array_like
        right ascention in degree.
    dec: array_like
        declination in degree.
    nest : bool, optional
        if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    

    returns
    -------
    weight_hp: array_like
            
    """
    hpix = radec2hpix(nside, ra, dec, nest=nest, nest2ring=nest2ring)
    npix = hp.nside2npix(nside)
    weight_hp = np.bincount(hpix, weights=weights, minlength=npix)
    return weight_hp

def radec2hpix(nside, ra, dec, nest=False, nest2ring=True):
    """ 
    Function transforms RA,DEC to HEALPix index
    
    parameters
    ----------
    nside : int
    
    ra : array_like
        right ascention in deg
    
    dec : array_like
        declination in deg
        
    nest : bool, optional
        if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    
    returns
    -------
    hpix : array_like
        HEALPix indices
    
    """
    hpix = hp.ang2pix(nside, np.radians(90 - dec), np.radians(ra), nest=nest)
    if (nest2ring & nest):
        print('converting hpix from nest to ring ordering')
        hpix = hp.nest2ring(nside, hpix)
    return hpix

def hpix2radec(nside, hpix, nest=False):
    """
    Function transforms HEALPix index (ring) to RA, DEC
    
    parameters 
    ----------
    nside : int
        
    hpix : array_like
        HEALPix indices
        
    nest : bool, optional
        if True, assume NESTED pixel ordering, otherwise, RING pixel ordering
    
    returns
    -------
    ra : array_like
        right ascention in deg
        
    dec : array_like
        declination in deg
        
    """
    theta, phi = hp.pixelfunc.pix2ang(nside, hpix, nest=nest)
    return np.degrees(phi), 90-np.degrees(theta)

def setup_logger():
    import logging 
    # CREATE LOGGER
    logger = logging.getLogger()
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