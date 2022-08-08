#!/usr/bin/env python
# coding: utf-8
#Initially copied from https://raw.githubusercontent.com/echaussidon/regressis/main/desi/generate_weight_DA02.py
#To be edited to allow general determination of DESI weights for target density fluctuations
import os
import shutil
import logging

import numpy as np
import pandas as pd
import healpy as hp

from regressis import PhotometricDataFrame, Regression, DR9Footprint, setup_logging
from regressis.utils import mkdir, setup_mplstyle, read_fits_to_pandas, build_healpix_map

from LSS import ssr_tools


def save_desi_data(LSS, survey, tracer, nside, dir_out, z_lim,regl=['_N','_S'],nran=18):
    """
    
    From clustering and randoms catalog build and save the healpix distribution of considered observed objects and the corresponding fracarea. 
    
    Parameters
    ----------
    LSS : str
        Path to locate where are the catalogs of galaxies / quasars.
    survey : str
        Which survey is used. Only relevant for the filename of the outputs. e.g., DA02 
    tracer : str
        Which tracer is used. e.g, BGS_ANY / LRG / ELG / QSO.
    nside : int
        Resolution of the healpix distribution map of the objects.
    dir_out : str
        Path where the ouputs will be saved.
    """
    #logger.info(f"Collect "+survey+" data for {tracer}:")

    dfs = []
    for reg in regl:
        dr = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+reg+'_clustering.dat.fits'))
        dfs.append(dr)
    data = pd.concat(dfs, ignore_index=True)
    data = data[(data['Z'] > z_lim[0]) & (data['Z'] < z_lim[1])]
    wts = data['WEIGHT_COMP'].values*data['WEIGHT_ZFAIL'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)

    #load photometric regions:
    north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=True, cut_desi=False).get_imaging_surveys()
    #logger.info("Number of pixels observed in each region:")
    #logger.info(f"        * North: {np.sum(map_data[north] > 0)} ({np.sum(map_data[north] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * South: {np.sum(map_data[south] > 0)} ({np.sum(map_data[south] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * Des:   {np.sum(map_data[des] > 0)}  ({np.sum(map_data[des] > 0)/np.sum(map_data > 0):2.2%})")

    ranl = []
    for reg in regl:
        #ran = pd.concat([read_fits_to_pandas(os.path.join(LSS, f'{tracer}zdone'+reg+'_{i}_clustering.ran.fits'), columns=['RA', 'DEC','Z']) for i in range(10)], ignore_index=True)
        for i in range(0,nran):
            ran = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+reg+'_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','Z']) 
            ranl.append(ran)
    randoms = pd.concat(ranl, ignore_index=True)
    print(len(data),len(randoms))
    # load in deg2 since we know the density of generated randoms in deg2
    map_randoms = build_healpix_map(nside, randoms['RA'].values, randoms['DEC'].values, in_deg2=True)
    # a random file is 2500 randoms per deg2
    mean = nran*2500
    #TO DO IN THE NEXT: or divide by the correct value in each pixel ! /global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits
    fracarea = map_randoms / mean
    fracarea[fracarea == 0] = np.NaN
    # remove pixels with too small fracarea
    sel = 1/fracarea > 5.0
    fracarea[sel] = np.NaN
    #logger.info(f"{np.sum(sel)} pixels are outlier on {np.sum(fracarea>0)}")

    ## savedata (without fracarea and not in degree !! --> we want just the number of object per pixel):
    filename_data = os.path.join(dir_out, f'{survey}_{tracer}_{nside}.npy')
    #logger.info(f'Save data: {filename_data}')
    np.save(filename_data, map_data)
    filename_fracarea = os.path.join(dir_out, f'{survey}_{tracer}_fracarea_{nside}.npy')
    #logger.info(f'Save corresponding fracarea: {filename_fracarea}\n')
    np.save(filename_fracarea, fracarea)

def save_desi_data_full(LSS, survey, tracer, nside, dir_out, z_lim,nran=18):
    """
    
    From clustering and randoms catalog build and save the healpix distribution of considered observed objects and the corresponding fracarea. 
    As apposed to above, use the "full" catalogs as inputs
    
    Parameters
    ----------
    LSS : str
        Path to locate where are the catalogs of galaxies / quasars.
    survey : str
        Which survey is used. Only relevant for the filename of the outputs. e.g., DA02 
    tracer : str
        Which tracer is used. e.g, BGS_ANY / LRG / ELG / QSO.
    nside : int
        Resolution of the healpix distribution map of the objects.
    dir_out : str
        Path where the ouputs will be saved.
    """
    #logger.info(f"Collect "+survey+" data for {tracer}:")

    zcol = 'Z_not4clus'
    
    cols = ['RA','DEC',zcol,'ZWARN','FRACZ_TILELOCID','DELTACHI2']
    if tracer[:3] == 'ELG':
        cols.append('o2c')
        cols.append('LOCATION_ASSIGNED')
    #data = fitsio.read(os.path.join(LSS, f'{tracer}'+'_full.dat.fits'))
    data = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+'_full.dat.fits'),columns=cols)
    if tracer == 'QSO':
        #good redshifts are currently just the ones that should have been defined in the QSO file when merged in full
        wz = data[zcol]*0 == 0
        wz &= data[zcol] != 999999
        wz &= data[zcol] != 1.e20
        wz &= data['ZWARN'] != 999999

    if tracer[:3] == 'ELG':
        wz = data['o2c'] > 0.9
        wz &= data['ZWARN']*0 == 0
        wz &= data['ZWARN'] != 999999
        print('length after oII cut '+str(len(data[wz])))
        wz &= data['LOCATION_ASSIGNED'] == 1
        print('length after also making sure location assigned '+str(len(data[wz])))

    if tracer == 'LRG':
        print('applying extra cut for LRGs')
        # Custom DELTACHI2 vs z cut from Rongpu
        wz = data['ZWARN'] == 0
        wz &= data['ZWARN']*0 == 0
        wz &= data['ZWARN'] != 999999

        selg = ssr_tools.LRG_goodz(data,zcol)
        wz &= selg

        #wz &= ff['DELTACHI2'] > dchi2
        print('length after Rongpu cut '+str(len(data[wz])))

    if tracer[:3] == 'BGS':
        wz = data['ZWARN'] == 0
        wz &= data['ZWARN']*0 == 0
        wz &= data['ZWARN'] != 999999

        print('applying extra cut for BGS')
        wz &= data['DELTACHI2'] > 40
        print('length after dchi2 cut '+str(len(data[wz])))

    wz &= data[zcol] > z_lim[0]
    wz &= data[zcol] < z_lim[1]

    data = data[wz]
    wts = 1./data['FRACZ_TILELOCID'].values#*data['WEIGHT_ZFAIL'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)

    #load photometric regions:
    north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=True, cut_desi=False).get_imaging_surveys()
    #logger.info("Number of pixels observed in each region:")
    #logger.info(f"        * North: {np.sum(map_data[north] > 0)} ({np.sum(map_data[north] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * South: {np.sum(map_data[south] > 0)} ({np.sum(map_data[south] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * Des:   {np.sum(map_data[des] > 0)}  ({np.sum(map_data[des] > 0)/np.sum(map_data > 0):2.2%})")

    ranl = []
    for i in range(0,nran):
        ran = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+'_'+str(i)+'_full.ran.fits'), columns=['RA', 'DEC']) 
        ranl.append(ran)
    randoms = pd.concat(ranl, ignore_index=True)
    print(len(data),len(randoms))
    # load in deg2 since we know the density of generated randoms in deg2
    map_randoms = build_healpix_map(nside, randoms['RA'].values, randoms['DEC'].values, in_deg2=True)
    # a random file is 2500 randoms per deg2
    mean = nran*2500
    #TO DO IN THE NEXT: or divide by the correct value in each pixel ! /global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits
    fracarea = map_randoms / mean
    fracarea[fracarea == 0] = np.NaN
    # remove pixels with too small fracarea
    sel = 1/fracarea > 5.0
    fracarea[sel] = np.NaN
    #logger.info(f"{np.sum(sel)} pixels are outlier on {np.sum(fracarea>0)}")

    ## savedata (without fracarea and not in degree !! --> we want just the number of object per pixel):
    filename_data = os.path.join(dir_out, f'{survey}_{tracer}_{nside}.npy')
    #logger.info(f'Save data: {filename_data}')
    np.save(filename_data, map_data)
    filename_fracarea = os.path.join(dir_out, f'{survey}_{tracer}_fracarea_{nside}.npy')
    #logger.info(f'Save corresponding fracarea: {filename_fracarea}\n')
    np.save(filename_fracarea, fracarea)


def _compute_weight(survey, tracer, footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, dataframe_params, max_plot_cart,pixweight_path=None, sgr_stream_path=None,feature_names=None):
    """

    Compute weight for a given tracer with a given parametrization

    Parameters
    ----------
    survey: str
        Which survey you want to use as SV3 or MAIN (for SV3 / MAIN targets) or DA02 / Y1 / etc. ...
        Useful only to load default map saved in data_dir and for the output name of the directory or filename.
    tracer: str
        Which tracer you want to use. Usefull only to load default map saved in data_dir and for
        the output name of the directory or file name.
    footprint: Footprint
        Contain all the footprint informations needed to extract the specific regions from an healpix map.
    suffix_tracer: str
        Additional suffix for tracer. Usefull only to load default map saved in data_dir and for
        the output name of the directory or filename.
    suffix_regressor : str
        Additional suffix to build regressor output directory. Useful to test on the same data different hyperparameters.
    cut_fracarea: bool
        If True create the dataframe with a selection on fracarea. In DA02, fracarea is already selected (set as nan where we don't want to use it) in the corresponding fracarea file.
    seed: int
        Fix the seed in ML algorithm for reproductibility
    param: dict
        dictionary with additional parameters to initialize :class:`PhotometricDataFrame`
    max_plot_cart: float
        Maximum value when plot map with plot_moll
    feature_names: list of str
        If not None use this list of feature during the regression otherwise use the default one.
    """
    dataframe = PhotometricDataFrame(survey, tracer, footprint, suffix_tracer, **dataframe_params)
    dataframe.set_features(pixmap=pixweight_path,sgr_stream=sgr_stream_path)
    dataframe.set_targets()
    dataframe.build(cut_fracarea=cut_fracarea)
    regression = Regression(dataframe, regressor='RF', suffix_regressor=suffix_regressor, n_jobs=40, use_kfold=True, feature_names=feature_names, compute_permutation_importance=True, overwrite=True, seed=seed, save_regressor=False)
    _ = regression.get_weight(save=True)
    regression.plot_maps_and_systematics(max_plot_cart=max_plot_cart, cut_fracarea=cut_fracarea)