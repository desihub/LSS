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


def save_desi_data(LSS, survey, tracer, nside, dir_out, z_lim):
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

    data = read_fits_to_pandas(os.path.join(LSS, f'{tracer}zdone_clustering.dat.fits'))
    data = data[(data['Z'] > z_lim[0]) & (data['Z'] < z_lim[1])]
    wts = data['WEIGHT_COMP'].values*data['WEIGHT_ZFAIL'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)

    #load photometric regions:
    north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=True, cut_desi=False).get_imaging_surveys()
    #logger.info("Number of pixels observed in each region:")
    #logger.info(f"        * North: {np.sum(map_data[north] > 0)} ({np.sum(map_data[north] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * South: {np.sum(map_data[south] > 0)} ({np.sum(map_data[south] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * Des:   {np.sum(map_data[des] > 0)}  ({np.sum(map_data[des] > 0)/np.sum(map_data > 0):2.2%})")

    randoms = pd.concat([read_fits_to_pandas(os.path.join(LSS, f'{tracer}zdone_{i}_clustering.ran.fits'), columns=['RA', 'DEC','Z']) for i in range(10)], ignore_index=True)
    # load in deg2 since we know the density of generated randoms in deg2
    map_randoms = build_healpix_map(nside, randoms['RA'].values, randoms['DEC'].values, in_deg2=True)
    # a random file is 2500 randoms per deg2
    mean = 10*2500
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


def _compute_weight(survey, tracer, footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, dataframe_params, max_plot_cart,pixweight_path=None, feature_names=None):
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
    dataframe.set_features(pixmap=pixweight_path)
    dataframe.set_targets()
    dataframe.build(cut_fracarea=cut_fracarea)
    regression = Regression(dataframe, regressor='RF', suffix_regressor=suffix_regressor, n_jobs=40, use_kfold=True, feature_names=feature_names, compute_permutation_importance=True, overwrite=True, seed=seed, save_regressor=False)
    _ = regression.get_weight(save=True)
    regression.plot_maps_and_systematics(max_plot_cart=max_plot_cart, cut_fracarea=cut_fracarea)