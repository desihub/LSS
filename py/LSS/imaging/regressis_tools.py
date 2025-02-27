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
import fitsio

from regressis import PhotometricDataFrame, Regression, footprint, setup_logging
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
    north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False).get_imaging_surveys()
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

def save_desi_data_full(LSS, survey, tracer, nside, dir_out, z_lim,nran=18,fracthresh=5.,foot=None):
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
    zlim : list
    	contains minimum,maximum redshifts to apply
    nran : int
        number of random files to use
    fracthresh : float
    	inverse of fraction coverage of healpix pixel maximum to be considered for building weights    
    """
    #logger.info(f"Collect "+survey+" data for {tracer}:")

    zcol = 'Z_not4clus'
    
    cols = ['RA','DEC',zcol,'ZWARN','FRACZ_TILELOCID','DELTACHI2','FRAC_TLOBS_TILES','WEIGHT_ZFAIL']
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
    wts = 1./data['FRACZ_TILELOCID'].values*1./data['FRAC_TLOBS_TILES'].values*data['WEIGHT_ZFAIL'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)

    #load photometric regions:
    #north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False).get_imaging_surveys()
    if foot is None:
        foot = footprint.DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    if tracer == 'QSO':
        north, south, des = foot.get_imaging_surveys()
    else:
        north, south_ngc, south_sgc = foot.update_map(foot.data['ISNORTH']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISNGC']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISSGC'])
    #logger.info("Number of pixels observed in each region:")
    #logger.info(f"        * North: {np.sum(map_data[north] > 0)} ({np.sum(map_data[north] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * South: {np.sum(map_data[south] > 0)} ({np.sum(map_data[south] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * Des:   {np.sum(map_data[des] > 0)}  ({np.sum(map_data[des] > 0)/np.sum(map_data > 0):2.2%})")

    ranl = []
    tran = tracer
    if tracer == 'BGS_BRIGHT-21.5':
        tran = 'BGS_BRIGHT'
    for i in range(0,nran):
        ran = read_fits_to_pandas(os.path.join(LSS, f'{tran}'+'_'+str(i)+'_full.ran.fits'), columns=['RA', 'DEC']) 
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
    zw = str(z_lim[0])+'_'+str(z_lim[1])
    tracer += zw
    filename_data = os.path.join(dir_out, f'{survey}_{tracer}_{nside}.npy')
    print('saved data to '+filename_data )
    #logger.info(f'Save data: {filename_data}')
    np.save(filename_data, map_data)
    filename_fracarea = os.path.join(dir_out, f'{survey}_{tracer}_fracarea_{nside}.npy')
    #logger.info(f'Save corresponding fracarea: {filename_fracarea}\n')
    np.save(filename_fracarea, fracarea)
    print('saved fracarea to '+filename_fracarea )


def _compute_weight(survey, tracer, footprint, suffix_tracer, suffix_regressor, cut_fracarea, seed, dataframe_params, max_plot_cart,pixweight_path=None, sgr_stream_path=None,feature_names=None,pixmap_external=None,feature_names_ext=None,use_sgr=False):
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
    print('about to make dataframe')
    dataframe = PhotometricDataFrame(survey, tracer, footprint, suffix_tracer, **dataframe_params)
    print('about to set feature')
    dataframe.set_features(pixmap=pixweight_path,pixmap_external=pixmap_external,sgr_stream=sgr_stream_path,sel_columns=feature_names,sel_columns_external=feature_names_ext, use_sgr_stream=use_sgr)
    print('about to set targets')
    dataframe.set_targets()
    print('about to build')
    output_dir = dataframe.output_dataframe_dir 
    dataframe.output_dataframe_dir = None
    dataframe.build(cut_fracarea=cut_fracarea)
    dataframe.output_dataframe_dir = output_dir
    print('about to do regression')
    regression = Regression(dataframe, regressor='RF', suffix_regressor=suffix_regressor, n_jobs=40, use_kfold=True, feature_names=feature_names, compute_permutation_importance=True, overwrite=True, seed=seed, save_regressor=False)
    print('about to get weight')
    _ = regression.get_weight(save=True)
    #regression.plot_maps_and_systematics(max_plot_cart=max_plot_cart, cut_fracarea=cut_fracarea)

def get_desi_data_clus_compute_weight(LSS, survey, tracer, nside, dir_out, z_lim,dataframe_params, nran=18,data_type='data',\
fracthresh=5.,foot=None, suffix_tracer=None, suffix_regressor=None, cut_fracarea=False, seed=42, \
max_plot_cart=300,pixweight_path=None, sgr_stream_path=None,\
feature_names=None,pixmap_external=None,feature_names_ext=None,use_sgr=False,use_map_veto=''):
    """
    
    Combine steps above with no in between write-out
    
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
    zlim : list
    	contains minimum,maximum redshifts to apply
    nran : int
        number of random files to use
    fracthresh : float
    	inverse of fraction coverage of healpix pixel maximum to be considered for building weights    
    """
    #logger.info(f"Collect "+survey+" data for {tracer}:")

    zcol = 'Z'
    
    cols = ['RA','DEC',zcol,'WEIGHT_COMP','WEIGHT','WEIGHT_FKP','FRAC_TLOBS_TILES']
    LSS = LSS.replace('global','dvs_ro')
    datan = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+'_NGC_clustering.dat.fits'),columns=cols)
    datas = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+'_SGC_clustering.dat.fits'),columns=cols)
    data = pd.concat([datan,datas])

    wz = data[zcol] > z_lim[0]
    wz &= data[zcol] < z_lim[1]

    data = data[wz]
    wts =  data['WEIGHT_COMP'].values/data['FRAC_TLOBS_TILES'].values#data['WEIGHT'].values*data['WEIGHT_FKP'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)

    #load photometric regions:
    if foot is None:
        foot = footprint.DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    if tracer == 'QSO':
        north, south, des = foot.get_imaging_surveys()
    else:
        north, south_ngc, south_sgc = foot.update_map(foot.data['ISNORTH']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISNGC']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISSGC'])

    ranl = []
    tran = tracer
    if tracer == 'BGS_BRIGHT-21.5':
        tran = 'BGS_BRIGHT'
    for i in range(0,nran):
        rann = read_fits_to_pandas(os.path.join(LSS, tran+'_NGC_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT','WEIGHT_FKP'])
        rans = read_fits_to_pandas(os.path.join(LSS, tran+'_SGC_'+str(i)+'_clustering.ran.fits'), columns=['RA', 'DEC','WEIGHT','WEIGHT_FKP']) 
        ran = pd.concat([rann,rans])
        ranl.append(ran)
    randoms = pd.concat(ranl, ignore_index=True)
    print(len(data),len(randoms))
    # load in deg2 since we know the density of generated randoms in deg2, but weights mess this up
    #wts =randoms['WEIGHT'].values*randoms['WEIGHT_FKP'].values  #np.ones(len(randoms))#
    #wts /= np.mean(wts)
    #map_randoms = build_healpix_map(nside, randoms['RA'].values, randoms['DEC'].values,weights=wts, in_deg2=True)
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
    zw = str(z_lim[0])+'_'+str(z_lim[1])
    tracer += zw

    print('about to make dataframe')
    dataframe = PhotometricDataFrame(survey, tracer, foot, suffix_tracer, **dataframe_params)
#dataframe.set_features(pixmap={'North': '/global/homes/e/edmondc/Software/regressis/data/pixweight-dr9-256.fits', 'South_mid_ngc': '/global/homes/e/edmondc/Software/regressis/data/pixweight-dr9-256.fits'}, pixmap_external='/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight_external.fits')
    print('about to set feature')
    dataframe.set_features(pixmap={'North':pixweight_path+'N.fits','South':pixweight_path+'S.fits','South_mid_ngc':pixweight_path+'S.fits','Des':pixweight_path+'S.fits','South_mid_sgc':pixweight_path+'S.fits'},pixmap_external=pixmap_external,sgr_stream=sgr_stream_path,sel_columns=feature_names,sel_columns_external=feature_names_ext, use_sgr_stream=use_sgr)
    print('about to set targets')
    dataframe.set_targets(targets=map_data,fracarea=fracarea)
    print('about to build')
    output_dir = dataframe.output_dataframe_dir 
    dataframe.output_dataframe_dir = None
    dataframe.build(cut_fracarea=cut_fracarea)
    dataframe.output_dataframe_dir = output_dir
    print('about to do regression')
    regression = Regression(dataframe, regressor='RF', suffix_regressor=suffix_regressor, n_jobs=40, use_kfold=True, feature_names=feature_names, compute_permutation_importance=True, overwrite=True, seed=seed, save_regressor=False)
    print('about to get weight')
    _ = regression.get_weight(save=True)


def get_desi_data_full_compute_weight(LSS, survey, tracer, nside, dir_out, z_lim,dataframe_params, nran=18,data_type='data',\
fracthresh=5.,foot=None, suffix_tracer=None, suffix_regressor=None, cut_fracarea=False, seed=42, \
max_plot_cart=300,pixweight_path=None, sgr_stream_path=None,\
feature_names=None,pixmap_external=None,feature_names_ext=None,use_sgr=False,use_map_veto='',regressor='RF'):
    """
    
    Combine steps above with no in between write-out
    
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
    zlim : list
    	contains minimum,maximum redshifts to apply
    nran : int
        number of random files to use
    fracthresh : float
    	inverse of fraction coverage of healpix pixel maximum to be considered for building weights    
    """
    #logger.info(f"Collect "+survey+" data for {tracer}:")

    print('in weight function')
    zcol = 'Z_not4clus'
    
    cols = ['RA','DEC',zcol,'ZWARN','FRACZ_TILELOCID','FRAC_TLOBS_TILES']#,'WEIGHT_ZFAIL']
    if data_type=='data':
        cols.append('WEIGHT_ZFAIL')
    else:
        cols.append('TARGETID')
    if tracer[:3] != 'QSO':
        cols.append('DELTACHI2')
    if tracer[:3] == 'ELG':
        cols.append('o2c')
        cols.append('LOCATION_ASSIGNED')
    #data = fitsio.read(os.path.join(LSS, f'{tracer}'+'_full.dat.fits'))
    LSS = LSS.replace('global','dvs_ro')
    data = read_fits_to_pandas(os.path.join(LSS, f'{tracer}'+'_full'+use_map_veto+'.dat.fits'),columns=cols)
    print('read data')
    if data_type=='mock':
        # cut full cat based on NGC+SGC clustering cats TARGETID
        dt_n = fitsio.read(os.path.join(LSS, tracer+'_NGC_clustering.dat.fits'),columns=['TARGETID'])
        dt_s = fitsio.read(os.path.join(LSS, tracer+'_SGC_clustering.dat.fits'),columns=['TARGETID'])
        df_concat = np.concatenate([dt_n,dt_s])

        mask = np.isin(data['TARGETID'],df_concat['TARGETID'])
        dsize_before = len(data)
        data = data[mask]
        print(dsize_before,len(data),dsize_before-len(data))
        
        data['WEIGHT_ZFAIL'] = np.ones_like(data['RA'])
    
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
    wts = 1./data['FRACZ_TILELOCID'].values*1./data['FRAC_TLOBS_TILES'].values*data['WEIGHT_ZFAIL'].values
    map_data = build_healpix_map(nside, data['RA'].values, data['DEC'].values, weights=wts, in_deg2=False)
    print('made data map')
    #load photometric regions:
    #north, south, des = DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False).get_imaging_surveys()
    if foot is None:
        foot = footprint.DR9Footprint(nside, mask_lmc=False, clear_south=True, mask_around_des=False, cut_desi=False)
    if tracer == 'QSO':
        north, south, des = foot.get_imaging_surveys()
    else:
        north, south_ngc, south_sgc = foot.update_map(foot.data['ISNORTH']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISNGC']), foot.update_map(foot.data['ISSOUTH'] & foot.data['ISSGC'])
    #logger.info("Number of pixels observed in each region:")
    #logger.info(f"        * North: {np.sum(map_data[north] > 0)} ({np.sum(map_data[north] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * South: {np.sum(map_data[south] > 0)} ({np.sum(map_data[south] > 0)/np.sum(map_data > 0):2.2%})")
    #logger.info(f"        * Des:   {np.sum(map_data[des] > 0)}  ({np.sum(map_data[des] > 0)/np.sum(map_data > 0):2.2%})")

    ranl = []
    tran = tracer
    if tracer == 'BGS_BRIGHT-21.5':
        tran = 'BGS_BRIGHT'
    for i in range(0,nran):
        print('reading random '+str(i))
        ran = fitsio.read(os.path.join(LSS, f'{tran}'+'_'+str(i)+'_full'+use_map_veto+'.ran.fits'), columns=['RA', 'DEC']) #read_fits_to_pandas(os.path.join(LSS, f'{tran}'+'_'+str(i)+'_full'+use_map_veto+'.ran.fits'), columns=['RA', 'DEC']) 
        ranl.append(ran)
    randoms = pd.DataFrame(np.concatenate(ranl))#pd.concat(ranl, ignore_index=True)
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
    zw = str(z_lim[0])+'_'+str(z_lim[1])
    tracer += zw

    print('about to make dataframe')
    dataframe = PhotometricDataFrame(survey, tracer, foot, suffix_tracer, **dataframe_params)
#dataframe.set_features(pixmap={'North': '/global/homes/e/edmondc/Software/regressis/data/pixweight-dr9-256.fits', 'South_mid_ngc': '/global/homes/e/edmondc/Software/regressis/data/pixweight-dr9-256.fits'}, pixmap_external='/global/cfs/cdirs/desi/survey/catalogs/pixweight_maps_all/pixweight_external.fits')
    print('about to set feature')
    dataframe.set_features(pixmap={'North':pixweight_path+'N.fits','South':pixweight_path+'S.fits','South_mid_ngc':pixweight_path+'S.fits','Des':pixweight_path+'S.fits','South_mid_sgc':pixweight_path+'S.fits'},pixmap_external=pixmap_external,sgr_stream=sgr_stream_path,sel_columns=feature_names,sel_columns_external=feature_names_ext, use_sgr_stream=use_sgr)
    print('about to set targets')
    dataframe.set_targets(targets=map_data,fracarea=fracarea)
    print('about to build')
    output_dir = dataframe.output_dataframe_dir 
    dataframe.output_dataframe_dir = None
    dataframe.build(cut_fracarea=cut_fracarea)
    dataframe.output_dataframe_dir = output_dir
    print('about to do regression')
    use_kfold = True
    if regressor == 'Linear':
        use_kfold = False
    regression = Regression(dataframe, regressor=regressor, suffix_regressor=suffix_regressor, n_jobs=40, use_kfold=use_kfold, feature_names=feature_names, compute_permutation_importance=True, overwrite=True, seed=seed, save_regressor=False)
    print('about to get weight')
    _ = regression.get_weight(save=True)
