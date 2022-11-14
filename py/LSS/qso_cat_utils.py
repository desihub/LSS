#!/usr/bin/env python
# coding: utf-8

"""
author:  edmond chaussidon (CEA saclay)
contact: edmond.chaussidon@cea.fr

Remarks:
    * 1) log:

         If you want to desactivate the log (ie) information display in your terminal.
         Add these two lines in your script once the module is loaded.
         # import logging
         # logging.getlog("QSO_CAT_UTILS").setLevel(logging.ERROR)

    * 2) Data:

        The QSO catalog will be (for the moment) available here:
                `/global/cfs/cdirs/desi/users/edmondc/QSO_catalog/`.
        For additional information, please read the README.md file.
        Any requests or comments are welcome.

    * 3) Quality cut:

        We apply here the following quality cuts (based on VI and TS paper):
            * NO cut on ZWARN !!
            * for release <= everest: fiber_ok = (cat['COADD_FIBERSTATUS']==0)
            * for release >= fuji: fiber_ok = (cat['COADD_FIBERSTATUS']==0) | (cat['COADD_FIBERSTATUS']==8388608) | (cat['COADD_FIBERSTATUS']==16777216)
            * the two last bits appeared in fuji, can add it for previous release without any impacts.
            * definition of maskbits: https://github.com/desihub/desispec/blob/master/py/desispec/maskbits.py
"""

import sys
import os
import glob
import logging

import fitsio
import numpy as np
import pandas as pd


log = logging.getLogger("QSO_CAT_UTILS")


def desi_target_from_survey(survey):
    """ Return the survey of DESI_TARGET as a function of survey used (cmx, sv1, sv2, sv3, main)."""
    if survey == 'special':
        # to avoid error, return one of the column, SV2_DESI_TARGET should be full of 0.
        return 'SV2_DESI_TARGET'
    if survey == 'cmx':
        return 'CMX_TARGET'
    elif survey == 'sv1':
        return 'SV1_DESI_TARGET'
    elif survey == 'sv2':
        return 'SV2_DESI_TARGET'
    elif survey == 'sv3':
        return 'SV3_DESI_TARGET'
    elif survey == 'main':
        return 'DESI_TARGET'


def read_fits_to_pandas(filename, ext=1, columns=None):
    """
    Read a .fits file and convert it into a :class:`pandas.DataFrame`.
    Warning: it does not work if a column contains a list or an array.
    Parameters
    ----------
    filename : str
        Path where the .fits file is saved.
    ext : int or str
        Extension to read.
    columns : list of str
        List of columns to read. Useful to avoid to use too much memory.
    Returns :
    ---------
    data_frame : pandas.DataFrame
        Data frame containing data in the fits file.
    """
    log.info(f'Read ext: {ext} from {filename}')
    file = fitsio.FITS(filename)[ext]
    if columns is not None:
        file = file[columns]
    return pd.DataFrame(file.read().byteswap().newbyteorder())


def save_dataframe_to_fits(dataframe, filename, extname="QSO_CAT", clobber=True):
    """
    Save info from pandas dataframe in a fits file.

    Remark: Here we do not expect complex structure into dataframe (ie) only int/float/bool are expected in columns.
            We can use df.to_records().
    Args:
        dataframe (pandas dataframe): dataframe containg the all the necessary QSO info
        filename (str):  name of the fits file
        extname (str): name of the hdu in which the dataframe will be written
        clobber (bool):  overwrite the fits file defined by filename ? default=True
    Returns:
        None
    """
    if dataframe.shape[0] == 0:
        log.warning("No info to save...")
    else:
        # No complex structure, to_records() is sufficient.
        fits = fitsio.FITS(filename, 'rw', clobber=clobber)
        if clobber:
            log.warning(f'OVERWRITE the file : {filename}')
        else:
            log.warning(f'EXPAND the file : {filename}')
        fits.write(dataframe.to_records(index=False), extname=extname)
        fits.close()


def compute_RF_TS_proba(dataframe):
    """
    Compute the probabilty to be selected with the Random Forest of the Target Selection algorithm.
    It add the MW_TRANSMISSION for each band and the PROBA_RF.

    Args:
        * dataframe (pandas DataFrame): dataframe containing at least the flux, ebv, target_RA, target_DEC

    """

    def compute_MW_transmission(dataframe):
        """ TODO """
        from desiutil.dust import ext_odonnell
        from desitarget.io import desitarget_resolve_dec
        from speclite import filters  # to correct the photometry

        decamwise = filters.load_filters('decam2014-g', 'decam2014-r', 'decam2014-z', 'wise2010-W1', 'wise2010-W2')
        bassmzlswise = filters.load_filters('BASS-g', 'BASS-r', 'MzLS-z', 'wise2010-W1', 'wise2010-W2')

        north = (dataframe['TARGET_RA'] > 80) & (dataframe['TARGET_RA'] < 300) & (dataframe['TARGET_DEC'] > desitarget_resolve_dec())

        RV = 3.1
        EBV = dataframe['EBV']

        mw_transmission = np.array([10**(-0.4 * EBV[i] * RV * ext_odonnell(bassmzlswise.effective_wavelengths.value, Rv=RV)) if north[i]
                                    else 10**(-0.4 * EBV[i] * RV * ext_odonnell(decamwise.effective_wavelengths.value, Rv=RV)) for i in range(EBV.size)])

        dataframe.insert(20, 'MW_TRANSMISSION_G', mw_transmission[:, 0])
        dataframe.insert(21, 'MW_TRANSMISSION_R', mw_transmission[:, 0])
        dataframe.insert(22, 'MW_TRANSMISSION_Z', mw_transmission[:, 0])
        dataframe.insert(23, 'MW_TRANSMISSION_W1', mw_transmission[:, 0])
        dataframe.insert(24, 'MW_TRANSMISSION_W2', mw_transmission[:, 0])

    def compute_colors(dataframe):
        """ TO DO"""
        from desitarget.cuts import shift_photo_north
        from desitarget.io import desitarget_resolve_dec

        gflux = dataframe['FLUX_G'].values / dataframe['MW_TRANSMISSION_G'].values
        rflux = dataframe['FLUX_R'].values / dataframe['MW_TRANSMISSION_R'].values
        zflux = dataframe['FLUX_Z'].values / dataframe['MW_TRANSMISSION_Z'].values
        W1flux = dataframe['FLUX_W1'].values / dataframe['MW_TRANSMISSION_W1'].values
        W2flux = dataframe['FLUX_W2'].values / dataframe['MW_TRANSMISSION_W2'].values

        gflux[np.isnan(gflux) | np.isinf(gflux)] = 0.
        rflux[np.isnan(rflux) | np.isinf(rflux)] = 0.
        zflux[np.isnan(zflux) | np.isinf(zflux)] = 0.
        W1flux[np.isnan(W1flux) | np.isinf(W1flux)] = 0.
        W2flux[np.isnan(W2flux) | np.isinf(W2flux)] = 0.

        # Shift the North photometry to match the South:
        north = (dataframe['TARGET_RA'] > 80) & (dataframe['TARGET_RA'] < 300) & (dataframe['TARGET_DEC'] > desitarget_resolve_dec())
        log.info(f'shift photometry for {north.sum()} objects')
        gflux[north], rflux[north], zflux[north] = shift_photo_north(gflux[north], rflux[north], zflux[north])

        # invalid value to avoid warning with log estimation --> deal with nan
        with np.errstate(divide='ignore', invalid='ignore'):
            g = np.where(gflux > 0, 22.5 - 2.5 * np.log10(gflux), 0.)
            r = np.where(rflux > 0, 22.5 - 2.5 * np.log10(rflux), 0.)
            z = np.where(zflux > 0, 22.5 - 2.5 * np.log10(zflux), 0.)
            W1 = np.where(W1flux > 0, 22.5 - 2.5 * np.log10(W1flux), 0.)
            W2 = np.where(W2flux > 0, 22.5 - 2.5 * np.log10(W2flux), 0.)

        g[np.isnan(g) | np.isinf(g)] = 0.
        r[np.isnan(r) | np.isinf(r)] = 0.
        z[np.isnan(z) | np.isinf(z)] = 0.
        W1[np.isnan(W1) | np.isinf(W1)] = 0.
        W2[np.isnan(W2) | np.isinf(W2)] = 0.

        # Compute the colors:
        colors = np.zeros((r.size, 11))
        colors[:, 0] = g - r
        colors[:, 1] = r - z
        colors[:, 2] = g - z
        colors[:, 3] = g - W1
        colors[:, 4] = r - W1
        colors[:, 5] = z - W1
        colors[:, 6] = g - W2
        colors[:, 7] = r - W2
        colors[:, 8] = z - W2
        colors[:, 9] = W1 - W2
        colors[:, 10] = r

        return colors

    def compute_proba(dataframe):
        """ TO DO """
        import desitarget.myRF as myRF
        rf_fileName = os.path.join(os.path.dirname(myRF.__file__), 'data/rf_model_dr9_final.npz')

        attributes = compute_colors(dataframe)

        log.info('Load Random Forest: ')
        log.info('    * ' + rf_fileName)
        log.info(f'Random Forest over: {len(attributes)} objects')
        log.info('    * start RF calculation...')
        myrf = myRF.myRF(attributes, '', numberOfTrees=500, version=2)
        myrf.loadForest(rf_fileName)
        proba_rf = myrf.predict_proba()

        dataframe.insert(25, 'PROBA_RF', proba_rf)

    # add the MW_TRANSMISSION columns
    compute_MW_transmission(dataframe)
    # Add the PROBA_RF column
    compute_proba(dataframe)


def qso_catalog_maker(redrock, mgii, qn, use_old_extname_for_redrock=False, use_old_extname_for_fitsio=False, keep_all=False):
    """
    Compile the different QSO identifications to build the QSO catalog from a RR, mgII, Qn file.
    Args:
        redrock (str): redrock file with redshifts (formerly zbest)
        mgii (str): mgii file containing the mgii afterburner output
        qn (str): qn file containing the qn afterburner (with new run of RR) output
        use_old_extname_for_redrock (bool); default=False, If true use ZBEST instead REDSHIFTS for extname in redrock file?
        use_old_extname_for_fitsio (bool): default=False, For FUJI extname QN+RR is remplaced by QN_RR to avoid error with newer version of fitsio (>= 1.1.3).
                                           To use desi_qso_qn_afterburner for everest and older files please activate this flag and use ONLY fitsio = 1.1.2.
                                           For daily production, this modification was done in: 18/01/2022.
        keep_all (bool): if True return all the targets. if False return only targets which are selected as QSO.
    Returns:
        QSO_cat (pandas dataframe): Dataframe containing all the information
    """
    from functools import reduce

    # selection of which column will be in the final QSO_cat:
    columns_zbest = ['TARGETID', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE']  # , 'SUBTYPE', 'DELTACHI2', 'CHI2']
    # remark: check if the name exist before selecting them
    columns_fibermap = ['TARGETID', 'TARGET_RA', 'TARGET_DEC', 'LOCATION', 'MORPHTYPE', 'COADD_FIBERSTATUS', 'COADD_NUMEXP', 'COADD_EXPTIME',
                        'EBV', 'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
                        'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2', 'MASKBITS',
                        'CMX_TARGET', 'SV1_DESI_TARGET', 'SV2_DESI_TARGET', 'SV3_DESI_TARGET', 'DESI_TARGET',
                        'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET', 'SCND_TARGET']

    columns_tsnr2 = ['TARGETID', 'TSNR2_QSO', 'TSNR2_LYA']
    # for david
    # columns_tsnr2 = ['TARGETID', 'TSNR2_GPBDARK_B', 'TSNR2_ELG_B', 'TSNR2_GPBBRIGHT_B', 'TSNR2_LYA_B', 'TSNR2_BGS_B', 'TSNR2_GPBBACKUP_B', 'TSNR2_QSO_B', 'TSNR2_LRG_B', 'TSNR2_GPBDARK_R', 'TSNR2_ELG_R', 'TSNR2_GPBBRIGHT_R', 'TSNR2_LYA_R', 'TSNR2_BGS_R', 'TSNR2_GPBBACKUP_R', 'TSNR2_QSO_R', 'TSNR2_LRG_R', 'TSNR2_GPBDARK_Z', 'TSNR2_ELG_Z', 'TSNR2_GPBBRIGHT_Z', 'TSNR2_LYA_Z', 'TSNR2_BGS_Z', 'TSNR2_GPBBACKUP_Z', 'TSNR2_QSO_Z', 'TSNR2_LRG_Z','TSNR2_GPBDARK', 'TSNR2_ELG', 'TSNR2_GPBBRIGHT', 'TSNR2_LYA', 'TSNR2_BGS', 'TSNR2_GPBBACKUP', 'TSNR2_QSO', 'TSNR2_LRG']

    columns_mgii = ['TARGETID', 'IS_QSO_MGII', 'DELTA_CHI2', 'A', 'SIGMA', 'B', 'VAR_A', 'VAR_SIGMA', 'VAR_B']
    columns_mgii_rename = {"DELTA_CHI2": "DELTA_CHI2_MGII", "A": "A_MGII", "SIGMA": "SIGMA_MGII", "B": "B_MGII", "VAR_A": "VAR_A_MGII", "VAR_SIGMA": "VAR_SIGMA_MGII", "VAR_B": "VAR_B_MGII"}

    columns_qn = ['TARGETID', 'Z_NEW', 'ZERR_NEW', 'Z_RR', 'Z_QN', 'IS_QSO_QN_NEW_RR',
                  'C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha',
                  'Z_LYA', 'Z_CIV', 'Z_CIII', 'Z_MgII', 'Z_Hbeta', 'Z_Halpha']

    # load data:
    zbest = read_fits_to_pandas(redrock, ext='ZBEST' if use_old_extname_for_redrock else 'REDSHIFTS', columns=columns_zbest)
    fibermap = read_fits_to_pandas(redrock, ext='FIBERMAP', columns=[name for name in columns_fibermap if name in fitsio.read(redrock, ext='FIBERMAP', rows=[0]).dtype.names])
    tsnr2 = read_fits_to_pandas(redrock, ext='TSNR2', columns=columns_tsnr2)
    mgii = read_fits_to_pandas(mgii, ext='MGII', columns=columns_mgii).rename(columns=columns_mgii_rename)
    qn = read_fits_to_pandas(qn, ext='QN+RR' if use_old_extname_for_fitsio else 'QN_RR', columns=columns_qn)

    # add DESI_TARGET column to avoid error of conversion when concatenate the different files with pd.concat() which fills with NaN columns that do not exit in a DataFrame.
    # convert int64 to float 64 --> destructive tranformation !!
    for DESI_TARGET in ['CMX_TARGET', 'SV1_DESI_TARGET', 'SV2_DESI_TARGET', 'SV3_DESI_TARGET', 'DESI_TARGET', 'SV1_SCND_TARGET', 'SV2_SCND_TARGET', 'SV3_SCND_TARGET', 'SCND_TARGET']:
        if not(DESI_TARGET in fibermap.columns):
            fibermap[DESI_TARGET] = np.zeros(fibermap['TARGETID'].size, dtype=np.int64)

    # QN afterburner is run with a threshold 0.5. With VI, we choose 0.95 as final threshold.
    # &= since IS_QSO_QN_NEW_RR contains only QSO for QN which are not QSO for RR.
    log.info('Increase the QN threshold selection from 0.5 to 0.95.')
    qn['IS_QSO_QN'] = np.max(np.array([qn[name] for name in ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']]), axis=0) > 0.95
    qn['IS_QSO_QN_NEW_RR'] &= qn['IS_QSO_QN']

    log.info('Merge on TARGETID all the info into a singe dataframe.')
    QSO_cat = reduce(lambda left, right: pd.merge(left, right, on=['TARGETID'], how='outer'), [zbest, fibermap, tsnr2, mgii, qn])

    # Add BITMASK:
    QSO_cat['QSO_MASKBITS'] = np.zeros(QSO_cat.shape[0], dtype='i')
    log.info('Selection with SPECTYPE.')
    QSO_cat.loc[QSO_cat['SPECTYPE'] == 'QSO', 'QSO_MASKBITS'] += 2**1
    log.info('Selection with MgII.')
    QSO_cat.loc[QSO_cat['IS_QSO_MGII'], 'QSO_MASKBITS'] += 2**2
    log.info('Selection with QN (add new z from Redrock with QN prior where it is relevant).')
    QSO_cat.loc[QSO_cat['IS_QSO_QN'], 'QSO_MASKBITS'] += 2**3
    QSO_cat.loc[QSO_cat['IS_QSO_QN_NEW_RR'], 'QSO_MASKBITS'] += 2**4
    QSO_cat.loc[QSO_cat['IS_QSO_QN_NEW_RR'], 'Z'] = QSO_cat['Z_NEW'][QSO_cat['IS_QSO_QN_NEW_RR']].values
    QSO_cat.loc[QSO_cat['IS_QSO_QN_NEW_RR'], 'ZERR'] = QSO_cat['ZERR_NEW'][QSO_cat['IS_QSO_QN_NEW_RR']].values

    # Add quality cuts: no cut on zwarn, cut on fiberstatus
    QSO_cat.loc[~((QSO_cat['COADD_FIBERSTATUS'] == 0) | (QSO_cat['COADD_FIBERSTATUS'] == 8388608) | (QSO_cat['COADD_FIBERSTATUS'] == 16777216)), 'QSO_MASKBITS'] = 0

    # remove useless columns:
    QSO_cat.drop(columns=['IS_QSO_MGII', 'IS_QSO_QN', 'IS_QSO_QN_NEW_RR', 'Z_NEW', 'ZERR_NEW'], inplace=True)

    # Correct bump at z~3.7
    sel_pb_redshift = (((QSO_cat['Z'] > 3.65) & (QSO_cat['Z'] < 3.9)) | ((QSO_cat['Z'] > 5.15) & (QSO_cat['Z'] < 5.35))) & ((QSO_cat['C_LYA'] < 0.95) | (QSO_cat['C_CIV'] < 0.95))
    log.info(f'Remove bump at z~3.7: exclude {sel_pb_redshift.sum()} QSOs.')
    QSO_cat.loc[sel_pb_redshift, 'QSO_MASKBITS'] = 0

    if keep_all:
        log.info('Return all the targets without any cut on QSO selection.')
        return QSO_cat
    else:
        QSO_cat = QSO_cat[QSO_cat['QSO_MASKBITS'] > 0]
        if QSO_cat.shape[0] == 0:
            log.info('No QSO found...')
        else:
            log.info(f"Final selection gives: {QSO_cat.shape[0]} QSO !")
        return QSO_cat


def qso_catalog_for_a_tile(path_to_tile, tile, last_night, survey, program):
    """
    Build the QSO catalog for the tile using the last_night. It is relevant for cumulative directory.
    This function is usefull to be called in pool.starmap under multiprocessing.

    Args:
        path_to_tile (str): Where the tiles are.
        tile (str): which tile do you want to treat.
        last_night (str): corresponding last night to tile
        survey (str): sv3/main ... only to add information to the catalog
        program (str): dark/bright/backup only to add information to the catalog

    Return:
        QSO_cat (DataFrame): pandas DataFrame containing the concatenation of run_catalog_maker in each available petal
    """

    def run_catalog_maker(path_to_tile, tile, night, petal, survey, program):
        """Run qso_catalog_maker in the considered tile-last_night-petal. If one file does not exist it return a void DataFrame."""
        redrock = os.path.join(path_to_tile, tile, night, f"redrock-{petal}-{tile}-thru{night}.fits")
        mgii_afterburner = os.path.join(path_to_tile, tile, night, f"qso_mgii-{petal}-{tile}-thru{night}.fits")
        qn_afterburner = os.path.join(path_to_tile, tile, night, f"qso_qn-{petal}-{tile}-thru{night}.fits")

        if os.path.isfile(redrock):
            if os.path.isfile(mgii_afterburner) & os.path.isfile(qn_afterburner):
                qso_cat = qso_catalog_maker(redrock, mgii_afterburner, qn_afterburner)
                qso_cat['TILEID'] = int(tile)
                qso_cat['LASTNIGHT'] = int(night)
                qso_cat['PETAL_LOC'] = int(petal)
                qso_cat['SURVEY'] = survey
                qso_cat['PROGRAM'] = program
            else:
                log.error(f'There is a problem with: {mgii_afterburner} | {qn_afterburner}')
                qso_cat = pd.DataFrame()
        else:
            # this can happen, it miss some petal.
            log.info(f'Redrock file does not exist: {redrock}')
            qso_cat = pd.DataFrame()
        return qso_cat

    return pd.concat([run_catalog_maker(path_to_tile, tile, last_night, petal, survey, program) for petal in range(10)], ignore_index=True)


def build_qso_catalog_from_tiles(redux='/global/cfs/cdirs/desi/spectro/redux/', release='fuji', dir_output='', npool=20, tiles_to_use=None):
    """
    Build the QSO catalog from the healpix directory.

    Warning: no retro compatibility for release <= everest (extname has changed --> the option can be added since it exists in qso_catalog_maker)

    Args:
        * redux (str): path where is saved the spectroscopic data.
        * release (str): which release do you want to use (everest, fuji, guadalupe, ect...).
        * dir_output (str): directory where the QSO catalog will be saved.
        * npool (int): nbr of workers used for the parallelisation.
        * tiles_to_use (list of str): Build the catalog only on this list of tiles. Default=None, use all the tiles collected from tiles-{release}.fits file.
    """
    import multiprocessing
    from itertools import repeat

    # remove desimodule log
    os.environ["DESI_LOGLEVEL"] = "ERROR"

    # Data directory
    DIR = os.path.join(redux, release, 'tiles', 'cumulative')

    # load tiles info:
    tile_info = fitsio.FITS(os.path.join(redux, release, f'tiles-{release}.fits'))[1][['TILEID', 'LASTNIGHT', 'SURVEY', 'PROGRAM']]

    tiles = np.array(tile_info['TILEID'][:], dtype='str')
    last_night = np.array(tile_info['LASTNIGHT'][:], dtype='str')
    survey = np.array(tile_info['SURVEY'][:], dtype='str')
    program = np.array(tile_info['PROGRAM'][:], dtype='str')

    if tiles_to_use is not None:
        sel = np.isin(tiles, tiles_to_use)
        tiles, last_night, survey, program = tiles[sel], last_night[sel], survey[sel], program[sel]

    log.info(f'There are {tiles.size} tiles to treat with npool={npool}')
    logging.getLogger("QSO_CAT_UTILS").setLevel(logging.ERROR)
    with multiprocessing.Pool(npool) as pool:
        arguments = zip(repeat(DIR), tiles, last_night, survey, program)
        QSO_cat = pd.concat(pool.starmap(qso_catalog_for_a_tile, arguments), ignore_index=True)
    logging.getLogger("QSO_CAT_UTILS").setLevel(logging.INFO)

    log.info('Compute the TS probas...')
    compute_RF_TS_proba(QSO_cat)

    save_dataframe_to_fits(QSO_cat, os.path.join(dir_output, f'QSO_cat_{release}_cumulative.fits'))


def qso_catalog_for_a_pixel(path_to_pix, pre_pix, pixel, survey, program, keep_all=False):
    """
    Build the QSO catalog for the tile using the last_night. It is relevant for cumulative directory.
    This function is usefull to be called in pool.starmap under multiprocessing.

    Args:
        * path_to_pix (str): Where the pixels are.
        * pre_pix (str): which pre_pix in healpix directory do you want to use.
        * pixel (str): which pixel do you want to use.
        * survey (str): which TS do you want to use (sv1/sv3/main)
        * program (str): either dark / bright / backup
        * keep_all (bool): if True return all the targets. if False return only targets which are selected as QSO.

    Return:
        QSO_cat (DataFrame): pandas DataFrame containing the QSO_catalog for the considered pixel.
    """
    redrock = os.path.join(path_to_pix, str(pre_pix), str(pixel), f"redrock-{survey}-{program}-{pixel}.fits")
    mgii_afterburner = os.path.join(path_to_pix, str(pre_pix), str(pixel), f"qso_mgii-{survey}-{program}-{pixel}.fits")
    qn_afterburner = os.path.join(path_to_pix, str(pre_pix), str(pixel), f"qso_qn-{survey}-{program}-{pixel}.fits")

    if os.path.isfile(redrock):
        if os.path.isfile(mgii_afterburner) & os.path.isfile(qn_afterburner):
            qso_cat = qso_catalog_maker(redrock, mgii_afterburner, qn_afterburner, keep_all=keep_all)
            qso_cat['HPXPIXEL'] = int(pixel)
            qso_cat['SURVEY'] = survey
            qso_cat['PROGRAM'] = program
        else:
            log.error(f'There is a problem with: {mgii_afterburner} | {qn_afterburner}')
            qso_cat = pd.DataFrame()
    else:
        # It is not expected, the pixel should not be created if no targets are inside ?
        log.error(f'Redrock file does not exist: {redrock}')
        qso_cat = pd.DataFrame()
    return qso_cat


def build_qso_catalog_from_healpix(redux='/global/cfs/cdirs/desi/spectro/redux/', release='fuji', survey='sv3', program='dark', dir_output='', npool=20, keep_qso_targets=True, keep_all=False):
    """
    Build the QSO catalog from the healpix directory.

    Warning: no retro compatibility for release <= everest (extname has changed --> the option can be added since it exists in qso_catalog_maker)

    Args:
        * redux (str): path where is saved the spectroscopic data.
        * release (str): which release do you want to use (everest, fuji, guadalupe, ect...).
        * survey (str): which survey of the target selection (sv1, sv3, main).
        * program (str) : either dark / bright or backup program.
        * dir_output (str): directory where the QSO catalog will be saved.
        * npool (int): nbr of workers used for the parallelisation.
        * keep_qso_targets (bool): if True save only QSO targets. default=True
        * keep_all (bool): if True return all the targets. if False return only targets which are selected as QSO. default=False
    """
    import multiprocessing
    from itertools import repeat

    # remove desimodule log
    os.environ["DESI_LOGLEVEL"] = "ERROR"

    # Data directory
    DIR = os.path.join(redux, release, 'healpix', survey, program)

    # Collect the pre-pixel and pixel number
    pre_pix_list = np.sort([os.path.basename(path) for path in glob.glob(os.path.join(DIR, "*"))])
    pre_pix_list_long, pixel_list = [], []
    for pre_pix in pre_pix_list:
        pixel_list_tmp = [os.path.basename(path) for path in glob.glob(os.path.join(DIR, pre_pix, "*"))]
        pre_pix_list_long += [pre_pix] * len(pixel_list_tmp)
        pixel_list += pixel_list_tmp

    log.info(f'There are {len(pixel_list)} pixels to treat with npool={npool}')
    logging.getLogger("QSO_CAT_UTILS").setLevel(logging.ERROR)
    with multiprocessing.Pool(npool) as pool:
        arguments = zip(repeat(DIR), pre_pix_list_long, pixel_list, repeat(survey), repeat(program), repeat(keep_all))
        QSO_cat = pd.concat(pool.starmap(qso_catalog_for_a_pixel, arguments), ignore_index=True)
    logging.getLogger("QSO_CAT_UTILS").setLevel(logging.INFO)

    if not keep_all:
        # to save computational time
        log.info('Compute the TS probas...')
        compute_RF_TS_proba(QSO_cat)

    if keep_qso_targets:
        log.info('Keep only qso targets...')
        save_dataframe_to_fits(QSO_cat.iloc[QSO_cat[desi_target_from_survey(survey)].values & 2**2 != 0], os.path.join(dir_output, f'QSO_cat_{release}_{survey}_{program}_healpix_only_qso_targets.fits'))

    suffix = ''
    if keep_all:
        suffix = '_all_targets'
    save_dataframe_to_fits(QSO_cat, os.path.join(dir_output, f'QSO_cat_{release}_{survey}_{program}_healpix{suffix}.fits'))


def afterburner_is_missing_in_tiles(redux='/global/cfs/cdirs/desi/spectro/redux/', release='fuji', outdir=''):
    """
    Goes throught all the directory of tiles and check if afterburner files exist when the associated redrock file exist also.
    If files are missing, they are saved in .txt file.
    Args:
        * redux (str): path where is saved the spectroscopic data.
        * release (str): which release do you want to check.
        * outdir (str): path where the .txt output will be saved in case if it lacks some afterburner files.
    """
    import tqdm

    dir_list = ['pernight', 'perexp', 'cumulative']
    suff_dir_list = ['', 'exp', 'thru']

    for dir, suff_dir in zip(dir_list, suff_dir_list):
        DIR = os.path.join(redux, release, 'tiles', dir)
        tiles = np.sort([os.path.basename(path) for path in glob.glob(os.path.join(DIR, '*'))])
        log.info(f'Inspection of {tiles.size} tiles in {DIR}...')
        pb_qn, pb_mgII = [], []

        for tile in tqdm.tqdm(tiles):
            nights = np.sort([os.path.basename(path) for path in glob.glob(os.path.join(DIR, tile, '*'))])
            for night in nights:
                for petal in range(10):
                    if os.path.isfile(os.path.join(DIR, tile, night, f"redrock-{petal}-{tile}-{suff_dir}{night}.fits")):
                        if not (os.path.isfile(os.path.join(DIR, tile, night, f"qso_qn-{petal}-{tile}-{suff_dir}{night}.fits")) or
                                os.path.isfile(os.path.join(DIR, tile, night, f"qso_qn-{petal}-{tile}-{suff_dir}{night}.notargets")) or
                                os.path.isfile(os.path.join(DIR, tile, night, f"qso_qn-{petal}-{tile}-{suff_dir}{night}.misscamera"))):
                            pb_qn += [[int(tile), int(night), int(petal)]]
                        if not (os.path.isfile(os.path.join(DIR, tile, night, f"qso_mgii-{petal}-{tile}-{suff_dir}{night}.fits")) or
                                os.path.isfile(os.path.join(DIR, tile, night, f"qso_mgii-{petal}-{tile}-{suff_dir}{night}.notargets")) or
                                os.path.isfile(os.path.join(DIR, tile, night, f"qso_mgii-{petal}-{tile}-{suff_dir}{night}.misscamera"))):
                            pb_mgII += [[int(tile), int(night), int(petal)]]

        log.info(f'Under the directory {DIR} it lacks:')
        log.info(f'    * {len(pb_qn)} QN files')
        log.info(f'    * {len(pb_mgII)} MgII files')
        if len(pb_qn) > 0:
            np.savetxt(os.path.join(outdir, f'pb_qn_{release}_{dir}.txt'), pb_qn, fmt='%d')
        if len(pb_mgII) > 0:
            np.savetxt(os.path.join(outdir, f'pb_mgII_{release}_{dir}.txt'), pb_mgII, fmt='%d')


def afterburner_is_missing_in_healpix(redux='/global/cfs/cdirs/desi/spectro/redux/', release='fuji', outdir=''):
    """
    Goes throught all the directory of healpix and check if afterburner files exist when the associated redrock file exist also.
    If files are missing, they are saved in .txt file.
    Args:
        * redux (str): path where is saved the spectroscopic data.
        * release (str): which release do you want to check.
        * outdir (str): path where the .txt output will be saved in case if it lacks some afterburner files.
    """
    import tqdm

    DIR = os.path.join(redux, release, 'healpix')

    # sv1 / sv3 / main
    survey_list = [os.path.basename(path) for path in glob.glob(os.path.join(DIR, '*'))]
    for survey in survey_list:

        # dark / bright / backup
        program_list = [os.path.basename(path) for path in glob.glob(os.path.join(DIR, survey, '*'))]
        for program in program_list:
            log.info(f'Inspection of {os.path.join(DIR, survey, program)}...')
            pb_qn, pb_mgII = [], []

            # collect the huge pixels directory
            healpix_huge_pixels = np.sort([os.path.basename(path) for path in glob.glob(os.path.join(DIR, survey, program, '*'))])
            for num in tqdm.tqdm(healpix_huge_pixels):
                pix_numbers = np.sort([os.path.basename(path) for path in glob.glob(os.path.join(DIR, survey, program, num, '*'))])
                for pix in pix_numbers:
                    if os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"redrock-{survey}-{program}-{pix}.fits")):
                        if not (os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_qn-{survey}-{program}-{pix}.fits")) or
                                os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_qn-{survey}-{program}-{pix}.notargets")) or
                                os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_qn-{survey}-{program}-{pix}.misscamera"))):
                            pb_qn += [[int(num), int(pix)]]
                        if not (os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_mgii-{survey}-{program}-{pix}.fits")) or
                                os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_mgii-{survey}-{program}-{pix}.notargets")) or
                                os.path.isfile(os.path.join(DIR, survey, program, num, pix, f"qso_mgii-{survey}-{program}-{pix}.misscamera"))):
                            pb_mgII += [[int(num), int(pix)]]

            log.info(f'Under the directory {os.path.join(DIR, survey, program)} it lacks:')
            log.info(f'    * {len(pb_qn)} QN files')
            log.info(f'    * {len(pb_mgII)} MgII files')
            if len(pb_qn) > 0:
                np.savetxt(os.path.join(outdir, f'pb_qn_healpix_{survey}_{program}.txt'), pb_qn, fmt='%d')
            if len(pb_mgII) > 0:
                np.savetxt(os.path.join(outdir, f'pb_mgII_healpix_{survey}_{program}.txt'), pb_mgII, fmt='%d')


if __name__ == '__main__':
    """ Please do not execute: This just examples of how you can use the following function in your proper code. """

    if sys.argv[1] == 'inspect_afterburners':
        """ Simple inspection of qso afterburner in fuji and guadalupe release. """
        log.info('Test the existence of QSO afterburners in Fuji and in Guadalupe.\nWe check is mgII and qn files were produced if the corresponding redrock file exits.')

        log.info('Inspect Fuji...')
        afterburner_is_missing_in_tiles(release='fuji')
        afterburner_is_missing_in_healpix(release='fuji')

        log.info('Inspect Guadalupe...')
        afterburner_is_missing_in_tiles(release='guadalupe')
        afterburner_is_missing_in_healpix(release='guadalupe')

    if sys.argv[1] == 'qso_catalog_from_files':
        log.info('Build QSO catalog for cumulative fuji tile 107 petal 3')
        redrock = '/global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/107/20210428/redrock-3-107-thru20210428.fits'
        mgII = '/global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/107/20210428/qso_mgII-3-107-thru20210428.fits'
        qn = '/global/cfs/cdirs/desi/spectro/redux/fuji/tiles/cumulative/107/20210428/qso_qn-3-107-thru20210428.fits'

        QSO_cat = qso_catalog_maker(redrock, mgII, qn)

    if sys.argv[1] == 'build_qso_catalog':
        """ Simple example of how to build the QSO catalog from healpix or cumulative directory"""

        redux = '/global/cfs/cdirs/desi/spectro/redux/'

        log.info('Build QSO catalog from cumulative directory for guadalupe release:')
        build_qso_catalog_from_tiles(redux=redux, release='guadalupe', dir_output='')

        log.info('Build QSO catalog from healpix directory for guadalupe release:')
        build_qso_catalog_from_healpix(redux=redux, release='guadalupe', survey='main', program='dark', dir_output='')
