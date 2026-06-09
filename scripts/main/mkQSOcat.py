#!/usr/bin/env python
# coding: utf-8

"""
Build the two main-survey dark-time QSO catalogs (healpix-based) directly from the
zcatalog ``zpix-main-dark.fits`` and ``zpix-main-dark-extra.fits`` files.

Two catalogs are produced:
    * "QSO targets"     : QSO-targeted objects confirmed as QSOs, i.e.
                          (DESI_TARGET QSO bit) & GOOD_Z_QSO & ~bad_qso. This matches
                          the "only_qso_targets" file of LSS.qso_cat_utils (which keeps
                          QSO targets with QSO_MASKBITS > 0, after the bad_qso cut).
        --> QSO_cat_{release}_main_dark_healpix_only_qso_targets_{version}.fits
    * "all QSOs"        : every object identified as a QSO, defined as the union of
                          four mutually-exclusive subsamples (QSO / ELG / WISE-VAR /
                          BGS targets that pass the relevant QSO identification).
        --> QSO_cat_{release}_main_dark_healpix_{version}.fits

Redshift convention (matching QSOcat_dev.ipynb):
    * QSO-target subsample : Z = Z_QSO (from the extra catalog), used where GOOD_Z_QSO.
    * ELG / WISE-VAR / BGS : Z = "zelg", i.e. the redrock Z with Z_NEW substituted in
                             wherever IS_QSO_QN_NEW_RR is set.
    * ZERR / ZWARN are taken straight from the extra (redrock) catalog, no substitution.

Output columns: only those present (or trivially derivable) in the two input
catalogs are written. 
"""

import os
import sys
import argparse
import logging

import numpy as np
import fitsio
import healpy as hp
from astropy.table import Table

from desitarget import targetmask
from LSS import common_tools as common

# ----------------------------------------------------------------------------- #
# logging
logname = 'mkQSOcat'
logger = logging.getLogger(logname)
logger.setLevel(logging.INFO)
ch = logging.StreamHandler()
ch.setLevel(logging.INFO)
ch.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
logger.addHandler(ch)

# ----------------------------------------------------------------------------- #
# constants

# columns to read from zpix-main-dark.fits (only what is needed for selection + output)
ZCAT_COLS = ['TARGETID',
             'COADD_FIBERSTATUS',
             'TARGET_RA',
             'TARGET_DEC',
             'OBJTYPE',
             'DESI_TARGET',
             'SCND_TARGET',
             'COADD_NUMEXP',
             'COADD_EXPTIME']

# columns to read from zpix-main-dark-extra.fits
EXTRA_COLS = ['TARGETID',
              'Z',
              'ZERR',
              'ZWARN',
              'SPECTYPE',
              'TSNR2_LYA',
              'TSNR2_QSO',
              'TSNR2_ELG',
              'TSNR2_LRG',
              'C_LYA',
              'C_CIV',
              'C_CIII',
              'C_MgII',
              'C_Hbeta',
              'C_Halpha',
              'Z_NEW',
              'GOOD_Z_QSO',
              'Z_QSO',
              'IS_QSO_MGII',
              'IS_QSO_QN_NEW_RR']

# afterburner confidence columns used to build the QN selections
C_COLS = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

# final output column order (intersection of the requested schema with what the two
# input catalogs can provide)
OUTPUT_COLS = ['TARGETID', 'Z', 'ZERR', 'ZWARN', 'SPECTYPE', 'COADD_FIBERSTATUS',
               'TARGET_RA', 'TARGET_DEC', 'OBJTYPE', 'DESI_TARGET', 'SCND_TARGET',
               'COADD_NUMEXP', 'COADD_EXPTIME', 'TSNR2_LYA', 'TSNR2_QSO',
               'C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha',
               'QSO_MASKBITS', 'HPXPIXEL', 'SURVEY', 'PROGRAM', 'TSNR2_ELG', 'TSNR2_LRG']

# healpix resolution used for the main-survey healpix grouping
HPX_NSIDE = 64

SURVEY = 'main'
PROGRAM = 'dark'

ZCAT_FN = 'zpix-main-dark.fits'
EXTRA_FN = 'zpix-main-dark-extra.fits'

EXTNAME = 'QSO_CAT'


# ----------------------------------------------------------------------------- #
# helpers

def dvs_ro(path):
    """Return the dvs_ro (read-only) NERSC mount equivalent of a /global path."""
    return path.replace('/global/', '/dvs_ro/', 1)


def get_catdir(release):
    """Return the directory holding the input zcatalog files for a given release."""
    if release == 'loa':
        return '/global/cfs/cdirs/desicollab/users/rongpu/data/redux/loa/zcatalog/v2/main/'
    # matterhorn / nevis
    return f'/global/cfs/cdirs/desi/spectro/redux/{release}/zcatalog/v2/main/'


def as_str(arr):
    """Return a unicode view of a (possibly byte-string) array for safe comparisons."""
    arr = np.asarray(arr)
    if arr.dtype.kind == 'S':
        return np.char.decode(arr)
    if arr.dtype.kind != 'U':
        return arr.astype(str)
    return arr


def get_colnames(filename):
    """Return the list of column names available in the first HDU of a fits file."""
    with fitsio.FITS(filename) as f:
        return f[1].get_colnames()


def align_to(reference_tids, tids, data):
    """
    Reorder ``data`` (a structured array indexed like ``tids``) so that its rows
    match the order of ``reference_tids``. TARGETID is assumed unique and present
    in both. Falls back to a no-op if the two TARGETID arrays are already identical.
    """
    if np.array_equal(reference_tids, tids):
        return data
    logger.info('TARGETID order differs between the two inputs; matching on TARGETID')
    sorter = np.argsort(tids)
    pos = np.searchsorted(tids, reference_tids, sorter=sorter)
    idx = sorter[pos]
    if not np.array_equal(tids[idx], reference_tids):
        raise ValueError('TARGETID sets differ between the two input catalogs; cannot align rows.')
    return data[idx]


def determine_version(outdir, release, explicit=None):
    """
    Determine the catalog version string. If ``explicit`` is given, use it as-is.
    Otherwise default to 'v0' and increment ('v1', 'v2', ...) until neither of the
    two output files already exists for that version.
    """
    if explicit is not None:
        return explicit
    n = 0
    while True:
        ver = f'v{n}'
        f_targets = os.path.join(outdir, f'QSO_cat_{release}_main_dark_healpix_only_qso_targets_{ver}.fits')
        f_qsos = os.path.join(outdir, f'QSO_cat_{release}_main_dark_healpix_{ver}.fits')
        if not os.path.exists(f_targets) and not os.path.exists(f_qsos):
            return ver
        n += 1


# ----------------------------------------------------------------------------- #
# main

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--release', required=True, choices=['loa', 'matterhorn', 'nevis'],
                        help='spectroscopic release used to locate the input zcatalogs')
    parser.add_argument('--outdir', default='.',
                        help='directory to write the output QSO catalogs (default: current directory)')
    parser.add_argument('--version', default=None,
                        help="catalog version string (e.g. 'v0'); default auto-increments from 'v0'")
    args = parser.parse_args()
    logger.info(f'arguments: {args}')

    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        logger.info(f'creating output directory {outdir}')
        os.makedirs(outdir, exist_ok=True)

    catdir = get_catdir(args.release)
    zcat_path = dvs_ro(os.path.join(catdir, ZCAT_FN))
    extra_path = dvs_ro(os.path.join(catdir, EXTRA_FN))
    logger.info(f'reading inputs from {catdir}')

    # ----- read zcat (basic catalog) ----------------------------------------- #
    zcat_avail = get_colnames(zcat_path)
    has_healpix = 'HEALPIX' in zcat_avail
    read_cols = list(ZCAT_COLS)
    if has_healpix:
        read_cols.append('HEALPIX')

    logger.info(f'reading {len(read_cols)} columns from {ZCAT_FN}')
    zcat = fitsio.read(zcat_path, columns=read_cols)
    logger.info(f'read {len(zcat)} rows from {ZCAT_FN}')

    # ----- read zextra (extra catalog) --------------------------------------- #
    logger.info(f'reading {len(EXTRA_COLS)} columns from {EXTRA_FN}')
    zextra = fitsio.read(extra_path, columns=EXTRA_COLS)
    logger.info(f'read {len(zextra)} rows from {EXTRA_FN}')

    # ----- align the two catalogs on TARGETID -------------------------------- #
    zextra = align_to(zcat['TARGETID'], zextra['TARGETID'], zextra)
    n = len(zcat)

    # ----- healpix (HPXPIXEL) ------------------------------------------------ #
    if has_healpix:
        logger.info('using HEALPIX column already present in the input catalog')
        hpxpixel = np.asarray(zcat['HEALPIX'])
    else:
        logger.info(f'computing HEALPIX (nside={HPX_NSIDE}, nested) from RA/DEC')
        hpxpixel = hp.ang2pix(HPX_NSIDE, zcat['TARGET_RA'], zcat['TARGET_DEC'],
                              lonlat=True, nest=True)

    # ----- QuasarNet confidence selections ----------------------------------- #
    max_c = np.max(np.array([zextra[name] for name in C_COLS]), axis=0)
    QN6 = max_c > 0.6
    QN99 = max_c > 0.99
    del max_c

    # ----- target bits ------------------------------------------------------- #
    varbit = targetmask.scnd_mask['WISE_VAR_QSO']
    elgbit = targetmask.desi_mask['ELG']
    bgsbit = targetmask.desi_mask['BGS_ANY']
    qsobit = targetmask.desi_mask['QSO']

    desi_target = np.asarray(zcat['DESI_TARGET'])
    scnd_target = np.asarray(zcat['SCND_TARGET'])
    is_QSO = (desi_target & qsobit) > 0
    is_ELG = (desi_target & elgbit) > 0
    is_BGS = (desi_target & bgsbit) > 0
    is_VAR = (scnd_target & varbit) > 0

    # ----- spectroscopic QSO identifications --------------------------------- #
    spectype_qso = as_str(zextra['SPECTYPE']) == 'QSO'
    is_mgii = np.asarray(zextra['IS_QSO_MGII']).astype(bool)
    is_qn_new_rr = np.asarray(zextra['IS_QSO_QN_NEW_RR']).astype(bool)
    good_z_qso = np.asarray(zextra['GOOD_Z_QSO']).astype(bool)

    is_OK_for_ELG = spectype_qso & QN6
    is_OK_for_BGS = spectype_qso & (QN6 | is_mgii)
    is_OK_for_VAR = spectype_qso | is_mgii | QN99

    # ----- quality cut ------------------------------------------------------- #
    # OBJTYPE != TGT removes e.g. sky fibers (excess around z~3.7 and at low z).
    bad_qso = as_str(zcat['OBJTYPE']) != 'TGT'
    fiberstatus = np.asarray(zcat['COADD_FIBERSTATUS'])
    # keep COADD_FIBERSTATUS == 0 or == 2**3 (see desispec maskbits)
    bad_qso |= ~((fiberstatus == 0) | (fiberstatus == 2**3))

    # ----- "all QSOs": union of four mutually-exclusive subsamples ----------- #
    # 1) QSO targets
    selqso = is_QSO
    member_qso = selqso & good_z_qso & ~bad_qso

    # 2) ELG targets (not QSO) identified as QSO via SPECTYPE & QN>0.6
    selelg = is_ELG & ~is_QSO
    member_elg = selelg & is_OK_for_ELG & ~bad_qso

    # 3) WISE variability secondary targets (not QSO, not ELG, not BGS)
    selvar = is_VAR & ~is_QSO & ~is_ELG & ~is_BGS
    member_var = selvar & is_OK_for_VAR & ~bad_qso

    # 4) BGS targets (not QSO, not ELG, not WISE-VAR)
    selbgs = is_BGS & ~is_QSO & ~is_ELG & ~is_VAR
    member_bgs = selbgs & is_OK_for_BGS & ~bad_qso

    member_all = member_qso | member_elg | member_var | member_bgs
    logger.info(f'identified QSOs: {int(member_qso.sum())} QSO, {int(member_elg.sum())} ELG, '
                f'{int(member_var.sum())} VAR, {int(member_bgs.sum())} BGS '
                f'-> {int(member_all.sum())} total')
    logger.info(f'QSO targets (all): {int(selqso.sum())}')

    # ----- redshift assembly ------------------------------------------------- #
    # zelg: redrock Z with Z_NEW substituted where IS_QSO_QN_NEW_RR
    z_out = np.where(is_qn_new_rr, zextra['Z_NEW'], zextra['Z']).astype(zextra['Z'].dtype)
    # QSO-target subsample uses the dedicated QSO redshift where it is good
    z_out[member_qso] = zextra['Z_QSO'][member_qso]

    # ----- QSO_MASKBITS (canonical bit definition) --------------------------- #
    # decision order: BGS < ELG < QSO ; WISE_VAR_QSO treated like QSO.
    qso_maskbits = np.zeros(n, dtype=np.int32)
    qso_maskbits[is_QSO & spectype_qso] += 2**1
    qso_maskbits[is_QSO & is_mgii] += 2**2
    qso_maskbits[is_QSO & QN99] += 2**3
    qso_maskbits[is_QSO & is_qn_new_rr & QN99] += 2**4
    qso_maskbits[is_BGS & (~is_ELG) & (~is_QSO) & is_OK_for_BGS] += 2**5
    qso_maskbits[is_ELG & (~is_QSO) & is_OK_for_ELG] += 2**6
    qso_maskbits[is_VAR & is_OK_for_VAR] += 2**7
    qso_maskbits[bad_qso] = 0

    # ----- assemble the output table ----------------------------------------- #
    out = Table()
    out['TARGETID'] = zcat['TARGETID']
    out['Z'] = z_out
    out['ZERR'] = zextra['ZERR']
    out['ZWARN'] = zextra['ZWARN']
    out['SPECTYPE'] = zextra['SPECTYPE']
    out['COADD_FIBERSTATUS'] = zcat['COADD_FIBERSTATUS']
    out['TARGET_RA'] = zcat['TARGET_RA']
    out['TARGET_DEC'] = zcat['TARGET_DEC']
    out['OBJTYPE'] = zcat['OBJTYPE']
    out['DESI_TARGET'] = zcat['DESI_TARGET']
    out['SCND_TARGET'] = zcat['SCND_TARGET']
    out['COADD_NUMEXP'] = zcat['COADD_NUMEXP']
    out['COADD_EXPTIME'] = zcat['COADD_EXPTIME']
    out['TSNR2_LYA'] = zextra['TSNR2_LYA']
    out['TSNR2_QSO'] = zextra['TSNR2_QSO']
    for name in C_COLS:
        out[name] = zextra[name]
    out['QSO_MASKBITS'] = qso_maskbits
    out['HPXPIXEL'] = hpxpixel
    out['SURVEY'] = np.full(n, SURVEY)
    out['PROGRAM'] = np.full(n, PROGRAM)
    out['TSNR2_ELG'] = zextra['TSNR2_ELG']
    out['TSNR2_LRG'] = zextra['TSNR2_LRG']

    # keep only the requested columns, in the requested order
    out = out[OUTPUT_COLS]

    # ----- write outputs ----------------------------------------------------- #
    version = determine_version(outdir, args.release, explicit=args.version)
    logger.info(f'using catalog version {version}')

    targets_fn = os.path.join(
        outdir, f'QSO_cat_{args.release}_main_dark_healpix_only_qso_targets_{version}.fits')
    qsos_fn = os.path.join(
        outdir, f'QSO_cat_{args.release}_main_dark_healpix_{version}.fits')

    logger.info(f'writing "QSO targets" catalog ({int(member_qso.sum())} rows) to {targets_fn}')
    common.write_LSS_scratchcp(out[member_qso], targets_fn, extname=EXTNAME, logger=logger)

    logger.info(f'writing "all QSOs" catalog ({int(member_all.sum())} rows) to {qsos_fn}')
    common.write_LSS_scratchcp(out[member_all], qsos_fn, extname=EXTNAME, logger=logger)

    logger.info('done')


if __name__ == '__main__':
    main()
