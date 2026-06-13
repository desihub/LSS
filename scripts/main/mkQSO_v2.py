#!/usr/bin/env python
# coding: utf-8

"""
Build the two main-survey dark-time QSO catalogs (healpix-based) directly from the
zcatalog ``zpix-main-dark.fits``, ``zpix-main-dark-extra.fits``,
``zpix-main-dark-imaging.fits``, and (for v2) ``exp_fibermap/zpix-main-dark-expfibermap.fits``.

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
                          Note that Z_QSO is the QuasarNET corrected redshift
    * ELG / WISE-VAR / BGS : Z = "zelg", i.e. the redrock Z with Z_NEW substituted in
                             wherever IS_QSO_QN_NEW_RR is set; ZERR, ZWARN, and SPECTYPE
                             are likewise substituted with their _NEW counterparts.

Coadd MJDs (``COADD_FIRSTMJD``, ``COADD_LASTMJD``, ``COADD_MEANMJD``) are taken from
``MIN_MJD``, ``MAX_MJD``, and ``MEAN_MJD`` in the main zpix catalog. Coadd nights
(``COADD_FIRSTNIGHT``, ``COADD_LASTNIGHT``) are derived from EXP_FIBERMAP for the
identified QSOs only.
"""

import os
import sys
import getpass
import argparse
import datetime
import logging

import numpy as np
import fitsio
import healpy as hp
from astropy.table import Table, join, unique

from desitarget import targetmask
from LSS import common_tools as common

# ----------------------------------------------------------------------------- #
# logging
SCRIPTNAME = os.path.splitext(os.path.basename(__file__))[0]
logger = logging.getLogger(SCRIPTNAME)
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
             'COADD_EXPTIME',
             'EFFTIME_SPEC',
             'MIN_MJD',
             'MAX_MJD',
             'MEAN_MJD']

# columns to read from zpix-main-dark-extra.fits
EXTRA_COLS = ['TARGETID',
              'Z',
              'ZERR',
              'ZWARN',
              'SPECTYPE',
              'TSNR2_LYA',
              'TSNR2_QSO',
              'TSNR2_ELG',
              'C_LYA',
              'C_CIV',
              'C_CIII',
              'C_MgII',
              'C_Hbeta',
              'C_Halpha',
              'Z_NEW',
              'ZERR_NEW',
              'ZWARN_NEW',
              'SPECTYPE_NEW',
              'GOOD_Z_QSO',
              'Z_QSO',
              'IS_QSO_MGII',
              'IS_QSO_QN_NEW_RR']

# columns to read from zpix-main-dark-imaging.fits
IMAGE_COLS = ['TARGETID',
              'EBV',
              'FLUX_G',
              'FLUX_R',
              'FLUX_Z',
              'FLUX_W1',
              'FLUX_W2',
              'FLUX_IVAR_G',
              'FLUX_IVAR_R',
              'FLUX_IVAR_Z',
              'FLUX_IVAR_W1',
              'FLUX_IVAR_W2']

# imaging columns written to the output (TARGETID used only for matching)
IMAGE_OUT_COLS = [col for col in IMAGE_COLS if col != 'TARGETID']

# afterburner confidence columns used to build the QN selections
C_COLS = ['C_LYA', 'C_CIV', 'C_CIII', 'C_MgII', 'C_Hbeta', 'C_Halpha']

# final output column order (intersection of the requested schema with what the two
# input catalogs can provide)
OUTPUT_COLS = ['TARGETID', 'HPXPIXEL', 'SURVEY', 'PROGRAM',
               'Z', 'ZERR', 'ZWARN', 'SPECTYPE',
               'TARGET_RA', 'TARGET_DEC', 'OBJTYPE', 'EBV',
               'FLUX_G', 'FLUX_R', 'FLUX_Z', 'FLUX_W1', 'FLUX_W2',
               'FLUX_IVAR_G', 'FLUX_IVAR_R', 'FLUX_IVAR_Z', 'FLUX_IVAR_W1', 'FLUX_IVAR_W2',
               'DESI_TARGET', 'SCND_TARGET',
               'COADD_FIBERSTATUS', 'COADD_NUMEXP', 'COADD_EXPTIME', 'EFFTIME_SPEC',
               'TSNR2_LYA', 'TSNR2_QSO', 'TSNR2_ELG',
               'QSO_MASKBITS', 'LASTNIGHT',
               'COADD_FIRSTNIGHT', 'COADD_FIRSTMJD', 'COADD_LASTNIGHT', 'COADD_LASTMJD',
               'COADD_MEANMJD']

# healpix resolution used for the main-survey healpix grouping
HPX_NSIDE = 64

SURVEY = 'main'
PROGRAM = 'dark'

ZCAT_FN = 'zpix-main-dark.fits'
EXTRA_FN = 'zpix-main-dark-extra.fits'
IMAG_FN = 'zpix-main-dark-imaging.fits'
ZTILE_FN = 'ztile-main-dark-cumulative.fits'
EXPFIBERMAP_FN = 'exp_fibermap/zpix-main-dark-expfibermap.fits'

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


def get_expfibermap_path(catdir, zcat_path):
    """Return path to EXP_FIBERMAP data (v2 separate file, or v1 extension in zpix)."""
    exp_path = dvs_ro(os.path.join(catdir, EXPFIBERMAP_FN))
    if os.path.isfile(exp_path):
        return exp_path, 'EXP_FIBERMAP'
    zcat_path = dvs_ro(zcat_path)
    with fitsio.FITS(zcat_path) as f:
        for hdu in f:
            if hdu.get_extname().upper() == 'EXP_FIBERMAP':
                return zcat_path, 'EXP_FIBERMAP'
    raise FileNotFoundError(
        f'EXP_FIBERMAP not found in {exp_path} or as an extension of {zcat_path}')


def add_lastnight(qf, catdir):
    """Add LASTNIGHT from the tiles zcatalog (left join on TARGETID)."""
    tilefn = dvs_ro(os.path.join(catdir, ZTILE_FN))
    logger.info(f'reading LASTNIGHT from {tilefn}')
    t = Table(fitsio.read(tilefn, columns=['TARGETID', 'LASTNIGHT']))
    t.sort('TARGETID')
    qf_tids = qf['TARGETID']
    pos = np.searchsorted(t['TARGETID'], qf_tids)
    pos = np.clip(pos, 0, len(t) - 1)
    match = t['TARGETID'][pos] == qf_tids
    lastnight = np.zeros(len(qf), dtype=t['LASTNIGHT'].dtype)
    lastnight[match] = t['LASTNIGHT'][pos[match]]
    qf['LASTNIGHT'] = lastnight
    n_missing = int(np.sum(lastnight == 0))
    logger.info(f'{n_missing} entries without LASTNIGHT info')
    return qf


def read_exp_nights_for_targets(exp_path, extname, targetids, chunk_size=5_000_000):
    """Read NIGHT from EXP_FIBERMAP, keeping only rows matching ``targetids``."""
    targetids = np.unique(targetids)
    chunks = []
    with fitsio.FITS(exp_path) as f:
        hdu = f[extname]
        nrows = hdu.get_info()['nrows']
        logger.info(f'scanning {nrows} EXP_FIBERMAP rows in chunks of {chunk_size}')
        for start in range(0, nrows, chunk_size):
            end = min(start + chunk_size, nrows)
            block = hdu.read(rows=range(start, end), columns=['TARGETID', 'NIGHT'])
            sel = np.isin(block['TARGETID'], targetids)
            if np.any(sel):
                chunks.append(block[sel])
    if not chunks:
        return Table(np.zeros(0, dtype=[('TARGETID', 'i8'), ('NIGHT', 'i4')]))
    return Table(np.concatenate(chunks))


def add_coadd_nights(qf, catdir, zcat_path):
    """Add COADD_FIRSTNIGHT and COADD_LASTNIGHT from EXP_FIBERMAP for output targets."""
    exp_path, extname = get_expfibermap_path(catdir, zcat_path)
    logger.info(f'reading EXP_FIBERMAP nights for {len(qf)} targets from {exp_path}')
    expinfo = read_exp_nights_for_targets(exp_path, extname, qf['TARGETID'])
    logger.info(f'found {len(expinfo)} exposure rows for output targets')
    if len(expinfo) == 0:
        qf['COADD_FIRSTNIGHT'] = np.zeros(len(qf), dtype=int)
        qf['COADD_LASTNIGHT'] = np.zeros(len(qf), dtype=int)
        return qf
    logger.info('getting FIRST night info')
    expinfo.sort('NIGHT')
    expinfo_first = unique(expinfo, keys=['TARGETID'])
    expinfo_first['NIGHT'].name = 'COADD_FIRSTNIGHT'
    expinfo_first.keep_columns(['TARGETID', 'COADD_FIRSTNIGHT'])
    qf = join(qf, expinfo_first, keys=['TARGETID'], join_type='left')
    del expinfo_first
    logger.info('getting LAST night info')
    expinfo_last = unique(expinfo, keys=['TARGETID'], keep='last')
    expinfo_last['NIGHT'].name = 'COADD_LASTNIGHT'
    expinfo_last.keep_columns(['TARGETID', 'COADD_LASTNIGHT'])
    qf = join(qf, expinfo_last, keys=['TARGETID'], join_type='left')
    return qf


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

    outdir = os.path.abspath(args.outdir)
    if not os.path.isdir(outdir):
        logger.info(f'creating output directory {outdir}')
        os.makedirs(outdir, exist_ok=True)

    version = determine_version(outdir, args.release, explicit=args.version)
    logfile = os.path.join(outdir, f'{SCRIPTNAME}_{args.release}_{version}.log')
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.INFO)
    fh.setFormatter(ch.formatter)
    logger.addHandler(fh)
    logger.info(f'logging to {logfile}')
    logger.info(f'arguments: {args}')
    logger.info(f'using catalog version {version}')

    catdir = get_catdir(args.release)
    zcat_path = dvs_ro(os.path.join(catdir, ZCAT_FN))
    extra_path = dvs_ro(os.path.join(catdir, EXTRA_FN))
    imag_path = dvs_ro(os.path.join(catdir, IMAG_FN))
    logger.info(f'reading inputs from {catdir}')

    # ----- read zcat (basic catalog) ----------------------------------------- #
    zcat_avail = get_colnames(zcat_path)
    read_cols = list(ZCAT_COLS)
    has_healpix = 'HEALPIX' in zcat_avail
    if has_healpix:
        read_cols.append('HEALPIX')
    has_uniqpix = 'UNIQPIX' in zcat_avail
    if has_uniqpix:
        read_cols.append('UNIQPIX')

    logger.info(f'reading {len(read_cols)} columns from {ZCAT_FN}')
    zcat = fitsio.read(zcat_path, columns=read_cols)
    logger.info(f'read {len(zcat)} rows from {ZCAT_FN}')

    # ----- read zextra (extra catalog) --------------------------------------- #
    logger.info(f'reading {len(EXTRA_COLS)} columns from {EXTRA_FN}')
    zextra = fitsio.read(extra_path, columns=EXTRA_COLS)
    logger.info(f'read {len(zextra)} rows from {EXTRA_FN}')

    # ----- read zimag (imaging catalog) -------------------------------------- #
    logger.info(f'reading {len(IMAGE_COLS)} columns from {IMAG_FN}')
    zimag = fitsio.read(imag_path, columns=IMAGE_COLS)
    logger.info(f'read {len(zimag)} rows from {IMAG_FN}')

    # ----- align the catalogs on TARGETID ------------------------------------ #
    zextra = align_to(zcat['TARGETID'], zextra['TARGETID'], zextra)
    zimag = align_to(zcat['TARGETID'], zimag['TARGETID'], zimag)
    n = len(zcat)

    # ----- healpix (HPXPIXEL) ------------------------------------------------ #
    if has_healpix:
        logger.info('using HEALPIX column already present in the input catalog')
        hpxpixel = np.asarray(zcat['HEALPIX'])
    else:
        logger.info(f'computing HEALPIX (nside={HPX_NSIDE}, nested) from RA/DEC')
        hpxpixel = hp.ang2pix(HPX_NSIDE, zcat['TARGET_RA'], zcat['TARGET_DEC'],
                              lonlat=True, nest=True)

    # ----- uniqpix ----------------------------------------------------------- #
    if has_uniqpix:
        logger.info('UNIQPIX column present in the input catalog')
        uniqpix = np.asarray(zcat['UNIQPIX'])

    # ----- QuasarNet confidence selections ----------------------------------- #
    # QN threshold changed from 0.95 to 0.99 starting with loa
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
    # Substitute QN-afterburner rerun values where IS_QSO_QN_NEW_RR is set.
    z_out = np.where(is_qn_new_rr, zextra['Z_NEW'], zextra['Z']).astype(zextra['Z'].dtype)
    zerr_out = np.where(is_qn_new_rr, zextra['ZERR_NEW'], zextra['ZERR'])
    zwarn_out = np.where(is_qn_new_rr, zextra['ZWARN_NEW'], zextra['ZWARN'])
    spectype_out = zextra['SPECTYPE'].copy()
    spectype_out[is_qn_new_rr] = zextra['SPECTYPE_NEW'][is_qn_new_rr]
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

    # ----- subset to identified QSOs before building the output table ---------- #
    idx = np.where(member_all)[0]
    n_out = len(idx)
    member_qso_out = member_qso[idx]
    logger.info(f'building output table for {n_out} identified QSOs')

    # ----- assemble the output table ----------------------------------------- #
    out = Table()
    out['TARGETID'] = zcat['TARGETID'][idx]
    out['Z'] = z_out[idx]
    out['ZERR'] = zerr_out[idx]
    out['ZWARN'] = zwarn_out[idx]
    out['SPECTYPE'] = spectype_out[idx]
    out['COADD_FIBERSTATUS'] = zcat['COADD_FIBERSTATUS'][idx]
    out['TARGET_RA'] = zcat['TARGET_RA'][idx]
    out['TARGET_DEC'] = zcat['TARGET_DEC'][idx]
    out['OBJTYPE'] = zcat['OBJTYPE'][idx]
    for name in IMAGE_OUT_COLS:
        out[name] = zimag[name][idx]
    out['DESI_TARGET'] = zcat['DESI_TARGET'][idx]
    out['SCND_TARGET'] = zcat['SCND_TARGET'][idx]
    out['COADD_NUMEXP'] = zcat['COADD_NUMEXP'][idx]
    out['COADD_EXPTIME'] = zcat['COADD_EXPTIME'][idx]
    out['EFFTIME_SPEC'] = zcat['EFFTIME_SPEC'][idx]
    out['COADD_FIRSTMJD'] = zcat['MIN_MJD'][idx]
    out['COADD_LASTMJD'] = zcat['MAX_MJD'][idx]
    out['COADD_MEANMJD'] = zcat['MEAN_MJD'][idx]
    out['TSNR2_LYA'] = zextra['TSNR2_LYA'][idx]
    out['TSNR2_QSO'] = zextra['TSNR2_QSO'][idx]
    out['TSNR2_ELG'] = zextra['TSNR2_ELG'][idx]
    out['QSO_MASKBITS'] = qso_maskbits[idx]
    out['HPXPIXEL'] = hpxpixel[idx]
    output_cols = list(OUTPUT_COLS)
    if has_uniqpix:
        out['UNIQPIX'] = uniqpix[idx]
        output_cols.insert(output_cols.index('HPXPIXEL'), 'UNIQPIX')
    out['SURVEY'] = np.full(n_out, SURVEY)
    out['PROGRAM'] = np.full(n_out, PROGRAM)

    out = add_lastnight(out, catdir)
    out = add_coadd_nights(out, catdir, zcat_path)

    # keep only the requested columns, in the requested order
    out = out[output_cols]

    # ----- write outputs ----------------------------------------------------- #
    targets_fn = os.path.join(
        outdir, f'QSO_cat_{args.release}_main_dark_healpix_only_qso_targets_{version}.fits')
    qsos_fn = os.path.join(
        outdir, f'QSO_cat_{args.release}_main_dark_healpix_{version}.fits')

    header_comments = [
        f"Created on {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
        f"Created by {getpass.getuser()} with {os.path.basename(__file__)}",
    ]

    logger.info(f'writing "QSO targets" catalog ({int(member_qso_out.sum())} rows) to {targets_fn}')
    common.write_LSS_scratchcp(out[member_qso_out], targets_fn, extname=EXTNAME,
                               comments=header_comments, logger=logger)

    logger.info(f'writing "all QSOs" catalog ({n_out} rows) to {qsos_fn}')
    common.write_LSS_scratchcp(out, qsos_fn, extname=EXTNAME,
                               comments=header_comments, logger=logger)

    logger.info('done')


if __name__ == '__main__':
    main()
