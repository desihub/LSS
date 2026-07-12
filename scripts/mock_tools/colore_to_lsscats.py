#!/usr/bin/env python
# coding: utf-8

"""
Convert CoLoRe-2LPT mock data + randoms FITS catalogs into the LSScats layout
expected by scripts/xirunpc.py.

Input formats (CoLoRe convention):
    * Mock catalog: BINTABLE with at least TARGET_RA, TARGET_DEC, Z, ZWARN,
      TARGETID. Often in extension ``ZCATALOG`` (e.g. ``zcat.fits``); some
      variants use an unnamed table at HDU 1 (e.g. ``zcat_gauss_400.fits``).
      The mock extension is auto-detected unless ``--mock-ext`` is set.
    * Randoms catalog: BINTABLE extension ``CATALOG`` with at least
      RA, DEC, Z, MOCKID (and IS_QSO_TARGET if using ``--qsos-only``).

Output files (default tracer ``QSO``)::
    {outdir}/QSO_NGC_clustering.dat.fits
    {outdir}/QSO_SGC_clustering.dat.fits
    {outdir}/QSO_NGC_{0..nran-1}_clustering.ran.fits
    {outdir}/QSO_SGC_{0..nran-1}_clustering.ran.fits

Random clustering files are written once; if they already exist in ``outdir``,
they are not regenerated (use ``--force-randoms`` to overwrite). By default,
randoms are subsampled to ``--randoms-factor`` times the mock data count
(default 10). Data files are always regenerated from the mock catalog.

Example::

    python scripts/mock_tools/colore_to_lsscats.py \\
        --mock /path/to/mock_zcatalog.fits \\
        --randoms /path/to/randoms.fits \\
        --outdir /path/to/LSScats \\
        --nran 4

Then run xirunpc.py with ``--basedir /path/to/LSScats``.
"""

import argparse
import os

import fitsio
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u

from LSS import common_tools as common


MOCK_REQUIRED_COLS = {'TARGET_RA', 'TARGET_DEC', 'Z', 'TARGETID'}
MOCK_PREFERRED_EXTS = ('ZCATALOG',)


def parse_fits_ext(value):
    """Parse a FITS extension argument as an HDU index or EXTNAME."""
    try:
        return int(value)
    except (TypeError, ValueError):
        return value


def find_table_extension(path, preferred_exts=(), required_cols=()):
    """Return EXTNAME or HDU index for the first matching BINTABLE."""
    with fitsio.FITS(path) as fits:
        for ext in preferred_exts:
            try:
                fits[ext]
            except (KeyError, OSError, ValueError):
                continue
            return ext

        for i, hdu in enumerate(fits):
            hdr = hdu.read_header()
            if hdr.get('XTENSION') != 'BINTABLE':
                continue
            cols = {c.upper() for c in hdu.get_colnames()}
            if required_cols and not required_cols.issubset(cols):
                continue
            extname = hdr.get('EXTNAME', '')
            return extname if extname else i

    raise ValueError(f'no matching BINTABLE found in {path}')


def resolve_mock_extension(path, mock_ext=None):
    """Choose the mock catalog HDU, auto-detecting when ``mock_ext`` is None."""
    if mock_ext is not None:
        return mock_ext
    ext = find_table_extension(
        path, preferred_exts=MOCK_PREFERRED_EXTS, required_cols=MOCK_REQUIRED_COLS)
    _log(f'mock extension auto-detected: {ext!r}')
    return ext


def resolve_randoms_extension(path, randoms_ext=None):
    """Choose the randoms catalog HDU, auto-detecting when ``randoms_ext`` is None."""
    if randoms_ext is not None:
        return randoms_ext
    ext = find_table_extension(path, preferred_exts=('CATALOG',), required_cols={'RA', 'DEC', 'Z', 'MOCKID'})
    _log(f'randoms extension auto-detected: {ext!r}')
    return ext


def _log(msg, logger=None):
    if logger is None:
        print(msg)
    else:
        logger.info(msg)


def read_mock_table(path, ext=None, columns=None):
    """Read the CoLoRe mock catalog."""
    if columns is None:
        columns = ['TARGET_RA', 'TARGET_DEC', 'Z', 'ZERR', 'ZWARN', 'TARGETID']
    ext = resolve_mock_extension(path, ext)
    _log(f'reading mock catalog from {path} [{ext}]')
    return Table(fitsio.read(path, ext=ext, columns=columns))


def read_randoms_table(path, ext=None, qsos_only=False):
    """Read the CoLoRe randoms catalog."""
    columns = ['RA', 'DEC', 'Z', 'MOCKID']
    if qsos_only:
        columns.append('IS_QSO_TARGET')
    ext = resolve_randoms_extension(path, ext)
    _log(f'reading randoms catalog from {path} [{ext}]')
    return Table(fitsio.read(path, ext=ext, columns=columns))


def mock_selection(tab, zmin=None, zmax=None, ismock=True):
    """Apply QSO mock redshift / TARGETID quality cuts (mkclusdat-style)."""
    sel = np.isfinite(tab['Z'])
    sel &= tab['Z'] != 999999
    sel &= tab['Z'] != 1.e20
    if 'ZWARN' in tab.colnames:
        sel &= tab['ZWARN'] != 999999
    if ismock:
        sel &= tab['TARGETID'] < 419430400000000
    if zmin is not None:
        sel &= tab['Z'] >= zmin
    if zmax is not None:
        sel &= tab['Z'] < zmax
    return sel


def randoms_selection(tab, zmin=None, zmax=None, qsos_only=False):
    """Keep randoms with valid redshifts; optionally require IS_QSO_TARGET."""
    sel = np.isfinite(tab['Z'])
    if qsos_only:
        sel &= tab['IS_QSO_TARGET'].astype(bool)
    if zmin is not None:
        sel &= tab['Z'] >= zmin
    if zmax is not None:
        sel &= tab['Z'] < zmax
    return sel


def to_clustering_table(tab, ra_col='RA', dec_col='DEC', id_col='TARGETID'):
    """Build a minimal clustering catalog table for xirunpc (--weight_type default)."""
    out = Table()
    out['RA'] = np.asarray(tab[ra_col], dtype=np.float64)
    out['DEC'] = np.asarray(tab[dec_col], dtype=np.float64)
    out['Z'] = np.asarray(tab['Z'], dtype=np.float64)
    out['TARGETID'] = np.asarray(tab[id_col], dtype=np.int64)
    out['WEIGHT'] = np.ones(len(out), dtype=np.float64)
    return out


def galactic_caps(ra, dec):
    """Return boolean mask for NGC (galactic b > 0), matching splitclusGC."""
    c = SkyCoord(ra * u.deg, dec * u.deg, frame='icrs')
    gc = c.transform_to('galactic')
    return gc.b.value > 0


def downsample_randoms(tab, n_data, factor, seed):
    """Uniformly subsample randoms to ``factor`` times the data count."""
    if factor is None or factor <= 0:
        _log(f'randoms: keeping all {len(tab)} rows (downsampling disabled)')
        return tab
    target_n = int(factor * n_data)
    if target_n >= len(tab):
        _log(f'randoms: keeping all {len(tab)} rows '
             f'(target {target_n} = {factor:g}x{n_data} data)')
        return tab
    rng = np.random.default_rng(seed)
    keep = np.sort(rng.choice(len(tab), size=target_n, replace=False))
    _log(f'randoms: subsampled {len(tab)} -> {target_n} rows ({factor:g}x {n_data} data)')
    return tab[keep]


def assign_random_realizations(n, nran, seed=42):
    """Randomly assign each random to a realization index in [0, nran)."""
    rng = np.random.default_rng(seed)
    return rng.integers(0, nran, size=n).astype(np.int32)


def random_output_paths(tracer, outdir, nran):
    """Return paths to all random clustering FITS files for this tracer."""
    paths = []
    for cap in ('NGC', 'SGC'):
        for iran in range(nran):
            paths.append(os.path.join(outdir, f'{tracer}_{cap}_{iran}_clustering.ran.fits'))
    return paths


def randoms_exist(tracer, outdir, nran):
    """True if every expected random clustering file is already present."""
    return all(os.path.isfile(p) for p in random_output_paths(tracer, outdir, nran))


def prepare_randoms(randoms_path, randoms_ext, qsos_only, zmin, zmax, nran,
                    seed, n_data, randoms_factor):
    """Read, select, downsample, split, and return random clustering tables."""
    randoms = read_randoms_table(randoms_path, ext=randoms_ext, qsos_only=qsos_only)
    rsel = randoms_selection(randoms, zmin=zmin, zmax=zmax, qsos_only=qsos_only)
    _log(f'randoms: kept {int(rsel.sum())} / {len(randoms)} rows after selection')
    randoms = randoms[rsel]
    randoms = downsample_randoms(randoms, n_data, randoms_factor, seed)

    ran = to_clustering_table(randoms, ra_col='RA', dec_col='DEC', id_col='MOCKID')
    iran = assign_random_realizations(len(randoms), nran, seed=seed + 1)
    ran_by_realization = []
    ngc_mask = galactic_caps(ran['RA'], ran['DEC'])
    for i in range(nran):
        in_real = iran == i
        ran_i = ran[in_real]
        ran_by_realization.append((ran_i[ngc_mask[in_real]], ran_i[~ngc_mask[in_real]]))
        _log(f'random realization {i}: {int(in_real.sum())} total, '
             f'{len(ran_by_realization[-1][0])} NGC, {len(ran_by_realization[-1][1])} SGC')
    return ran_by_realization


def write_data_files(data_ngc, data_sgc, tracer, outdir, logger=None):
    """Write mock data clustering FITS files."""
    os.makedirs(outdir, exist_ok=True)
    comments = ['Clustering catalog produced by colore_to_lsscats.py']
    for cap, tab in [('NGC', data_ngc), ('SGC', data_sgc)]:
        outf = os.path.join(outdir, f'{tracer}_{cap}_clustering.dat.fits')
        _log(f'writing {len(tab)} data rows to {outf}', logger)
        common.write_LSS(tab, outf, comments=comments)


def write_random_files(ran_by_realization, tracer, outdir, logger=None):
    """Write random clustering FITS files."""
    os.makedirs(outdir, exist_ok=True)
    comments = ['Clustering catalog produced by colore_to_lsscats.py']
    for iran, (ran_ngc, ran_sgc) in enumerate(ran_by_realization):
        for cap, tab in [('NGC', ran_ngc), ('SGC', ran_sgc)]:
            outf = os.path.join(outdir, f'{tracer}_{cap}_{iran}_clustering.ran.fits')
            _log(f'writing {len(tab)} random rows (realization {iran}, {cap}) to {outf}', logger)
            common.write_LSS(tab, outf, comments=comments)


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--mock', required=True,
                        help='path to CoLoRe mock FITS (with ZCATALOG or --mock-ext extension)')
    parser.add_argument('--randoms', default=None,
                        help='path to CoLoRe randoms FITS (required unless random outputs exist)')
    parser.add_argument('--outdir', required=True,
                        help='output directory for LSScats clustering files')
    parser.add_argument('--tracer', default='QSO',
                        help='tracer name used in output filenames (default: QSO)')
    parser.add_argument('--nran', type=int, default=4,
                        help='number of random catalog realizations (default: 4)')
    parser.add_argument('--mock-ext', default=None, type=parse_fits_ext,
                        help='FITS extension name or HDU index for mock catalog '
                             '(default: auto-detect ZCATALOG or first table with mock columns)')
    parser.add_argument('--randoms-ext', default=None, type=parse_fits_ext,
                        help='FITS extension name or HDU index for randoms catalog '
                             '(default: auto-detect CATALOG or first matching table)')
    parser.add_argument('--zmin', type=float, default=None,
                        help='optional minimum redshift cut applied at prep time')
    parser.add_argument('--zmax', type=float, default=None,
                        help='optional maximum redshift cut applied at prep time')
    parser.add_argument('--no-mock-targetid-cut', action='store_true',
                        help='do not apply mock TARGETID < 419430400000000 cut')
    parser.add_argument('--qsos-only', action='store_true',
                        help='keep only randoms with IS_QSO_TARGET set')
    parser.add_argument('--seed', type=int, default=1216,
                        help='RNG seed for randoms subsampling and nran split (default: 42)')
    parser.add_argument('--randoms-factor', type=float, default=10,
                        help='subsample randoms to this many times the mock data count '
                             '(default: 10; <=0 keeps all randoms after selection)')
    parser.add_argument('--force-randoms', action='store_true',
                        help='regenerate random clustering files even if they already exist')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    skip_randoms = randoms_exist(args.tracer, args.outdir, args.nran) and not args.force_randoms

    if skip_randoms:
        _log('random clustering files already exist; skipping random catalog processing')
    else:
        if args.randoms is None:
            parser.error('--randoms is required when random clustering files are missing '
                         '(or pass --force-randoms with --randoms to overwrite)')

    mock = read_mock_table(args.mock, extname=args.mock_ext)
    sel = mock_selection(mock, zmin=args.zmin, zmax=args.zmax,
                         ismock=not args.no_mock_targetid_cut)
    _log(f'mock: kept {int(sel.sum())} / {len(mock)} rows after selection')
    mock = mock[sel]

    data = to_clustering_table(mock, ra_col='TARGET_RA', dec_col='TARGET_DEC',
                               id_col='TARGETID')
    data_ngc = galactic_caps(data['RA'], data['DEC'])
    ngc_data = data[data_ngc]
    sgc_data = data[~data_ngc]
    _log(f'data split: {len(ngc_data)} NGC, {len(sgc_data)} SGC')

    write_data_files(ngc_data, sgc_data, args.tracer, args.outdir)

    if not skip_randoms:
        ran_by_realization = prepare_randoms(
            args.randoms, args.randoms_ext, args.qsos_only,
            args.zmin, args.zmax, args.nran, args.seed,
            n_data=len(data), randoms_factor=args.randoms_factor)
        write_random_files(ran_by_realization, args.tracer, args.outdir)

    _log('done')
    _log(f'run xirunpc.py with --basedir {os.path.abspath(args.outdir)} '
         f'--tracer {args.tracer} --region NGC SGC')


if __name__ == '__main__':
    main()
